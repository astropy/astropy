# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""

import os
import re
import shutil
from collections.abc import Callable, Generator
from contextlib import contextmanager, nullcontext
from copy import deepcopy
from dataclasses import KW_ONLY, dataclass, field, replace
from enum import Enum, auto
from functools import wraps
from inspect import cleandoc
from pathlib import Path
from threading import RLock
from types import TracebackType
from typing import (
    Literal,
    NewType,
    ParamSpec,
    TypeAlias,
    TypedDict,
    assert_never,
    final,
)
from warnings import warn

from astropy.config._tempfile_shim import TemporaryDirectory
from astropy.utils.exceptions import AstropyUserWarning

__all__ = [
    "get_cache_dir",
    "get_cache_dir_path",
    "get_config_dir",
    "get_config_dir_path",
    "set_temp_cache",
    "set_temp_config",
    "temporary_cache_dir_path",
    "temporary_config_dir_path",
]

P = ParamSpec("P")
PathGetter: TypeAlias = Callable[[str], Path]


class _DirType(Enum):
    CACHE = auto()
    CONFIG = auto()

    @property
    def path_getter(self) -> PathGetter:
        match self:
            case _DirType.CACHE:
                return get_cache_dir_path
            case _DirType.CONFIG:
                return get_config_dir_path
            case _ as unreachable:
                assert_never(unreachable)

    @property
    def legacy_context_manager(self) -> "type[_SetTempPath]":
        match self:
            case _DirType.CACHE:
                return set_temp_cache
            case _DirType.CONFIG:
                return set_temp_config
            case _ as unreachable:
                assert_never(unreachable)

    def as_str(self) -> str:
        return self.name.lower()


@final
@dataclass(slots=True, frozen=True, kw_only=True)
class _Namespace:
    root: str
    fragments: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not all(re.fullmatch(r"[\w-]+", e) for e in self.elements):
            raise ValueError("Found invalid namespace elements.")

    @classmethod
    def from_str(cls, s: str, /) -> "_Namespace":
        elements = s.split(".")
        return _Namespace(root=elements[0], fragments=tuple(elements[1:]))

    @property
    def elements(self) -> list[str]:
        return [self.root, *self.fragments]

    def join(self) -> str:
        # should round-trip with from_str
        return ".".join(self.elements)


@final
@dataclass(slots=True, frozen=True, kw_only=True)
class _DirectoryElements:
    base_node: Path | None = None
    sub_nodes: list[str] = field(default_factory=list)

    def join(self) -> Path:
        if self.base_node is None:
            raise AssertionError
        return self.base_node.joinpath(*self.sub_nodes)


# Any number of threads may call any combination of set_temp_cache and set_temp_config,
# in any order. Hence they need to share a single lock (otherwise deadlocks would be possible).
# In particular, the possibility of nesting more than one instance of each context imposes
# the use of a re-entrant lock (RLock).
_PATHS_LOCK = RLock()


Error = NewType("Error", str)


def _resolve(var: str) -> Path | Error:
    # if a Path is returned, it is guaranteed to:
    # - exist
    # - be absolute (in conformance with the XDG standard)
    if (s := os.getenv(var)) is None:
        raise AssertionError

    if not (p := Path(s)).exists():
        return Error(f"{var} is set to {s}, but no such file or directory was found.")

    if p.is_file():
        return Error(f"{var} is set to {s}, which is a file (expected a directory).")

    assert p.is_dir()  # redundancy, fine if skipped in prod

    if not p.is_absolute():
        return Error(
            f"{var} is set to {s}, which is relative (expected an absolute path)."
        )

    return p


@final
@dataclass(slots=True, frozen=True)
class _AstropyEnvvar:
    dirtype: _DirType

    @property
    def name(self) -> str:
        match self.dirtype:
            case _DirType.CACHE:
                return "ASTROPY_CACHE_DIR"
            case _DirType.CONFIG:
                return "ASTROPY_CONFIG_DIR"
            case _ as unreachable:
                assert_never(unreachable)

    def is_defined(self) -> bool:
        return self.name in os.environ

    def resolve(self) -> Path | Error:
        assert self.is_defined()
        return _resolve(self.name)


@final
@dataclass(slots=True, frozen=True)
class _XDGEnvvar:
    dirtype: _DirType

    @property
    def name(self) -> str:
        match self.dirtype:
            case _DirType.CACHE:
                return "XDG_CACHE_HOME"
            case _DirType.CONFIG:
                return "XDG_CONFIG_HOME"
            case _ as unreachable:
                assert_never(unreachable)

    @property
    def value(self) -> str | None:
        return os.getenv(self.name)

    def is_defined(self) -> bool:
        return self.value is not None

    def resolve(self) -> Path | Error:
        assert self.is_defined()
        match _resolve(self.name):
            case str() as err:
                return Error(f"{err} This environment variable will be ignored.")
            case Path() as p:
                return p
            case _ as unreachable:
                assert_never(unreachable)


@final
@dataclass(slots=True, frozen=True, kw_only=True)
class _DirectoryEnvvarSet:
    astropy: _AstropyEnvvar
    xdg: _XDGEnvvar


@final
@dataclass(slots=True, frozen=True)
class _DirectoryFinder:
    dirtype: _DirType
    _: KW_ONLY
    # in Python 3.15+, this should be a frozendict
    # (i.e., don't rely on this attribute's being mutable)
    overrides: dict[_Namespace, Path] = field(default_factory=dict)

    @property
    def envvars(self) -> _DirectoryEnvvarSet:
        return _DirectoryEnvvarSet(
            astropy=_AstropyEnvvar(self.dirtype),
            xdg=_XDGEnvvar(self.dirtype),
        )

    def default_base_node(self) -> Path:
        # default to a XDG-compliant scheme
        return Path.home() / f".{self.dirtype.as_str()}"

    def find_directory_elements(self, namespace: str) -> _DirectoryElements:
        ns = _Namespace.from_str(namespace)
        de = _DirectoryElements(sub_nodes=ns.elements)

        if ns in self.overrides:
            return replace(de, base_node=self.overrides[ns])

        if self.envvars.astropy.is_defined() and ns.root == "astropy":
            match self.envvars.astropy.resolve():
                case Path() as base_node:
                    return replace(de, base_node=base_node, sub_nodes=ns.fragments)
                case str() as err:
                    warn(err, AstropyUserWarning, stacklevel=2)
                case _ as unreachable:
                    assert_never(unreachable)

        if self.envvars.xdg.is_defined():
            match self.envvars.xdg.resolve():
                case Path() as base_node:
                    return replace(de, base_node=base_node)
                case str() as err:
                    warn(err, AstropyUserWarning, stacklevel=2)
                case _ as unreachable:
                    assert_never(unreachable)

        return replace(de, base_node=self.default_base_node())

    def find_namespaced_node(self, namespace: str) -> Path:
        return self.find_directory_elements(namespace).join()


class _TempDirKwargs(TypedDict):
    dir: os.PathLike[str] | str | None
    suffix: str | None
    prefix: str | None
    delete: bool


@contextmanager
def _temporary_dir_ctx(
    *,
    directory_finder_name: Literal["_CacheFinder", "_ConfigFinder"],
    ctx_gen,
    namespace: str,
    dir_: os.PathLike[str] | str | None = None,
    suffix: str | None = None,
    prefix: str | None = None,
    delete: bool = True,
) -> Generator[Path, None, None]:
    # the common implementation for temporary_cache_dir_path and temporary_config_dir_path

    kwargs: _TempDirKwargs = {
        "dir": dir_,
        "suffix": suffix,
        "prefix": prefix,
        "delete": delete,
    }
    with _PATHS_LOCK, TemporaryDirectory(**kwargs) as tmp_dir, ctx_gen():
        tmp_path = Path(tmp_dir)
        initial_df: _DirectoryFinder = globals()[directory_finder_name]
        df: _DirectoryFinder = replace(
            initial_df,
            overrides=initial_df.overrides | {_Namespace.from_str(namespace): tmp_path},
        )
        globals()[directory_finder_name] = df

        yield_val = df.find_namespaced_node(namespace)
        if yield_val != tmp_path:
            assert yield_val.is_relative_to(tmp_path)
            yield_val.mkdir(parents=True)

        try:
            yield yield_val
        finally:
            globals()[directory_finder_name] = initial_df


@contextmanager
def _clear_cfgobjs():
    from .configuration import _cfgobjs

    with _PATHS_LOCK:
        initial_cfgobjs = deepcopy(_cfgobjs)
        _cfgobjs.clear()

        yield

        _cfgobjs.clear()
        _cfgobjs.update(initial_cfgobjs)


_CacheFinder = _DirectoryFinder(_DirType.CACHE)


@contextmanager
def temporary_cache_dir_path(
    dir_: os.PathLike[str] | str | None = None,
    /,
    *,
    namespace: str,
    suffix: str | None = None,
    prefix: str | None = None,
    delete: bool = True,
) -> Generator[Path, None, None]:
    """
    A context manager to create a temporary directory and define it as the cache
    directory associated with some namespace.

    The directory returned will shadow bot ASTROPY_CACHE_DIR and XDG_CACHE_HOME
    environment variable if defined.

    This may also be used as a decorator on a function to set the cache path
    just within that function.

    Thread safety is guaranteed, but concurrency isn't:
    only a single thread at a time may execute code within this context.

    .. versionadded:: 8.0.0

    Parameters
    ----------
    namespace: str, keyword-only, mandatory
      a unique identifier for the namespace associated to a temporary directory,
      which will used to name the directory itself.
      This string must be none-empty, and can only contain alphanumeric characters,
      ``_``, ``-`` or ``.``.
      ``.`` is special cased to represent a path separator (see ``os.sep``) in the
      output directory.

    All arguments from :ref:`tempfile.TemporaryDirectory` are optionally supported,
    with a couple differences:
    - ``dir`` is positional-only and must come first
    - all other arguments are keyword-only
    - ``delete`` is supported even on Python 3.11

    See Also
    --------
    temporary_config_dir_path: a similar function for configuration directories
    set_temp_cache: a legacy function with similar goals but a much less predictable behavior
    """
    with _temporary_dir_ctx(
        namespace=namespace,
        dir_=dir_,
        suffix=suffix,
        prefix=prefix,
        delete=delete,
        directory_finder_name="_CacheFinder",
        ctx_gen=nullcontext,
    ) as tmp_path:
        yield tmp_path


_ConfigFinder = _DirectoryFinder(_DirType.CONFIG)


@contextmanager
def temporary_config_dir_path(
    dir_: os.PathLike[str] | str | None = None,
    /,
    *,
    namespace: str,
    suffix: str | None = None,
    prefix: str | None = None,
    delete: bool = True,
) -> Generator[Path, None, None]:
    """
    A context manager to create a temporary directory and define it as the configuration
    directory associated with some namespace.

    The directory returned will shadow bot ASTROPY_CONFIG_DIR and XDG_CONFIG_HOME
    environment variable if defined.

    This may also be used as a decorator on a function to set the config path
    just within that function.

    Thread safety is guaranteed, but concurrency isn't:
    only a single thread at a time may execute code within this context.

    .. versionadded:: 8.0.0

    Parameters
    ----------
    namespace: str, keyword-only, mandatory
      a unique identifier for the namespace associated to a temporary directory,
      which will used to name the directory itself.
      This string must be none-empty, and can only contain alphanumeric characters,
      ``_``, ``-`` or ``.``.
      ``.`` is special cased to represent a path separator (see ``os.sep``) in the
      output directory.

    All arguments from :ref:`tempfile.TemporaryDirectory` are optionally supported,
    with a couple differences:
    - ``dir`` is positional-only and must come first
    - all other arguments are keyword-only
    - ``delete`` is supported even on Python 3.11

    See Also
    --------
    temporary_cache_dir_path: a similar function for cache directories
    set_temp_config: a legacy function with similar goals but a much less predictable behavior
    """
    with _temporary_dir_ctx(
        namespace=namespace,
        dir_=dir_,
        suffix=suffix,
        prefix=prefix,
        delete=delete,
        directory_finder_name="_ConfigFinder",
        ctx_gen=_clear_cfgobjs,
    ) as tmp_path:
        yield tmp_path


def get_config_dir_path(rootname: str = "astropy") -> Path:
    """
    Determines the configuration directory associated with a namespace and creates the
    directory if it doesn't exist.

    This directory is typically ``$XDG_CONFIG_HOME/<namespace>``, but can be overwritten
    with the ``ASTROPY_CONFIG_DIR`` environment variable, or with
    ``temporary_config_dir_path``.

    .. versionchanged: 8.0.0
      in previous versions, the return value pointed to ``$HOME/.astropy/config`` by default

    .. versionchanged: 8.0.0
      added support for ``ASTROPY_CONFIG_DIR``

    .. versionchanged: 8.0.0
      no symlinks are ever created

    Parameters
    ----------
    rootname : str, optional (default: 'astropy')
        Namespace of the associated with the directory. For example, for ``'pkgname'``,
        the directory would be ``$XDG_CONFIG_HOME/pkgname``

    Returns
    -------
    configdir : Path
        The absolute path to the configuration directory.

    """
    node = _ConfigFinder.find_namespaced_node(rootname)
    node.mkdir(parents=True, exist_ok=True)
    return node


def get_config_dir(rootname: str = "astropy") -> str:
    return str(get_config_dir_path(rootname))


if get_config_dir_path.__doc__ is not None:
    # guard against PYTHONOPTIMIZE mode
    get_config_dir.__doc__ = cleandoc(
        get_config_dir_path.__doc__
        + """
    See Also
    --------
    get_config_dir_path : same as this function, except that the return value is a pathlib.Path

    """
    )


def get_cache_dir_path(rootname: str = "astropy") -> Path:
    """
    Determines the cache directory associated with a namespace and creates the
    directory if it doesn't exist.

    This directory is typically ``$XDG_CACHE_HOME/<namespace>``, but can be overwritten
    with the ``ASTROPY_CACHE_DIR`` environment variable, or with
    ``temporary_cache_dir_path``.

    .. versionchanged: 8.0.0
      in previous versions, the return value pointed to ``$HOME/.astropy/cache`` by default

    .. versionchanged: 8.0.0
      added support for ``ASTROPY_CACHE_DIR``

    .. versionchanged: 8.0.0
      no symlinks are ever created

    Parameters
    ----------
    rootname : str, optional (default: 'astropy')
        Namespace of the associated with the directory. For example, for ``'pkgname'``,
        the directory would be ``$XDG_CACHE_HOME/pkgname``

    Returns
    -------
    cachedir : Path
        The absolute path to the cache directory.

    """
    node = _CacheFinder.find_namespaced_node(rootname)
    node.mkdir(parents=True, exist_ok=True)
    return node


def get_cache_dir(rootname: str = "astropy") -> str:
    return str(get_cache_dir_path(rootname))


if get_cache_dir_path.__doc__ is not None:
    # guard against PYTHONOPTIMIZE mode
    get_cache_dir.__doc__ = cleandoc(
        get_cache_dir_path.__doc__
        + """
    See Also
    --------
    get_cache_dir_path : same as this function, except that the return value is a pathlib.Path

    """
    )


class _SetTempPath:
    # This base class serves as a deduplication layer for its only two intended
    # children (set_temp_cache and set_temp_config)
    _dirtype: _DirType

    def __init__(
        self, path: os.PathLike[str] | str | None = None, delete: bool = False
    ) -> None:
        if path is None:
            self._path = None
        else:
            if (path := Path(path)).is_file():
                raise FileExistsError(
                    f"{path} is a file and cannot be used as a directory"
                )
            self._path = path.resolve()

        self._delete = delete

    @classmethod
    def _get_directory_finder_name(cls) -> str:
        match cls._dirtype:
            case _DirType.CACHE:
                return "_CacheFinder"
            case _DirType.CONFIG:
                return "_ConfigFinder"
            case _ as unreachable:
                assert_never(unreachable)

    @classmethod
    def _get_directory_finder(cls) -> _DirectoryFinder:
        return globals()[cls._get_directory_finder_name()]

    @classmethod
    def _get_current_override(cls) -> Path | None:
        return cls._get_directory_finder().overrides.get(_Namespace(root="astropy"))

    def __enter__(self) -> str:
        _PATHS_LOCK.acquire()
        self._prev_path = self._get_current_override()
        try:
            finder_name = self._get_directory_finder_name()
            initial_df: _DirectoryFinder = globals()[finder_name]
            self._prev_df = initial_df
            ns = _Namespace(root="astropy")
            if self._path is None:
                overrides = initial_df.overrides.copy()
                if ns in overrides:
                    del overrides[ns]
            else:
                overrides = initial_df.overrides | {ns: self._path}

            df: _DirectoryFinder = replace(initial_df, overrides=overrides)
        except Exception:
            _PATHS_LOCK.release()
            raise
        try:
            globals()[finder_name] = df
            return str(df.find_namespaced_node("astropy"))
        except Exception:
            globals()[finder_name] = initial_df
            _PATHS_LOCK.release()
            raise

    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        finder_name = self._get_directory_finder_name()

        try:
            if self._delete and self._path is not None:
                shutil.rmtree(self._path)
        finally:
            globals()[finder_name] = self._prev_df
            _PATHS_LOCK.release()

    def __call__(self, func: Callable[P, object]) -> Callable[P, None]:
        """Implements use as a decorator."""

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> None:
            with self:
                func(*args, **kwargs)

        return wrapper


class set_temp_config(_SetTempPath):
    """
    Context manager to set a temporary path for the Astropy config, primarily
    for use with testing.

    If the path set by this context manager does not already exist it will be
    created, if possible.

    This may also be used as a decorator on a function to set the config path
    just within that function.

    Thread safety is guaranteed since astropy 7.2.1, but concurrency isn't:
    only a single thread at a time may execute code within this context.

    .. versionchanged: 8.0.0
       This function is soft-deprecated. It won't emit deprecation warnings but its
       use is discouraged for new code, as the exact behavior is hard to predict.
       Prefer ref:`temporary_config_dir_path` where available.

    Parameters
    ----------
    path : str, optional
        The directory (which must exist) in which to find the Astropy config
        files, or create them if they do not already exist.  If None, this
        restores the config path to the user's default config path as returned
        by `get_config_dir` as though this context manager were not in effect
        (this is useful for testing).  In this case the ``delete`` argument is
        always ignored.

    delete : bool, optional
        If True, cleans up the temporary directory after exiting the temp
        context (default: False).

    See Also
    --------
    temporary_config_dir_path: a function with similar goals but a much more predictable behavior
    """

    _dirtype = _DirType.CONFIG

    def __enter__(self) -> str:
        # Special case for the config case, where we need to reset all the
        # cached config objects.  We do keep the cache, since some of it
        # may have been set programmatically rather than be stored in the
        # config file (e.g., iers.conf.auto_download=False for our tests).
        from .configuration import _cfgobjs

        _PATHS_LOCK.acquire()

        try:
            _cfgobjs_copy = _cfgobjs.copy()
        except Exception:
            _PATHS_LOCK.release()
            raise

        self._cfgobjs_copy = _cfgobjs_copy

        try:
            _cfgobjs.clear()
            return super().__enter__()
        except Exception:
            _cfgobjs.update(self._cfgobjs_copy)
            del self._cfgobjs_copy
            _PATHS_LOCK.release()
            raise

    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        try:
            from .configuration import _cfgobjs

            _cfgobjs.clear()
            _cfgobjs.update(self._cfgobjs_copy)
            del self._cfgobjs_copy
            super().__exit__(type, value, tb)
        finally:
            _PATHS_LOCK.release()


class set_temp_cache(_SetTempPath):
    """
    Context manager to set a temporary path for the Astropy download cache,
    primarily for use with testing (though there may be other applications
    for setting a different cache directory, for example to switch to a cache
    dedicated to large files).

    If the path set by this context manager does not already exist it will be
    created, if possible.

    This may also be used as a decorator on a function to set the cache path
    just within that function.

    Thread safety is guaranteed since astropy 7.2.1, but concurrency isn't:
    only a single thread at a time may execute code within this context.

    .. versionchanged: 8.0.0
       This function is soft-deprecated. It won't emit deprecation warnings but its
       use is discouraged for new code, as the exact behavior is hard to predict.
       Prefer :ref:`temporary_cache_dir_path` where available.

    Parameters
    ----------
    path : str
        The directory (which must exist) in which to find the Astropy cache
        files, or create them if they do not already exist.  If None, this
        restores the cache path to the user's default cache path as returned
        by `get_cache_dir` as though this context manager were not in effect
        (this is useful for testing).  In this case the ``delete`` argument is
        always ignored.

    delete : bool, optional
        If True, cleans up the temporary directory after exiting the temp
        context (default: False).

    See Also
    --------
    temporary_cache_dir_path: a function with similar goals but a much more predictable behavior
    """

    _dirtype = _DirType.CACHE
