# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""

import os
import shutil
import sys
from collections.abc import Callable, Generator
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import KW_ONLY, dataclass, replace
from enum import Enum, auto
from functools import partial, wraps
from inspect import cleandoc
from pathlib import Path
from threading import RLock
from types import TracebackType
from typing import (
    Literal,
    NewType,
    ParamSpec,
    Protocol,
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
    "temporary_cache_path",
    "temporary_config_path",
]

P = ParamSpec("P")
PathGetter: TypeAlias = Callable[[], Path]


class _DirType(Enum):
    CACHE = auto()
    CONFIG = auto()

    def as_str(self) -> str:
        return self.name.lower()

    def specialize_node(self, node: Path) -> Path:
        match self:
            case _DirType.CACHE:
                return node / "cache"
            case _DirType.CONFIG:
                return node
            case _ as unreachable:
                assert_never(unreachable)

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


# Any number of threads may call any combination of set_temp_cache and set_temp_config,
# in any order. Hence they need to share a single lock (otherwise deadlocks would be possible).
# In particular, the possibility of nesting more than one instance of each context imposes
# the use of a re-entrant lock (RLock).
_PATHS_MUTEX = RLock()


Error = NewType("Error", str)


def _resolve(var: str) -> Path | Error:
    # if a Path is returned, it is guaranteed to:
    # - exist
    # - be absolute (in conformance with the XDG standard)
    if (s := os.getenv(var)) is None:
        raise AssertionError

    if not (p := Path(s)).exists():
        return Error(f"{var} is set to {s!r}, but no such file or directory was found.")

    if p.is_file():
        return Error(f"{var} is set to {s!r}, which is a file (expected a directory).")

    assert p.is_dir()  # redundancy, fine if skipped in prod

    if not p.is_absolute():
        return Error(
            f"{var} is set to {s!r}, which is relative (expected an absolute path)."
        )

    return p


@final
class _DirectoryEnvvar(Protocol):
    dirtype: _DirType

    @property
    def name(self) -> str: ...
    def resolve(self) -> Path | Error:
        # if a Path is returned, the following conditions are guaranteed:
        # - it's absolute
        # - its parent exists
        # - no file with the same name exist
        # importantly, this method is not responsible for ensuring the return path
        # exists, just that it is available
        ...


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

    def resolve(self) -> Path | Error:
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

    def resolve(self) -> Path | Error:
        match _resolve(self.name):
            case str() as err:
                return Error(f"{err} This environment variable will be ignored.")
            case Path() as p:
                pass
            case _ as unreachable:
                assert_never(unreachable)

        if (subp := p / "astropy").is_file():
            return Error(
                f"{self.name} is set to {os.getenv(self.name)!r}, but this directory "
                "already contains a file under 'astropy', where a directory was expected. "
                "This environment variable will be ignored."
            )

        return subp


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
    home_override: Path | None = None

    def __post_init__(self) -> None:
        if self.home_override is None:
            return
        if not self.home_override.exists():
            raise FileNotFoundError(f"No such file or directory {self.home_override}")
        if not self.home_override.is_absolute():
            raise ValueError(
                f"override path must be absolute, but {self.home_override} is not."
            )
        if self.home_override.is_file():
            raise ValueError(
                f"override path must be a directory, but {self.home_override} is a file."
            )

        assert self.home_override.is_dir()  # redundancy, fine if skipped

    @property
    def envvars(self) -> _DirectoryEnvvarSet:
        return _DirectoryEnvvarSet(
            astropy=_AstropyEnvvar(self.dirtype),
            xdg=_XDGEnvvar(self.dirtype),
        )

    def default_base_node(self, base_name: str = "astropy", /) -> Path:
        return Path.home() / f".{base_name}"

    def find_base_node(self, base_name: str = "astropy", /) -> Path:
        if self.home_override is not None:
            return self.home_override

        envvars: list[_DirectoryEnvvar] = [self.envvars.astropy, self.envvars.xdg]
        for var in filter(lambda v: v.name in os.environ, envvars):
            match var.resolve():
                case str() as err:
                    warn(err, AstropyUserWarning, stacklevel=2)
                    continue
                case Path() as res:
                    return res
                case _ as unreachable:
                    assert_never(unreachable)

        return self.default_base_node(base_name)

    def find_specialized_node(self, base_name: str = "astropy", /) -> Path:
        base_node = self.find_base_node(base_name)
        return self.dirtype.specialize_node(base_node).absolute()


class _TempDirKwargs(TypedDict):
    dir: os.PathLike[str] | str | None
    suffix: str | None
    prefix: str | None
    delete: bool


@contextmanager
def _temporary_dir_ctx(
    dir_: os.PathLike[str] | str | None = None,
    /,
    *,
    suffix: str | None = None,
    prefix: str | None = None,
    delete: bool = True,
    base_name: str = "astropy",
    directory_finder_name: Literal["_CacheFinder", "_ConfigFinder"],
) -> Generator[Path, None, None]:
    # the common implementation for temporary_cache_path and temporary_config_path

    kwargs: _TempDirKwargs = {
        "dir": dir_,
        "suffix": suffix,
        "prefix": prefix,
        "delete": delete,
    }
    with _PATHS_MUTEX, TemporaryDirectory(**kwargs) as tmp_dir:
        tmp_path = Path(tmp_dir)
        initial_df: _DirectoryFinder = globals()[directory_finder_name]
        df: _DirectoryFinder = replace(initial_df, home_override=tmp_path)
        globals()[directory_finder_name] = df

        yield_val = df.find_specialized_node(base_name)
        if yield_val != tmp_path:
            assert yield_val.is_relative_to(tmp_path)
            yield_val.mkdir(parents=True)
        yield yield_val

        globals()[directory_finder_name] = initial_df


@contextmanager
def _clear_cfgobjs():
    from .configuration import _cfgobjs

    with _PATHS_MUTEX:
        initial_cfgobjs = deepcopy(_cfgobjs)
        _cfgobjs.clear()

        yield

        _cfgobjs.clear()
        _cfgobjs.update(initial_cfgobjs)


_CacheFinder = _DirectoryFinder(_DirType.CACHE)
temporary_cache_path = partial(
    _temporary_dir_ctx,
    directory_finder_name="_CacheFinder",
)

_ConfigFinder = _DirectoryFinder(_DirType.CONFIG)
temporary_config_path = partial(
    _clear_cfgobjs()(_temporary_dir_ctx),
    directory_finder_name="_ConfigFinder",
)


def get_config_dir_path(rootname: str = "astropy") -> Path:
    """
    Determines the package configuration directory name and creates the
    directory if it doesn't exist.

    This directory is typically ``$HOME/.astropy/config``, but if the
    XDG_CONFIG_HOME environment variable is set and the
    ``$XDG_CONFIG_HOME/astropy`` directory exists, it will be that directory.
    If neither exists, the former will be created and symlinked to the latter.

    Parameters
    ----------
    rootname : str
        Name of the root configuration directory. For example, if ``rootname =
        'pkgname'``, the configuration directory would be ``<home>/.pkgname/``
        rather than ``<home>/.astropy`` (depending on platform).

    Returns
    -------
    configdir : Path
        The absolute path to the configuration directory.

    """
    return set_temp_config._get_dir_path(rootname)


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
    Determines the Astropy cache directory name and creates the directory if it
    doesn't exist.

    This directory is typically ``$HOME/.astropy/cache``, but if the
    XDG_CACHE_HOME environment variable is set and the
    ``$XDG_CACHE_HOME/astropy`` directory exists, it will be that directory.
    If neither exists, the former will be created and symlinked to the latter.

    Parameters
    ----------
    rootname : str
        Name of the root cache directory. For example, if
        ``rootname = 'pkgname'``, the cache directory will be
        ``<cache>/.pkgname/``.

    Returns
    -------
    cachedir : Path
        The absolute path to the cache directory.

    """
    return set_temp_cache._get_dir_path(rootname)


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
        if path is not None:
            path = Path(path).resolve()

        self._path = path
        self._delete = delete
        self._prev_path = self._get_current_override()

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
        return cls._get_directory_finder().home_override

    @classmethod
    def _get_xdg_envvar_name(cls) -> str:
        xdg_envvar = cls._get_directory_finder().envvars.xdg
        return xdg_envvar.name

    def __enter__(self) -> str:
        _PATHS_MUTEX.acquire()
        try:
            finder_name = self._get_directory_finder_name()
            initial_df: _DirectoryFinder = globals()[finder_name]
            df: _DirectoryFinder = replace(initial_df, home_override=self._path)
            self._prev_df = initial_df
        except Exception:
            _PATHS_MUTEX.release()
            raise
        try:
            globals()[finder_name] = df
            return str(self.__class__._get_dir_path(rootname="astropy"))
        except Exception:
            globals()[finder_name] = initial_df
            _PATHS_MUTEX.release()
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
            _PATHS_MUTEX.release()

    def __call__(self, func: Callable[P, object]) -> Callable[P, None]:
        """Implements use as a decorator."""

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> None:
            with self:
                func(*args, **kwargs)

        return wrapper

    @staticmethod
    def _warn_env_var_is_ignored(
        reason: str, *, env_var_name: str, stacklevel: int = 4
    ) -> None:
        env_var_value = os.getenv(env_var_name)
        warning_msg_template = (
            f"{env_var_name} is set to '{env_var_value}', "
            "but {reason}. This environment variable will be ignored."
        )
        warn(
            warning_msg_template.format(reason=reason),
            category=AstropyUserWarning,
            stacklevel=stacklevel,
        )

    @classmethod
    def _get_dir_path(cls, rootname: str) -> Path:
        if (xch := cls._get_current_override()) is not None:
            path = xch / rootname
            if not path.is_file():
                path.mkdir(exist_ok=True)
            return path.resolve()

        envvar_name = cls._get_xdg_envvar_name()

        if (env_dir_str := os.getenv(envvar_name)) is None:
            return cls._find_or_create_root_dir(linkto=None, pkgname=rootname)

        warn = partial(
            cls._warn_env_var_is_ignored,
            env_var_name=envvar_name,
        )
        if not (env_dir_path := Path(env_dir_str)).is_dir():
            warn("no such directory was found")
            return cls._find_or_create_root_dir(linkto=None, pkgname=rootname)

        path = env_dir_path / rootname
        if path.is_symlink():
            return cls._find_or_create_root_dir(linkto=None, pkgname=rootname)
        elif path.is_dir():
            return path.resolve()

        linkto: Path | None
        if path.is_file():
            warn(
                f"this directory already contains a file under '{path.name}', "
                "where a directory was expected"
            )
            # renounce linking early to avoid redundant warnings
            linkto = None
        else:
            # symlink will be set to this if the directory is created
            linkto = path

        return cls._find_or_create_root_dir(linkto=linkto, pkgname=rootname)

    @classmethod
    def _find_or_create_root_dir(
        cls,
        linkto: Path | None,
        pkgname: str = "astropy",
    ) -> Path:
        innerdir = Path.home() / f".{pkgname}"
        maindir = innerdir / cls._dirtype.as_str()

        if maindir.is_file():
            raise OSError(
                f"Intended {pkgname} {cls._dirtype.as_str()} directory {maindir} is actually a file."
            )

        warn = partial(
            cls._warn_env_var_is_ignored,
            env_var_name=cls._get_xdg_envvar_name(),
        )

        # At this point, linkto is either None or a path that doesn't exist yet
        # Other linkto conditions already filtered in _get_dir_path() above.
        if linkto is not None:
            if sys.platform.startswith("win"):
                warn("is not supported on Windows")
                linkto = None
            elif maindir.is_dir():
                warn(
                    f"the default location, {maindir}, already exists, and takes precedence"
                )
                linkto = None  # redundant, but shouldn't make a difference

        if not maindir.is_dir():
            # first create .astropy dir if needed
            if innerdir.is_file():
                raise OSError(
                    f"Intended {pkgname} {cls._dirtype.as_str()} directory {maindir} is actually a file."
                )
            maindir.mkdir(parents=True, exist_ok=True)
            if linkto is not None:
                linkto.symlink_to(maindir)

        return maindir.resolve()


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
    """

    _dirtype = _DirType.CONFIG

    def __enter__(self) -> str:
        # Special case for the config case, where we need to reset all the
        # cached config objects.  We do keep the cache, since some of it
        # may have been set programmatically rather than be stored in the
        # config file (e.g., iers.conf.auto_download=False for our tests).
        _PATHS_MUTEX.acquire()

        try:
            from .configuration import _cfgobjs

            self._cfgobjs_copy = _cfgobjs.copy()
            _cfgobjs.clear()
            return super().__enter__()
        except Exception:
            _PATHS_MUTEX.release()
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
            _PATHS_MUTEX.release()


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
    """

    _dirtype = _DirType.CACHE
