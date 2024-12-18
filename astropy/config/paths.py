# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""

from __future__ import annotations

import os
import shutil
import sys
from functools import wraps
from inspect import cleandoc
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable
    from types import TracebackType
    from typing import Literal, ParamSpec

    P = ParamSpec("P")

__all__ = [
    "get_cache_dir",
    "get_cache_dir_path",
    "get_config_dir",
    "get_config_dir_path",
    "set_temp_cache",
    "set_temp_config",
]


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
    _temp_path: Path | None = None

    # This base class serves as a deduplication layer for its only two intended
    # children (set_temp_cache and set_temp_config)
    _directory_type: Literal["cache", "config"]
    _directory_env_var: Literal["XDG_CACHE_HOME", "XDG_CONFIG_HOME"]

    def __init__(
        self, path: os.PathLike[str] | str | None = None, delete: bool = False
    ) -> None:
        if path is not None:
            path = Path(path).resolve()

        self._path = path
        self._delete = delete
        self._prev_path = self.__class__._temp_path

    def __enter__(self) -> str:
        self.__class__._temp_path = self._path
        try:
            return str(self.__class__._get_dir_path(rootname="astropy"))
        except Exception:
            self.__class__._temp_path = self._prev_path
            raise

    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        self.__class__._temp_path = self._prev_path

        if self._delete and self._path is not None:
            shutil.rmtree(self._path)

    def __call__(self, func: Callable[P, object]) -> Callable[P, None]:
        """Implements use as a decorator."""

        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> None:
            with self:
                func(*args, **kwargs)

        return wrapper

    @classmethod
    def _get_dir_path(cls, rootname: str) -> Path:
        if (xch := cls._temp_path) is not None:
            path = xch / rootname
            if not path.is_file():
                path.mkdir(exist_ok=True)
            return path.resolve()

        if (
            (dir_ := os.getenv(cls._directory_env_var)) is not None
            and (xch := Path(dir_)).exists()
            and not (xchpth := xch / rootname).is_symlink()
        ):
            if xchpth.exists():
                return xchpth.resolve()

            # symlink will be set to this if the directory is created
            linkto = xchpth
        else:
            linkto = None

        return cls._find_or_create_root_dir(linkto, rootname)

    @classmethod
    def _find_or_create_root_dir(
        cls,
        linkto: Path | None,
        pkgname: str = "astropy",
    ) -> Path:
        innerdir = Path.home() / f".{pkgname}"
        maindir = innerdir / cls._directory_type

        if maindir.is_file():
            raise OSError(
                f"Intended {pkgname} {cls._directory_type} directory {maindir} is actually a file."
            )
        if not maindir.is_dir():
            # first create .astropy dir if needed
            if innerdir.is_file():
                raise OSError(
                    f"Intended {pkgname} {cls._directory_type} directory {maindir} is actually a file."
                )
            maindir.mkdir(parents=True, exist_ok=True)
            if (
                not sys.platform.startswith("win")
                and linkto is not None
                and not linkto.exists()
            ):
                os.symlink(maindir, linkto)

        return maindir.resolve()


class set_temp_config(_SetTempPath):
    """
    Context manager to set a temporary path for the Astropy config, primarily
    for use with testing.

    If the path set by this context manager does not already exist it will be
    created, if possible.

    This may also be used as a decorator on a function to set the config path
    just within that function.

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

    _directory_type = "config"
    _directory_env_var = "XDG_CONFIG_HOME"

    def __enter__(self) -> str:
        # Special case for the config case, where we need to reset all the
        # cached config objects.  We do keep the cache, since some of it
        # may have been set programmatically rather than be stored in the
        # config file (e.g., iers.conf.auto_download=False for our tests).
        from .configuration import _cfgobjs

        self._cfgobjs_copy = _cfgobjs.copy()
        _cfgobjs.clear()
        return super().__enter__()

    def __exit__(
        self,
        type: type[BaseException] | None,
        value: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        from .configuration import _cfgobjs

        _cfgobjs.clear()
        _cfgobjs.update(self._cfgobjs_copy)
        del self._cfgobjs_copy
        super().__exit__(type, value, tb)


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

    _directory_type = "cache"
    _directory_env_var = "XDG_CACHE_HOME"
