# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains functions to determine where configuration and
data/cache files used by Astropy should be placed.
"""

import os
import shutil
import sys
from functools import wraps
from pathlib import Path

__all__ = ["get_config_dir", "get_cache_dir", "set_temp_config", "set_temp_cache"]


def get_config_dir(rootname: str = "astropy") -> str:
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
    configdir : str
        The absolute path to the configuration directory.

    """
    # symlink will be set to this if the directory is created
    linkto = None

    # If using set_temp_config, that overrides all
    if set_temp_config._temp_path is not None:
        xch = set_temp_config._temp_path
        config_path = xch / rootname
        config_path.mkdir(exist_ok=True)
        return str(config_path.resolve())

    # first look for XDG_CONFIG_HOME
    if (
        (xdg_config_home := os.getenv("XDG_CONFIG_HOME")) is not None
        and (xch := Path(xdg_config_home)).exists()
        and not (xchpth := xch / rootname).is_symlink()
    ):
        if xchpth.exists():
            return str(xchpth.resolve())
        else:
            linkto = xchpth

    return str(_find_or_create_root_dir("config", linkto, rootname))


def get_cache_dir(rootname: str = "astropy") -> str:
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
    cachedir : str
        The absolute path to the cache directory.

    """
    # symlink will be set to this if the directory is created
    linkto = None

    # If using set_temp_cache, that overrides all
    if set_temp_cache._temp_path is not None:
        xch = set_temp_cache._temp_path
        cache_path = xch / rootname
        if not cache_path.is_file():
            cache_path.mkdir(exist_ok=True)
        return str(cache_path.resolve())

    # first look for XDG_CACHE_HOME
    if (
        (xdg_config_home := os.getenv("XDG_CONFIG_HOME")) is not None
        and (xch := Path(xdg_config_home)).exists()
        and not (xchpth := xch / rootname).is_symlink()
    ):
        if xchpth.exists():
            return str(xchpth.resolve())
        else:
            linkto = xchpth

    return str(_find_or_create_root_dir("cache", linkto, rootname))


class _SetTempPath:
    _temp_path = None
    _default_path_getter = None

    def __init__(self, path=None, delete=False):
        if path is not None:
            path = Path(path).resolve()

        self._path = path
        self._delete = delete
        self._prev_path = self.__class__._temp_path

    def __enter__(self):
        self.__class__._temp_path = self._path
        try:
            return self._default_path_getter("astropy")
        except Exception:
            self.__class__._temp_path = self._prev_path
            raise

    def __exit__(self, *args):
        self.__class__._temp_path = self._prev_path

        if self._delete and self._path is not None:
            shutil.rmtree(self._path)

    def __call__(self, func):
        """Implements use as a decorator."""

        @wraps(func)
        def wrapper(*args, **kwargs):
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

    _default_path_getter = staticmethod(get_config_dir)

    def __enter__(self):
        # Special case for the config case, where we need to reset all the
        # cached config objects.  We do keep the cache, since some of it
        # may have been set programmatically rather than be stored in the
        # config file (e.g., iers.conf.auto_download=False for our tests).
        from .configuration import _cfgobjs

        self._cfgobjs_copy = _cfgobjs.copy()
        _cfgobjs.clear()
        return super().__enter__()

    def __exit__(self, *args):
        from .configuration import _cfgobjs

        _cfgobjs.clear()
        _cfgobjs.update(self._cfgobjs_copy)
        del self._cfgobjs_copy
        super().__exit__(*args)


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

    _default_path_getter = staticmethod(get_cache_dir)


def _find_or_create_root_dir(
    dirnm: str,
    linkto: Path | None,
    pkgname: str = "astropy",
) -> Path:
    innerdir = Path.home() / f".{pkgname}"
    maindir = innerdir / dirnm

    if maindir.is_file():
        raise OSError(
            f"Intended {pkgname} {dirnm} directory {maindir} is actually a file."
        )
    if not maindir.is_dir():
        # first create .astropy dir if needed
        if innerdir.is_file():
            raise OSError(
                f"Intended {pkgname} {dirnm} directory {maindir} is actually a file."
            )
        maindir.mkdir(parents=True, exist_ok=True)
        if (
            not sys.platform.startswith("win")
            and linkto is not None
            and not linkto.exists()
        ):
            os.symlink(maindir, linkto)

    return maindir.resolve()
