# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import sys
from contextlib import nullcontext
from pathlib import Path
from threading import Lock

import pytest

from astropy.config import configuration, paths
from astropy.utils.exceptions import AstropyUserWarning

OLD_CONFIG = {}

_IGNORE_CONFIG_PATHS_GLOBAL_STATE_LOCK = Lock()


@pytest.fixture
def ignore_config_paths_global_state(monkeypatch, tmp_path_factory):
    with _IGNORE_CONFIG_PATHS_GLOBAL_STATE_LOCK:
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)

        monkeypatch.setattr(paths.set_temp_cache, "_temp_path", None)
        monkeypatch.setattr(paths.set_temp_config, "_temp_path", None)

        mock_home_dir = tmp_path_factory.mktemp("MOCK_HOME")

        def mock_home():
            return mock_home_dir

        monkeypatch.setattr(Path, "home", mock_home)

        yield


def setup_module():
    OLD_CONFIG.clear()
    OLD_CONFIG.update(configuration._cfgobjs)


def teardown_module():
    configuration._cfgobjs.clear()
    configuration._cfgobjs.update(OLD_CONFIG)


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_paths():
    assert "astropy" in paths.get_config_dir()
    assert "astropy" in paths.get_cache_dir()
    assert str(paths.get_config_dir_path()) == paths.get_config_dir()
    assert str(paths.get_cache_dir_path()) == paths.get_cache_dir()
    assert paths.get_config_dir_path().is_absolute()
    assert paths.get_cache_dir_path().is_absolute()

    assert "testpkg" in paths.get_config_dir(rootname="testpkg")
    assert "testpkg" in paths.get_cache_dir(rootname="testpkg")


ENV_VAR_AND_FUNC = [
    pytest.param(
        "XDG_CACHE_HOME",
        paths.get_cache_dir_path,
        id="xdg-cache",
    ),
    pytest.param(
        "XDG_CONFIG_HOME",
        paths.get_config_dir_path,
        id="xdg-config",
    ),
]


@pytest.mark.parametrize("env_var, func", ENV_VAR_AND_FUNC)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_xdg_variables(monkeypatch, tmp_path, env_var, func):
    target_dir = tmp_path / "astropy"
    target_dir.mkdir()
    monkeypatch.setenv(env_var, str(tmp_path))
    assert func() == target_dir


@pytest.mark.parametrize("env_var, func", ENV_VAR_AND_FUNC)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_dir(monkeypatch, tmp_path, env_var, func):
    default_path = func()
    target_dir = tmp_path / "nonexistent"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_msg = (
        rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
        rf"but no such directory was found\. "
        r"This environment variable will be ignored\.$"
    )
    with pytest.warns(AstropyUserWarning, match=expected_msg):
        new_path = func()

    assert new_path == default_path


def get_directory_type(func) -> str:
    for t in ["cache", "config"]:
        if t in func.__name__:
            return t
    raise RuntimeError


@pytest.mark.parametrize("env_var, func", ENV_VAR_AND_FUNC)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_subdir_and_default(monkeypatch, tmp_path, env_var, func):
    default_parent_dir = Path.home() / ".astropy"
    assert not default_parent_dir.exists()

    target_dir = tmp_path
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir / get_directory_type(func)

    assert not expected_location.exists()

    if sys.platform.startswith("win"):
        expected_msg = (
            rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
            r"but is not supported on Windows\. "
            r"This environment variable will be ignored\.$"
        )
        ctx = pytest.warns(AstropyUserWarning, match=expected_msg)
    else:
        ctx = nullcontext()

    with ctx:
        ret = func()

    assert ret.is_dir() or sys.platform.startswith("win")
    if not sys.platform.startswith("win"):
        assert expected_location.is_dir()
        assert default_location.is_dir()
        assert expected_location.is_symlink()
        assert expected_location.resolve() == default_location


# --- remaining tests stay unchanged ---
# You can keep the rest of the tests as you had them before,
# they do not require modification for the Windows warning fix.
