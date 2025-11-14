# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import re
import subprocess
import sys
from contextlib import nullcontext
from inspect import cleandoc
from pathlib import Path
from threading import Lock

import pytest

from astropy.config import configuration, create_config_file, paths, set_temp_config
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

OLD_CONFIG = {}

_IGNORE_CONFIG_PATHS_GLOBAL_STATE_LOCK = Lock()


@pytest.fixture
def ignore_config_paths_global_state(monkeypatch, tmp_path_factory):
    # ignore global state of the test session
    # and preserve thread safety across all users of this fixture
    with _IGNORE_CONFIG_PATHS_GLOBAL_STATE_LOCK:
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)

        monkeypatch.setattr(paths.set_temp_cache, "_temp_path", None)
        monkeypatch.setattr(paths.set_temp_config, "_temp_path", None)

        # also mock $HOME as it's part of the global state taken into account
        # for path detection
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
    # Regression test for #17514 - XDG_CACHE_HOME had no effect
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
    # have the env var point to a writable location,
    # where an 'astropy' subdir is missing, but can be created silently
    default_parent_dir = Path.home() / ".astropy"
    assert not default_parent_dir.exists()

    target_dir = tmp_path
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir / get_directory_type(func)

    assert not expected_location.exists()

    if sys.platform.startswith("win"):
        # FIXME: legacy behavior is to never symlink on windows...
        # should we never symlink anywhere instead for consistency
        # (not to mention simplicity) ?
        expected_msg = (
            rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
            rf"but is not supported on Windows\. "
            r"This environment variable will be ignored\.$"
        )
        ctx = pytest.warns(AstropyUserWarning, match=expected_msg)
    else:
        ctx = nullcontext()

    with ctx:
        ret = func()

    assert ret.is_dir()

    # FIXME: should the environment variable *completely* take precedence here,
    # and should we *avoid* silently creating (and linking) the default location ?
    # see https://github.com/astropy/astropy/issues/18791
    assert ret == default_location

    if sys.platform.startswith("win"):
        assert not expected_location.exists()
        return

    assert expected_location.is_dir()
    assert default_location.is_dir()
    assert expected_location.exists()
    assert expected_location.is_dir()
    assert expected_location.is_symlink()
    assert expected_location.resolve() == default_location


@pytest.mark.parametrize("env_var, func", ENV_VAR_AND_FUNC)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_subdir_but_default_exists(
    monkeypatch, tmp_path, env_var, func
):
    # have the env var point to a writable location,
    # where an 'astropy' subdir is missing, but can be created silently
    default_parent_dir = Path.home() / ".astropy"

    target_dir = tmp_path
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir / get_directory_type(func)

    # this is the crucial variation in setup
    # as compared to test_env_variables_missing_subdir_and_default
    default_location.mkdir(parents=True)

    assert not expected_location.exists()

    expected_msg_template = (
        rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
        r"but {reason}. "
        r"This environment variable will be ignored\.$"
    )
    if sys.platform.startswith("win"):
        # FIXME: legacy behavior is to never symlink on windows...
        # should we never symlink anywhere instead for consistency
        # (not to mention simplicity) ?
        reason = "is not supported on Windows"
    else:
        reason = (
            rf"the default location, {re.escape(str(default_location))}, "
            "already exists, and takes precedence"
        )

    expected_msg = expected_msg_template.format(reason=reason)

    with pytest.warns(AstropyUserWarning, match=expected_msg):
        ret = func()

    assert ret.is_dir()

    # FIXME: should the environment variable *completely* take precedent here,
    # ignoring the existence of default_parent_dir ?
    # see https://github.com/astropy/astropy/issues/18791
    assert not expected_location.exists()


@pytest.mark.parametrize("env_var, func", ENV_VAR_AND_FUNC)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_setup_file_exists(monkeypatch, tmp_path, env_var, func):
    # check what happens if we request a location that's already
    # taken, but is a file
    default_path = func()
    target_dir = tmp_path / "subdir"
    target_dir.mkdir()
    expected_path = target_dir / "astropy"
    expected_path.touch()  # create a file

    monkeypatch.setenv(env_var, str(target_dir))

    expected_msg = (
        rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
        rf"but this directory already contains a file under '{expected_path.name}', "
        r"where a directory was expected\. "
        r"This environment variable will be ignored\.$"
    )
    with pytest.warns(AstropyUserWarning, match=expected_msg):
        new_path = func()
    assert new_path == default_path


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_set_temp_config(tmp_path):
    # Check that we start in an understood state.
    assert configuration._cfgobjs == OLD_CONFIG

    orig_config_dir = paths.get_config_dir(rootname="astropy")
    (temp_config_dir := tmp_path / "config").mkdir()
    temp_astropy_config = temp_config_dir / "astropy"

    # Test decorator mode
    @paths.set_temp_config(temp_config_dir)
    def test_func():
        assert paths.get_config_dir(rootname="astropy") == str(temp_astropy_config)
        assert paths.get_config_dir_path(rootname="astropy") == temp_astropy_config

        # Test temporary restoration of original default
        with paths.set_temp_config() as d:
            assert d == orig_config_dir == paths.get_config_dir(rootname="astropy")

    test_func()

    # Test context manager mode (with cleanup)
    with paths.set_temp_config(temp_config_dir, delete=True):
        assert paths.get_config_dir(rootname="astropy") == str(temp_astropy_config)
        assert paths.get_config_dir_path(rootname="astropy") == temp_astropy_config

    assert not temp_config_dir.exists()
    # Check that we have returned to our old configuration.
    assert configuration._cfgobjs == OLD_CONFIG


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_set_temp_cache(tmp_path):
    orig_cache_dir = paths.get_cache_dir(rootname="astropy")
    (temp_cache_dir := tmp_path / "cache").mkdir()
    temp_astropy_cache = temp_cache_dir / "astropy"

    # Test decorator mode
    @paths.set_temp_cache(temp_cache_dir)
    def test_func():
        assert paths.get_cache_dir(rootname="astropy") == str(temp_astropy_cache)

        # Test temporary restoration of original default
        with paths.set_temp_cache() as d:
            assert d == orig_cache_dir == paths.get_cache_dir(rootname="astropy")

    test_func()

    # Test context manager mode (with cleanup)
    with paths.set_temp_cache(temp_cache_dir, delete=True):
        assert paths.get_cache_dir(rootname="astropy") == str(temp_astropy_cache)

    assert not temp_cache_dir.exists()


def test_set_temp_cache_resets_on_exception(tmp_path):
    """Test for regression of  bug #9704"""
    t = paths.get_cache_dir()
    (a := tmp_path / "a").write_text("not a good cache\n")
    with pytest.raises(OSError), paths.set_temp_cache(a):
        pass
    assert t == paths.get_cache_dir()


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_config_file():
    from astropy.config.configuration import get_config, reload_config

    apycfg = get_config("astropy")
    assert apycfg.filename.endswith("astropy.cfg")

    cfgsec = get_config("astropy.config")
    assert cfgsec.depth == 1
    assert cfgsec.name == "config"
    assert cfgsec.parent.filename.endswith("astropy.cfg")

    # try with a different package name, still inside astropy config dir:
    testcfg = get_config("testpkg", rootname="astropy")
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert ".astropy" in parts or "astropy" in parts
    assert parts[-1] == "testpkg.cfg"
    configuration._cfgobjs["testpkg"] = None  # HACK

    # try with a different package name, no specified root name (should
    #   default to astropy):
    testcfg = get_config("testpkg")
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert ".astropy" in parts or "astropy" in parts
    assert parts[-1] == "testpkg.cfg"
    configuration._cfgobjs["testpkg"] = None  # HACK

    # try with a different package name, specified root name:
    testcfg = get_config("testpkg", rootname="testpkg")
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert ".testpkg" in parts or "testpkg" in parts
    assert parts[-1] == "testpkg.cfg"
    configuration._cfgobjs["testpkg"] = None  # HACK

    # try with a subpackage with specified root name:
    testcfg_sec = get_config("testpkg.somemodule", rootname="testpkg")
    parts = os.path.normpath(testcfg_sec.parent.filename).split(os.sep)
    assert ".testpkg" in parts or "testpkg" in parts
    assert parts[-1] == "testpkg.cfg"
    configuration._cfgobjs["testpkg"] = None  # HACK

    reload_config("astropy")


def check_config(conf):
    # test that the output contains some lines that we expect
    assert "# unicode_output = False" in conf
    assert "[io.fits]" in conf
    assert "[table]" in conf
    assert "# replace_warnings = ," in conf
    assert "[table.jsviewer]" in conf
    assert "# css_urls = https://cdn.datatables.net/2.1.8/css/dataTables.dataTables.min.css" in conf  # fmt: skip
    assert "[visualization.wcsaxes]" in conf
    assert "## Whether to log exceptions before raising them." in conf
    assert "# log_exceptions = False" in conf


def test_generate_config(tmp_path):
    from astropy.config.configuration import generate_config

    out = io.StringIO()
    generate_config("astropy", out)
    conf = out.getvalue()

    outfile = tmp_path / "astropy.cfg"
    generate_config("astropy", outfile)
    with open(outfile) as fp:
        conf2 = fp.read()

    for c in (conf, conf2):
        check_config(c)


def test_generate_config2(tmp_path):
    """Test that generate_config works with the default filename."""

    with set_temp_config(tmp_path):
        from astropy.config.configuration import generate_config

        generate_config("astropy")

    assert os.path.exists(tmp_path / "astropy" / "astropy.cfg")

    with open(tmp_path / "astropy" / "astropy.cfg") as fp:
        conf = fp.read()

    check_config(conf)


def test_generate_config_subclasses(tmp_path):
    """Test that generate_config works with subclasses of ConfigNamespace."""
    from astropy.config.configuration import ConfigItem, ConfigNamespace

    class MyPackageNamespace(ConfigNamespace):
        pass

    class RecursiveTestConf(MyPackageNamespace):
        ti = ConfigItem(5, "this is a Description")

    with set_temp_config(tmp_path):
        from astropy.config.configuration import generate_config

        generate_config("astropy")

    assert os.path.exists(tmp_path / "astropy" / "astropy.cfg")

    with open(tmp_path / "astropy" / "astropy.cfg") as fp:
        conf = fp.read()

    assert "# ti = 5" in conf


def test_create_config_file(tmp_path, caplog):
    with set_temp_config(tmp_path):
        create_config_file("astropy")

    # check that the config file has been created
    assert (
        "The configuration file has been successfully written"
        in caplog.records[0].message
    )
    assert os.path.exists(tmp_path / "astropy" / "astropy.cfg")

    with open(tmp_path / "astropy" / "astropy.cfg") as fp:
        conf = fp.read()
    check_config(conf)

    caplog.clear()

    # now modify the config file
    conf = conf.replace("# unicode_output = False", "unicode_output = True")
    with open(tmp_path / "astropy" / "astropy.cfg", mode="w") as fp:
        fp.write(conf)

    with set_temp_config(tmp_path):
        create_config_file("astropy")

    # check that the config file has not been overwritten since it was modified
    assert (
        "The configuration file already exists and seems to have been customized"
        in caplog.records[0].message
    )

    caplog.clear()

    with set_temp_config(tmp_path):
        create_config_file("astropy", overwrite=True)

    # check that the config file has been overwritten
    assert (
        "The configuration file has been successfully written"
        in caplog.records[0].message
    )


def test_configitem():
    from astropy.config.configuration import ConfigItem, ConfigNamespace, get_config

    ci = ConfigItem(34, "this is a Description")

    class Conf(ConfigNamespace):
        tstnm = ci

    conf = Conf()

    assert ci.module == "astropy.config.tests.test_configs"
    assert ci() == 34
    assert ci.description == "this is a Description"

    assert conf.tstnm == 34

    assert str(conf) == cleandoc(
        """
        Configuration parameters for `astropy.config.tests.test_configs`

        ConfigItem: tstnm
          cfgtype='integer'
          defaultvalue=34
          description='this is a Description'
          module=astropy.config.tests.test_configs
          value=34
        """
    )

    sec = get_config(ci.module)
    assert sec["tstnm"] == 34

    ci.description = "updated Descr"
    ci.set(32)
    assert ci() == 32

    # It's useful to go back to the default to allow other test functions to
    # call this one and still be in the default configuration.
    ci.description = "this is a Description"
    ci.set(34)
    assert ci() == 34

    # Test iterator for one-item namespace
    result = list(conf)
    assert result == ["tstnm"]
    result = list(conf.keys())
    assert result == ["tstnm"]
    result = list(conf.values())
    assert result == [ci]
    result = list(conf.items())
    assert result == [("tstnm", ci)]


def test_configitem_types():
    from astropy.config.configuration import ConfigItem, ConfigNamespace

    ci1 = ConfigItem(34)
    ci2 = ConfigItem(34.3)
    ci3 = ConfigItem(True)
    ci4 = ConfigItem("astring")

    class Conf(ConfigNamespace):
        tstnm1 = ci1
        tstnm2 = ci2
        tstnm3 = ci3
        tstnm4 = ci4

    conf = Conf()

    assert isinstance(conf.tstnm1, int)
    assert isinstance(conf.tstnm2, float)
    assert isinstance(conf.tstnm3, bool)
    assert isinstance(conf.tstnm4, str)

    with pytest.raises(TypeError):
        conf.tstnm1 = 34.3
    conf.tstnm2 = 12  # this would should succeed as up-casting
    with pytest.raises(TypeError):
        conf.tstnm3 = "fasd"
    with pytest.raises(TypeError):
        conf.tstnm4 = 546.245

    # Test iterator for multi-item namespace. Assume ordered by insertion order.
    item_names = list(conf)
    assert item_names == ["tstnm1", "tstnm2", "tstnm3", "tstnm4"]
    result = list(conf.keys())
    assert result == item_names
    result = list(conf.values())
    assert result == [ci1, ci2, ci3, ci4]
    result = list(conf.items())
    assert result == [
        ("tstnm1", ci1),
        ("tstnm2", ci2),
        ("tstnm3", ci3),
        ("tstnm4", ci4),
    ]


def test_configitem_options(tmp_path):
    from astropy.config.configuration import ConfigItem, ConfigNamespace, get_config

    cio = ConfigItem(["op1", "op2", "op3"])

    class Conf(ConfigNamespace):
        tstnmo = cio

    sec = get_config(cio.module)

    assert isinstance(cio(), str)
    assert cio() == "op1"
    assert sec["tstnmo"] == "op1"

    cio.set("op2")
    with pytest.raises(TypeError):
        cio.set("op5")
    assert sec["tstnmo"] == "op2"

    # now try saving
    apycfg = sec
    while apycfg.parent is not apycfg:
        apycfg = apycfg.parent
    f = tmp_path / "astropy.cfg"
    with open(f, "wb") as fd:
        apycfg.write(fd)

    assert "tstnmo = op2" in f.read_text().splitlines()


def test_help(capsys):
    from astropy import conf

    use_color_msg = cleandoc(
        """
        ConfigItem: use_color
          cfgtype='boolean'
          defaultvalue={is_not_windows}
          description='When True, use ANSI color escape sequences when writing to the console.'
          module=astropy
          value={is_not_windows}
        """
    ).format(is_not_windows=sys.platform != "win32")
    conf.help("use_color")
    assert capsys.readouterr().out == use_color_msg + "\n"

    conf.help()
    help_text = capsys.readouterr().out
    assert help_text.startswith(
        "Configuration parameters for `astropy`.\n\nConfigItem: unicode_output"
    )
    assert use_color_msg in help_text


def test_help_invalid_config_item():
    from astropy import conf

    with pytest.raises(
        KeyError,
        match=(
            "'bad_name' is not among configuration items "
            r"\('unicode_output', 'use_color', 'max_lines', 'max_width'\)"
        ),
    ):
        conf.help("bad_name")


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_config_noastropy_fallback(monkeypatch):
    """
    Tests to make sure configuration items fall back to their defaults when
    there's a problem accessing the astropy directory
    """

    # make sure the _find_or_create_root_dir function fails as though the
    # astropy dir could not be accessed
    @classmethod
    def osraiser(cls, linkto, pkgname=None):
        raise OSError

    monkeypatch.setattr(paths._SetTempPath, "_find_or_create_root_dir", osraiser)

    # also have to make sure the stored configuration objects are cleared
    monkeypatch.setattr(configuration, "_cfgobjs", {})

    # make sure the config dir search fails
    with pytest.raises(OSError):
        paths.get_config_dir(rootname="astropy")

    with pytest.raises(OSError):
        paths.get_config_dir_path(rootname="astropy")

    # now run the basic tests, and make sure the warning about no astropy
    # is present
    test_configitem()


def test_configitem_setters():
    from astropy.config.configuration import ConfigItem, ConfigNamespace

    class Conf(ConfigNamespace):
        tstnm12 = ConfigItem(42, "this is another Description")

    conf = Conf()

    assert conf.tstnm12 == 42
    with conf.set_temp("tstnm12", 45):
        assert conf.tstnm12 == 45
    assert conf.tstnm12 == 42

    conf.tstnm12 = 43
    assert conf.tstnm12 == 43

    with conf.set_temp("tstnm12", 46):
        assert conf.tstnm12 == 46

    # Make sure it is reset even with Exception
    try:
        with conf.set_temp("tstnm12", 47):
            raise Exception
    except Exception:
        pass

    assert conf.tstnm12 == 43


def test_empty_config_file():
    from astropy.config.configuration import is_unedited_config_file

    def get_content(fn):
        with open(get_pkg_data_filename(fn), encoding="latin-1") as fd:
            return fd.read()

    content = get_content("data/empty.cfg")
    assert is_unedited_config_file(content)

    content = get_content("data/not_empty.cfg")
    assert not is_unedited_config_file(content)


class TestAliasRead:
    def setup_class(self):
        configuration._override_config_file = get_pkg_data_filename("data/alias.cfg")

    def test_alias_read(self):
        from astropy.utils.data import conf

        with pytest.warns(
            AstropyDeprecationWarning,
            match=r"Config parameter 'name_resolve_timeout' in section "
            r"\[coordinates.name_resolve\].*",
        ) as w:
            conf.reload()
            assert conf.remote_timeout == 42

        assert len(w) == 1

    def teardown_class(self):
        from astropy.utils.data import conf

        configuration._override_config_file = None
        conf.reload()


def test_configitem_unicode():
    from astropy.config.configuration import ConfigItem, ConfigNamespace, get_config

    cio = ConfigItem("ასტრონომიის")

    class Conf(ConfigNamespace):
        tstunicode = cio

    sec = get_config(cio.module)

    assert isinstance(cio(), str)
    assert cio() == "ასტრონომიის"
    assert sec["tstunicode"] == "ასტრონომიის"


def test_warning_move_to_top_level():
    # Check that the warning about deprecation config items in the
    # file works.  See #2514
    from astropy import conf

    configuration._override_config_file = get_pkg_data_filename("data/deprecated.cfg")

    try:
        with pytest.warns(AstropyDeprecationWarning) as w:
            conf.reload()
            conf.max_lines
        assert len(w) == 1
    finally:
        configuration._override_config_file = None
        conf.reload()


def test_no_home():
    # "import astropy" fails when neither $HOME or $XDG_CONFIG_HOME
    # are set.  To test, we unset those environment variables for a
    # subprocess and try to import astropy.

    test_path = os.path.dirname(__file__)
    astropy_path = os.path.abspath(os.path.join(test_path, "..", "..", ".."))

    env = os.environ.copy()
    paths = [astropy_path]
    if env.get("PYTHONPATH"):
        paths.append(env.get("PYTHONPATH"))
    env["PYTHONPATH"] = os.pathsep.join(paths)

    for val in ["HOME", "XDG_CONFIG_HOME"]:
        if val in env:
            del env[val]

    retcode = subprocess.check_call([sys.executable, "-c", "import astropy"], env=env)

    assert retcode == 0
