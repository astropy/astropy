# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import re
import subprocess
import sys
from contextlib import chdir, nullcontext
from inspect import cleandoc
from pathlib import Path
from typing import Literal
from uuid import uuid4

import pytest

from astropy.config import (
    configuration,
    create_config_file,
    paths,
    set_temp_config,
    temporary_cache_path,
    temporary_config_path,
)
from astropy.config.paths import (
    _AstropyEnvvar,
    _DirectoryFinder,
    _DirType,
    _resolve,
    _XDGEnvvar,
)
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning

OLD_CONFIG = {}


def setup_module():
    OLD_CONFIG.clear()
    OLD_CONFIG.update(configuration._cfgobjs)


def teardown_module():
    configuration._cfgobjs.clear()
    configuration._cfgobjs.update(OLD_CONFIG)


@pytest.mark.parametrize(
    "cls, regexp",
    [
        pytest.param(
            _AstropyEnvvar, r"^ASTROPY_(?P<dirtype>[A-Z]+)_DIR$", id="ASTROPY_"
        ),
        pytest.param(_XDGEnvvar, r"^XDG_(?P<dirtype>[A-Z]+)_HOME$", id="XDG_"),
    ],
)
@pytest.mark.parametrize("dirtype", _DirType)
def test_direnvvar_impl(cls, regexp, dirtype):
    ev = cls(dirtype)
    assert (m := re.match(regexp, ev.name)) is not None
    assert m.group("dirtype") == dirtype.name


# TODO: split into smaller tests
def test_resolve_envvar(monkeypatch, tmp_path: Path):
    uid = str(uuid4())
    d1 = tmp_path / uid
    value = str(d1)
    monkeypatch.setenv("TEST_ENV_VAR", value)
    err_template = "TEST_ENV_VAR is set to {value!r}, {reason}."

    err1 = _resolve("TEST_ENV_VAR")
    assert isinstance(err1, str)
    assert err1 == err_template.format(
        value=value,
        reason="but no such file or directory was found",
    )

    d1.touch()  # create a file
    err2 = _resolve("TEST_ENV_VAR")
    assert isinstance(err2, str)
    assert err2 == err_template.format(
        value=value,
        reason="which is a file (expected a directory)",
    )
    d1.unlink()

    d1.mkdir()
    p = _resolve("TEST_ENV_VAR")
    assert isinstance(p, Path) and p.is_absolute()
    assert p == d1

    with chdir(tmp_path):
        d2 = d1.relative_to(Path.cwd())
        monkeypatch.setenv("TEST_ENV_VAR", str(d2))
        err3 = _resolve("TEST_ENV_VAR")
    assert isinstance(err3, str)
    assert err3 == err_template.format(
        value=str(d2),
        reason="which is relative (expected an absolute path)",
    )


@pytest.mark.parametrize("dirtype", _DirType)
def test_xdgenvvar_resolve_not_found(dirtype, monkeypatch, tmp_path):
    envvar = _XDGEnvvar(dirtype)
    uid = str(uuid4())
    missing_path = tmp_path / uid
    monkeypatch.setenv(envvar.name, str(missing_path))

    err = envvar.resolve()
    assert isinstance(err, str)
    assert err == (
        f"{envvar.name} is set to {str(missing_path)!r}, "
        "but no such file or directory was found. "
        "This environment variable will be ignored."
    )


@pytest.mark.parametrize("dirtype", _DirType)
def test_xdgenvvar_resolve_invalid_content(dirtype, monkeypatch, tmp_path):
    envvar = _XDGEnvvar(dirtype)
    monkeypatch.setenv(envvar.name, str(tmp_path))
    tmp_path.joinpath("astropy").touch()

    err = envvar.resolve()
    assert isinstance(err, str)
    assert err == (
        f"{envvar.name} is set to {str(tmp_path)!r}, "
        "but this directory already contains a file under 'astropy', "
        "where a directory was expected. "
        "This environment variable will be ignored."
    )


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.parametrize(
    "setup",
    [
        pytest.param(lambda d: None, id="missing-dir"),
        pytest.param(lambda d: d.mkdir(), id="existing-dir"),
    ],
)
def test_xdgenvvar_resolve_success(dirtype, setup, monkeypatch, tmp_path):
    envvar = _XDGEnvvar(dirtype)
    monkeypatch.setenv(envvar.name, str(tmp_path))

    expected = tmp_path / "astropy"
    setup(expected)
    p = envvar.resolve()
    assert isinstance(p, Path)
    assert p.is_absolute()
    assert p == expected


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_find_default(dirtype, tmp_path: Path):
    df = _DirectoryFinder(dirtype)
    assert df.find_base_node("mybase") == df.default_base_node("mybase")


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_override_not_found(dirtype, tmp_path: Path):
    with pytest.raises(
        FileNotFoundError,
        match=rf"^No such file or directory {re.escape(str(tmp_path / 'spam'))}",
    ):
        _DirectoryFinder(dirtype, home_override=tmp_path / "spam")


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_override_file(dirtype, tmp_path: Path):
    egg = tmp_path / "egg"
    egg.touch()
    with pytest.raises(
        ValueError,
        match=rf"^override path must be a directory, but {re.escape(str(egg))} is a file\.$",
    ):
        _DirectoryFinder(dirtype, home_override=egg)


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_override_relpath(dirtype, tmp_path: Path):
    egg = tmp_path.joinpath("egg")
    egg.touch()
    with chdir(tmp_path):
        eggs = egg.relative_to(Path.cwd())
        with pytest.raises(
            ValueError,
            match=rf"^override path must be absolute, but {re.escape(str(eggs))} is not\.$",
        ):
            _DirectoryFinder(dirtype, home_override=eggs)


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_override_success(dirtype, tmp_path: Path):
    bacon = tmp_path / "bacon"
    bacon.mkdir()
    df = _DirectoryFinder(dirtype, home_override=bacon)
    assert df.find_base_node() == bacon


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.parametrize("var_kind", ["astropy", "xdg"])
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_find_from_envvars_single_set(
    var_kind: Literal["astropy", "xdg"], dirtype: _DirType, tmp_path: Path, monkeypatch
):
    df = _DirectoryFinder(dirtype)
    p = tmp_path / "spam"

    envvar_name = getattr(df.envvars, var_kind).name
    monkeypatch.setenv(envvar_name, str(p))

    err_template = "{var} is set to {val!r}, which is a file (expected a directory)."

    # file
    p.touch()
    with pytest.warns(
        AstropyUserWarning,
        match=re.escape(
            err_template.format(
                var=envvar_name,
                val=str(p),
            )
        ),
    ):
        res = df.find_base_node()
    assert res == df.default_base_node()
    p.unlink()

    # dir
    p.mkdir()
    res = df.find_base_node()
    match var_kind:
        case "astropy":
            expected = p
        case "xdg":
            expected = p / "astropy"
        case _:
            raise AssertionError
    assert res == expected


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_find_from_envvars_double_set(
    dirtype: _DirType, tmp_path: Path, monkeypatch
):
    df = _DirectoryFinder(dirtype)
    p0 = tmp_path / "spam"
    p1 = tmp_path / "bacon"

    monkeypatch.setenv(df.envvars.astropy.name, str(p0))
    monkeypatch.setenv(df.envvars.xdg.name, str(p1))

    err_template = "{var} is set to {val!r}, which is a file (expected a directory)."
    # file + file
    p0.touch()
    p1.touch()
    with (
        pytest.warns(
            AstropyUserWarning,
            match=re.escape(
                err_template.format(
                    var=df.envvars.astropy.name,
                    val=str(p0),
                )
            ),
        ),
        pytest.warns(
            AstropyUserWarning,
            match=re.escape(
                err_template.format(
                    var=df.envvars.xdg.name,
                    val=str(p1),
                )
            ),
        ),
    ):
        res = df.find_base_node()
    assert res == df.default_base_node()
    p0.unlink()
    p1.unlink()

    # file + dir
    p0.touch()
    p1.mkdir()
    with pytest.warns(
        AstropyUserWarning,
        match=re.escape(
            err_template.format(
                var=df.envvars.astropy.name,
                val=str(p0),
            )
        ),
    ):
        res = df.find_base_node()
    assert res == p1 / "astropy"
    p0.unlink()
    p1.rmdir()

    # dir + file
    p0.mkdir()
    p1.touch()
    res = df.find_base_node()
    assert res == p0
    p0.rmdir()
    p1.unlink()

    # dir + dir
    p0.mkdir()
    p1.mkdir()
    res = df.find_base_node()
    assert res == p0


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


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_xdg_variables(monkeypatch, tmp_path, dirtype):
    # Regression test for #17514 - XDG_CACHE_HOME had no effect
    target_dir = tmp_path / "astropy"
    target_dir.mkdir()
    monkeypatch.setenv(f"XDG_{dirtype.name}_HOME", str(tmp_path))

    assert dirtype.path_getter() == target_dir


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_dir(monkeypatch, tmp_path, dirtype):
    default_path = dirtype.path_getter()
    target_dir = tmp_path / "nonexistent"
    env_var = f"XDG_{dirtype.name}_HOME"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_msg = (
        rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
        rf"but no such directory was found\. "
        r"This environment variable will be ignored\.$"
    )
    with pytest.warns(AstropyUserWarning, match=expected_msg):
        new_path = dirtype.path_getter()

    assert new_path == default_path


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_subdir_and_default(monkeypatch, tmp_path, dirtype):
    # have the env var point to a writable location,
    # where an 'astropy' subdir is missing, but can be created silently
    default_parent_dir = Path.home() / ".astropy"
    assert not default_parent_dir.exists()

    target_dir = tmp_path
    env_var = f"XDG_{dirtype.name}_HOME"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir / dirtype.as_str()

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
        ret = dirtype.path_getter()

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


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_subdir_but_default_exists(
    monkeypatch, tmp_path, dirtype
):
    # have the env var point to a writable location,
    # where an 'astropy' subdir is missing, but can be created silently
    default_parent_dir = Path.home() / ".astropy"

    target_dir = tmp_path
    env_var = f"XDG_{dirtype.name}_HOME"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir / dirtype.as_str()

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
        ret = dirtype.path_getter()

    assert ret.is_dir()

    # FIXME: should the environment variable *completely* take precedence here,
    # ignoring the existence of default_parent_dir ?
    # see https://github.com/astropy/astropy/issues/18791
    assert not expected_location.exists()


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_setup_file_exists(monkeypatch, tmp_path, dirtype):
    # check what happens if we request a location that's already
    # taken, but is a file
    default_path = dirtype.path_getter()
    target_dir = tmp_path / "subdir"
    target_dir.mkdir()
    expected_path = target_dir / "astropy"
    expected_path.touch()  # create a file
    env_var = f"XDG_{dirtype.name}_HOME"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_msg = (
        rf"^{env_var} is set to '{re.escape(str(target_dir))}', "
        rf"but this directory already contains a file under '{expected_path.name}', "
        r"where a directory was expected\. "
        r"This environment variable will be ignored\.$"
    )
    with pytest.warns(AstropyUserWarning, match=expected_msg):
        new_path = dirtype.path_getter()
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
    with pytest.raises(Exception) as _:
        # we except the context manager itself to raise an exception
        with paths.set_temp_cache(a):
            # this line should never run. If it does,
            # it'll raise an unexpected exception and fail the test
            1 / 0
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


def test_alias_read():
    from astropy.utils.data import conf

    with (
        set_temp_config(Path(__file__).with_name("data") / "alias"),
        pytest.warns(
            AstropyDeprecationWarning,
            match=(
                "^Config parameter 'name_resolve_timeout' in section "
                r"\[coordinates\.name_resolve\].*"
            ),
        ),
    ):
        assert conf.remote_timeout == 42


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

    with (
        set_temp_config(Path(__file__).with_name("data") / "deprecated"),
        pytest.warns(AstropyDeprecationWarning),
    ):
        conf.max_lines


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


@pytest.mark.parametrize("ctx_manager", [temporary_cache_path, temporary_config_path])
@pytest.mark.parametrize(
    "kwargs",
    [
        pytest.param({}, id="delete-implicit"),
        pytest.param({"delete": True}, id="delete-explicit"),
        pytest.param({"delete": False}, id="no-delete"),
    ],
)
def test_set_temp_dir_delete(ctx_manager, kwargs):
    # reason: old impl allows using existing directories, which creates all sorts of quirks:
    # - pre-existing content may be deleted
    # - multiple threads may use the same dir, which may get pulled from under them
    with ctx_manager(**kwargs) as tmp_path:
        target = tmp_path / "test"
        target.touch()

    expect_deletion = kwargs.get("delete", True)
    assert target.exists() is not expect_deletion
    assert tmp_path.exists() is not expect_deletion
