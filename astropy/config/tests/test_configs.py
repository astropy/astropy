# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import re
import subprocess
import sys
from contextlib import chdir, suppress
from inspect import cleandoc
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

import pytest

from astropy import conf
from astropy.config import (
    configuration,
    create_config_file,
    paths,
    set_temp_config,
    temporary_cache_dir_path,
    temporary_config_dir_path,
)
from astropy.config.configuration import ConfigItem, ConfigNamespace
from astropy.config.paths import (
    _DirectoryFinder,
    _DirType,
    _Envvar,
    _Namespace,
    _resolve,
    _SpecSource,
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
    "s, expected",
    [
        ("foo", _Namespace(root="foo", fragments=())),
        ("foo.bar", _Namespace(root="foo", fragments=("bar",))),
        ("foo.bar.baz", _Namespace(root="foo", fragments=("bar", "baz"))),
        ("foo_bar.baz", _Namespace(root="foo_bar", fragments=("baz",))),
        ("foo.bar_baz", _Namespace(root="foo", fragments=("bar_baz",))),
        ("foo-bar.baz", _Namespace(root="foo-bar", fragments=("baz",))),
        ("foo.bar-baz", _Namespace(root="foo", fragments=("bar-baz",))),
        ("0.1.2", _Namespace(root="0", fragments=("1", "2"))),
    ],
)
def test_namespace_from_str(s, expected):
    ns = _Namespace.from_str(s)
    assert ns == expected

    # check roundtrips
    assert ns.join() == s


@pytest.mark.parametrize("s", ["", "()", "foo/", "foo\\", "foo/bar", "foo\\bar"])
def test_namespace_from_invalid_str(s):
    with pytest.raises(ValueError, match=r"^Found invalid namespace elements\.$"):
        _Namespace.from_str(s)


@pytest.mark.parametrize(
    "spec, regexp",
    [
        pytest.param(
            _SpecSource.ASTROPY, r"^ASTROPY_(?P<dirtype>[A-Z]+)_DIR$", id="ASTROPY_"
        ),
        pytest.param(_SpecSource.XDG, r"^XDG_(?P<dirtype>[A-Z]+)_HOME$", id="XDG_"),
    ],
)
@pytest.mark.parametrize("dirtype", _DirType)
def test_direnvvar_impl(spec, regexp, dirtype):
    ev = _Envvar(spec=spec, dirtype=dirtype)
    assert (m := re.match(regexp, ev.name)) is not None
    assert m.group("dirtype") == dirtype.name


# TODO: split into smaller tests
def test_resolve_envvar(monkeypatch, tmp_path: Path):
    uid = str(uuid4())
    d1 = tmp_path / uid
    value = str(d1)
    monkeypatch.setenv("TEST_ENV_VAR", value)
    err_template = "TEST_ENV_VAR is set to {value}, {reason}. This environment variable will be ignored."

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
    envvar = _Envvar(spec=_SpecSource.XDG, dirtype=dirtype)
    uid = str(uuid4())
    missing_path = tmp_path / uid
    monkeypatch.setenv(envvar.name, str(missing_path))

    err = envvar.resolve()
    assert isinstance(err, str)
    assert err == (
        f"{envvar.name} is set to {str(missing_path)}, "
        "but no such file or directory was found. "
        "This environment variable will be ignored."
    )


@pytest.mark.parametrize("dirtype", _DirType)
def test_xdgenvvar_resolve_success(dirtype, monkeypatch, tmp_path):
    envvar = _Envvar(spec=_SpecSource.XDG, dirtype=dirtype)
    expected = tmp_path
    monkeypatch.setenv(envvar.name, str(expected))

    p = envvar.resolve()
    assert isinstance(p, Path)
    assert p.is_absolute()
    assert p == expected


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_find_default(dirtype):
    df = _DirectoryFinder(dirtype)
    de = df.find_directory_elements("mynamespace")
    assert de.base_node == df.default_base_node()


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_directoryfinder_override_success(dirtype, tmp_path: Path):
    namespace = "bacon"
    df = _DirectoryFinder(dirtype, overrides={_Namespace.from_str(namespace): tmp_path})
    assert df.find_namespaced_node(namespace) == tmp_path / namespace


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

    err_template = "{var} is set to {val}, which is a file (expected a directory)."

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
        res = df.find_directory_elements("astropy")
    assert res.base_node == df.default_base_node()
    p.unlink()

    # dir
    p.mkdir()
    res = df.find_directory_elements("astropy")
    assert res.base_node == p


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

    err_template = "{var} is set to {val}, which is a file (expected a directory)."
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
        res = df.find_directory_elements("astropy")
    assert res.base_node == df.default_base_node()
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
        res = df.find_directory_elements("astropy")
    assert res.base_node == p1
    p0.unlink()
    p1.rmdir()

    # dir + file
    p0.mkdir()
    p1.touch()
    res = df.find_directory_elements("astropy")
    assert res.base_node == p0
    p0.rmdir()
    p1.unlink()

    # dir + dir
    p0.mkdir()
    p1.mkdir()
    res = df.find_directory_elements("astropy")
    assert res.base_node == p0


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
        rf"^{env_var} is set to {re.escape(str(target_dir))}, "
        rf"but no such file or directory was found\. "
        r"This environment variable will be ignored\.$"
    )
    with pytest.warns(AstropyUserWarning, match=expected_msg):
        new_path = dirtype.path_getter()

    assert new_path == default_path


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.parametrize(
    "setup",
    [
        pytest.param(lambda defloc: None, id="default-missing"),
        pytest.param(lambda defloc: defloc.mkdir(parents=True), id="default-exists"),
    ],
)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_env_variables_missing_subdir_and_default(
    monkeypatch, tmp_path, dirtype: _DirType, setup
):
    # have the env var point to a writable location,
    # where an 'astropy' subdir is missing, but can be created silently
    default_parent_dir = Path.home() / ".astropy"
    assert not default_parent_dir.exists()

    target_dir = tmp_path
    env_var = f"XDG_{dirtype.name}_HOME"
    monkeypatch.setenv(env_var, str(target_dir))

    expected_location = target_dir / "astropy"
    default_location = default_parent_dir

    setup(default_location)
    previous_state = default_location.exists()

    assert not expected_location.exists()
    ret = dirtype.path_getter()

    assert ret.is_dir()
    assert ret.is_absolute()
    assert not ret.is_symlink()
    assert ret == expected_location
    assert default_location.exists() is previous_state


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

    with pytest.raises(FileExistsError):
        dirtype.path_getter()


@pytest.mark.usefixtures("ignore_config_paths_global_state")
@pytest.mark.parametrize("dirtype", _DirType)
def test_set_temp_config(tmp_path, dirtype):
    # Check that we start in an understood state.
    assert configuration._cfgobjs == OLD_CONFIG

    orig_dir = dirtype.path_getter(rootname="astropy")
    (temp_dir := tmp_path / "test").mkdir()
    temp_astropy_subdir = temp_dir / "astropy"

    # Test decorator mode
    @dirtype.legacy_context_manager(temp_dir)
    def test_func():
        assert dirtype.path_getter(rootname="astropy") == temp_astropy_subdir

        # Test temporary restoration of original default
        with dirtype.legacy_context_manager() as d:
            assert str(dirtype.path_getter(rootname="astropy")) == str(orig_dir)
            assert Path(d) == orig_dir

    test_func()

    # Test context manager mode (with cleanup)
    with dirtype.legacy_context_manager(temp_dir, delete=True):
        assert dirtype.path_getter(rootname="astropy") == temp_astropy_subdir

    assert not temp_dir.exists()
    # Check that we have returned to our old configuration.
    assert configuration._cfgobjs == OLD_CONFIG


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.parametrize(
    "create_kwargs",
    [
        pytest.param({}, id="implicit"),
        pytest.param({"ensure_exists": True}, id="explicit"),
    ],
)
def test_get_path_without_creation(dirtype: _DirType, create_kwargs, tmp_path):
    with dirtype.legacy_context_manager(tmp_path) as tmp_dir:
        assert not os.path.exists(tmp_dir)
        dirtype.path_getter("astropy", ensure_exists=False)
        assert not os.path.exists(tmp_dir)
        dirtype.path_getter("astropy", **create_kwargs)
        assert os.path.exists(tmp_dir)


@pytest.mark.parametrize("dirtype", _DirType)
@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_fallback_to_legacy_dir_if_found(dirtype):
    df = _DirectoryFinder(dirtype)
    namespace = str(uuid4())
    node = df.find_namespaced_node(namespace)
    legacy_node = df.legacy_default_base_node(namespace)
    assert node != legacy_node

    legacy_node.mkdir(parents=True)
    assert df.find_namespaced_node(namespace) == legacy_node


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

    conf_file = tmp_path / "astropy" / "astropy.cfg"
    with set_temp_config(tmp_path) as foo:
        from astropy.config.configuration import generate_config

        generate_config("astropy")

    assert conf_file.is_file()
    check_config(conf_file.read_text())


class _MyPackageNamespace(ConfigNamespace):
    pass


class _RecursiveTestConf(_MyPackageNamespace):
    ti = ConfigItem(5, "this is a Description")


def test_generate_config_subclasses(tmp_path):
    """Test that generate_config works with subclasses of ConfigNamespace."""

    conf_file = tmp_path / "astropy" / "astropy.cfg"
    with set_temp_config(tmp_path):
        from astropy.config.configuration import generate_config

        generate_config("astropy")

    assert conf_file.is_file()
    assert "# ti = 5" in conf_file.read_text()


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_create_config_file_1(tmp_path, caplog):
    conf_file = tmp_path / "astropy" / "astropy.cfg"
    with set_temp_config(tmp_path):
        create_config_file("astropy")

    # check that the config file has been created
    assert conf_file.is_file()
    check_config(conf_file.read_text())
    assert (
        "The configuration file has been successfully written"
        in caplog.records[0].message
    )


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_create_config_file_no_overwrite(tmp_path, caplog):
    conf_file = tmp_path / "astropy" / "astropy.cfg"
    with set_temp_config(tmp_path):
        create_config_file("astropy")
    caplog.clear()

    # check that the config file has been created
    assert conf_file.is_file()
    check_config(conf1 := conf_file.read_text())

    # now modify the config file
    conf2 = conf1.replace("# unicode_output = False", "unicode_output = True")
    assert conf2 != conf1
    conf_file.write_text(conf2)

    with set_temp_config(tmp_path):
        create_config_file("astropy")

    # check that the config file has not been overwritten since it was modified
    assert conf_file.read_text() == conf2
    assert (
        "The configuration file already exists and seems to have been customized"
        in caplog.records[0].message
    )


@pytest.mark.usefixtures("ignore_config_paths_global_state")
def test_create_config_file_overwrite(tmp_path, caplog):
    conf_file = tmp_path / "astropy" / "astropy.cfg"
    with set_temp_config(tmp_path):
        create_config_file("astropy")
    caplog.clear()

    # check that the config file has been created
    assert conf_file.is_file()
    check_config(conf1 := conf_file.read_text())

    # now modify the config file
    conf2 = conf1.replace("# unicode_output = False", "unicode_output = True")
    assert conf2 != conf1
    conf_file.write_text(conf2)

    with set_temp_config(tmp_path):
        create_config_file("astropy", overwrite=True)

    # check that the config file has been overwritten
    assert conf_file.read_text() == conf1
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


@pytest.mark.no_optimized_interpreter
def test_help(capsys):

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


@pytest.mark.no_optimized_interpreter
def test_help_invalid_config_item():

    with pytest.raises(
        KeyError,
        match=(
            "'bad_name' is not among configuration items "
            r"\('unicode_output', 'use_color', 'max_lines', 'max_width'\)"
        ),
    ):
        conf.help("bad_name")


@pytest.mark.only_optimized_interpreter
def test_help_in_optimized_python():
    from astropy import conf

    with pytest.raises(
        RuntimeError,
        match=("The help method is not available under Python's optimized mode."),
    ):
        conf.help()


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


@pytest.mark.parametrize(
    "ctx_manager", [temporary_cache_dir_path, temporary_config_dir_path]
)
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
    with ctx_manager(namespace="my.namespace", **kwargs) as tmp_path:
        target = tmp_path / "test"
        target.touch()

    expect_deletion = kwargs.get("delete", True)
    assert target.exists() is not expect_deletion
    assert tmp_path.exists() is not expect_deletion


@pytest.mark.usefixtures("ignore_config_paths_global_state")
class TestEnvvarRegressions:
    def test_astropy_cache_defined_but_invalid(self, monkeypatch, tmp_path):
        tmp_file = tmp_path / "file"
        tmp_file.touch()
        monkeypatch.setenv("ASTROPY_CACHE_DIR", str(tmp_file))
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))

        with pytest.warns(
            AstropyUserWarning,
            match=(
                rf"^ASTROPY_CACHE_DIR is set to {re.escape(str(tmp_file))}, "
                r"which is a file \(expected a directory\)\. "
                r"This environment variable will be ignored\.$"
            ),
        ):
            res = paths.get_cache_dir_path()

        assert res == tmp_path / "astropy"

    @pytest.mark.usefixtures("ignore_config_paths_global_state")
    def test_config_objs_leak(self, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

        ref = 0x4D3D3D3
        conf.max_width = ref

        with TemporaryDirectory() as td, set_temp_config(td):
            assert conf.max_width == None

        with temporary_config_dir_path(namespace="astropy"):
            assert conf.max_width == None

    @pytest.mark.usefixtures("ignore_config_paths_global_state")
    def test_resources_cleanup_overrides(self, monkeypatch, tmp_path):
        monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

        ref = paths.get_cache_dir()

        class UniqueException(Exception): ...

        with suppress(UniqueException), temporary_cache_dir_path(namespace="astropy"):
            raise UniqueException

        assert paths.get_cache_dir() == ref

    @pytest.mark.usefixtures("ignore_config_paths_global_state")
    def test_resources_cleanup_cfgobjs(self, monkeypatch, tmp_path):

        _ = conf.max_width  # populate _cfgobjs once
        ref = sorted(configuration._cfgobjs)

        class UniqueException(Exception):
            pass

        with suppress(UniqueException), temporary_config_dir_path(namespace="astropy"):
            raise UniqueException("simulated failure inside the with-block")

        res = sorted(configuration._cfgobjs)
        assert res == ref
