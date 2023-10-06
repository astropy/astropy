# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import os
import subprocess
import sys
from inspect import cleandoc

import pytest

from astropy.config import configuration, create_config_file, paths, set_temp_config
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyDeprecationWarning

OLD_CONFIG = {}


def setup_module():
    OLD_CONFIG.clear()
    OLD_CONFIG.update(configuration._cfgobjs)


def teardown_module():
    configuration._cfgobjs.clear()
    configuration._cfgobjs.update(OLD_CONFIG)


def test_paths():
    assert "astropy" in paths.get_config_dir()
    assert "astropy" in paths.get_cache_dir()

    assert "testpkg" in paths.get_config_dir(rootname="testpkg")
    assert "testpkg" in paths.get_cache_dir(rootname="testpkg")


def test_set_temp_config(tmp_path, monkeypatch):
    # Check that we start in an understood state.
    assert configuration._cfgobjs == OLD_CONFIG
    # Temporarily remove any temporary overrides of the configuration dir.
    monkeypatch.setattr(paths.set_temp_config, "_temp_path", None)

    orig_config_dir = paths.get_config_dir(rootname="astropy")
    (temp_config_dir := tmp_path / "config").mkdir()
    temp_astropy_config = temp_config_dir / "astropy"

    # Test decorator mode
    @paths.set_temp_config(temp_config_dir)
    def test_func():
        assert paths.get_config_dir(rootname="astropy") == str(temp_astropy_config)

        # Test temporary restoration of original default
        with paths.set_temp_config() as d:
            assert d == orig_config_dir == paths.get_config_dir(rootname="astropy")

    test_func()

    # Test context manager mode (with cleanup)
    with paths.set_temp_config(temp_config_dir, delete=True):
        assert paths.get_config_dir(rootname="astropy") == str(temp_astropy_config)

    assert not temp_config_dir.exists()
    # Check that we have returned to our old configuration.
    assert configuration._cfgobjs == OLD_CONFIG


def test_set_temp_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(paths.set_temp_cache, "_temp_path", None)

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
    assert "# css_urls = https://cdn.datatables.net/1.10.12/css/jquery.dataTables.css," in conf  # fmt: skip
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
    with open(f, encoding="utf-8") as fd:
        lns = [x.strip() for x in fd.readlines()]

    assert "tstnmo = op2" in lns


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


def test_config_noastropy_fallback(monkeypatch):
    """
    Tests to make sure configuration items fall back to their defaults when
    there's a problem accessing the astropy directory
    """

    # make sure the config directory is not searched
    monkeypatch.setenv("XDG_CONFIG_HOME", "foo")
    monkeypatch.delenv("XDG_CONFIG_HOME")
    monkeypatch.setattr(paths.set_temp_config, "_temp_path", None)

    # make sure the _find_or_create_root_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto, pkgname=None):
        raise OSError

    monkeypatch.setattr(paths, "_find_or_create_root_dir", osraiser)

    # also have to make sure the stored configuration objects are cleared
    monkeypatch.setattr(configuration, "_cfgobjs", {})

    with pytest.raises(OSError):
        # make sure the config dir search fails
        paths.get_config_dir(rootname="astropy")

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
