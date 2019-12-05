# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os
import sys
import subprocess

import pytest

from astropy.tests.helper import catch_warnings

from astropy.utils.data import get_pkg_data_filename
from astropy.config import configuration
from astropy.config import paths
from astropy.utils.exceptions import AstropyDeprecationWarning


def test_paths():
    assert 'astropy' in paths.get_config_dir()
    assert 'astropy' in paths.get_cache_dir()

    assert 'testpkg' in paths.get_config_dir(rootname='testpkg')
    assert 'testpkg' in paths.get_cache_dir(rootname='testpkg')


def test_set_temp_config(tmpdir, monkeypatch):
    monkeypatch.setattr(paths.set_temp_config, '_temp_path', None)

    orig_config_dir = paths.get_config_dir(rootname='astropy')
    temp_config_dir = str(tmpdir.mkdir('config'))
    temp_astropy_config = os.path.join(temp_config_dir, 'astropy')

    # Test decorator mode
    @paths.set_temp_config(temp_config_dir)
    def test_func():
        assert paths.get_config_dir(rootname='astropy') == temp_astropy_config

        # Test temporary restoration of original default
        with paths.set_temp_config() as d:
            assert d == orig_config_dir == paths.get_config_dir(rootname='astropy')

    test_func()

    # Test context manager mode (with cleanup)
    with paths.set_temp_config(temp_config_dir, delete=True):
        assert paths.get_config_dir(rootname='astropy') == temp_astropy_config

    assert not os.path.exists(temp_config_dir)


def test_set_temp_cache(tmpdir, monkeypatch):
    monkeypatch.setattr(paths.set_temp_cache, '_temp_path', None)

    orig_cache_dir = paths.get_cache_dir(rootname='astropy')
    temp_cache_dir = str(tmpdir.mkdir('cache'))
    temp_astropy_cache = os.path.join(temp_cache_dir, 'astropy')

    # Test decorator mode
    @paths.set_temp_cache(temp_cache_dir)
    def test_func():
        assert paths.get_cache_dir(rootname='astropy') == temp_astropy_cache

        # Test temporary restoration of original default
        with paths.set_temp_cache() as d:
            assert d == orig_cache_dir == paths.get_cache_dir(rootname='astropy')

    test_func()

    # Test context manager mode (with cleanup)
    with paths.set_temp_cache(temp_cache_dir, delete=True):
        assert paths.get_cache_dir(rootname='astropy') == temp_astropy_cache

    assert not os.path.exists(temp_cache_dir)


def test_set_temp_cache_resets_on_exception(tmpdir):
    """Test for regression of  bug #9704"""
    t = paths.get_cache_dir()
    a = tmpdir / 'a'
    with open(a, 'wt') as f:
        f.write("not a good cache\n")
    with pytest.raises(OSError):
        with paths.set_temp_cache(a):
            pass
    assert t == paths.get_cache_dir()


def test_config_file():
    from astropy.config.configuration import get_config, reload_config

    apycfg = get_config('astropy')
    assert apycfg.filename.endswith('astropy.cfg')

    cfgsec = get_config('astropy.config')
    assert cfgsec.depth == 1
    assert cfgsec.name == 'config'
    assert cfgsec.parent.filename.endswith('astropy.cfg')

    # try with a different package name, still inside astropy config dir:
    testcfg = get_config('testpkg', rootname='astropy')
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert '.astropy' in parts or 'astropy' in parts
    assert parts[-1] == 'testpkg.cfg'
    configuration._cfgobjs['testpkg'] = None # HACK

    # try with a different package name, no specified root name (should
    #   default to astropy):
    testcfg = get_config('testpkg')
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert '.astropy' in parts or 'astropy' in parts
    assert parts[-1] == 'testpkg.cfg'
    configuration._cfgobjs['testpkg'] = None # HACK

    # try with a different package name, specified root name:
    testcfg = get_config('testpkg', rootname='testpkg')
    parts = os.path.normpath(testcfg.filename).split(os.sep)
    assert '.testpkg' in parts or 'testpkg' in  parts
    assert parts[-1] == 'testpkg.cfg'
    configuration._cfgobjs['testpkg'] = None # HACK

    # try with a subpackage with specified root name:
    testcfg_sec = get_config('testpkg.somemodule', rootname='testpkg')
    parts = os.path.normpath(testcfg_sec.parent.filename).split(os.sep)
    assert '.testpkg' in parts or 'testpkg' in  parts
    assert parts[-1] == 'testpkg.cfg'
    configuration._cfgobjs['testpkg'] = None # HACK

    reload_config('astropy')


def test_configitem():

    from astropy.config.configuration import ConfigNamespace, ConfigItem, get_config

    ci = ConfigItem(34, 'this is a Description')

    class Conf(ConfigNamespace):
        tstnm = ci

    conf = Conf()

    assert ci.module == 'astropy.config.tests.test_configs'
    assert ci() == 34
    assert ci.description == 'this is a Description'

    assert conf.tstnm == 34

    sec = get_config(ci.module)
    assert sec['tstnm'] == 34

    ci.description = 'updated Descr'
    ci.set(32)
    assert ci() == 32

    # It's useful to go back to the default to allow other test functions to
    # call this one and still be in the default configuration.
    ci.description = 'this is a Description'
    ci.set(34)
    assert ci() == 34


def test_configitem_types():

    from astropy.config.configuration import ConfigNamespace, ConfigItem

    cio = ConfigItem(['op1', 'op2', 'op3'])

    class Conf(ConfigNamespace):
        tstnm1 = ConfigItem(34)
        tstnm2 = ConfigItem(34.3)
        tstnm3 = ConfigItem(True)
        tstnm4 = ConfigItem('astring')

    conf = Conf()

    assert isinstance(conf.tstnm1, int)
    assert isinstance(conf.tstnm2, float)
    assert isinstance(conf.tstnm3, bool)
    assert isinstance(conf.tstnm4, str)

    with pytest.raises(TypeError):
        conf.tstnm1 = 34.3
    conf.tstnm2 = 12  # this would should succeed as up-casting
    with pytest.raises(TypeError):
        conf.tstnm3 = 'fasd'
    with pytest.raises(TypeError):
        conf.tstnm4 = 546.245


def test_configitem_options(tmpdir):

    from astropy.config.configuration import ConfigNamespace, ConfigItem, get_config

    cio = ConfigItem(['op1', 'op2', 'op3'])

    class Conf(ConfigNamespace):
        tstnmo = cio

    conf = Conf()

    sec = get_config(cio.module)

    assert isinstance(cio(), str)
    assert cio() == 'op1'
    assert sec['tstnmo'] == 'op1'

    cio.set('op2')
    with pytest.raises(TypeError):
        cio.set('op5')
    assert sec['tstnmo'] == 'op2'

    # now try saving
    apycfg = sec
    while apycfg.parent is not apycfg:
        apycfg = apycfg.parent
    f = tmpdir.join('astropy.cfg')
    with open(f.strpath, 'wb') as fd:
        apycfg.write(fd)
    with open(f.strpath, 'r', encoding='utf-8') as fd:
        lns = [x.strip() for x in f.readlines()]

    assert 'tstnmo = op2' in lns


def test_config_noastropy_fallback(monkeypatch):
    """
    Tests to make sure configuration items fall back to their defaults when
    there's a problem accessing the astropy directory
    """

    # make sure the config directory is not searched
    monkeypatch.setenv('XDG_CONFIG_HOME', 'foo')
    monkeypatch.delenv('XDG_CONFIG_HOME')
    monkeypatch.setattr(paths.set_temp_config, '_temp_path', None)

    # make sure the _find_or_create_root_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto, pkgname=None):
        raise OSError
    monkeypatch.setattr(paths, '_find_or_create_root_dir', osraiser)

    # also have to make sure the stored configuration objects are cleared
    monkeypatch.setattr(configuration, '_cfgobjs', {})

    with pytest.raises(OSError):
        # make sure the config dir search fails
        paths.get_config_dir(rootname='astropy')

    # now run the basic tests, and make sure the warning about no astropy
    # is present
    with catch_warnings(configuration.ConfigurationMissingWarning) as w:
        test_configitem()
    assert len(w) == 1
    w = w[0]
    assert 'Configuration defaults will be used' in str(w.message)


def test_configitem_setters():

    from astropy.config.configuration import ConfigNamespace, ConfigItem

    class Conf(ConfigNamespace):
        tstnm12 = ConfigItem(42, 'this is another Description')

    conf = Conf()

    assert conf.tstnm12 == 42
    with conf.set_temp('tstnm12', 45):
        assert conf.tstnm12 == 45
    assert conf.tstnm12 == 42

    conf.tstnm12 = 43
    assert conf.tstnm12 == 43

    with conf.set_temp('tstnm12', 46):
        assert conf.tstnm12 == 46

    # Make sure it is reset even with Exception
    try:
        with conf.set_temp('tstnm12', 47):
            raise Exception
    except Exception:
        pass

    assert conf.tstnm12 == 43


def test_empty_config_file():
    from astropy.config.configuration import is_unedited_config_file

    def get_content(fn):
        with open(get_pkg_data_filename(fn), 'rt', encoding='latin-1') as fd:
            return fd.read()

    content = get_content('data/empty.cfg')
    assert is_unedited_config_file(content)

    content = get_content('data/not_empty.cfg')
    assert not is_unedited_config_file(content)

    content = get_content('data/astropy.0.3.cfg')
    assert is_unedited_config_file(content)

    content = get_content('data/astropy.0.3.windows.cfg')
    assert is_unedited_config_file(content)


class TestAliasRead:

    def setup_class(self):
        configuration._override_config_file = get_pkg_data_filename('data/alias.cfg')

    def test_alias_read(self):
        from astropy.utils.data import conf

        with catch_warnings() as w:
            conf.reload()
            assert conf.remote_timeout == 42

        assert len(w) == 1
        assert str(w[0].message).startswith(
            "Config parameter 'name_resolve_timeout' in section "
            "[coordinates.name_resolve]")

    def teardown_class(self):
        from astropy.utils.data import conf

        configuration._override_config_file = None
        conf.reload()


def test_configitem_unicode(tmpdir):

    from astropy.config.configuration import ConfigNamespace, ConfigItem, get_config

    cio = ConfigItem('ასტრონომიის')

    class Conf(ConfigNamespace):
        tstunicode = cio

    conf = Conf()

    sec = get_config(cio.module)

    assert isinstance(cio(), str)
    assert cio() == 'ასტრონომიის'
    assert sec['tstunicode'] == 'ასტრონომიის'


def test_warning_move_to_top_level():
    # Check that the warning about deprecation config items in the
    # file works.  See #2514
    from astropy import conf

    configuration._override_config_file = get_pkg_data_filename('data/deprecated.cfg')

    try:
        with catch_warnings(AstropyDeprecationWarning) as w:
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
    astropy_path = os.path.abspath(
        os.path.join(test_path, '..', '..', '..'))

    env = os.environ.copy()
    paths = [astropy_path]
    if env.get('PYTHONPATH'):
        paths.append(env.get('PYTHONPATH'))
    env['PYTHONPATH'] = os.pathsep.join(paths)

    for val in ['HOME', 'XDG_CONFIG_HOME']:
        if val in env:
            del env[val]

    retcode = subprocess.check_call(
        [sys.executable, '-c', 'import astropy'],
        env=env)

    assert retcode == 0


def test_unedited_template():
    # Test that the config file is written at most once
    config_dir = os.path.join(os.path.dirname(__file__), '..', '..')
    configuration.update_default_config('astropy', config_dir)
    assert configuration.update_default_config('astropy', config_dir) is False
