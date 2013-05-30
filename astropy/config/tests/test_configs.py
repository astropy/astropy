# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

try:
    # used by test_get_config_items
    from ..configuration import ConfigurationItem

    TESTCONF1 = ConfigurationItem('test1', 1, 'descr')
    TESTCONF2 = ConfigurationItem('test2', 2, 'descr')
except:
    # if this fails on import, don't worry - the tests will catch it.
    pass


def test_paths():
    from ..paths import get_config_dir, get_cache_dir

    assert 'astropy' in get_config_dir()
    assert 'astropy' in get_cache_dir()


def test_config_file():
    from ..configuration import get_config, reload_config, save_config
    from os.path import exists

    apycfg = get_config('astropy')
    assert apycfg.filename.endswith('astropy.cfg')

    cfgsec = get_config('astropy.config')
    assert cfgsec.depth == 1
    assert cfgsec.name == 'config'
    assert cfgsec.parent.filename.endswith('astropy.cfg')

    reload_config('astropy')

    # saving shouldn't change the file, because reload should have made sure it
    # is based on the current file.  But don't do it if there's no file
    if exists(apycfg.filename):
        save_config('astropy')


def test_configitem():
    from ..configuration import ConfigurationItem, get_config

    ci = ConfigurationItem('tstnm', 34, 'this is a Description')

    assert ci.module == 'astropy.config.tests.test_configs'
    assert ci() == 34
    assert ci.description == 'this is a Description'

    sec = get_config(ci.module)
    assert sec['tstnm'] == 34
    assert sec.comments['tstnm'][0] == ''
    assert sec.comments['tstnm'][1] == 'this is a Description'

    ci.description = 'updated Descr'
    ci.set(32)
    assert ci() == 32
    assert sec.comments['tstnm'][1] == 'updated Descr'

    # It's useful to go back to the default to allow other test functions to
    # call this one and still be in the default configuration.
    ci.description = 'this is a Description'
    ci.set(34)
    assert ci() == 34
    assert sec.comments['tstnm'][1] == 'this is a Description'


def test_configitem_save(tmpdir):
    from ..configuration import ConfigurationItem, get_config
    from shutil import copy

    ci = ConfigurationItem('tstnm2', 42, 'this is another Description')
    apycfg = get_config(ci.module)

    # now try saving

    while apycfg.parent is not apycfg:
        apycfg = apycfg.parent
    f = tmpdir.join('astropy.cfg')
    with open(f.strpath, 'w') as fd:
        apycfg.write(fd)
    with io.open(f.strpath, 'rU') as fd:
        lns = [x.strip() for x in fd.readlines()]

    assert 'tstnm2 = 42' in lns
    assert '# this is another Description' in lns

    oldfn = apycfg.filename
    try:
        # We had used LocalPath's `copy` method here, but it copies
        # the file in TEXT mode, destroying newlines on Windows.  We
        # need to do a binary copy.
        copy(tmpdir.join('astropy.cfg').strpath,
             tmpdir.join('astropy.cfg2').strpath)
        # tmpdir.join('astropy.cfg').copy(tmpdir.join('astropy.cfg2'))
        apycfg.filename = tmpdir.join('astropy.cfg2').realpath().strpath

        ci.set(30)
        ci.save()

        with io.open(apycfg.filename, 'rU') as f:
            lns = [x.strip() for x in f.readlines()]
            assert '[config.tests.test_configs]' in lns
            assert 'tstnm2 = 30' in lns

        ci.save(31)

        with io.open(apycfg.filename, 'rU') as f:
            lns = [x.strip() for x in f.readlines()]
            assert '[config.tests.test_configs]' in lns
            assert 'tstnm2 = 31' in lns

        # also try to save one that doesn't yet exist
        apycfg.filename = tmpdir.join('astropy.cfg3').realpath().strpath
        ci.save()

        with io.open(apycfg.filename, 'rU') as f:
            lns = [x.strip() for x in f.readlines()]
            assert '[config.tests.test_configs]' in lns
            assert 'tstnm2 = 30' in lns

    finally:
        apycfg.filename = oldfn


def test_configitem_types():
    from ..configuration import ConfigurationItem
    from ...tests.helper import pytest

    ci1 = ConfigurationItem('tstnm1', 34)
    assert isinstance(ci1(), int)

    ci2 = ConfigurationItem('tstnm2', 34.3)
    assert isinstance(ci2(), float)

    ci3 = ConfigurationItem('tstnm3', True)
    assert isinstance(ci3(), bool)

    ci4 = ConfigurationItem('tstnm4', 'astring')
    assert isinstance(ci4(), str)

    with pytest.raises(TypeError):
        ci1.set(34.3)
    ci2.set(12)  # this would should succeed as up-casting
    with pytest.raises(TypeError):
        ci3.set('fasd')
    with pytest.raises(TypeError):
        ci4.set(546.245)


def test_configitem_options(tmpdir):
    from ..configuration import ConfigurationItem, get_config
    from ...tests.helper import pytest

    cio = ConfigurationItem('tstnmo', ['op1', 'op2', 'op3'])
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
    with open(f.strpath, 'w') as fd:
        apycfg.write(fd)
    with io.open(f.strpath, 'rU') as fd:
        lns = [x.strip() for x in f.readlines()]

    assert '# Options: op1, op2, op3' in lns
    assert 'tstnmo = op2' in lns


def test_config_noastropy_fallback(monkeypatch, recwarn):
    """
    Tests to make sure configuration items fall back to their defaults when
    there's a problem accessing the astropy directory
    """
    from ...tests.helper import pytest
    from .. import paths, configuration

    # make sure the config directory is not searched
    monkeypatch.setenv('XDG_CONFIG_HOME', 'foo')
    monkeypatch.delenv('XDG_CONFIG_HOME')

    # make sure the _find_or_create_astropy_dir function fails as though the
    # astropy dir could not be accessed
    def osraiser(dirnm, linkto):
        raise OSError
    monkeypatch.setattr(paths, '_find_or_create_astropy_dir', osraiser)

    # also have to make sure the stored configuration objects are cleared
    monkeypatch.setattr(configuration, '_cfgobjs', {})

    with pytest.raises(OSError):
        # make sure the config dir search fails
        paths.get_config_dir()

    # now run the basic tests, and make sure the warning about no astropy
    # is present
    test_configitem()
    assert len(recwarn.list) > 0
    w = recwarn.pop()
    assert w.category == configuration.ConfigurationMissingWarning
    assert 'Configuration defaults will be used' in str(w.message)
    assert 'and configuration cannot be saved due to' in str(w.message)


def test_get_config_items():
    """ Checks if the get_config_items function is working correctly, using
    `ConfigurationItem` objects from this module.
    """
    import sys

    from ..configuration import get_config_items

    itemslocal = get_config_items(sys.modules['astropy.config.tests.test_configs'])
    itemslocalnone = get_config_items(None)
    itemsname = get_config_items('astropy.config.tests.test_configs')

    assert itemslocal == itemsname
    assert itemslocal == itemslocalnone

    assert 'TESTCONF1' in itemslocal.keys()
    assert 'TESTCONF2' in itemslocal.keys()


def test_configitem_setters():
    from ..configuration import ConfigurationItem

    ci = ConfigurationItem('tstnm12', 42, 'this is another Description')

    assert ci() == 42
    with ci.set_temp(45):
        assert ci() == 45
    assert ci() == 42

    ci.set(43)
    assert ci() == 43

    with ci.set_temp(46):
        assert ci() == 46

    assert ci() == 43
