# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import catch_warnings

import io
import os
import shutil

from ...utils.data import get_pkg_data_filename


def test_paths():
    from ..paths import get_config_dir, get_cache_dir

    assert 'astropy' in get_config_dir()
    assert 'astropy' in get_cache_dir()


def test_config_file():
    from ..configuration import get_config, reload_config, save_config

    apycfg = get_config('astropy')
    assert apycfg.filename.endswith('astropy.cfg')

    cfgsec = get_config('astropy.config')
    assert cfgsec.depth == 1
    assert cfgsec.name == 'config'
    assert cfgsec.parent.filename.endswith('astropy.cfg')

    reload_config('astropy')


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


def test_config_noastropy_fallback(monkeypatch):
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
    with catch_warnings() as w:
        test_configitem()
    assert len(w) > 0
    w = w[0]
    assert w.category == configuration.ConfigurationMissingWarning
    assert 'Configuration defaults will be used' in str(w.message)


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

    # Make sure it is reset even with Exception
    try:
        with ci.set_temp(47):
            raise Exception
    except:
        pass

    assert ci() == 43


def test_empty_config_file():
    from ..configuration import is_unedited_config_file

    fn = get_pkg_data_filename('data/empty.cfg')
    assert is_unedited_config_file(fn)

    fn = get_pkg_data_filename('data/not_empty.cfg')
    assert not is_unedited_config_file(fn)

    fn = get_pkg_data_filename('data/astropy.0.3.cfg')
    assert is_unedited_config_file(fn)
