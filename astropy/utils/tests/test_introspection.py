# Licensed under a 3-clause BSD style license - see LICENSE.rst

# namedtuple is needed for find_mod_objs so it can have a non-local module
from collections import namedtuple

import pytest

from .. import introspection
from ..introspection import (find_current_module, find_mod_objs,
                             isinstancemethod, minversion)


def test_pkg_finder():
    """
    Tests that the `find_current_module` function works. Note that
    this also implicitly tests compat.misc._patched_getmodule
    """
    mod1 = 'astropy.utils.introspection'
    mod2 = 'astropy.utils.tests.test_introspection'
    mod3 = 'astropy.utils.tests.test_introspection'
    assert find_current_module(0).__name__ == mod1
    assert find_current_module(1).__name__ == mod2
    assert find_current_module(0, True).__name__ == mod3


def test_find_current_mod():
    from sys import getrecursionlimit

    thismodnm = __name__

    assert find_current_module(0) is introspection
    assert find_current_module(1).__name__ == thismodnm
    assert find_current_module(getrecursionlimit() + 1) is None

    assert find_current_module(0, True).__name__ == thismodnm
    assert find_current_module(0, [introspection]).__name__ == thismodnm
    assert find_current_module(0, ['astropy.utils.introspection']).__name__ == thismodnm

    with pytest.raises(ImportError):
        find_current_module(0, ['faddfdsasewrweriopunjlfiurrhujnkflgwhu'])


def test_find_mod_objs():
    lnms, fqns, objs = find_mod_objs('astropy')

    # this import  is after the above call intentionally to make sure
    # find_mod_objs properly imports astropy on its own
    import astropy

    # just check for astropy.test ... other things might be added, so we
    # shouldn't check that it's the only thing
    assert 'test' in lnms
    assert astropy.test in objs

    lnms, fqns, objs = find_mod_objs(__name__, onlylocals=False)
    assert 'namedtuple' in lnms
    assert 'collections.namedtuple' in fqns
    assert namedtuple in objs

    lnms, fqns, objs = find_mod_objs(__name__, onlylocals=True)
    assert 'namedtuple' not in lnms
    assert 'collections.namedtuple' not in fqns
    assert namedtuple not in objs


def test_minversion():
    from types import ModuleType
    test_module = ModuleType(str("test_module"))
    test_module.__version__ = '0.12.2'
    good_versions = ['0.12', '0.12.1', '0.12.0.dev']
    bad_versions = ['1', '1.2rc1']
    for version in good_versions:
        assert minversion(test_module, version)
    for version in bad_versions:
        assert not minversion(test_module, version)


def test_minversion_with_dev_versions():
    """Test minversion with development version strings.

    Regression test for issue where LooseVersion fails to compare
    versions like '1.14.3' with '1.14dev'.
    See: https://bugs.python.org/issue30272
    """
    from types import ModuleType

    # Test case 1: Version with dev suffix
    test_module = ModuleType(str("test_module"))
    test_module.__version__ = '1.14.3'
    assert minversion(test_module, '1.14dev')
    assert minversion(test_module, '1.14')
    assert minversion(test_module, '1.14.0')
    assert minversion(test_module, '1.14.2')
    assert not minversion(test_module, '1.14.4')
    assert not minversion(test_module, '1.15')

    # Test case 2: Module version is dev
    test_module.__version__ = '1.14dev'
    assert minversion(test_module, '1.13')
    assert not minversion(test_module, '1.14')
    assert not minversion(test_module, '1.14.0')

    # Test case 3: Various rc/alpha/beta versions
    test_module.__version__ = '2.0.1'
    assert minversion(test_module, '2.0rc1')
    assert minversion(test_module, '2.0a1')
    assert minversion(test_module, '2.0b1')
    assert minversion(test_module, '2.0.0')
    assert not minversion(test_module, '2.1rc1')

    # Test case 4: Inclusive vs non-inclusive comparison
    test_module.__version__ = '1.14'
    assert minversion(test_module, '1.14', inclusive=True)
    assert not minversion(test_module, '1.14', inclusive=False)
    assert minversion(test_module, '1.13', inclusive=False)
