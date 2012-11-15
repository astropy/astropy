# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import misc

#namedtuple is needed for find_mod_objs so it can have a non-local module
from collections import namedtuple
import warnings


def test_pkg_finder():
    """
    Tests that the `utils.misc.find_current_module` function works. Note that
    this also implicitly tests compat.misc._patched_getmodule
    """
    mod1 = 'astropy.utils.misc'
    mod2 = 'astropy.utils.tests.test_misc'
    mod3 = 'astropy.utils.tests.test_misc'
    assert misc.find_current_module(0).__name__ == mod1
    assert misc.find_current_module(1).__name__ == mod2
    assert misc.find_current_module(0, True).__name__ == mod3


def test_find_mod_objs():
    lnms, fqns, objs = misc.find_mod_objs('astropy')

    # this import  is after the above call intentionally to make sure
    # find_mod_objs properly imports astropy on its own
    import astropy

    # just check for astropy.test ... other things might be added, so we
    # shouldn't check that it's the only thing
    assert 'test' in lnms
    assert astropy.test in objs

    lnms, fqns, objs = misc.find_mod_objs('astropy.utils.tests.test_misc',
                                          onlylocals=False)
    assert 'namedtuple' in lnms
    assert 'collections.namedtuple' in fqns
    assert namedtuple in objs

    lnms, fqns, objs = misc.find_mod_objs('astropy.utils.tests.test_misc',
                                          onlylocals=True)
    assert 'namedtuple' not in lnms
    assert 'collections.namedtuple' not in fqns
    assert namedtuple not in objs


def test_find_current_mod():
    from sys import getrecursionlimit
    from ...tests.helper import pytest

    thismodnm = __name__

    assert misc.find_current_module(0) is misc
    assert misc.find_current_module(1).__name__ == thismodnm
    assert misc.find_current_module(getrecursionlimit() + 1) is None

    assert misc.find_current_module(0, True).__name__ == thismodnm
    assert misc.find_current_module(0, [misc]).__name__ == thismodnm
    assert misc.find_current_module(0, ['astropy.utils.misc']).__name__ == thismodnm

    with pytest.raises(ImportError):
        misc.find_current_module(0, ['faddfdsasewrweriopunjlfiurrhujnkflgwhu'])


def test_isiterable():
    from numpy import array

    assert misc.isiterable(2) is False
    assert misc.isiterable([2]) is True
    assert misc.isiterable([1, 2, 3]) is True
    assert misc.isiterable(array(2)) is False
    assert misc.isiterable(array([1, 2, 3])) is True


def test_deprecated_attribute():
    class DummyClass:
        def __init__(self):
            self._foo = 42

        def set_private(self):
            self._foo = 100

        foo = misc.deprecated_attribute('foo', '0.2')

    dummy = DummyClass()

    with warnings.catch_warnings(record=True) as w:
        warnings.resetwarnings()
        warnings.simplefilter('always')
        x = dummy.foo

    assert len(w) == 1
    assert str(w[0].message) == ("The foo attribute is deprecated and may be "
                                 "removed in a future version.")

    with warnings.catch_warnings(record=True) as w:
        warnings.resetwarnings()
        warnings.simplefilter('always')
        dummy.set_private()

    assert len(w) == 0
