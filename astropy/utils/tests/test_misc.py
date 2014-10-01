# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
import os
import pickle

#namedtuple is needed for find_mod_objs so it can have a non-local module
from collections import namedtuple

import numpy as np

from .. import data, misc
from ..exceptions import AstropyDeprecationWarning
from ...tests.helper import remote_data, catch_warnings


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
    assert misc.isiterable(2) is False
    assert misc.isiterable([2]) is True
    assert misc.isiterable([1, 2, 3]) is True
    assert misc.isiterable(np.array(2)) is False
    assert misc.isiterable(np.array([1, 2, 3])) is True


def test_deprecated_attribute():
    class DummyClass:
        def __init__(self):
            self._foo = 42

        def set_private(self):
            self._foo = 100

        foo = misc.deprecated_attribute('foo', '0.2')

    dummy = DummyClass()

    with catch_warnings(AstropyDeprecationWarning) as w:
        x = dummy.foo

    assert len(w) == 1
    assert str(w[0].message) == ("The foo attribute is deprecated and may be "
                                 "removed in a future version.")

    with catch_warnings() as w:
        dummy.set_private()

    assert len(w) == 0


# This needs to be defined outside of the test function, because we
# want to try to pickle it.
@misc.deprecated('100.0')
class TestA(object):
    """
    This is the class docstring.
    """
    def __init__(self):
        """
        This is the __init__ docstring
        """
        pass


def test_deprecated_class():
    orig_A = TestA.__bases__[0]

    # The only thing that should be different about the new class
    # is __doc__, __init__, __bases__ and __subclasshook__.
    for x in dir(orig_A):
        if x not in ('__doc__', '__init__', '__bases__', '__dict__',
                     '__subclasshook__'):
            assert getattr(TestA, x) == getattr(orig_A, x)

    with catch_warnings(AstropyDeprecationWarning) as w:
        TestA()

    assert len(w) == 1
    assert 'function' not in TestA.__doc__
    assert 'deprecated' in TestA.__doc__
    assert 'function' not in TestA.__init__.__doc__
    assert 'deprecated' in TestA.__init__.__doc__

    # Make sure the object is picklable
    pickle.dumps(TestA)


def test_deprecated_class_with_super():
    """
    Regression test for an issue where classes that used `super()` in their
    ``__init__`` did not actually call the correct class's ``__init__`` in the
    MRO.
    """

    @misc.deprecated('100.0')
    class TestB(object):
        def __init__(self, a, b):
            super(TestB, self).__init__()

    with catch_warnings(AstropyDeprecationWarning) as w:
        TestB(1, 2)

    assert len(w) == 1
    assert 'function' not in TestB.__doc__
    assert 'deprecated' in TestB.__doc__
    assert 'function' not in TestB.__init__.__doc__
    assert 'deprecated' in TestB.__init__.__doc__


def test_deprecated_static_and_classmethod():
    """
    Regression test for issue introduced by
    https://github.com/astropy/astropy/pull/2811 and mentioned also here:
    https://github.com/astropy/astropy/pull/2580#issuecomment-51049969
    where it appears that deprecated staticmethods didn't work on Python 2.6.
    """

    class A(object):
        @misc.deprecated('1.0')
        @staticmethod
        def B():
            pass

        @misc.deprecated('1.0')
        @classmethod
        def C(cls):
            pass

    with catch_warnings(AstropyDeprecationWarning) as w:
        A.B()

    assert len(w) == 1
    assert 'deprecated' in A.B.__doc__

    with catch_warnings(AstropyDeprecationWarning) as w:
        A.C()

    assert len(w) == 1
    assert 'deprecated' in A.C.__doc__


@remote_data
def test_api_lookup():
    strurl = misc.find_api_page('astropy.utils.misc', 'dev', False, timeout=3)
    objurl = misc.find_api_page(misc, 'dev', False, timeout=3)

    assert strurl == objurl
    assert strurl == 'http://devdocs.astropy.org/utils/index.html#module-astropy.utils.misc'


def test_skip_hidden():
    path = data._find_pkg_data_path('data')
    for root, dirs, files in os.walk(path):
        assert '.hidden_file.txt' in files
        assert 'local.dat' in files
        # break after the first level since the data dir contains some other
        # subdirectories that don't have these files
        break

    for root, dirs, files in misc.walk_skip_hidden(path):
        assert '.hidden_file.txt' not in files
        assert 'local.dat' in files
        break


def test_JsonCustomEncoder():
    assert json.dumps(np.arange(3), cls=misc.JsonCustomEncoder) == '[0, 1, 2]'
    assert json.dumps(1+2j, cls=misc.JsonCustomEncoder) == '[1.0, 2.0]'
    assert json.dumps(set([1, 2, 1]), cls=misc.JsonCustomEncoder) == '[1, 2]'
    assert json.dumps(b'hello world \xc3\x85',
                      cls=misc.JsonCustomEncoder) == '"hello world \\u00c5"'
    assert json.dumps({1: 2},
                      cls=misc.JsonCustomEncoder) == '{"1": 2}'  # default
