# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

#namedtuple is needed for find_mod_objs so it can have a non-local module
from collections import namedtuple

from ...extern import six
from ...tests.helper import pytest
from .. import introspection
from ..introspection import (find_current_module, find_mod_objs,
                             isinstancemethod)


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


def test_isinstancemethod():
    """
    Note, this is an exact copy of the doctest in `isinstancemethod`'s
    docstring.

    It is included here as well so that it can be tested on Python 2 and 3
    which require very different implementations.  Once we enable running
    doctests on Python 3 this extra test can be dropped.
    """

    class MetaClass(type):
        def a_classmethod(cls): pass

    @six.add_metaclass(MetaClass)
    class MyClass(object):
        def an_instancemethod(self): pass

        @classmethod
        def another_classmethod(cls): pass

        @staticmethod
        def a_staticmethod(): pass

    assert not isinstancemethod(MyClass, MyClass.a_classmethod)
    assert not isinstancemethod(MyClass, MyClass.another_classmethod)
    assert not isinstancemethod(MyClass, MyClass.a_staticmethod)
    assert isinstancemethod(MyClass, MyClass.an_instancemethod)
