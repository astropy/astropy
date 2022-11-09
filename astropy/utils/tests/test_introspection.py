# Licensed under a 3-clause BSD style license - see LICENSE.rst

# namedtuple is needed for find_mod_objs so it can have a non-local module
from collections import namedtuple
from unittest import mock

import pytest
import yaml

from astropy.utils import introspection
from astropy.utils.introspection import find_current_module, find_mod_objs, minversion


def test_pkg_finder():
    """
    Tests that the `find_current_module` function works. Note that
    this also implicitly tests compat.misc._patched_getmodule
    """
    mod1 = "astropy.utils.introspection"
    mod2 = "astropy.utils.tests.test_introspection"
    mod3 = "astropy.utils.tests.test_introspection"
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
    assert find_current_module(0, ["astropy.utils.introspection"]).__name__ == thismodnm

    with pytest.raises(ImportError):
        find_current_module(0, ["faddfdsasewrweriopunjlfiurrhujnkflgwhu"])


def test_find_mod_objs():
    lnms, fqns, objs = find_mod_objs("astropy")

    # this import  is after the above call intentionally to make sure
    # find_mod_objs properly imports astropy on its own
    import astropy

    # just check for astropy.test ... other things might be added, so we
    # shouldn't check that it's the only thing
    assert "test" in lnms
    assert astropy.test in objs

    lnms, fqns, objs = find_mod_objs(__name__, onlylocals=False)
    assert "namedtuple" in lnms
    assert "collections.namedtuple" in fqns
    assert namedtuple in objs

    lnms, fqns, objs = find_mod_objs(__name__, onlylocals=True)
    assert "namedtuple" not in lnms
    assert "collections.namedtuple" not in fqns
    assert namedtuple not in objs


def test_minversion():
    import numpy

    good_versions = ["1.16", "1.16.1", "1.16.0.dev", "1.16dev"]
    bad_versions = ["100000", "100000.2rc1"]
    for version in good_versions:
        assert minversion(numpy, version)
        assert minversion("numpy", version)
    for version in bad_versions:
        assert not minversion(numpy, version)
        assert not minversion("numpy", version)

    assert minversion(yaml, "3.1")
    assert minversion("yaml", "3.1")


def test_find_current_module_bundle():
    """
    Tests that the `find_current_module` function would work if used inside
    an application bundle. Since we can't test this directly, we test what
    would happen if inspect.getmodule returned `None`, which is what happens
    inside PyInstaller and py2app bundles.
    """
    with mock.patch("inspect.getmodule", return_value=None):
        mod1 = "astropy.utils.introspection"
        mod2 = "astropy.utils.tests.test_introspection"
        mod3 = "astropy.utils.tests.test_introspection"
        assert find_current_module(0).__name__ == mod1
        assert find_current_module(1).__name__ == mod2
        assert find_current_module(0, True).__name__ == mod3
