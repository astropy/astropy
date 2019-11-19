# -*- coding: utf-8 -*-

# Licensed under a 3-clause BSD style license - see LICENSE.rst

# test helper.run_tests function
from astropy import test as run_tests

from astropy.tests import helper


# run_tests should raise ValueError when asked to run on a module it can't find
def test_module_not_found():
    with helper.pytest.raises(ValueError):
        run_tests(package='fake.module')


# run_tests should raise ValueError when passed an invalid pastebin= option
def test_pastebin_keyword():
    with helper.pytest.raises(ValueError):
        run_tests(pastebin='not_an_option')


def test_unicode_literal_conversion():
    assert isinstance('ångström', str)
