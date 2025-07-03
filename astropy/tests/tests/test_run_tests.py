# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

# test helper.run_tests function
from astropy import test as run_tests


# run_tests should raise ValueError when asked to run on a module it can't find
def test_module_not_found():
    with pytest.raises(ValueError), pytest.deprecated_call(match="The test runner"):
        run_tests(package="fake.module")


# run_tests should raise ValueError when passed an invalid pastebin= option
def test_pastebin_keyword():
    with pytest.raises(ValueError), pytest.deprecated_call(match="The test runner"):
        run_tests(pastebin="not_an_option")
