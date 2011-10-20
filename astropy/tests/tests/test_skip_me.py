# this test doesn't actually use any online data, it should just be skipped
# by run_tests because it has the helper.big_data decorator.

from astropy.tests import helper

@helper.big_data
def test_skip_me():
    assert True
