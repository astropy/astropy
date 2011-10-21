# this test doesn't actually use any online data, it should just be skipped
# by run_tests because it has the remote_data decorator.

from astropy.tests.helper import remote_data

@remote_data
def test_skip_me():
    assert True
