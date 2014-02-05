import urllib
from ...tests.helper import pytest

def test_outgoing_fails():
    with pytest.raises(IOError):
        urllib.urlopen('http://www.astropy.org')

def test_localconnect_succeeds():
    """
    Ensure that connections to localhost are allowed, since these are genuinely
    not remotedata.
    """
    urllib.urlopen('http://localhost')
    urllib.urlopen('http://127.0.0.1')

