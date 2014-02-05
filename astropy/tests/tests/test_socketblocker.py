import urllib
from ..helper import pytest


def test_outgoing_fails():
    with pytest.raises(IOError):
        urllib.urlopen('http://www.astropy.org')


def test_localconnect_succeeds():
    """
    Ensure that connections to localhost are allowed, since these are genuinely
    not remotedata.
    """
    urllib.urlopen('http://localhost').close()
    urllib.urlopen('http://127.0.0.1').close()

