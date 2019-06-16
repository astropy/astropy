import pytest

from astropy import physical_constants, astronomical_constants
import astropy.constants as const


def test_version_match():
    pversion = physical_constants.get()
    refpversion = const.h.__class__.__name__.lower()
    assert pversion == refpversion
    aversion = astronomical_constants.get()
    refaversion = const.M_sun.__class__.__name__.lower()
    assert aversion == refaversion


def test_previously_imported():
    with pytest.raises(RuntimeError):
        physical_constants.set('codata2018')

    with pytest.raises(RuntimeError):
        astronomical_constants.set('iau2015')
