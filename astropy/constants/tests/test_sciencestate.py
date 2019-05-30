import pytest


def test_version_match():
    from astropy import physical_constants, astronomical_constants
    pversion = physical_constants.get()
    import astropy.constants as const
    refpversion = const.h.__class__.__name__.lower()
    assert pversion == refpversion
    aversion = astronomical_constants.get()
    refaversion = const.M_sun.__class__.__name__.lower()
    assert aversion == refaversion


def test_previously_imported():
    from astropy import physical_constants, astronomical_constants
    import astropy.units
    import astropy.constants

    with pytest.raises(RuntimeError):
        physical_constants.set('codata2018')

    with pytest.raises(RuntimeError):
        astronomical_constants.set('iau2015')
