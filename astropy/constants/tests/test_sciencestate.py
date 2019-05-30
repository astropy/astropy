import sys
import pytest


@pytest.fixture
def clean_sys_modules():
    try:
        del sys.modules['astropy.units']
    except KeyError:
        pass
    assert 'astropy.units' not in sys.modules  # Sanity check.

    try:
        del sys.modules['astropy.constants']
    except KeyError:
        pass
    assert 'astropy.constants' not in sys.modules  # Sanity check.
    yield
    import astropy.units
    import astropy.constants


@pytest.mark.usefixtures('clean_sys_modules')
def test_invalid_config():
    """Test invalid config items"""
    from astropy import physical_constants, astronomical_constants
    if 'astropy.units' in sys.modules:
        del sys.modules['astropy.units']
    if 'astropy.constants' in sys.modules:
        del sys.modules['astropy.constants']
    with pytest.raises(ValueError):
        physical_constants.set('cooldata2014')

    with pytest.raises(ValueError):
        astronomical_constants.set('iau2006')


@pytest.mark.usefixtures('clean_sys_modules')
def test_previously_imported():
    from astropy import physical_constants, astronomical_constants
    import astropy.units
    import astropy.constants

    with pytest.raises(RuntimeError):
        physical_constants.set('codata2018')

    with pytest.raises(RuntimeError):
        astronomical_constants.set('iau2015')


@pytest.mark.usefixtures('clean_sys_modules')
def test_version_match():
    from astropy import physical_constants, astronomical_constants
    pversion = physical_constants.get()
    import astropy.constants as const
    refpversion = const.h.__class__.__name__.lower()
    assert pversion == refpversion
    aversion = astronomical_constants.get()
    refaversion = const.M_sun.__class__.__name__.lower()
    assert aversion == refaversion
