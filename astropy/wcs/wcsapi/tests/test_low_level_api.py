from pytest import raises

from astropy.wcs.wcsapi.low_level_api import validate_physical_types


def test_validate_physical_types():

    # Check valid cases
    validate_physical_types(['pos.eq.ra', 'pos.eq.ra'])
    validate_physical_types(['spect.dopplerVeloc.radio', 'custom:spam'])
    validate_physical_types(['time', None])

    # Make sure validation is case sensitive
    with raises(ValueError) as exc:
        validate_physical_types(['pos.eq.ra', 'Pos.eq.dec'])
    assert exc.value.args[0] == 'Invalid physical type: Pos.eq.dec'

    # Make sure nonsense types are picked up
    with raises(ValueError) as exc:
        validate_physical_types(['spam'])
    assert exc.value.args[0] == 'Invalid physical type: spam'
