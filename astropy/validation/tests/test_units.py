import pytest

from astropy import units as u
from astropy.validation.errors import PhysicalInconsistencyError
from astropy.validation.units import validate_units_consistency


def test_valid_flux_distance():
    assert validate_units_consistency(1 * u.erg / u.s / u.cm**2, 10 * u.pc, strict=True)


def test_invalid_flux_distance():
    with pytest.raises(PhysicalInconsistencyError):
        validate_units_consistency(1 * u.erg, 10 * u.s, strict=True)
