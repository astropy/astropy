import pytest
import numpy as np
import astropy.units as u

_m = [
    1,
    1 * u.one,
    np.array([-1, 0, 1]),
]

_wavelength = [
    500 * u.nm,
    [300, 400, 500] * u.AA,
    2 * u.eV,
]

_d = [
    1 * u.um,
    [0.01, 0.02] * u.mm,
]


def _grating_equation(
    m: float | np.ndarray | u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
) -> u.Quantity:
    """Simple version of the grating equation which uses quantities."""
    return np.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation(m: u.Quantity, wavelength: u.Quantity, d: u.Quantity):

    @u.quantity_ufunc_overload
    def grating_equation_ufunc(
        m: float | np.ndarray | u.Quantity[u.one],
        wavelength: u.Quantity[u.um],
        d: u.Quantity[u.um],
    ) -> u.Quantity[u.rad]:
        return np.arcsin(m * wavelength / d)

    try:
        result_expected = _grating_equation(m, wavelength, d)
    except u.UnitConversionError:
        with pytest.raises(u.UnitConversionError):
            grating_equation_ufunc(m, wavelength, d)

    result = grating_equation_ufunc(m, wavelength, d)

    assert np.all(result == result_expected)
