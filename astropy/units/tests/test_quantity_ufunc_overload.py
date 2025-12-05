import numba
import pytest
import numpy as np
import astropy.units as u
from astropy.units.quantity_helper import helpers


@pytest.fixture(scope="module")
def ufunc_helpers():

    _ufunc_helpers = helpers.UFUNC_HELPERS

    yield helpers.UFUNC_HELPERS

    helpers.UFUNC_HELPERS = _ufunc_helpers


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
    [0.01, 0.02, 0.03] * u.mm,
]

_equivalencies = u.spectral()


def grating_equation(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    """Simple version of the grating equation which uses quantities."""
    return np.arcsin(m * wavelength / d)


@u.quantity_ufunc_overload
@numba.vectorize
def grating_equation_ufunc(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    return np.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation(
    ufunc_helpers,
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):

    try:
        result_expected = grating_equation(m, wavelength, d)
    except u.UnitTypeError:
        with pytest.raises(u.UnitConversionError):
            grating_equation_ufunc(m, wavelength, d)
        return

    result = grating_equation_ufunc(m, wavelength, d)

    assert np.allclose(result, result_expected)


@u.quantity_ufunc_overload(equivalencies=_equivalencies)
@numba.vectorize
def grating_equation_ufunc_equiv(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    return np.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation_equivalencies(
    ufunc_helpers,
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):

    wavelength_ = wavelength.to(u.nm, equivalencies=_equivalencies)

    result_expected = grating_equation(m, wavelength_, d)
    result = grating_equation_ufunc_equiv(m, wavelength, d)

    assert np.allclose(result, result_expected)


@u.quantity_ufunc_overload
@numba.vectorize
def grating_equation_ufunc_missing_units(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity,
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    return np.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation_missing_units(
    ufunc_helpers,
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):
    with pytest.raises(u.UnitConversionError):
        grating_equation_ufunc_missing_units(m, wavelength, d)
