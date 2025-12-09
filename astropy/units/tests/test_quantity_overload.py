import copy
import math

import astropy.units as u
import numba
import numpy as np
import pytest
from astropy.units.quantity_helper import converters

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


@u.quantity_overload
def grating_equation_overload(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which uses np.emath.arcsin(),
    which is not implemented by astropy.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation(
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):

    try:
        result_expected = grating_equation(m, wavelength, d)
    except u.UnitTypeError:
        with pytest.raises(u.UnitConversionError):
            grating_equation_overload(m, wavelength, d)
        return

    result = grating_equation_overload(m, wavelength, d)

    assert np.allclose(result, result_expected)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation_ndarray_inputs(
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):
    with pytest.raises(u.UnitConversionError):
        grating_equation_overload(m, wavelength.value, d)


@u.quantity_overload(equivalencies=u.spectral())
def grating_equation_overload_equiv(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which uses np.emath.arcsin(),
    and allows the wavelength to be specified in terms on photon energy.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation_equivalencies(
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):

    wavelength_ = wavelength.to(u.nm, equivalencies=_equivalencies)

    result_expected = grating_equation(m, wavelength_, d)
    result = grating_equation_overload_equiv(m, wavelength, d)

    assert np.allclose(result, result_expected)


@u.quantity_overload
def grating_equation_overload_equiv_missing_unit(
    m: float | np.ndarray | u.Quantity[u.one],
    wavelength: u.Quantity,
    d: u.Quantity[u.um],
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which has a missing type annotation.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("m", _m)
@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
def test_grating_equation_missing_unit(
    m: u.Quantity,
    wavelength: u.Quantity,
    d: u.Quantity,
):
    with pytest.raises(u.UnitConversionError):
        grating_equation_overload_equiv_missing_unit(m, wavelength, d)
