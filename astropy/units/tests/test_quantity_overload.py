import numpy as np
import pytest

import astropy.units as u

_wavelength = [
    500 * u.nm,
    [300, 400, 500] * u.eV,
]

_d = [
    1 * u.um,
    [0.01, 0.02, 0.03] * u.mm,
]

_m = [
    1,
    [-1, 0, 1] * u.one,
]

_equivalencies = u.spectral()


def grating_equation(
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray | u.Quantity[u.one] = 1,
) -> u.Quantity[u.rad]:
    """Simple version of the grating equation which uses quantities."""
    return np.arcsin(m * wavelength / d)


@u.quantity_overload
def grating_equation_overload(
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray = 1,
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which uses np.emath.arcsin(),
    which is not implemented by astropy.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    try:
        result_expected = grating_equation(wavelength, d, m=m)
    except u.UnitTypeError:
        with pytest.raises(u.UnitConversionError):
            grating_equation_overload(wavelength, d, m=m)
        return

    result = grating_equation_overload(wavelength, d, m=m)

    assert np.allclose(result, result_expected)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_ndarray_inputs(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    with pytest.raises(u.UnitConversionError):
        grating_equation_overload(wavelength.value, d, m=m)


@u.quantity_overload(equivalencies=u.spectral())
def grating_equation_overload_equiv(
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray = 1,
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which uses np.emath.arcsin(),
    and allows the wavelength to be specified in terms on photon energy.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_equivalencies(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    wavelength_ = wavelength.to(u.nm, equivalencies=_equivalencies)

    result_expected = grating_equation(wavelength_, d, m=m)
    result = grating_equation_overload_equiv(wavelength, d, m=m)

    assert np.allclose(result, result_expected)


@u.quantity_overload
def grating_equation_overload_missing_unit(
    wavelength: u.Quantity,
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray = 1,
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which has a missing type annotation.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_missing_unit(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    with pytest.raises(u.UnitConversionError):
        grating_equation_overload_missing_unit(wavelength, d, m=m)


@u.quantity_overload
def grating_equation_overload_multiple_units(
    wavelength: u.Quantity[u.um] | u.Quantity[u.m],
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray = 1,
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which has multiple units specified.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_multiple_units(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    with pytest.raises(
        expected_exception=TypeError,
        match="Only one unit specification allowed, got .*",
    ):
        grating_equation_overload_multiple_units(wavelength, d, m=m)


@u.quantity_overload(equivalencies=_equivalencies)
def grating_equation_overload_variational(
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
    *args: None,
    m: int | np.ndarray = 1,
    **kwargs: None,
) -> u.Quantity[u.rad]:
    """
    A version of the grating equation which has variational arguments.
    """
    return np.emath.arcsin(m * wavelength / d)


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_variational(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    wavelength_ = wavelength.to(u.nm, equivalencies=_equivalencies)

    result_expected = grating_equation(wavelength_, d, m=m)
    result = grating_equation_overload_variational(wavelength, d, None, m=m, foo=None)

    assert np.allclose(result, result_expected)


@u.quantity_overload(equivalencies=_equivalencies)
def grating_equation_overload_multiple_return(
    wavelength: u.Quantity[u.um],
    d: u.Quantity[u.um],
    *,
    m: int | np.ndarray | u.Quantity[u.one] = 1,
) -> tuple[u.Quantity[u.rad], float | np.ndarray]:
    """
    A version of the grating equation which has multiple return values.
    """
    return np.emath.arcsin(m * wavelength / d), 0


@pytest.mark.parametrize("wavelength", _wavelength)
@pytest.mark.parametrize("d", _d)
@pytest.mark.parametrize("m", _m)
def test_grating_equation_multiple_return(
    wavelength: u.Quantity,
    d: u.Quantity,
    m: u.Quantity,
):
    wavelength_ = wavelength.to(u.nm, equivalencies=_equivalencies)

    result_expected = grating_equation(wavelength_, d, m=m)
    result_1, result_2 = grating_equation_overload_multiple_return(wavelength, d, m=m)

    assert np.allclose(result_1, result_expected)
    assert np.allclose(result_2, 0)
