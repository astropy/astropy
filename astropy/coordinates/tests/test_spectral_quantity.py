import pytest
import numpy as np

from numpy.testing import assert_allclose

from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose

from astropy.coordinates.spectral_quantity import SpectralQuantity

SPECTRAL_UNITS = (u.GHz, u.micron, u.keV, (1 / u.nm).unit, u.km / u.s)


class TestSpectralQuantity:

    @pytest.mark.parametrize('unit', SPECTRAL_UNITS)
    def test_init_value(self, unit):
        SpectralQuantity(1, unit=unit)

    @pytest.mark.parametrize('unit', SPECTRAL_UNITS)
    def test_init_quantity(self, unit):
        SpectralQuantity(1 * unit)

    @pytest.mark.parametrize('unit', SPECTRAL_UNITS)
    def test_init_spectralquantity(self, unit):
        SpectralQuantity(SpectralQuantity(1, unit=unit))

    @pytest.mark.parametrize('unit', (u.kg, u.byte))
    def test_init_invalid(self, unit):
        with pytest.raises(u.UnitsError, match='SpectralQuantity instances require units'):
            SpectralQuantity(1, unit=unit)
        with pytest.raises(u.UnitsError, match='SpectralQuantity instances require units'):
            SpectralQuantity(1 * unit)

    @pytest.mark.parametrize(('unit1', 'unit2'), zip(SPECTRAL_UNITS, SPECTRAL_UNITS))
    def test_spectral_conversion(self, unit1, unit2):
        sq1 = SpectralQuantity(1 * unit1)
        sq2 = sq1.to(unit2)
        sq3 = sq2.to(str(unit1))  # check that string units work
        assert isinstance(sq2, SpectralQuantity)
        assert isinstance(sq3, SpectralQuantity)
        assert_quantity_allclose(sq1, sq3)

    def test_doppler_conversion(self):

        sq1 = SpectralQuantity(1 * u.km / u.s, doppler_convention='optical', doppler_rest=500 * u.nm)

        sq2 = sq1.to(u.m / u.s)
        assert_allclose(sq2.value, 1000)

        sq3 = sq1.to(u.m / u.s, doppler_convention='radio')
        assert_allclose(sq3.value, 999.996664)

        sq4 = sq1.to(u.m / u.s, doppler_convention='relativistic')
        assert_allclose(sq4.value, 999.998332)

        sq5 = sq1.to(u.m / u.s, doppler_rest=499.9 * u.nm)
        assert_allclose(sq5.value, 60970.685737)

        val5 = sq1.to_value(u.m / u.s, doppler_rest=499.9 * u.nm)
        assert_allclose(val5, 60970.685737)

    def test_doppler_conversion_validation(self):

        sq1 = SpectralQuantity(1 * u.GHz)
        sq2 = SpectralQuantity(1 * u.km / u.s)

        with pytest.raises(ValueError, match='doppler_convention not set, cannot convert to/from velocities'):
            sq1.to(u.km / u.s)

        with pytest.raises(ValueError, match='doppler_convention not set, cannot convert to/from velocities'):
            sq2.to(u.GHz)

        with pytest.raises(ValueError, match='doppler_rest not set, cannot convert to/from velocities'):
            sq1.to(u.km / u.s, doppler_convention='radio')

        with pytest.raises(ValueError, match='doppler_rest not set, cannot convert to/from velocities'):
            sq2.to(u.GHz, doppler_convention='radio')

        with pytest.raises(u.UnitsError, match="Argument 'doppler_rest' to function 'to' must be in units"):
            sq1.to(u.km / u.s, doppler_convention='radio', doppler_rest=5 * u.kg)

        with pytest.raises(u.UnitsError, match="Argument 'doppler_rest' to function 'to' must be in units"):
            sq2.to(u.GHz, doppler_convention='radio', doppler_rest=5 * u.kg)

        with pytest.raises(ValueError, match="doppler_convention should be one of optical/radio/relativistic"):
            sq1.to(u.km / u.s, doppler_convention='banana', doppler_rest=5 * u.GHz)

        with pytest.raises(ValueError, match="doppler_convention should be one of optical/radio/relativistic"):
            sq2.to(u.GHz, doppler_convention='banana', doppler_rest=5 * u.GHz)

        with pytest.raises(ValueError, match='Original doppler_convention not set'):
            sq2.to(u.km / u.s, doppler_convention='radio')

        with pytest.raises(ValueError, match='Original doppler_rest not set'):
            sq2.to(u.km / u.s, doppler_rest=5 * u.GHz)

    def test_doppler_set_parameters(self):

        sq1 = SpectralQuantity(1 * u.km / u.s)

        with pytest.raises(ValueError, match="doppler_convention should be one of optical/radio/relativistic"):
            sq1.doppler_convention = 'banana'

        assert sq1.doppler_convention is None

        sq1.doppler_convention = 'radio'

        assert sq1.doppler_convention == 'radio'

        with pytest.raises(AttributeError, match="doppler_convention has already been set, and cannot be changed"):
            sq1.doppler_convention = 'optical'

        assert sq1.doppler_convention == 'radio'

        with pytest.raises(u.UnitsError, match="Argument 'value' to function 'doppler_rest' must be in units"):
            sq1.doppler_rest = 5 * u.kg

        sq1.doppler_rest = 5 * u.GHz

        assert_quantity_allclose(sq1.doppler_rest, 5 * u.GHz)

        with pytest.raises(AttributeError, match="doppler_rest has already been set, and cannot be changed"):
            sq1.doppler_rest = 4 * u.GHz

        assert_quantity_allclose(sq1.doppler_rest, 5 * u.GHz)

    def test_arithmetic(self):

        # Checks for arithmetic - some operations should return SpectralQuantity,
        # while some should just return plain Quantity

        # First, operations that should return SpectralQuantity

        sq1 = SpectralQuantity(10 * u.AA)
        sq2 = sq1 * 2
        assert isinstance(sq2, SpectralQuantity)
        assert sq2.value == 20
        assert sq2.unit == u.AA

        sq2 = sq1 / 2
        assert isinstance(sq2, SpectralQuantity)
        assert sq2.value == 5
        assert sq2.unit == u.AA

        sq3 = SpectralQuantity(10 * u.AA)
        sq3 *= 2
        assert isinstance(sq3, SpectralQuantity)
        assert sq3.value == 20
        assert sq3.unit == u.AA

        sq4 = SpectralQuantity(10 * u.AA)
        sq4 /= 2
        assert isinstance(sq4, SpectralQuantity)
        assert sq4.value == 5
        assert sq4.unit == u.AA

        sq5 = SpectralQuantity(10 * u.AA)
        with pytest.raises(TypeError, match='Cannot store the result of this operation in SpectralQuantity'):
            sq5 += 10 * u.AA

        # Note different order to sq2
        sq6 = SpectralQuantity(10 * u.AA)
        sq6 = 2 * sq1
        assert isinstance(sq6, SpectralQuantity)
        assert sq6.value == 20
        assert sq6.unit == u.AA

        # Next, operations that should return Quantity

        q1 = sq1 / u.s
        assert isinstance(q1, u.Quantity) and not isinstance(q1, SpectralQuantity)
        assert q1.value == 10
        assert q1.unit.is_equivalent(u.AA / u.s)

        q2 = sq1 / u.kg
        assert isinstance(q2, u.Quantity) and not isinstance(q2, SpectralQuantity)
        assert q2.value == 10
        assert q2.unit.is_equivalent(u.AA / u.kg)

        q3 = sq1 + 10 * u.AA
        assert isinstance(q3, u.Quantity) and not isinstance(q3, SpectralQuantity)
        assert q3.value == 20
        assert q3.unit == u.AA

        q4 = sq1 / SpectralQuantity(5 * u.AA)
        assert isinstance(q4, u.Quantity) and not isinstance(q4, SpectralQuantity)
        assert q4.value == 2
        assert q4.unit == u.one

    def test_ufuncs(self):

        # Checks for ufuncs - some operations should return SpectralQuantity,
        # while some should just return plain Quantity

        # First, operations that should return SpectralQuantity

        sq1 = SpectralQuantity([10, 20, 30] * u.AA)
        for ufunc in (np.min, np.max):
            sq2 = ufunc(sq1)
            assert isinstance(sq2, SpectralQuantity)
            assert sq2.value == ufunc(sq1.value)
            assert sq2.unit == u.AA

    def test_functions(self):

        # Checks for other functions - some operations should return SpectralQuantity,
        # while some should just return plain Quantity

        # First, operations that should return SpectralQuantity

        sq1 = SpectralQuantity([10, 20, 30] * u.AA)
        for func in (np.nanmin, np.nanmax):
            sq2 = func(sq1)
            assert isinstance(sq2, SpectralQuantity)
            assert sq2.value == func(sq1.value)
            assert sq2.unit == u.AA

        # Next, operations that should return Quantity

        for func in (np.sum,):
            q3 = func(sq1)
            assert isinstance(q3, u.Quantity) and not isinstance(q3, SpectralQuantity)
            assert q3.value == func(sq1.value)
            assert q3.unit == u.AA

    @pytest.mark.xfail
    def test_functions_std(self):

        # np.std should return a Quantity but it returns a SpectralQuantity. We
        # make this a separate xfailed test for now, but once this passes,
        # np.std could also just be added to the main test_functions test.
        # See https://github.com/astropy/astropy/issues/10245 for more details.

        # Checks for other functions - some operations should return SpectralQuantity,
        # while some should just return plain Quantity

        # First, operations that should return SpectralQuantity

        sq1 = SpectralQuantity([10, 20, 30] * u.AA)
        q1 = np.std(sq1)
        assert isinstance(q1, u.Quantity) and not isinstance(q1, SpectralQuantity)
        assert q1.value == np.sum(sq1.value)
        assert q1.unit == u.AA
