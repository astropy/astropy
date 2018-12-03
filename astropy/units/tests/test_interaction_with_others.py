# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test that astropy's Unit class can deal with units from other projects."""

import pytest

from astropy.units import Unit, UnitBase, dimensionless_unscaled


pytestmark = pytest.mark.filterwarnings("ignore")


class OtherUnitModuleTests:
    # Subclasses must define a setup_class that sets
    # - oum for the module to the units.
    # - dimensionless to its dimensionless unit.

    @pytest.mark.parametrize('unit', ['m', 's', 'kg', 'g', 'cm'])
    def test_simple(self, unit):
        u_unit = getattr(self.oum, unit)
        assert str(u_unit == unit)
        u_m = Unit(u_unit)
        assert u_m == Unit(unit)

    def test_dimensionless(self):
        u_dimensionless = Unit(self.dimensionless)
        assert isinstance(u_dimensionless, UnitBase)
        assert u_dimensionless == dimensionless_unscaled

    def test_more_complicated(self):
        u_unit = self.oum.m / self.oum.s
        u_m_per_s = Unit(u_unit)
        assert isinstance(u_m_per_s, UnitBase)
        assert u_m_per_s == Unit('m/s')


class TestUnyt(OtherUnitModuleTests):
    def setup_class(cls):
        try:
            import unyt
        except ImportError:
            pytest.skip('unyt not available.')

        cls.oum = unyt
        cls.dimensionless = unyt.dimensionless


@pytest.mark.filterwarnings('ignore')
class TestAmuse(OtherUnitModuleTests):
    def setup_class(cls):
        try:
            from amuse.units import units as amuse_units
        except ImportError:
            pytest.skip('amuse not available.')

        cls.oum = amuse_units
        cls.dimensionless = amuse_units.none

    def teardown_class(cls):
        # Annoyingly, amuse thinks it is OK to just override
        # numpy operators by monkeypatching, which can obviously
        # screw up further tests.  So, we need to get rid of it.
        from amuse.units import quantities
        import atexit
        if quantities._previous_operators:
            quantities.unset_numpy_operators()
            atexit.unregister(quantities.unset_numpy_operators)
