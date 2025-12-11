"""
Tests for to_system() with custom unit systems.

This test module covers issue #19045 where to_system() fails for units
in custom unit systems whose bases are non-canonical (composite) units
like velocity (m/s) or the gravitational constant G (m^3/(kg*s^2)).

These tests can be added to astropy/units/tests/test_units.py
"""

import numpy as np
import pytest
import astropy.units as u
from astropy.units.core import CompositeUnit, UnitConversionError


# =============================================================================
# Helper function (the fix) - should be added to core.py
# =============================================================================

def _decompose_to_arbitrary_bases(unit, bases):
    """
    Decompose a unit into arbitrary base units using linear algebra.
    
    This is a more general version of decompose() that can handle composite
    base units (like velocity or the gravitational constant G).
    """
    SI_DIMS = [u.m, u.kg, u.s, u.A, u.K, u.mol, u.cd, u.rad]
    n_dims = len(SI_DIMS)
    
    def get_dimension_vector(un):
        decomposed = un.decompose()
        vec = np.zeros(n_dims)
        if not decomposed.bases:
            return vec
        for base, power in zip(decomposed.bases, decomposed.powers):
            for i, si_dim in enumerate(SI_DIMS):
                if base == si_dim or base.is_equivalent(si_dim):
                    vec[i] = power
                    break
        return vec
    
    bases = list(bases)
    n_bases = len(bases)
    
    A = np.zeros((n_dims, n_bases))
    base_scales = []
    
    for j, base in enumerate(bases):
        A[:, j] = get_dimension_vector(base)
        base_scales.append(base.decompose().scale)
    
    d = get_dimension_vector(unit)
    target_scale = unit.decompose().scale
    
    try:
        p, residuals, rank, s = np.linalg.lstsq(A, d, rcond=None)
        if not np.allclose(A @ p, d, rtol=1e-10, atol=1e-10):
            raise UnitConversionError(
                f"'{unit}' cannot be expressed in terms of the given bases"
            )
    except np.linalg.LinAlgError:
        raise UnitConversionError(
            f"Cannot solve for decomposition of '{unit}' into given bases"
        )
    
    result_scale = target_scale
    for j, power in enumerate(p):
        result_scale /= base_scales[j] ** power
    
    result_bases = []
    result_powers = []
    
    for j, power in enumerate(p):
        if not np.isclose(power, 0, atol=1e-10):
            result_bases.append(bases[j])
            if np.isclose(power, round(power), rtol=1e-10):
                result_powers.append(int(round(power)))
            else:
                result_powers.append(float(power))
    
    return CompositeUnit(
        result_scale if result_bases else target_scale,
        result_bases,
        result_powers,
        _error_check=False
    )


def to_system_fixed(unit, system):
    """Fixed version of to_system that works with composite base units."""
    try:
        decomposed = unit.decompose(bases=system.bases)
    except UnitConversionError:
        decomposed = _decompose_to_arbitrary_bases(unit, system.bases)
    
    try:
        composed = decomposed.compose(units=system)
        if composed:
            return sorted(
                composed,
                key=lambda x: len(set(x.bases).difference(system.bases)),
            )
    except Exception:
        pass
    
    return [decomposed]


# =============================================================================
# Test Classes (Custom Unit Systems)
# =============================================================================

class VelocitySystem:
    """Unit system with velocity as a base unit instead of time."""
    m = u.m
    velocity = u.m / u.s
    bases = {u.m, u.m / u.s}


class FrequencySystem:
    """Unit system using Hz (frequency) instead of time."""
    m = u.m
    Hz = u.Hz
    bases = {u.m, u.Hz}


class GravitationalSystem:
    """
    Unit system for gravitational simulations.
    
    Bases: [length, velocity, G]
    This is the system that originally failed in issue #19045.
    """
    from astropy.constants import G as G_const
    kpc = u.kpc
    km_s = u.km / u.s
    G_unit = G_const.unit
    bases = {u.kpc, u.km / u.s, G_const.unit}


class MasslessSystem:
    """System without mass dimension - used to test failures."""
    m = u.m
    s = u.s
    bases = {u.m, u.s}


# =============================================================================
# Tests
# =============================================================================

class TestToSystemWithCompositeBases:
    """Tests for to_system() with composite base units."""
    
    def test_velocity_base_acceleration(self):
        """
        Test expressing acceleration in a system with velocity as a base.
        
        acceleration [L/T^2] = velocity^2 / length = [L/T]^2 / [L]
        """
        accel = u.m / u.s**2
        result = to_system_fixed(accel, VelocitySystem)
        
        assert len(result) >= 1
        assert accel.is_equivalent(result[0])
        assert np.isclose(accel.to(result[0]), 1.0)
    
    def test_velocity_base_energy(self):
        """
        Test expressing energy in a system with velocity as a base.
        
        energy [M*L^2/T^2] = mass * velocity^2
        """
        energy = u.kg * u.m**2 / u.s**2
        
        # This system doesn't have mass, so should fail
        with pytest.raises(UnitConversionError):
            to_system_fixed(energy, VelocitySystem)
    
    def test_frequency_base_velocity(self):
        """
        Test expressing velocity in a system with Hz as a base.
        
        velocity [L/T] = length * Hz = [L] * [1/T]
        """
        velocity = u.m / u.s
        result = to_system_fixed(velocity, FrequencySystem)
        
        assert len(result) >= 1
        assert velocity.is_equivalent(result[0])
        assert np.isclose(velocity.to(result[0]), 1.0)
    
    def test_frequency_base_acceleration(self):
        """
        Test expressing acceleration in a system with Hz as a base.
        
        acceleration [L/T^2] = length * Hz^2
        """
        accel = u.m / u.s**2
        result = to_system_fixed(accel, FrequencySystem)
        
        assert len(result) >= 1
        assert accel.is_equivalent(result[0])
    
    def test_gravitational_system_density(self):
        """
        Test expressing density in the gravitational system (original issue).
        
        This is the exact case from issue #19045.
        density [M/L^3] = velocity^2 / (length^2 * G)
        
        Since G has dimensions [L^3 / (M*T^2)]:
        velocity^2 / (length^2 * G) = [L/T]^2 / ([L]^2 * [L^3/(M*T^2)])
                                    = [L^2/T^2] * [M*T^2/L^5]
                                    = [M/L^3] ✓
        """
        density = u.Msun / u.kpc**3
        result = to_system_fixed(density, GravitationalSystem)
        
        assert len(result) >= 1
        assert density.is_equivalent(result[0])
        # Allow for small numerical error
        assert np.isclose(density.to(result[0]), 1.0, rtol=1e-5)
    
    def test_gravitational_system_mass(self):
        """
        Test expressing mass in the gravitational system.
        
        mass = velocity^2 * length / G
        """
        mass = u.Msun
        result = to_system_fixed(mass, GravitationalSystem)
        
        assert len(result) >= 1
        assert mass.is_equivalent(result[0])
    
    def test_gravitational_system_time(self):
        """
        Test expressing time in the gravitational system.
        
        time = length / velocity
        """
        time = u.Gyr
        result = to_system_fixed(time, GravitationalSystem)
        
        assert len(result) >= 1
        assert time.is_equivalent(result[0])
    
    def test_impossible_decomposition(self):
        """Test that impossible decompositions fail with clear error."""
        # Can't express mass in a system without mass dimension
        mass = u.kg
        with pytest.raises(UnitConversionError):
            to_system_fixed(mass, MasslessSystem)
    
    def test_dimensionless_unit(self):
        """Test that dimensionless units work correctly."""
        dimless = u.dimensionless_unscaled
        result = to_system_fixed(dimless, VelocitySystem)
        
        assert len(result) >= 1
        assert dimless.is_equivalent(result[0])


class TestToSystemRegressions:
    """Regression tests to ensure standard behavior is preserved."""
    
    def test_si_system(self):
        """Test that standard SI conversion still works."""
        pressure = u.Pa
        result = to_system_fixed(pressure, u.si)
        
        assert len(result) >= 1
        assert u.Pa in result or any(r.is_equivalent(u.Pa) for r in result)
    
    def test_cgs_system(self):
        """Test that standard CGS conversion still works."""
        pressure = u.Pa
        result = to_system_fixed(pressure, u.cgs)
        
        assert len(result) >= 1
        # Pa = 10 Ba in CGS
        assert any(r.is_equivalent(u.Ba) for r in result)
    
    def test_compound_unit_si(self):
        """Test compound unit conversion to SI."""
        area_per_time = u.m**2 / u.s
        result = to_system_fixed(area_per_time, u.si)
        
        assert len(result) >= 1
        assert area_per_time.is_equivalent(result[0])


# =============================================================================
# Run tests if executed directly
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Running tests for to_system() with composite bases")
    print("=" * 70)
    
    test_class = TestToSystemWithCompositeBases()
    
    tests = [
        ("velocity_base_acceleration", test_class.test_velocity_base_acceleration),
        ("frequency_base_velocity", test_class.test_frequency_base_velocity),
        ("frequency_base_acceleration", test_class.test_frequency_base_acceleration),
        ("gravitational_system_density", test_class.test_gravitational_system_density),
        ("gravitational_system_mass", test_class.test_gravitational_system_mass),
        ("gravitational_system_time", test_class.test_gravitational_system_time),
        ("dimensionless_unit", test_class.test_dimensionless_unit),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            test_func()
            print(f"✓ {name}")
            passed += 1
        except Exception as e:
            print(f"✗ {name}: {e}")
            failed += 1
    
    # Test expected failures (without pytest)
    print("\nExpected failure tests:")
    
    # Test velocity_base_energy - should fail
    try:
        energy = u.kg * u.m**2 / u.s**2
        to_system_fixed(energy, VelocitySystem)
        print("✗ velocity_base_energy: Should have raised UnitConversionError")
        failed += 1
    except UnitConversionError:
        print("✓ velocity_base_energy (correctly raised error)")
        passed += 1
    except Exception as e:
        print(f"✗ velocity_base_energy: Wrong exception: {type(e).__name__}: {e}")
        failed += 1
    
    # Test impossible_decomposition - should fail
    try:
        mass = u.kg
        to_system_fixed(mass, MasslessSystem)
        print("✗ impossible_decomposition: Should have raised UnitConversionError")
        failed += 1
    except UnitConversionError:
        print("✓ impossible_decomposition (correctly raised error)")
        passed += 1
    except Exception as e:
        print(f"✗ impossible_decomposition: Wrong exception: {type(e).__name__}: {e}")
        failed += 1
    
    # Regression tests
    print("\nRegression tests:")
    
    regression_class = TestToSystemRegressions()
    regression_tests = [
        ("si_system", regression_class.test_si_system),
        ("cgs_system", regression_class.test_cgs_system),
        ("compound_unit_si", regression_class.test_compound_unit_si),
    ]
    
    for name, test_func in regression_tests:
        try:
            test_func()
            print(f"✓ {name}")
            passed += 1
        except Exception as e:
            print(f"✗ {name}: {e}")
            failed += 1
    
    print("\n" + "=" * 70)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 70)

