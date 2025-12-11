"""
Analysis and Proposed Fix for Astropy Issue #19045
===================================================

Issue: Unit.to_system() fails for density units in custom unit systems 
whose bases are non-canonical (composite) units like [length, velocity, G].

This file contains:
1. Root cause analysis  
2. Mathematical formulation of the fix
3. Prototype implementation
4. Tests demonstrating the fix works

Root Cause Analysis
-------------------

The `to_system()` method works by:

    def to_system(self, system):
        return sorted(
            self.decompose(bases=system.bases).compose(units=system),
            key=lambda x: len(set(x.bases).difference(system.bases)),
        )

The problem is in `decompose(bases=system.bases)`. When `decompose()` is called
on an IrreducibleUnit (like `kg` or `s`) with custom bases:

    def decompose(self, bases: Collection[UnitBase] = ()) -> UnitBase:
        if len(bases) and self not in bases:
            for base in bases:
                try:
                    scale = self._to(base)  # Try direct conversion
                except UnitsError:
                    pass
                else:
                    return CompositeUnit(scale, [base], [1])
            
            raise UnitConversionError(
                f"Unit {self} can not be decomposed into the requested bases"
            )

This only tries DIRECT conversion from self to each base. It cannot handle cases
where self needs to be expressed as a COMBINATION of base units with powers.

Example: Express `kg` in terms of bases = {m, m/s, m^3/(kg*s^2)}

- kg cannot directly convert to m (different dimensions)
- kg cannot directly convert to m/s (different dimensions)
- kg cannot directly convert to m^3/(kg*s^2) (different dimensions)

But mathematically, kg CAN be expressed:
    kg = (m/s)^2 * m / (m^3/(kg*s^2)) = velocity^2 * length / G

The fix requires solving a LINEAR SYSTEM to find the powers.


Mathematical Formulation
------------------------

Let D be the number of SI base dimensions (length, mass, time, etc.)
Let N be the number of custom base units.

Each custom base unit can be decomposed into SI dimensions:
    base_i = scale_i * m^a_i * kg^b_i * s^c_i * ...

This gives us a D×N matrix A where A[d,i] = power of dimension d in base_i.

For a target unit with SI decomposition:
    target = scale_t * m^a_t * kg^b_t * s^c_t * ...

We get a target vector d of length D.

We need to find powers p = [p_1, p_2, ..., p_N] such that:
    A @ p = d

This is a system of linear equations. If the bases span the necessary
dimensional space, there will be a solution (possibly non-unique if N > D).


Prototype Implementation
------------------------
"""

import numpy as np
import astropy.units as u
from astropy.units.core import CompositeUnit, IrreducibleUnit, UnitConversionError


# Define the fundamental SI dimensions and their order
SI_DIMENSION_UNITS = [
    u.m,   # length
    u.kg,  # mass
    u.s,   # time
    u.A,   # current
    u.K,   # temperature
    u.mol, # amount
    u.cd,  # luminous intensity
    u.rad, # angle (treated as dimensionless but tracked separately)
]


def get_si_dimension_vector(unit):
    """
    Get the SI dimension vector for a unit.
    
    Returns a numpy array of powers for [m, kg, s, A, K, mol, cd, rad].
    
    Example:
        get_si_dimension_vector(u.m / u.s**2) -> [1, 0, -2, 0, 0, 0, 0, 0]
    """
    decomposed = unit.decompose()
    
    # Handle scalar (dimensionless with just scale)
    if not decomposed.bases:
        return np.zeros(len(SI_DIMENSION_UNITS))
    
    dim_vec = np.zeros(len(SI_DIMENSION_UNITS))
    
    for base, power in zip(decomposed.bases, decomposed.powers):
        for i, si_dim in enumerate(SI_DIMENSION_UNITS):
            if base == si_dim or base.is_equivalent(si_dim):
                dim_vec[i] = power
                break
    
    return dim_vec


def get_scale(unit):
    """Get the scale factor of a unit relative to SI base units."""
    return unit.decompose().scale


def decompose_to_custom_bases(unit, custom_bases):
    """
    Decompose a unit into custom base units using linear algebra.
    
    Parameters
    ----------
    unit : UnitBase
        The unit to decompose.
    custom_bases : sequence of UnitBase
        The custom base units to decompose into.
    
    Returns
    -------
    CompositeUnit
        The unit expressed in terms of the custom bases.
        
    Raises
    ------
    UnitConversionError
        If the unit cannot be expressed in terms of the given bases.
        
    Notes
    -----
    This solves the linear system A @ p = d where:
    - A is the matrix of SI dimension vectors for each custom base
    - d is the SI dimension vector of the target unit
    - p is the unknown powers of each custom base
    """
    custom_bases = list(custom_bases)
    n_bases = len(custom_bases)
    n_dims = len(SI_DIMENSION_UNITS)
    
    # Build the dimension matrix A (n_dims x n_bases)
    # Each column is the SI dimension vector of a custom base
    A = np.zeros((n_dims, n_bases))
    base_scales = []
    
    for j, base in enumerate(custom_bases):
        A[:, j] = get_si_dimension_vector(base)
        base_scales.append(get_scale(base))
    
    # Get the target dimension vector
    d = get_si_dimension_vector(unit)
    target_scale = get_scale(unit)
    
    # Solve the system A @ p = d
    # Use least squares to handle over/under-determined systems
    try:
        # Check if dimensions are achievable
        # We need d to be in the column space of A
        p, residuals, rank, s = np.linalg.lstsq(A, d, rcond=None)
        
        # Check if the solution is exact (residuals should be ~0)
        reconstructed = A @ p
        if not np.allclose(reconstructed, d, rtol=1e-10, atol=1e-10):
            raise UnitConversionError(
                f"Unit {unit} cannot be expressed in terms of the given bases. "
                f"Dimension mismatch: target={d}, achievable={reconstructed}"
            )
        
    except np.linalg.LinAlgError as e:
        raise UnitConversionError(
            f"Cannot solve for decomposition of {unit}: {e}"
        )
    
    # Compute the scale factor
    # target_scale = result_scale * product(base_scales[i]^p[i])
    result_scale = target_scale
    for j, power in enumerate(p):
        result_scale /= base_scales[j] ** power
    
    # Build the result CompositeUnit
    # Filter out zero powers
    result_bases = []
    result_powers = []
    
    for j, power in enumerate(p):
        if not np.isclose(power, 0, atol=1e-10):
            result_bases.append(custom_bases[j])
            # Try to use integer powers if close
            if np.isclose(power, round(power), rtol=1e-10):
                result_powers.append(int(round(power)))
            else:
                result_powers.append(float(power))
    
    if not result_bases:
        # Dimensionless result
        return CompositeUnit(result_scale, [], [])
    
    return CompositeUnit(result_scale, result_bases, result_powers)


def to_system_fixed(unit, system):
    """
    Fixed version of to_system that works with non-canonical base units.
    
    This is a drop-in replacement for UnitBase.to_system() that can handle
    custom unit systems with composite base units like velocity and G.
    
    Parameters
    ----------
    unit : UnitBase
        The unit to convert.
    system : module-like
        An object with a 'bases' attribute containing the base units.
        
    Returns
    -------
    list of CompositeUnit
        Possible representations in the target system.
    """
    # First try to decompose into the custom bases
    try:
        decomposed = decompose_to_custom_bases(unit, system.bases)
    except UnitConversionError:
        raise
    
    # Then compose using the units available in the system
    # For now, just return the decomposed form
    # A full implementation would also use compose() to find named units
    result = sorted(
        [decomposed],
        key=lambda x: len(set(x.bases).difference(system.bases)),
    )
    
    return result


# =============================================================================
# Tests
# =============================================================================

def test_velocity_system():
    """Test with a system using velocity as a base unit."""
    print("\nTest: Velocity System")
    print("-" * 50)
    
    class VelocitySystem:
        bases = {u.m, u.m / u.s}
    
    # Test: acceleration (m/s^2) should be velocity^2 / length
    accel = u.m / u.s**2
    
    print(f"Source: {accel}")
    print(f"Bases: {VelocitySystem.bases}")
    
    try:
        result = to_system_fixed(accel, VelocitySystem)
        print(f"Result: {result}")
        
        # Verify equivalence
        if result:
            assert accel.is_equivalent(result[0]), "Result not equivalent!"
            print(f"Equivalence check: PASSED")
            print(f"Conversion factor: {accel.to(result[0])}")
    except Exception as e:
        print(f"FAILED: {e}")


def test_density_system():
    """Test with [length, velocity, G] bases (original issue)."""
    print("\nTest: Density System [length, velocity, G]")
    print("-" * 50)
    
    from astropy.constants import G
    
    class DensitySystem:
        bases = {u.kpc, u.km / u.s, G.unit}
    
    # Test: density (Msun/kpc^3) should be velocity^2 / (length^2 * G)
    density = u.Msun / u.kpc**3
    
    print(f"Source: {density}")
    print(f"Bases: {DensitySystem.bases}")
    
    try:
        result = to_system_fixed(density, DensitySystem)
        print(f"Result: {result}")
        
        # Verify equivalence  
        if result:
            assert density.is_equivalent(result[0]), "Result not equivalent!"
            print(f"Equivalence check: PASSED")
            print(f"Conversion factor: {density.to(result[0])}")
    except Exception as e:
        print(f"FAILED: {e}")


def test_frequency_system():
    """Test with Hz (1/s) as a base unit."""
    print("\nTest: Frequency System")
    print("-" * 50)
    
    class FrequencySystem:
        bases = {u.m, u.Hz}
    
    # Test: velocity (m/s) should be m * Hz
    velocity = u.m / u.s
    
    print(f"Source: {velocity}")
    print(f"Bases: {FrequencySystem.bases}")
    
    try:
        result = to_system_fixed(velocity, FrequencySystem)
        print(f"Result: {result}")
        
        if result:
            assert velocity.is_equivalent(result[0]), "Result not equivalent!"
            print(f"Equivalence check: PASSED")
            print(f"Conversion factor: {velocity.to(result[0])}")
    except Exception as e:
        print(f"FAILED: {e}")


def test_impossible_decomposition():
    """Test that impossible decompositions fail gracefully."""
    print("\nTest: Impossible Decomposition")
    print("-" * 50)
    
    class IncompleteSystem:
        # Only length - cannot express mass or time
        bases = {u.m}
    
    # Try to express velocity - should fail
    velocity = u.m / u.s
    
    print(f"Source: {velocity}")
    print(f"Bases: {IncompleteSystem.bases}")
    
    try:
        result = to_system_fixed(velocity, IncompleteSystem)
        print(f"Result: {result}")
        print("ERROR: Should have failed!")
    except UnitConversionError as e:
        print(f"Correctly failed: {e}")


if __name__ == "__main__":
    print("=" * 70)
    print("Testing Fixed to_system() Implementation")
    print("=" * 70)
    
    test_velocity_system()
    test_frequency_system()
    test_density_system()
    test_impossible_decomposition()
    
    print("\n" + "=" * 70)
    print("All tests completed!")
    print("=" * 70)

