"""
Proposed Patch for Astropy Issue #19045
========================================

This file contains the proposed changes to astropy/units/core.py to fix
the to_system() method for custom unit systems with non-canonical bases.

The key insight is that decompose(bases=...) currently only works with
SI-like irreducible base units. For composite bases like velocity (m/s)
or G (m^3/(kg*s^2)), we need to solve a linear system to find the powers.

PROPOSED CHANGES:
-----------------

1. Add a new helper function `_decompose_to_arbitrary_bases()` that uses
   linear algebra to decompose a unit into arbitrary (possibly composite)
   base units.

2. Modify `to_system()` to use this new function when the standard
   decompose() fails.

HOW TO APPLY THIS PATCH:
------------------------

The changes below should be added to astropy/units/core.py.
Search for the to_system method and add the helper function above it,
then modify to_system to use the new logic.

"""

# =============================================================================
# NEW CODE TO ADD TO astropy/units/core.py
# =============================================================================

# Add these imports at the top of core.py:
# import numpy as np  # (already present in some form)

# Add this helper function before the to_system method in UnitBase class:

PATCH_CODE = '''
def _decompose_to_arbitrary_bases(
    unit: "UnitBase",
    bases: Collection["UnitBase"],
) -> "CompositeUnit":
    """
    Decompose a unit into arbitrary base units using linear algebra.
    
    This is a more general version of decompose() that can handle composite
    base units (like velocity or the gravitational constant G).
    
    Parameters
    ----------
    unit : UnitBase
        The unit to decompose.
    bases : Collection of UnitBase
        The base units to decompose into. These can be composite units.
        
    Returns
    -------
    CompositeUnit
        The unit expressed in terms of the given bases.
        
    Raises
    ------
    UnitConversionError
        If the unit cannot be expressed in terms of the given bases
        (i.e., the dimensional space spanned by the bases doesn't include
        the dimensions of the target unit).
        
    Notes
    -----
    This works by:
    1. Decomposing both the target unit and all base units into SI dimensions
    2. Building a matrix A where column j is the dimension vector of base j
    3. Solving A @ p = d for the powers p, where d is the target dimensions
    4. Constructing the result as bases[0]^p[0] * bases[1]^p[1] * ...
    
    The SI dimensions used are: length, mass, time, current, temperature,
    amount, luminous intensity, and angle.
    """
    import numpy as np
    from . import si
    
    # Define the fundamental SI dimensions
    SI_DIMS = [si.m, si.kg, si.s, si.A, si.K, si.mol, si.cd, si.rad]
    n_dims = len(SI_DIMS)
    
    def get_dimension_vector(u):
        """Get SI dimension powers for a unit."""
        decomposed = u.decompose()
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
    
    # Build dimension matrix A (n_dims × n_bases)
    A = np.zeros((n_dims, n_bases))
    base_scales = []
    
    for j, base in enumerate(bases):
        A[:, j] = get_dimension_vector(base)
        base_scales.append(base.decompose().scale)
    
    # Get target dimensions
    d = get_dimension_vector(unit)
    target_scale = unit.decompose().scale
    
    # Solve A @ p = d using least squares
    try:
        p, residuals, rank, s = np.linalg.lstsq(A, d, rcond=None)
        
        # Verify solution is exact
        if not np.allclose(A @ p, d, rtol=1e-10, atol=1e-10):
            raise UnitConversionError(
                f"'{unit}' cannot be expressed in terms of the given bases"
            )
    except np.linalg.LinAlgError:
        raise UnitConversionError(
            f"Cannot solve for decomposition of '{unit}' into given bases"
        )
    
    # Compute scale factor
    result_scale = target_scale
    for j, power in enumerate(p):
        result_scale /= base_scales[j] ** power
    
    # Build result, filtering zero powers and rounding near-integers
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


# MODIFIED to_system method:

def to_system(self, system):
    """Convert this unit into ones belonging to the given system.

    Since more than one result may be possible, a list is always
    returned.

    Parameters
    ----------
    system : module
        The module that defines the unit system.  Commonly used
        ones include `astropy.units.si` and `astropy.units.cgs`.

        To use your own module it must contain unit objects and a
        sequence member named ``bases`` containing the base units of
        the system. The bases can be composite units (like velocity
        or the gravitational constant).

    Returns
    -------
    units : list of `CompositeUnit`
        With an attempt to sort simpler units to the start (see examples).

    Examples
    --------
    >>> import astropy.units as u
    >>> (u.N / u.m**2).to_system(u.si)  # preference for simpler units
    [Unit("Pa"), Unit("N / m2"), Unit("J / m3")]
    >>> u.Pa.to_system(u.cgs)
    [Unit("10 Ba"), Unit("10 P / s")]
    >>> u.Ba.to_system(u.si)
    [Unit("0.1 Pa"), Unit("0.1 N / m2"), Unit("0.1 J / m3")]
    >>> (u.AU/u.yr).to_system(u.cgs)  # preference for base units
    [Unit("474047 cm / s"), Unit("474047 Gal s")]
    >>> (u.m / u.s**2).to_system(u.cgs)
    [Unit("100 cm / s2"), Unit("100 Gal")]
    
    Custom unit systems with composite bases are also supported:
    
    >>> class VelocitySystem:
    ...     bases = {u.m, u.m / u.s}  # length and velocity as bases
    >>> (u.m / u.s**2).to_system(VelocitySystem)  # acceleration
    [Unit("m / s2")]

    """
    # Try standard decompose first (works for SI-like irreducible bases)
    try:
        decomposed = self.decompose(bases=system.bases)
    except UnitConversionError:
        # Fall back to linear algebra approach for composite bases
        decomposed = _decompose_to_arbitrary_bases(self, system.bases)
    
    # Try to compose into named units from the system
    try:
        composed = decomposed.compose(units=system)
        if composed:
            return sorted(
                composed,
                key=lambda x: len(set(x.bases).difference(system.bases)),
            )
    except UnitsError:
        pass
    
    # If compose fails, return the decomposed unit directly
    # This happens for custom systems without many named units
    return [decomposed]
'''

print(PATCH_CODE)

# =============================================================================
# TEST THE PATCHED BEHAVIOR
# =============================================================================

if __name__ == "__main__":
    import numpy as np
    import astropy.units as u
    from astropy.units.core import CompositeUnit, UnitConversionError
    from astropy.units import si
    
    # Implement the patch locally for testing
    def _decompose_to_arbitrary_bases(unit, bases):
        """Local implementation of the patch for testing."""
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
    
    def to_system_patched(unit, system):
        """Patched version of to_system for testing."""
        try:
            decomposed = unit.decompose(bases=system.bases)
        except UnitConversionError:
            decomposed = _decompose_to_arbitrary_bases(unit, system.bases)
        
        # Try to compose into named units from the system
        try:
            composed = decomposed.compose(units=system)
            if composed:
                return sorted(
                    composed,
                    key=lambda x: len(set(x.bases).difference(system.bases)),
                )
        except Exception:
            pass
        
        # If compose fails, return the decomposed unit directly
        # This happens for custom systems without many named units
        return [decomposed]
    
    print("\n" + "=" * 70)
    print("Testing Patched to_system() Behavior")
    print("=" * 70)
    
    # Test 1: Original failing case from issue
    print("\nTest 1: Original Issue - Density with [kpc, km/s, G] bases")
    print("-" * 50)
    
    from astropy.constants import G as G_const
    
    class CustomSystem:
        kpc = u.kpc
        km_s = u.km / u.s
        G_unit = G_const.unit
        bases = {u.kpc, u.km / u.s, G_const.unit}
    
    density = u.Msun / u.kpc**3
    print(f"Input: {density}")
    
    print("\nOld behavior (fails):")
    try:
        result = density.to_system(CustomSystem)
        print(f"  Result: {result}")
    except Exception as e:
        print(f"  Error: {e}")
    
    print("\nPatched behavior (works):")
    try:
        result = to_system_patched(density, CustomSystem)
        print(f"  Result: {result}")
        if result:
            print(f"  Equivalence check: {density.is_equivalent(result[0])}")
    except Exception as e:
        print(f"  Error: {e}")
    
    # Test 2: Standard SI still works
    print("\nTest 2: Standard SI System (regression check)")
    print("-" * 50)
    
    pressure = u.Pa
    print(f"Input: {pressure}")
    
    print("\nPatched behavior:")
    result = to_system_patched(pressure, u.si)
    print(f"  Result: {result[:3]}...")
    
    # Test 3: CGS still works
    print("\nTest 3: Standard CGS System (regression check)")
    print("-" * 50)
    
    print(f"Input: {pressure}")
    print("\nPatched behavior:")
    result = to_system_patched(pressure, u.cgs)
    print(f"  Result: {result[:3]}...")
    
    print("\n" + "=" * 70)
    print("Patch testing complete!")
    print("=" * 70)

