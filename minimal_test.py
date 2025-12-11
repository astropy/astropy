"""
Minimal reproducer for Astropy issue #19045.

This demonstrates the simplest possible case where to_system() fails
with a custom unit system that has non-canonical (composite) base units.

The core issue:
- to_system() calls decompose(bases=system.bases)
- decompose() expects bases to be SI-like irreducible units (m, s, kg, etc.)
- When bases contain composite units like velocity (m/s) or G (m^3/(kg*s^2)),
  decompose() cannot express other units in terms of them

The error occurs because decompose() works by breaking down units into
irreducible SI units first, then checking if each component is in the bases.
Composite bases don't match this logic.
"""

import astropy.units as u

print("=" * 70)
print("Minimal reproducer for Astropy issue #19045")
print("=" * 70)
print()

# =============================================================================
# CASE 1: Simplest failing case - velocity as a base unit
# =============================================================================
print("CASE 1: Velocity (m/s) as a base unit")
print("-" * 70)

class VelocitySystem:
    """System with velocity as a base instead of separate length and time."""
    m = u.m
    velocity = u.m / u.s  # Composite base unit
    bases = {u.m, u.m / u.s}  # length and velocity

# Try to convert acceleration (m/s^2) to this system
# Expected: acceleration = velocity^2 / m  (since [L/T^2] = [L/T]^2 / [L])
accel = u.m / u.s**2

print(f"Source unit: {accel} (acceleration)")
print(f"System bases: {VelocitySystem.bases}")
print()

# This should work: accel = velocity^2 / m
target = (u.m / u.s)**2 / u.m  # = m/s^2
print(f"Mathematical target: velocity^2 / length = {target}")
print(f"Equivalence check: {accel.is_equivalent(target)}")
print()

print("Test with to():")
try:
    result = accel.to(target)
    print(f"  SUCCESS: {accel}.to({target}) = {result}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

print("\nTest with to_system():")
try:
    result = accel.to_system(VelocitySystem)
    print(f"  SUCCESS: {result}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

print()

# =============================================================================
# CASE 2: Another simple case - using Hz (1/s) as a base
# =============================================================================
print("CASE 2: Hz (1/s) as a base unit")  
print("-" * 70)

class FrequencySystem:
    """System using Hz instead of s^-1."""
    m = u.m
    Hz = u.Hz  # This is 1/s, a composite unit
    bases = {u.m, u.Hz}

# Try to convert velocity (m/s) - should be m * Hz
velocity = u.m / u.s
print(f"Source unit: {velocity} (velocity)")
print(f"System bases: {FrequencySystem.bases}")
print()

target = u.m * u.Hz
print(f"Mathematical target: length * Hz = {target}")
print(f"Target decomposed: {target.decompose()}")
print()

print("Test with to():")
try:
    result = velocity.to(target)
    print(f"  SUCCESS: {velocity}.to({target}) = {result}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

print("\nTest with to_system():")
try:
    result = velocity.to_system(FrequencySystem)
    print(f"  SUCCESS: {result}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

print()

# =============================================================================
# CASE 3: Standard SI system (should work)
# =============================================================================
print("CASE 3: Standard SI system (control test)")
print("-" * 70)

print(f"Source unit: {velocity}")
print("Test with to_system(u.si):")
try:
    result = velocity.to_system(u.si)
    print(f"  SUCCESS: {result}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

print()
print("=" * 70)
print("Analysis")
print("=" * 70)
print("""
The root cause is in to_system():

    def to_system(self, system):
        return sorted(
            self.decompose(bases=system.bases).compose(units=system),  # <-- HERE
            key=lambda x: len(set(x.bases).difference(system.bases)),
        )

When decompose(bases=...) is called with composite units in bases,
it fails because decompose() expects bases to be irreducible SI units.

The decompose method works by breaking units down to their fundamental
SI components (m, s, kg, etc.) and checking if those are in the bases set.
Composite units like m/s or m^3/(kg*s^2) don't exist in that decomposed form.

For example, to express kg in terms of {m, m/s, m^3/(kg*s^2)}:
- We need to use G (m^3/(kg*s^2)) which CONTAINS kg
- So kg = (velocity^2 * length) / G = (m/s)^2 * m / (m^3/(kg*s^2))
- But decompose() doesn't know how to do this algebraic manipulation

A fix would need to:
1. Solve a system of linear equations in the exponents
2. Find the powers of each base unit that produce the target unit
""")

