"""
Test script to reproduce Astropy issue #19045.

This demonstrates that Unit.to_system() fails for a density unit in a custom
unit system whose bases are [length, velocity, G], even though the same
conversion works with to().

Expected behavior:
- density = velocity^2 / (length^2 * G)
- This is dimensionally correct: [M/L^3] = [L/T]^2 / ([L]^2 * [L^3/(M*T^2)])
                                         = [L^2/T^2] / [L^5/(M*T^2)]
                                         = [M/L^3] ✓
"""

import astropy.units as u
from astropy.constants import G

print("=" * 70)
print("Testing Astropy issue #19045: to_system() fails with custom unit system")
print("=" * 70)
print()

# Define the density unit we want to convert
density_unit = u.Msun / u.kpc**3
print(f"Source unit: {density_unit}")
print(f"Physical type: {density_unit.physical_type}")
print()

# Define target units in the custom system
length = u.kpc
velocity = u.km / u.s
G_unit = G.unit  # m^3 / (kg * s^2)

print("Custom unit system base units:")
print(f"  length: {length}")
print(f"  velocity: {velocity}")
print(f"  G: {G_unit}")
print()

# The expected result in the custom system:
# density = velocity^2 / (length^2 * G)
target_unit = velocity**2 / (length**2 * G_unit)
print(f"Expected target unit: {target_unit}")
print()

# Test 1: Using to() - this should work
print("-" * 70)
print("Test 1: Using to() method")
print("-" * 70)
try:
    result = density_unit.to(target_unit)
    print(f"SUCCESS: {density_unit} = {result} {target_unit}")
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")
print()

# Test 2: Using to_system() - this fails with the bug
print("-" * 70)
print("Test 2: Using to_system() method")
print("-" * 70)

# Create a module-like object with bases attribute
class CustomSystem:
    """Custom unit system with [length, velocity, G] as bases."""
    kpc = u.kpc
    km_s = u.km / u.s
    G_unit = G.unit
    bases = {u.kpc, u.km / u.s, G.unit}

print(f"System bases: {CustomSystem.bases}")
print()

try:
    result = density_unit.to_system(CustomSystem)
    print(f"SUCCESS: {density_unit}.to_system(CustomSystem) = {result}")
except Exception as e:
    print(f"FAILED: {type(e).__name__}: {e}")
print()

# Additional diagnostics
print("-" * 70)
print("Diagnostic information")
print("-" * 70)

# What does decompose give us?
print(f"\ndensity_unit.decompose():")
print(f"  {density_unit.decompose()}")

print(f"\nTrying to decompose with custom bases:")
try:
    decomposed = density_unit.decompose(bases=CustomSystem.bases)
    print(f"  {decomposed}")
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")

# Show what the bases decompose to
print(f"\nBase unit decompositions:")
for base in CustomSystem.bases:
    print(f"  {base} -> {base.decompose()}")

print()
print("=" * 70)
print("End of test")
print("=" * 70)

