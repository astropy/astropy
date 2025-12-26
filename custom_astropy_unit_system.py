"""
Custom Astropy Unit System with [length, velocity, G] as base units.

This module demonstrates the issue #19045 where to_system() fails for 
density units in a custom unit system that uses non-canonical base units.

The expected behavior: density = velocity^2 / (length^2 * G)
"""

import astropy.units as u
from astropy.constants import G


# Define a custom unit system with [length, velocity, G] as bases
# This is a non-standard choice since velocity and G are composite units
# in terms of SI base units (length/time and length^3/(mass*time^2) respectively)

# Create custom base units
length = u.kpc  # kiloparsec for length
velocity = u.km / u.s  # km/s for velocity

# Gravitational constant as a unit
G_unit = G.to(u.kpc**3 / (u.Msun * u.Gyr**2)).unit

# Define the custom unit system as a module-like namespace
class CustomUnitSystem:
    """
    A custom unit system using [length, velocity, G] as base units.
    
    This is useful for gravitational N-body simulations where these
    are natural units.
    """
    # The units available in this system
    kpc = u.kpc
    km_s = u.km / u.s
    Msun = u.Msun
    Gyr = u.Gyr
    
    # The base units for this system (required by to_system)
    bases = {u.kpc, u.km / u.s, G.unit}


# Alternative simpler system for testing
class SimpleCustomSystem:
    """
    A simpler custom unit system for testing.
    """
    m = u.m
    s = u.s
    kg = u.kg
    bases = {u.m, u.s, u.kg}

