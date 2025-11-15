# Direct ITRS to Observed Transformations

## Overview

This module provides direct transformations between ITRS (International Terrestrial Reference System) and observed coordinate frames (AltAz and HADec). These transformations stay entirely within the ITRS coordinate system and are specifically designed for nearby objects like satellites, aircraft, mountains, and buildings.

## Motivation

The traditional transformation path for converting ITRS coordinates to observed frames goes through ICRS/GCRS/CIRS, which introduces complications related to geocentric versus topocentric aberration. This is appropriate for distant astronomical objects but causes confusion and apparent inaccuracies when observing nearby objects that are fixed to or moving near Earth's surface.

The issue is that ITRS positions are tied to the Earth, but the existing ITRS->ITRS transform at different times incorrectly references these positions to the solar system barycenter, potentially displacing them by millions of kilometers.

## Key Features

1. **Time-Invariant ITRS Positions**: The transformations treat ITRS coordinates as time-invariant, which is appropriate for Earth-fixed or near-Earth objects.

2. **No Aberration Corrections**: Since these are geometric transformations within the ITRS, no stellar aberration corrections are applied.

3. **Direct Rotations**: The transformations use simple rotation matrices to convert between geocentric ITRS Cartesian coordinates and topocentric observed coordinates.

## When to Use These Transformations

### Use direct ITRS transformations for:
- Satellites (LEO, MEO, GEO)
- Aircraft
- Ground-based features (mountains, buildings, antennas)
- Any object whose position is naturally expressed in ITRS coordinates
- Scenarios where you want to avoid aberration effects

### Use traditional ICRS/GCRS transformations for:
- Stars, galaxies, and other distant astronomical objects
- Solar system bodies (planets, asteroids, comets)
- Scenarios where you need full astrometric corrections including aberration

## Usage Examples

### Satellite Tracking

```python
from astropy.coordinates import ITRS, AltAz, EarthLocation
from astropy.time import Time
from astropy import units as u

# Ground station location
ground_station = EarthLocation(lon=-104*u.deg, lat=40*u.deg, height=1650*u.m)

# Satellite position in ITRS (from TLE propagation, for example)
sat_position = ITRS(x=5000*u.km, y=3000*u.km, z=3000*u.km,
                    obstime=Time('2020-06-01T12:00:00'))

# Get pointing angles for the antenna
pointing = sat_position.transform_to(AltAz(obstime=sat_position.obstime,
                                            location=ground_station))

print(f"Azimuth: {pointing.az}")
print(f"Altitude: {pointing.alt}")
print(f"Distance: {pointing.distance}")
```

### Mountain Peak Observation

```python
from astropy.coordinates import ITRS, AltAz, EarthLocation
from astropy.time import Time

# Observer location
observer = EarthLocation(lon=-105*u.deg, lat=39*u.deg, height=2000*u.m)

# Mountain peak in ITRS coordinates
peak = EarthLocation(lon=-105.5*u.deg, lat=39.5*u.deg, height=4000*u.m)
peak_itrs = peak.get_itrs(Time('2020-01-01'))

# Get viewing direction
view = peak_itrs.transform_to(AltAz(obstime=Time('2020-01-01'),
                                     location=observer))

print(f"Look toward azimuth {view.az} at altitude {view.alt}")
```

## Technical Details

### Transformation Matrices

#### ITRS to AltAz
The transformation from ITRS to AltAz involves:
1. Rotation by observer's longitude around the z-axis
2. Rotation by (90° - latitude) around the y-axis
3. Negation of the x-axis to create the left-handed AltAz system

#### ITRS to HADec
The transformation from ITRS to HADec involves:
1. Rotation by observer's longitude around the z-axis
2. Negation of the y-axis to create the left-handed HADec system

### Coordinate System Conventions

- **ITRS**: Right-handed, geocentric Cartesian system
- **AltAz**: Left-handed, topocentric system (azimuth increases eastward from north)
- **HADec**: Left-handed, topocentric equatorial system (hour angle increases westward)

## Implementation Notes

1. The transformations are registered in the frame transformation graph, so they work seamlessly with Astropy's coordinate transformation infrastructure.

2. The transformations automatically handle:
   - Conversion from geocentric to topocentric coordinates
   - Proper handling of observer location
   - Distance preservation

3. The `obstime` attribute:
   - For ITRS -> Observed: The output frame's obstime is adopted without transforming the ITRS position
   - For Observed -> ITRS: The ITRS position is computed at the observed frame's obstime

## Testing

Comprehensive tests are provided in `test_itrs_observed_transforms.py`, including:
- Overhead position tests (altitude should be 90°)
- Horizon tests
- Round-trip transformations
- Multiple observer locations
- Array handling
- Realistic satellite tracking scenarios

## References

- Issue #13319: Discussion of ITRS to AltAz transformation issues
- `test_intermediate_transformations.test_straight_overhead()`: Previous workaround approach
- ITRS Definition: IERS Technical Note 36

## Authors

- Implementation based on proposal in issue #13319
- Tests and documentation by Stitch SWE
