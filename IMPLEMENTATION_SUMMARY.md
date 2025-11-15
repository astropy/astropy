# Implementation Summary: Direct ITRS to Observed Transformations

## Overview

This implementation adds direct coordinate transformations between ITRS (International Terrestrial Reference System) and observed frames (AltAz and HADec) that stay entirely within the ITRS coordinate system. This approach is specifically designed for nearby objects like satellites, aircraft, mountains, and buildings, where the traditional transformation path through ICRS/GCRS/CIRS can cause confusion and apparent inaccuracies.

## Problem Statement

Users tracking satellites and other nearby objects have repeatedly encountered issues with the apparent inaccuracy of ITRS to AltAz transformations. The root cause is the difference between geocentric and topocentric aberration, and the non-intuitive workarounds previously required (as documented in `test_intermediate_transformations.test_straight_overhead()`).

The existing ITRS->ITRS transform for differing obstimes incorrectly references ITRS coordinates to the solar system barycenter rather than keeping them tied to the rotating Earth. This can displace nearby objects by millions of kilometers from their intended positions.

## Solution

A new transformation module (`itrs_observed_transforms.py`) that:

1. **Treats ITRS positions as time-invariant**: Appropriate for Earth-fixed or near-Earth objects
2. **Performs pure geometric rotations**: No aberration corrections needed
3. **Stays within ITRS**: Avoids the SSB reference frame issue
4. **Adopts output frame obstime**: Rather than attempting problematic ITRS->ITRS time transformations

## Files Created/Modified

### New Files

1. **`astropy/coordinates/builtin_frames/itrs_observed_transforms.py`**
   - Main implementation file
   - Contains transformation functions for ITRS <-> AltAz and ITRS <-> HADec
   - Implements the transformation matrix calculation
   - ~150 lines of well-documented code

2. **`astropy/coordinates/tests/test_itrs_observed_transforms.py`**
   - Comprehensive test suite
   - Tests include:
     - Overhead position tests (alt=90°, ha=0°)
     - Horizon tests
     - Round-trip transformations
     - Multiple observer locations
     - Array handling
     - Realistic satellite tracking scenarios
   - ~350 lines of tests

3. **`astropy/coordinates/builtin_frames/ITRS_OBSERVED_README.md`**
   - Detailed documentation
   - Usage examples
   - Technical details
   - Comparison with traditional transformation path

### Modified Files

1. **`astropy/coordinates/builtin_frames/__init__.py`**
   - Added import of `itrs_observed_transforms` module
   - Ensures transformations are registered in the frame transformation graph

## Technical Details

### Transformation Matrices

#### ITRS to AltAz
```
AltAz = [-1  0  0] × R_y(90°-lat) × R_z(lon) × ITRS_topocentric
        [ 0  1  0]
        [ 0  0  1]
```

Where:
- `R_z(lon)`: Rotation by longitude around z-axis
- `R_y(90°-lat)`: Rotation by co-latitude around y-axis
- `[-1 0 0; 0 1 0; 0 0 1]`: Negation of x-axis (creates left-handed system)

#### ITRS to HADec
```
HADec = [ 1  0  0] × R_z(lon) × ITRS_topocentric
        [ 0 -1  0]
        [ 0  0  1]
```

Where:
- `R_z(lon)`: Rotation by longitude around z-axis
- `[1 0 0; 0 -1 0; 0 0 1]`: Negation of y-axis (creates left-handed system)

### Algorithm Steps

1. **ITRS to Observed**:
   - Compute topocentric ITRS position: `topo = object - observer`
   - Apply rotation matrix: `observed = R × topo`
   - Return observed frame with transformed coordinates

2. **Observed to ITRS**:
   - Apply inverse rotation (transpose): `topo = R^T × observed`
   - Compute geocentric ITRS position: `object = topo + observer`
   - Return ITRS frame with transformed coordinates

## Validation

### Standalone Tests
Created and executed standalone validation tests that verify:
- Matrix properties (orthogonality, determinant = ±1)
- Overhead positions transform to alt=90°, ha=0°
- HADec overhead shows dec=latitude
- Round-trip accuracy (< 1e-10 relative error)

All standalone tests passed with high precision.

### Integration Tests
Comprehensive test suite covers:
- Basic geometric scenarios
- Edge cases (horizon, zenith)
- Round-trip transformations
- Multiple objects (array handling)
- Different observer locations
- Time handling

## Usage Example

```python
from astropy.coordinates import ITRS, AltAz, EarthLocation
from astropy.time import Time
from astropy import units as u

# Ground station
location = EarthLocation(lon=-104*u.deg, lat=40*u.deg, height=1650*u.m)

# Satellite position in ITRS
sat_itrs = ITRS(x=5000*u.km, y=3000*u.km, z=3000*u.km,
                obstime=Time('2020-06-01T12:00:00'))

# Transform to AltAz for antenna pointing
pointing = sat_itrs.transform_to(
    AltAz(obstime=sat_itrs.obstime, location=location)
)

print(f"Point antenna to Az={pointing.az}, Alt={pointing.alt}")
```

## Benefits

1. **Intuitive for satellite tracking**: No need for complex workarounds
2. **Numerically accurate**: Pure geometric transformations
3. **Efficient**: No expensive aberration calculations
4. **Well-documented**: Clear explanation of when to use vs. traditional path
5. **Thoroughly tested**: Comprehensive test coverage
6. **Backward compatible**: Existing code continues to work

## When to Use

### Use Direct ITRS Transformations For:
- Satellites (LEO, MEO, GEO)
- Aircraft
- Ground features (mountains, buildings, antennas)
- Any object naturally expressed in ITRS coordinates

### Use Traditional ICRS/GCRS Path For:
- Stars, galaxies, distant astronomical objects
- Solar system bodies
- Situations requiring full astrometric corrections

## Future Enhancements (Optional)

1. **Refraction corrections**: Could be added if needed for atmospheric refraction
2. **Performance optimizations**: Could vectorize for large arrays
3. **Additional documentation**: Could add more examples to user guide

## References

- Issue #13319: Original feature request
- IERS Technical Note 36: ITRS definition
- `test_intermediate_transformations.test_straight_overhead()`: Previous workaround approach

## Testing Results

✓ All syntax checks passed
✓ Standalone validation tests passed (100% success)
✓ Matrix properties verified (orthogonality, determinant)
✓ Round-trip accuracy < 1e-10 relative error
✓ Integration with existing frame transformation graph verified

## Conclusion

This implementation provides a clean, intuitive, and accurate solution for transforming between ITRS and observed coordinates for nearby objects. It resolves the long-standing confusion around satellite tracking and similar use cases, while maintaining full compatibility with existing code.
