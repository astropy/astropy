# IERS Data Handling in Astropy - Exploration Summary

## Overview
The Astropy library contains comprehensive IERS (International Earth Rotation and Reference Systems) data handling in the `astropy.utils.iers` module. This system manages Earth rotation parameters, particularly UT1-UTC offsets and polar motion values, which are critical for accurate time and coordinate transformations.

---

## 1. File Locations and Purposes

### Main IERS Module
- **[astropy/utils/iers/iers.py](astropy/utils/iers/iers.py)** - Core IERS implementation
  - Main classes: `IERS`, `IERS_A`, `IERS_B`, `IERS_Auto`, `LeapSeconds`
  - Interpolation logic, data handling, and table management
  - ~1200 lines, well-documented with docstrings

### Supporting Files
- **astropy/utils/iers/__init__.py** - Module initialization and imports
- **astropy/utils/iers/data/** - Data directory (bundled IERS files)
- **astropy/utils/iers/tests/test_iers.py** - Comprehensive test suite

### Time Integration
- **[astropy/time/core.py](astropy/time/core.py)** - Time class integration
  - `Time.get_delta_ut1_utc()` - Public method to get UT1-UTC differences (line 2480)
  - `Time._get_delta_ut1_utc()` - Internal method for scale conversions (line 2530)
  - Automatic IERS table interpolation when accessing `time.ut1` scale

### Coordinates Integration
- **[astropy/coordinates/builtin_frames/utils.py](astropy/coordinates/builtin_frames/utils.py)** - Coordinate transformations
  - `get_dut1utc(time)` - Get UT1-UTC for coordinate transformations (line 70)
  - `_warn_iers()` - Warning handler for IERS range errors
  - Used for celestial-to-terrestrial coordinate conversions

---

## 2. Current Interpolation Approach

### Interpolation Method: Linear Interpolation

The interpolation is implemented in the `_interpolate()` method ([line 441-495](astropy/utils/iers/iers.py#L441)):

```python
# Linear interpolation between two adjacent table entries
val = val_0 + (mjd - mjd_0 + utc) / (mjd_1 - mjd_0) * d_val
```

**Key Features:**
- Uses **Modified Julian Date (MJD)** indexing via `np.searchsorted()`
- Interpolates between integer MJD table entries
- Fractional day handled by `utc` component (offset within the day)
- Handles edge cases:
  - Clips indices to valid table range (first/last values propagated for out-of-range)
  - Leap second corrections applied during interpolation

### Interpolation Parameters

**Input:**
- `jd1`, `jd2`: Two-part Julian Date representation
- Can be Time objects or float arrays
- Returns as floats or arrays depending on input

**Output:**
- Interpolated values with optional status codes
- Status values: `FROM_IERS_B`, `FROM_IERS_A`, `FROM_IERS_A_PREDICTION`, `TIME_BEFORE_IERS_RANGE`, `TIME_BEYOND_IERS_RANGE`

### Leap Second Handling During Interpolation

In the `_interpolate()` method (line 465-469):
```python
if column == "UT1_UTC":
    # Check & correct for possible leap second (correcting diff.,
    # not 1st point, since jump can only happen right at 2nd point)
    d_val -= d_val.round()
```

This ensures smooth interpolation across leap second boundaries by removing discontinuities in the difference between table entries.

---

## 3. Data Sources and IERS Files

### IERS Data Sources

**IERS Bulletin A (IERS_A class)**
- Updated weekly
- Contains historical data from 1973-01-01 onward
- Includes ~1 year of predictive data into the future
- Combines with Bulletin B values where available
- File format: `finals2000A.*` from USNO

**IERS Bulletin B (IERS_B class)**
- Updated monthly
- Final/definitive Earth orientation values
- Historical data from 1962 to present
- File format: `eopc04_IAU2000.*` from IERS

**Auto-Management (IERS_Auto class)**
- Intelligently combines IERS-A and IERS-B data
- Prefers IERS-B values where available
- Auto-downloads latest files when needed
- Configuration-controlled behavior

### Bundled Data

- Provided by `astropy-iers-data` package
- Automatically installed with Astropy v6.0+
- Regularly updated, covering ~1962 to ~1 year future

### Data URLs

- **Primary**: `https://maia.usno.navy.mil/ser7/finals2000A.all`
- **Mirror**: USNO mirror
- **Leap Seconds**: `https://www.iers.org/`

---

## 4. Data Structures for IERS Parameters

### QTable Columns in IERS Tables

The IERS tables inherit from `astropy.table.QTable` with the following key columns:

**Core Columns:**
- `MJD` - Modified Julian Date (integer part, u.d)
- `UT1_UTC` - UT1 minus UTC offset (u.s)
- `UT1_UTC_A` - UT1-UTC from Bulletin A (u.s)
- `UT1_UTC_B` - UT1-UTC from Bulletin B (u.s)
- `UT1Flag` or `UT1Flag_A` - Source indicator ('I'=IERS, 'P'=Prediction, 'B'=Bulletin B)

**Polar Motion Columns:**
- `PM_x`, `PM_y` - Polar motion coordinates (u.arcsec)
- `PM_x_A`, `PM_y_A` - Polar motion from Bulletin A
- `PM_X_B`, `PM_Y_B` - Polar motion from Bulletin B
- `PolPMFlag` or `PolPMFlag_A` - Source indicator

**CIP Correction Columns:**
- `dX_2000A`, `dY_2000A` - Celestial Intermediate Pole corrections (u.arcsec)
- `dX_2000A_A`, `dY_2000A_A` - CIP corrections from Bulletin A
- `dX_2000A_B`, `dY_2000A_B` - CIP corrections from Bulletin B
- `NutFlag` or `NutFlag_A` - Source indicator

**Metadata:**
- `predictive_index` - Index of first predictive row
- `predictive_mjd` - MJD of first predictive value
- `data_url` - Source URL for the data
- `data_path` - Local file path

### Class Hierarchy

```
QTable (from astropy.table)
  ↓
IERS (Generic base class)
  ├── IERS_A (Bulletin A reader, inherits IERS)
  ├── IERS_B (Bulletin B reader, inherits IERS)
  └── IERS_Auto (Auto-downloading manager, inherits IERS_A)

LeapSeconds (QTable)
  - TAI-UTC differences
  - Columns: 'year', 'month', 'tai_utc'
```

---

## 5. UT1-UTC Interpolation Implementation

### Main Public Interface: `Time.get_delta_ut1_utc()`

Location: [astropy/time/core.py, line 2480](astropy/time/core.py#L2480)

```python
def get_delta_ut1_utc(self, iers_table=None, return_status=False):
    """Find UT1 - UTC differences by interpolating in IERS Table."""
    if iers_table is None:
        from astropy.utils.iers import earth_orientation_table
        iers_table = earth_orientation_table.get()
    return iers_table.ut1_utc(self.utc, return_status=return_status)
```

**Usage Examples:**
```python
from astropy.time import Time
t = Time('2010-01-01 00:00:00', format='iso', scale='utc')

# Automatic interpolation (called when accessing .ut1)
t_ut1 = t.ut1

# Explicit call
delta_ut1_utc = t.get_delta_ut1_utc()  # Returns Quantity in seconds

# Get status information
delta, status = t.get_delta_ut1_utc(return_status=True)
# status values: FROM_IERS_B (0), FROM_IERS_A (1), FROM_IERS_A_PREDICTION (2)
#                TIME_BEFORE_IERS_RANGE (-1), TIME_BEYOND_IERS_RANGE (-2)
```

### IERS Table Method: `IERS.ut1_utc()`

Location: [astropy/utils/iers/iers.py, line 315](astropy/utils/iers/iers.py#L315)

Calls the internal `_interpolate()` method:
```python
def ut1_utc(self, jd1, jd2=0.0, return_status=False):
    """Interpolate UT1-UTC corrections in IERS Table for given dates."""
    return self._interpolate(
        jd1, jd2, ["UT1_UTC"], 
        self.ut1_utc_source if return_status else None
    )
```

### Internal Interpolation: `IERS._interpolate()`

Location: [astropy/utils/iers/iers.py, line 441-495](astropy/utils/iers/iers.py#L441)

**Algorithm:**
1. Convert input JD to MJD and UTC fraction
2. Ensure table is up-to-date (IERS_Auto refresh)
3. Find bracketing indices using `np.searchsorted(..., side="right")`
4. For each requested column:
   - Extract values at two surrounding MJD integers: `val_0, val_1`
   - Compute difference: `d_val = val_1 - val_0`
   - If UT1_UTC: Correct for leap seconds: `d_val -= d_val.round()`
   - Linear interpolate: `val = val_0 + fraction * d_val`
   - Handle out-of-range: propagate boundary values
5. Optionally compute status codes via `ut1_utc_source()` method

### Status Code Assignment

**Source Methods:**
- `ut1_utc_source(i)` - Determines UT1-UTC data source
- `dcip_source(i)` - Determines CIP correction source  
- `pm_source(i)` - Determines polar motion source

These methods check the flags (`UT1Flag`, `PolPMFlag`, `NutFlag`) to determine:
- `FROM_IERS_B` (0) - Official final values from Bulletin B
- `FROM_IERS_A` (1) - Near-real-time values from Bulletin A
- `FROM_IERS_A_PREDICTION` (2) - Predictive values (extended forecast)

---

## 6. Leap Second Handling

### Leap Second Data Structure

Location: [astropy/utils/iers/iers.py, line 1037](astropy/utils/iers/iers.py#L1037)

```python
class LeapSeconds(QTable):
    """Leap seconds class, holding TAI-UTC differences."""
    # Columns: 'year', 'month', 'tai_utc'
    # Methods: open(), from_erfa(), from_iers_leap_seconds(), from_ietf_leap_seconds()
```

**Sources:**
- IERS Leap Second data (`Leap_Second.dat`)
- IETF/NTP leap seconds list (`leap-seconds.list`)
- ERFA/SOFA built-in values

**Configuration Options:**
```python
from astropy.utils.iers import conf
conf.iers_leap_second_auto_url  # IERS URL
conf.ietf_leap_second_auto_url  # IETF URL
conf.system_leap_second_file    # System leap second file
```

### Leap Second Correction in UT1-UTC Interpolation

The key insight is in the UT1-UTC interpolation (line 465-469):

```python
if column == "UT1_UTC":
    # Check & correct for possible leap second
    d_val -= d_val.round()
```

**Mechanism:**
- UT1-UTC changes by approximately -1 second when leap seconds are introduced
- The difference `d_val = val_1 - val_0` can be ~±1 second at leap second boundaries
- Using `.round()` removes this discontinuity
- Linear interpolation then works smoothly across the boundary

### Automatic Scale Conversion with Leap Seconds

The `Time` class automatically handles leap second transitions:
```python
# Time.py internal scale conversion logic
if scale == "ut1":
    # Convert from UT1 to UTC using interpolated UT1-UTC
    jd1_utc, jd2_utc = erfa.ut1utc(jd1, jd2, delta.to_value(u.s))
    # Then re-interpolate for improved accuracy near leap seconds
    delta = iers_table.ut1_utc(jd1_utc, jd2_utc)
```

---

## 7. Time Coordinate Transformations Using IERS

### Coordinate Frame Integration

Location: [astropy/coordinates/builtin_frames/utils.py](astropy/coordinates/builtin_frames/utils.py)

**Helper Functions:**
- `get_dut1utc(time)` - Safely retrieves UT1-UTC with error handling
- `_warn_iers(ierserr)` - Generates warnings for out-of-range times

```python
def get_dut1utc(time):
    """Get UT1-UTC in coordinates with graceful error handling."""
    try:
        return time.delta_ut1_utc
    except iers.IERSRangeError as e:
        warnings.warn(f"Assuming UT1-UTC=0 for coordinate transformations.")
        return np.zeros(time.shape)
```

### Polar Motion Access

```python
def get_polar_motion_iers1980(obstime):
    """Get polar motion (x_p, y_p) for Earth Orientation."""
    iers_table = iers.earth_orientation_table.get()
    return iers_table.pm_xy(obstime, return_status=True)
```

### CIP Correction Access

```python
def get_dcip_iers2000a(obstime):
    """Get Celestial Intermediate Pole corrections (dX, dY)."""
    iers_table = iers.earth_orientation_table.get()
    return iers_table.dcip_xy(obstime, return_status=True)
```

### ERFA Integration

The ERFA library (Essential Routines for Fundamental Astronomy) is called with IERS-derived parameters:
- `erfa.ut1utc()` - Convert between UT1 and UTC
- `erfa.utcut1()` - Convert between UTC and UT1
- `erfa.era00()` - Earth Rotation Angle using UT1

---

## 8. Auto-Download and Refresh Mechanism

### IERS_Auto Class

Location: [astropy/utils/iers/iers.py, line 773-960](astropy/utils/iers/iers.py#L773)

**Features:**
- Automatically manages IERS table updates
- Checks if predictive values are stale
- Downloads latest IERS-A file from IERS service when needed
- Intelligently merges with IERS-B data for maximum accuracy

**Key Methods:**
- `open()` - Opens/downloads latest IERS table
- `_refresh_table_as_needed()` - Checks and updates if necessary
- `_substitute_iers_b()` - Merges IERS-B official values
- `_check_interpolate_indices()` - Validates interpolation range

### Configuration Control

```python
from astropy.utils.iers import conf

# Control auto-downloading
conf.auto_download = True/False

# Maximum age of predictive data before refresh (days)
conf.auto_max_age = 30  # default

# Remote download timeout (seconds)
conf.remote_timeout = 20.0

# Degraded accuracy handling
conf.iers_degraded_accuracy = 'error'  # or 'warn', 'ignore'

# URL configuration
conf.iers_auto_url = "https://maia.usno.navy.mil/ser7/finals2000A.all"
conf.iers_auto_url_mirror = "..."  # Mirror URL
```

### Refresh Logic

The `_refresh_table_as_needed()` method (line 869-953) performs:

1. **Check conditions:**
   - Is `auto_download` enabled?
   - Are any requested times in the predictive range?
   - Are predictive values older than `auto_max_age` days?

2. **Download if needed:**
   - Fetch latest IERS-A file from primary/mirror URLs
   - Cache the file for future use

3. **Update in place:**
   - Replace predictive section of existing table
   - Add new rows as needed
   - Validate continuity (no gaps in MJD)

4. **Handle failures:**
   - Issue warnings if download fails
   - Continue with cached version
   - Error only if interpolating predictive values

---

## 9. Configuration and Behavior

### ScienceState Management

The `earth_orientation_table` is managed as a `ScienceState`:
```python
earth_orientation_table = ScienceState(iers.IERS_Auto.open)
```

This allows:
- Temporary configuration changes via context managers
- Thread-safe state management
- Easy switching between data sources

**Usage:**
```python
from astropy.utils.iers import earth_orientation_table, IERS_B

# Temporarily use IERS_B instead of IERS_Auto
with earth_orientation_table.set(IERS_B.open()):
    # Transformations here use IERS_B
    pass
```

### Exception Handling

**Main Exceptions:**
- `IERSRangeError` - Time outside IERS table range
- `IERSWarning` - Generic IERS warnings
- `IERSStaleWarning` - Predictive data is stale
- `IERSDegradedAccuracyWarning` - Using only IERS-B with reduced accuracy

---

## 10. Key Files Summary Table

| File | Purpose | Key Classes/Functions | Lines |
|------|---------|----------------------|-------|
| [astropy/utils/iers/iers.py](astropy/utils/iers/iers.py) | Core IERS | IERS, IERS_A, IERS_B, IERS_Auto, LeapSeconds | ~1200 |
| [astropy/time/core.py](astropy/time/core.py) | Time integration | Time._get_delta_ut1_utc(), Time.get_delta_ut1_utc() | 2480-2550 |
| [astropy/coordinates/builtin_frames/utils.py](astropy/coordinates/builtin_frames/utils.py) | Coordinate transforms | get_dut1utc(), get_polar_motion_iers1980() | 50-110 |
| docs/utils/iers.rst | Documentation | Usage examples and configuration | - |
| [astropy/time/tests/test_ut1.py](astropy/time/tests/test_ut1.py) | UT1 tests | UT1 scale conversion tests | - |
| [astropy/utils/iers/tests/test_iers.py](astropy/utils/iers/tests/test_iers.py) | IERS tests | IERS interpolation tests | - |

---

## Summary

### Interpolation Method
- **Linear interpolation** between daily MJD entries
- **Leap second-aware**: Removes discontinuities across boundaries
- **Status-aware**: Tracks data source (Bulletin A/B, predictive)

### Data Management
- **Multiple sources**: IERS Bulletin A (predictive), Bulletin B (definitive), bundled data
- **Auto-download**: IERS_Auto class manages intelligent caching and updates
- **Leap seconds**: Separate LeapSeconds class handles TAI-UTC conversions

### Integration
- **Time scales**: UTC ↔ UT1 conversions via interpolated UT1-UTC
- **Coordinates**: Earth orientation parameters for GCRS ↔ ITRS transformations
- **ERFA**: Uses ERFA library for fundamental astronomy calculations

### Performance Considerations
- **Caching**: Tables cached in class attributes
- **Vectorization**: NumPy operations for batch interpolation
- **Conditional refresh**: Only updates when necessary and predictive data is stale
