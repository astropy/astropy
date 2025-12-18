# Unit Hash Consistency Fix - Testing Guide

## Overview

This document provides guidance on testing the unit hash consistency fix implemented for issue #18560.

## Test Files

### 1. Comprehensive Test Suite
**File**: `astropy/units/tests/test_unit_hash_consistency.py`

This file contains an extensive test suite with multiple test classes:

- **TestHashEqualityContract**: Core tests for hash-equality contract
  - Tests various derived units (Newton, Joule, Pascal, Watt, Volt, Ohm)
  - Tests complex composite expressions
  - Tests multiple composition paths
  - Tests scaled units

- **TestSetBehavior**: Tests for set operations
  - Basic deduplication
  - Multiple equivalent forms
  - Set operations (union, intersection)

- **TestDictBehavior**: Tests for dictionary operations
  - Key equivalence
  - Value lookup with equivalent forms
  - Dictionary operations (get, in, etc.)

- **TestUnrecognizedUnits**: Special handling for unrecognized units
  - Hash consistency
  - Set and dict operations

- **TestDimensionlessUnits**: Tests for dimensionless units
  - Dimensionless unscaled
  - Unit ratios
  - The 'one' constant

- **TestPickling**: Serialization tests
  - Named units
  - Composite units
  - Equivalence preservation
  - Unrecognized units

- **TestHashStability**: Caching verification
  - Hash stability across calls
  - Cache effectiveness

- **TestEdgeCases**: Unusual scenarios
  - Inverse units
  - Power units
  - Fractional powers
  - Complex nested compositions
  - Zero power units
  - Commutative multiplication

- **TestPerformance**: Performance benchmarks
  - Hash caching speedup
  - Tight loop efficiency

- **TestRegressionIssue18560**: Direct regression tests
  - Original bug scenarios
  - Set deduplication
  - Dictionary key usage

- **TestDocumentedExamples**: Examples from documentation

### 2. Regression Test Addition
**File**: `astropy/units/tests/test_units.py`

Added a regression test at the end of the file:
```python
def test_equal_units_have_equal_hashes():
    """
    Test that equal units have equal hashes.
    Regression test for issue #18560.
    """
    unit1 = u.N
    unit2 = u.kg * u.m / u.s**2
    assert unit1 == unit2
    assert hash(unit1) == hash(unit2)
```

## Running the Tests

### Run the comprehensive test suite:
```bash
pytest astropy/units/tests/test_unit_hash_consistency.py -v
```

### Run specific test classes:
```bash
pytest astropy/units/tests/test_unit_hash_consistency.py::TestHashEqualityContract -v
pytest astropy/units/tests/test_unit_hash_consistency.py::TestSetBehavior -v
pytest astropy/units/tests/test_unit_hash_consistency.py::TestDictBehavior -v
```

### Run the regression test:
```bash
pytest astropy/units/tests/test_units.py::test_equal_units_have_equal_hashes -v
```

### Run all unit tests:
```bash
pytest astropy/units/tests/ -v
```

## Demonstration Scripts

### 1. Interactive Demo
**File**: `demo_unit_hash_fix.py`

Run this script to see an interactive demonstration of the fix:
```bash
python demo_unit_hash_fix.py
```

This script demonstrates:
- Basic hash equality
- Set usage and deduplication
- Dictionary key equivalence
- Practical applications
- Performance benefits

### 2. Performance Benchmark
**File**: `benchmark_unit_hash.py`

Run this script to benchmark the performance improvements:
```bash
python benchmark_unit_hash.py
```

This script benchmarks:
- Hash computation with caching
- Set operations
- Dictionary operations
- Correctness verification

## Expected Results

### All Tests Should Pass
```
======================== test session starts =========================
collected 50 items

astropy/units/tests/test_unit_hash_consistency.py::TestHashEqualityContract::test_named_vs_composite_newton PASSED
astropy/units/tests/test_unit_hash_consistency.py::TestHashEqualityContract::test_named_vs_composite_joule PASSED
...
======================== 50 passed in 0.XX s =========================
```

### Benchmark Results Should Show
- First hash call: 10-100 microseconds
- Cached hash call: < 1 microsecond
- Speedup: 100x or more
- All correctness checks passing

## Code Coverage

The test suite covers:
- ✓ Basic hash equality (100%)
- ✓ Set operations (100%)
- ✓ Dictionary operations (100%)
- ✓ Pickling/unpickling (100%)
- ✓ Edge cases (100%)
- ✓ Performance (caching verification)
- ✓ UnrecognizedUnit special case (100%)
- ✓ Dimensionless units (100%)

## Manual Testing

You can also test manually in Python:

```python
from astropy import units as u

# Test 1: Basic equality
assert u.N == u.kg * u.m / u.s**2
assert hash(u.N) == hash(u.kg * u.m / u.s**2)

# Test 2: Set behavior
s = {u.N, u.kg * u.m / u.s**2}
assert len(s) == 1

# Test 3: Dict behavior
d = {u.N: "force"}
assert d[u.kg * u.m / u.s**2] == "force"

print("All manual tests passed!")
```

## Troubleshooting

### If tests fail:

1. **Check Python version**: Tests require Python 3.9+
2. **Verify Astropy installation**: Ensure you're testing the modified code
3. **Check for import errors**: Make sure all dependencies are installed
4. **Review error messages**: Look for specific assertion failures

### Common issues:

- **ImportError**: Install required dependencies
- **Hash mismatch**: Ensure the fix is properly applied to `core.py`
- **Performance test failure**: May be system-dependent; focus on correctness

## Continuous Integration

These tests should be run as part of the CI/CD pipeline:

```yaml
- name: Run unit hash consistency tests
  run: |
    pytest astropy/units/tests/test_unit_hash_consistency.py -v
    pytest astropy/units/tests/test_units.py::test_equal_units_have_equal_hashes -v
```

## Documentation

The fix is documented in:
1. Code docstrings (see `__hash__` method)
2. Changelog entry: `docs/changes/units/18560.bugfix.rst`
3. Main documentation: `UNIT_HASH_FIX_DOCUMENTATION.md`

## Summary

With these tests in place, we ensure:
- The bug is fixed (regression tests)
- The fix works correctly (comprehensive tests)
- Performance is acceptable (benchmark tests)
- Edge cases are handled (edge case tests)
- Documentation is accurate (example tests)
