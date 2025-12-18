# Unit Hash Consistency Fix

## Problem Statement

The Astropy units system had a critical hash consistency issue that violated Python's hash-equality contract. This caused unexpected behavior when using units in sets and dictionaries.

### The Core Issue

Python's data model requires that:
> If two objects compare equal (via `__eq__`), they must have the same hash value.

However, in Astropy's unit system, equivalent units that compared as equal had different hash values:

```python
import astropy.units as u

# These two units are mathematically equivalent
unit1 = u.N                    # Newton (named unit)
unit2 = u.kg * u.m / u.s**2   # Newton in base SI units (composite unit)

# They correctly compare as equal
assert unit1 == unit2  # True ✓

# But they had DIFFERENT hash values (WRONG!)
assert hash(unit1) == hash(unit2)  # False ✗
```

### Why This Was a Problem

This hash inconsistency caused serious issues with Python's hash-based data structures:

1. **Sets were broken**: Adding equivalent units to a set kept both instead of treating them as duplicates
   ```python
   s = {u.N, u.kg * u.m / u.s**2}
   len(s)  # Returns 2, but should be 1!
   ```

2. **Dictionary keys didn't work correctly**: The same physical unit could be stored under different keys
   ```python
   d = {u.N: "Newton"}
   d[u.kg * u.m / u.s**2]  # KeyError! (even though the units are equal)
   ```

3. **Violated Python's fundamental contract**: This could lead to unpredictable behavior in any code that relies on hash-based lookups

### Root Cause

The issue stemmed from how the `__hash__()` method was implemented in the `UnitBase` class:

**Original Implementation** (astropy/units/core.py):
```python
def __hash__(self):
    return hash((self.scale, tuple(self.bases), tuple(self.powers)))
```

This implementation hashed the **raw representation** of a unit:
- Named units like `u.N` have specific bases and powers in their internal representation
- Composite units like `u.kg * u.m / u.s**2` have different bases and powers
- Even though they represent the same physical quantity, they had different raw structures
- Result: Different hashes for equal units ✗

## Solution

The fix ensures that equal units always produce the same hash by using their **decomposed form** for hashing.

### Key Changes

**Modified File**: `astropy/units/core.py`

#### 1. Main Hash Implementation (Line 364-377)

**Before:**
```python
def __hash__(self):
    return hash((self.scale, tuple(self.bases), tuple(self.powers)))
```

**After:**
```python
def __hash__(self) -> int:
    return self._hash

@cached_property
def _hash(self) -> int:
    # Use decomposed form to ensure equal units have equal hashes.
    # This is necessary because units like u.N and u.kg*u.m/u.s**2
    # are equal (per __eq__) but have different raw representations.
    decomposed = self.decompose()
    return hash((decomposed.scale, *[x.name for x in decomposed.bases], *map(str, decomposed.powers)))
```

**Why This Works:**
- **Decomposition**: Both `u.N` and `u.kg * u.m / u.s**2` decompose to the exact same base SI units
- **Canonical form**: The decomposed form is a canonical representation that's identical for equal units
- **Caching**: Using `@cached_property` ensures the hash is computed only once per unit instance
- **Performance**: Hash computation is expensive, so caching improves performance significantly

#### 2. Pickling Support (Line 380-385)

**Added:**
```python
def __getstate__(self) -> dict[str, object]:
    # If we get pickled, we should *not* store the memoized members since
    # hashes of strings vary between sessions.
    state = self.__dict__.copy()
    state.pop("_hash", None)
    state.pop("_physical_type_id", None)
    return state
```

**Why This Is Needed:**
- Python's hash values for strings can vary between different Python sessions (hash randomization)
- When unpickling, the cached hash might be invalid
- By excluding `_hash` from the pickled state, we force recomputation on unpickle
- This ensures hash consistency across serialization boundaries

#### 3. UnrecognizedUnit Special Case (Line 1978-1987)

**Modified:**
```python
@cached_property
def _hash(self) -> int:
    # For unrecognized units, hash based on the name only
    # since they can't be decomposed to a meaningful form
    return hash(self.name)

# Explicitly inherit __hash__ from UnitBase to override Python's
# default behavior of setting __hash__ = None when __eq__ is defined
__hash__ = UnitBase.__hash__
```

**Why This Is Different:**
- `UnrecognizedUnit` represents units that couldn't be parsed (e.g., `FOO`)
- These units can't be decomposed, so we hash based on their name instead
- Two `UnrecognizedUnit` instances with the same name should hash equally
- The explicit `__hash__ = UnitBase.__hash__` line is necessary because `UnrecognizedUnit` defines `__eq__`, which would normally set `__hash__ = None`

### Technical Details

#### Hash Computation Algorithm

For standard units:
1. Decompose the unit to its base SI form (e.g., meters, kilograms, seconds)
2. Extract the scale factor
3. Extract the base unit names (as strings for stability)
4. Extract the powers (as strings for consistency)
5. Compute hash from the tuple: `(scale, base_names..., power_strings...)`

#### Why Use Decomposed Form?

| Unit | Raw Representation | Decomposed Form |
|------|-------------------|-----------------|
| `u.N` | bases=[N], powers=[1] | bases=[kg, m, s], powers=[1, 1, -2] |
| `u.kg * u.m / u.s**2` | bases=[kg, m, s], powers=[1, 1, -2] | bases=[kg, m, s], powers=[1, 1, -2] |

The decomposed forms are **identical**, ensuring equal hashes!

#### Performance Considerations

- **Caching is critical**: Decomposition is computationally expensive
- `@cached_property` ensures decomposition happens only once per unit instance
- The cached hash is stored in `_hash` attribute
- Subsequent hash calls return the cached value in O(1) time
- Memory trade-off: Small additional memory per unit for cached hash

## Verification

The fix was thoroughly tested to ensure correctness:

### Test Results

```python
import astropy.units as u

# Test 1: Equal units have equal hashes
assert u.N == u.kg * u.m / u.s**2  # ✓ Equal
assert hash(u.N) == hash(u.kg * u.m / u.s**2)  # ✓ Same hash

# Test 2: Works with multiple derived units
assert u.J == u.N * u.m  # ✓ Equal
assert hash(u.J) == hash(u.N * u.m)  # ✓ Same hash

# Test 3: Set behavior is correct
s = {u.N, u.kg * u.m / u.s**2}
assert len(s) == 1  # ✓ Only one element (treats as duplicate)

# Test 4: Dictionary behavior is correct
d = {u.N: "Newton", u.kg * u.m / u.s**2: "kg*m/s^2"}
assert len(d) == 1  # ✓ Only one key
assert d[u.N] == "kg*m/s^2"  # ✓ Second assignment overwrote first

# Test 5: UnrecognizedUnit works
foo1 = u.Unit("FOO", parse_strict="silent")
foo2 = u.Unit("FOO", parse_strict="silent")
assert foo1 == foo2  # ✓ Equal
assert hash(foo1) == hash(foo2)  # ✓ Same hash
```

All tests pass! ✓

## Impact

### Benefits

1. **Correctness**: Fixes fundamental violation of Python's hash-equality contract
2. **Sets work properly**: Equivalent units are now correctly treated as duplicates
3. **Dictionaries work properly**: Can use any equivalent unit form as a dictionary key
4. **No breaking changes**: The fix is transparent to existing code
5. **Performance**: Caching ensures minimal performance overhead

### Compatibility

- **Backward compatible**: Existing code continues to work
- **No API changes**: All public interfaces remain the same
- **Improved correctness**: Code that previously had subtle bugs will now work correctly

### Edge Cases Handled

1. **UnrecognizedUnit**: Special hash implementation based on unit name
2. **Pickling**: Hash cache is excluded from pickled state to avoid stale hashes
3. **Complex composite units**: All arithmetic combinations of units work correctly
4. **Scaled units**: Units with different scales hash differently (as expected)

## Conclusion

This fix resolves a critical consistency issue in Astropy's unit system that violated Python's fundamental data model requirements. By ensuring that equal units always have equal hashes, the fix enables proper use of units in sets, dictionaries, and other hash-based data structures.

The implementation is efficient (using caching), correct (thoroughly tested), and maintainable (well-documented). It represents a significant improvement in the reliability and correctness of the Astropy units system.
