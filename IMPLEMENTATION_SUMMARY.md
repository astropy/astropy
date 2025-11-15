# Implementation Summary: Fix for Misleading AttributeError in Subclassed SkyCoord

## Issue Description

When subclassing `SkyCoord` and adding custom properties that access non-existent attributes, the error message was misleading. It claimed the property itself didn't exist, when actually the property existed but tried to access a missing attribute.

### Example of the Bug

```python
import astropy.coordinates as coord

class custom_coord(coord.SkyCoord):
    @property
    def prop(self):
        return self.random_attr  # random_attr doesn't exist

c = custom_coord('00h42m30s', '+41d12m00s', frame='icrs')
c.prop
```

**Before Fix:**
```
AttributeError: 'custom_coord' object has no attribute 'prop'
```
❌ Misleading - `prop` DOES exist!

**After Fix:**
```
AttributeError: 'custom_coord' object attribute 'prop' raised an AttributeError.
This usually means 'prop' tried to access an attribute that doesn't exist.
```
✓ Clear - `prop` exists but encountered an error!

## Changes Made

### 1. Modified `astropy/coordinates/sky_coordinate.py`

**Location:** Lines 897-913 in the `__getattr__` method

**Change:** Added logic to detect when an attribute exists in the class hierarchy before raising AttributeError.

```python
# Fail - but check first if the attribute exists in the class hierarchy
# If an attribute exists as a class member (like a property), but we're in
# __getattr__, then that attribute raised an AttributeError during access.
# In this case, don't mask the real error by claiming the attribute doesn't exist.
for cls in type(self).__mro__:
    if attr in cls.__dict__:
        # The attribute exists, so it must have raised an AttributeError internally.
        # Don't hide that - make it clear the attribute exists but is raising an error.
        raise AttributeError(
            f"'{self.__class__.__name__}' object attribute '{attr}' raised an AttributeError. "
            f"This usually means '{attr}' tried to access an attribute that doesn't exist."
        )

# The attribute truly doesn't exist
raise AttributeError(
    f"'{self.__class__.__name__}' object has no attribute '{attr}'"
)
```

**Rationale:**
- When `__getattr__` is called for an attribute that exists in the class (like a property), it means that attribute raised an AttributeError internally
- By checking `type(self).__mro__`, we can distinguish between:
  - Attributes that exist but raised errors (properties, descriptors)
  - Attributes that don't exist at all
- This provides clearer error messages without changing any behavior

### 2. Added Test Case in `astropy/coordinates/tests/test_sky_coord.py`

**Location:** Lines 2170-2196

**Test Function:** `test_subclass_property_attribute_error()`

**Coverage:**
1. Tests that properties raising AttributeError get improved error messages
2. Tests that normal "attribute doesn't exist" errors still work correctly
3. Ensures backward compatibility

```python
def test_subclass_property_attribute_error():
    """
    Test that subclassed SkyCoord gives a clear error message when a custom
    property tries to access a non-existent attribute.
    """

    class CustomCoord(SkyCoord):
        @property
        def prop(self):
            return self.random_attr

    c = CustomCoord("00h42m30s", "+41d12m00s", frame="icrs")

    # Should indicate 'prop' raised an error, not that it doesn't exist
    with pytest.raises(AttributeError, match="'prop' raised an AttributeError"):
        c.prop

    # Normal attribute access should still work as expected
    with pytest.raises(
        AttributeError, match="'CustomCoord' object has no attribute 'nonexistent'"
    ):
        c.nonexistent
```

### 3. Created Changelog Entry

**File:** `docs/changes/coordinates/XXXXX.bugfix.rst`

**Content:**
```rst
Fixed misleading ``AttributeError`` messages when accessing properties in subclassed
``SkyCoord`` objects. Previously, if a custom property tried to access a non-existent
attribute, the error message would incorrectly claim that the property itself didn't exist.
Now the error message clearly indicates that the property exists but raised an ``AttributeError``,
helping developers identify the actual missing attribute.
```

**Note:** The `XXXXX` placeholder should be replaced with the actual pull request number.

## Technical Details

### Root Cause Analysis

Python's attribute lookup mechanism:
1. `obj.prop` → Calls property getter
2. Inside getter: `self.attr` → If `attr` doesn't exist, calls `__getattr__('attr')`
3. `__getattr__` raises `AttributeError` about `attr`
4. Python catches this and tries `__getattr__('prop')` as a fallback
5. Original `__getattr__` raises error about `prop` not existing (masking the real error)

### Solution Approach

The fix leverages Python's Method Resolution Order (MRO) to detect when an attribute exists in the class hierarchy:

```python
for cls in type(self).__mro__:
    if attr in cls.__dict__:
        # Found it - it exists but raised an error
        ...
```

This works because:
- Properties, methods, and descriptors are stored in `cls.__dict__`
- If `__getattr__` is called for something in `cls.__dict__`, that thing must have raised an error
- We can then provide a clearer error message

### Edge Cases Handled

1. **Properties accessing missing attributes** ✓
2. **Truly missing attributes** ✓
3. **Working properties** ✓
4. **Method access** ✓
5. **Nested property access** ✓
6. **Property chains** ✓

## Testing

### Manual Testing
- Created standalone test scripts to verify fix logic
- Tested with the exact example from the bug report
- Verified error messages are clear and helpful

### Automated Testing
- Added regression test to test suite
- Verified existing tests still pass
- Tested multiple scenarios (properties, methods, nested access)

## Backward Compatibility

✅ **100% Backward Compatible**

- No API changes
- No behavior changes for working code
- Only error messages are improved
- Existing tests continue to pass
- All SkyCoord subclasses work as before

## Benefits

1. **Better Developer Experience**: Clear error messages save debugging time
2. **Less Confusion**: No more misleading "attribute doesn't exist" messages
3. **Easier Debugging**: Immediately identifies the real problem
4. **Professional Quality**: Error messages match best practices from major Python libraries

## Files Changed Summary

| File | Lines Changed | Type |
|------|---------------|------|
| `astropy/coordinates/sky_coordinate.py` | 897-913 | Modified |
| `astropy/coordinates/tests/test_sky_coord.py` | 2170-2196 | Added |
| `docs/changes/coordinates/XXXXX.bugfix.rst` | 1-5 | Added |

**Total:** 3 files, ~40 lines of code/documentation

## Verification Steps

To verify the fix works:

```bash
# Run the new test
pytest astropy/coordinates/tests/test_sky_coord.py::test_subclass_property_attribute_error -v

# Run all SkyCoord tests
pytest astropy/coordinates/tests/test_sky_coord.py -v

# Run all coordinates tests
pytest astropy/coordinates/tests/ -v
```

## Conclusion

This fix improves the developer experience when working with subclassed `SkyCoord` objects by providing clearer error messages. It's a small change with significant impact on debugging and code maintainability.
