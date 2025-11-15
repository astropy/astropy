# Fix for Misleading AttributeError in Subclassed SkyCoord

## Problem

When subclassing `SkyCoord` and adding custom properties that try to access non-existent attributes, the error message was misleading. For example:

```python
import astropy.coordinates as coord

class custom_coord(coord.SkyCoord):
    @property
    def prop(self):
        return self.random_attr  # random_attr doesn't exist

c = custom_coord('00h42m30s', '+41d12m00s', frame='icrs')
c.prop  # Tries to access the property
```

**Before the fix**, this would raise:
```
AttributeError: 'custom_coord' object has no attribute 'prop'
```

This was misleading because `prop` **does** exist as a property; the real problem is that `random_attr` doesn't exist.

## Root Cause

The issue occurs due to Python's attribute lookup mechanism:

1. When `c.prop` is accessed, Python calls the property getter
2. Inside the getter, `self.random_attr` triggers `SkyCoord.__getattr__('random_attr')`
3. `__getattr__` raises `AttributeError` about `random_attr` not existing
4. Python's attribute mechanism catches this exception and thinks the property itself failed
5. Python then calls `__getattr__('prop')` as a fallback
6. The `__getattr__` method doesn't know that `prop` exists as a property, so it raises a new error claiming `prop` doesn't exist

This masks the real error about `random_attr`.

## Solution

The fix modifies the `__getattr__` method in `SkyCoord` to check if the requested attribute exists in the class hierarchy (e.g., as a property or method). If it does, the method knows that the attribute exists but raised an error internally, so it provides a clearer error message:

```python
def __getattr__(self, attr):
    # ... existing logic ...

    # Check if attribute exists in class hierarchy
    for cls in type(self).__mro__:
        if attr in cls.__dict__:
            # Attribute exists but raised an error
            raise AttributeError(
                f"'{self.__class__.__name__}' object attribute '{attr}' raised an AttributeError. "
                f"This usually means '{attr}' tried to access an attribute that doesn't exist."
            )

    # Attribute truly doesn't exist
    raise AttributeError(
        f"'{self.__class__.__name__}' object has no attribute '{attr}'"
    )
```

**After the fix**, the same code raises:
```
AttributeError: 'custom_coord' object attribute 'prop' raised an AttributeError.
This usually means 'prop' tried to access an attribute that doesn't exist.
```

This makes it much clearer that:
1. The property `prop` exists
2. The problem is inside the property implementation
3. The property likely tried to access a non-existent attribute

## Files Modified

1. **`astropy/coordinates/sky_coordinate.py`**
   - Modified the `__getattr__` method to detect when an attribute exists in the class hierarchy
   - Added logic to provide clearer error messages for properties that raise `AttributeError`

2. **`astropy/coordinates/tests/test_sky_coord.py`**
   - Added `test_subclass_property_attribute_error()` to verify the fix works correctly
   - Tests both the improved error message and that normal attribute errors still work

3. **`docs/changes/coordinates/XXXXX.bugfix.rst`**
   - Added changelog entry describing the bug fix

## Testing

The fix includes comprehensive tests:

1. **Test that properties raising AttributeError get clear messages**: Verifies that when a property tries to access a non-existent attribute, the error message indicates the property exists but raised an error.

2. **Test that normal AttributeError still works**: Verifies that accessing truly non-existent attributes still produces the expected "has no attribute" message.

3. **Test that working properties still function**: Verifies that properties that work correctly are not affected by the fix.

## Backward Compatibility

This change improves error messages without breaking any existing functionality:

- **No API changes**: The fix only modifies error messages, not any public APIs
- **Existing code continues to work**: All existing `SkyCoord` subclasses will continue to function identically
- **Better debugging experience**: Developers get clearer error messages when debugging custom subclasses

## Benefits

1. **Clearer debugging**: Developers can immediately identify that a property exists but is accessing a non-existent attribute
2. **Less confusion**: No more misleading messages claiming a property doesn't exist when it actually does
3. **Better developer experience**: Saves time when debugging custom `SkyCoord` subclasses
