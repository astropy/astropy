#!/usr/bin/env python
"""
Demonstration of the fix for Quantity.__array_ufunc__() duck array interoperability.

This script shows how the fix enables duck-typed arrays to work properly with
Quantity when units need conversion.
"""

import dataclasses
import numpy as np
import astropy.units as u


@dataclasses.dataclass
class DuckArray(np.lib.mixins.NDArrayOperatorsMixin):
    """
    A minimal duck-typed array that wraps a Quantity.

    This is similar to the example from the issue:
    https://github.com/Kankelborg-Group/named_arrays
    """
    ndarray: u.Quantity

    @property
    def unit(self) -> u.UnitBase:
        return self.ndarray.unit

    def __array_ufunc__(self, function, method, *inputs, **kwargs):
        # Convert DuckArray inputs to their underlying Quantity
        inputs = [inp.ndarray if isinstance(inp, DuckArray) else inp for inp in inputs]

        # Try to delegate to the first ndarray-like input
        for inp in inputs:
            if isinstance(inp, np.ndarray):
                result = inp.__array_ufunc__(function, method, *inputs, **kwargs)
                if result is not NotImplemented:
                    return DuckArray(result)

        return NotImplemented


def main():
    print("=" * 70)
    print("Duck Array Interoperability Demo")
    print("=" * 70)

    # Example 1: Same units - this always worked
    print("\n1. Same units (always worked):")
    print("   DuckArray(1*mm) + (1*mm)")
    try:
        result = DuckArray(1 * u.mm) + (1 * u.mm)
        print(f"   ✓ Result: {result.ndarray}")
    except Exception as e:
        print(f"   ✗ Error: {type(e).__name__}: {e}")

    # Example 2: Different but equivalent units - the problematic case
    print("\n2. Different but equivalent units (the fix!):")
    print("   (1*m) + DuckArray(1*mm)")
    try:
        result = (1 * u.m) + DuckArray(1 * u.mm)
        print(f"   ✓ Result: {result.ndarray}")
        print(f"   ✓ Type: {type(result).__name__}")
    except Exception as e:
        print(f"   ✗ Error: {type(e).__name__}: {e}")

    # Example 3: Reverse order
    print("\n3. Reverse order:")
    print("   DuckArray(1*mm) + (1*m)")
    try:
        result = DuckArray(1 * u.mm) + (1 * u.m)
        print(f"   ✓ Result: {result.ndarray}")
        print(f"   ✓ Type: {type(result).__name__}")
    except Exception as e:
        print(f"   ✗ Error: {type(e).__name__}: {e}")

    # Example 4: Incompatible units - should still fail
    print("\n4. Incompatible units (should fail):")
    print("   (1*m) + DuckArray(1*s)")
    try:
        result = (1 * u.m) + DuckArray(1 * u.s)
        print(f"   ✗ Unexpected success: {result.ndarray}")
    except u.UnitConversionError as e:
        print(f"   ✓ Expected error: {type(e).__name__}")
    except Exception as e:
        print(f"   ? Unexpected error: {type(e).__name__}: {e}")

    # Example 5: Other operations
    print("\n5. Multiplication:")
    print("   (2*m) * DuckArray(3*s)")
    try:
        result = (2 * u.m) * DuckArray(3 * u.s)
        print(f"   ✓ Result: {result.ndarray}")
        print(f"   ✓ Type: {type(result).__name__}")
    except Exception as e:
        print(f"   ✗ Error: {type(e).__name__}: {e}")

    print("\n" + "=" * 70)
    print("Demo complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
