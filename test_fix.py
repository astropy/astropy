#!/usr/bin/env python
"""Simple test to verify the fix for compound model fit_deriv"""

import numpy as np

# Test basic imports
print("Testing imports...")
try:
    from astropy.modeling import models
    from astropy.modeling.fitting import LinearLSQFitter
    print("✓ Imports successful")
except Exception as e:
    print(f"✗ Import failed: {e}")
    exit(1)

# Test 1: Mapping | Chebyshev1D (no left params)
print("\nTest 1: Mapping | Chebyshev1D")
try:
    fitter = LinearLSQFitter()
    model = models.Mapping((0,)) | models.Chebyshev1D(2)
    x = np.arange(10, dtype=float)
    y = 2 * x + 3

    result = fitter(model, x, y)
    print(f"✓ Test 1 PASSED: fit completed")
    print(f"  Prediction close to y: {np.allclose(result(x), y, atol=0.1)}")
except Exception as e:
    print(f"✗ Test 1 FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Chebyshev1D | Chebyshev1D (both have params)
print("\nTest 2: Chebyshev1D | Chebyshev1D")
try:
    fitter = LinearLSQFitter()
    model = models.Chebyshev1D(2) | models.Chebyshev1D(2)
    x = np.arange(10, dtype=float)
    y = np.arange(10, dtype=float)

    result = fitter(model, x, y)
    print(f"✓ Test 2 PASSED: fit completed")
except Exception as e:
    print(f"✗ Test 2 FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Identity | Chebyshev2D (2D orthopolynomial in compound)
print("\nTest 3: Identity | Chebyshev2D")
try:
    fitter = LinearLSQFitter()
    y_coord, x_coord = np.mgrid[0:5, 0:5]
    z = x_coord + y_coord

    compound_model = models.Identity(2) | models.Chebyshev2D(2, 2)
    compound_model.fittable = True
    compound_model.linear = True

    compound_fit = fitter(compound_model, x_coord, y_coord, z)
    print(f"✓ Test 3 PASSED: fit completed")
    print(f"  Prediction close to z: {np.allclose(compound_fit(x_coord, y_coord), z, atol=0.1)}")
except Exception as e:
    print(f"✗ Test 3 FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\nAll tests completed!")
