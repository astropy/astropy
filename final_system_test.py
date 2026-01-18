#!/usr/bin/env python3
"""
FINAL COMPREHENSIVE SYSTEM TEST
NumPy 2.5 Shape Mutation Fixes - Complete Validation
"""
import warnings
import numpy as np

def test_imports():
    """Test all fixed modules import without warnings"""
    print("=== IMPORT TESTING ===")
    
    modules = [
        'astropy.coordinates.solar_system',
        'astropy.coordinates.representation.base',
        'astropy.time.core', 
        'astropy.io.ascii.ecsv',
        'astropy.io.fits.util',
        'astropy.io.fits.fitsrec',
        'astropy.uncertainty.core'
    ]
    
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        
        success = 0
        for module in modules:
            try:
                __import__(module)
                success += 1
                print(f"✓ {module}")
            except Exception as e:
                print(f"✗ {module}: {e}")
        
        shape_warnings = [warning for warning in w if 'shape' in str(warning.message).lower()]
        
    return success, len(modules), len(shape_warnings)

def test_functionality():
    """Test actual functionality works"""
    print("\n=== FUNCTIONALITY TESTING ===")
    
    tests_passed = 0
    total_tests = 0
    
    # Test 1: Time objects
    total_tests += 1
    try:
        from astropy.time import Time
        t = Time(['2000-01-01', '2000-01-02'])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            t_reshaped = t.reshape(2, 1)
            shape_warnings = [warning for warning in w if 'shape' in str(warning.message).lower()]
            
        if len(shape_warnings) == 0:
            print("✓ Time reshape functionality")
            tests_passed += 1
        else:
            print(f"✗ Time reshape: {len(shape_warnings)} warnings")
    except Exception as e:
        print(f"✗ Time test failed: {e}")
    
    # Test 2: Coordinates
    total_tests += 1
    try:
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        coords = SkyCoord([1, 2]*u.deg, [3, 4]*u.deg)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            coords_reshaped = coords.reshape(2, 1)
            shape_warnings = [warning for warning in w if 'shape' in str(warning.message).lower()]
            
        if len(shape_warnings) == 0:
            print("✓ Coordinates reshape functionality")
            tests_passed += 1
        else:
            print(f"✗ Coordinates reshape: {len(shape_warnings)} warnings")
    except Exception as e:
        print(f"✗ Coordinates test failed: {e}")
    
    # Test 3: Shape setting on representations
    total_tests += 1
    try:
        from astropy.coordinates.representation import SphericalRepresentation
        from astropy import units as u
        
        s = SphericalRepresentation(
            lon=[0, 1, 2, 3, 4, 5]*u.deg,
            lat=[0, 0, 0, 0, 0, 0]*u.deg, 
            distance=[1, 1, 1, 1, 1, 1]*u.kpc
        )
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            s.shape = (2, 3)  # Uses our fixed shape setter
            shape_warnings = [warning for warning in w if 'shape' in str(warning.message).lower()]
            
        if len(shape_warnings) == 0 and s.shape == (2, 3):
            print("✓ Representation shape setting")
            tests_passed += 1
        else:
            print(f"✗ Representation shape setting: {len(shape_warnings)} warnings")
    except Exception as e:
        print(f"✗ Representation test failed: {e}")
    
    return tests_passed, total_tests

def main():
    """Run complete system validation"""
    print("=" * 60)
    print("FINAL COMPREHENSIVE SYSTEM TEST")
    print("NumPy 2.5 Shape Mutation Fixes")
    print("=" * 60)
    print(f"NumPy version: {np.__version__}")
    
    # Test imports
    import_success, import_total, shape_warnings = test_imports()
    
    # Test functionality  
    func_success, func_total = test_functionality()
    
    # Summary
    print("\n" + "=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    
    print(f"Import Tests: {import_success}/{import_total} passed")
    print(f"Functionality Tests: {func_success}/{func_total} passed") 
    print(f"Shape Warnings: {shape_warnings}")
    
    if (import_success == import_total and 
        func_success == func_total and 
        shape_warnings == 0):
        print("\n🎉 ALL TESTS PASSED - FIXES ARE COMPLETE!")
        print("✅ No NumPy shape mutation warnings")
        print("✅ All functionality working correctly")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
