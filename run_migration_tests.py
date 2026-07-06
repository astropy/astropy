#!/usr/bin/env python
"""
Test runner for np.rollaxis() -> np.moveaxis() migration verification.
Runs tests for all modules affected by the migration.
"""

import subprocess
import sys

# Test modules affected by the migration
TEST_MODULES = [
    "astropy/modeling/tests/test_model_sets.py",
    "astropy/modeling/tests/test_parameters.py",
    "astropy/modeling/tests/test_fitters.py",
    "astropy/coordinates/tests/test_representation_methods.py",
    "astropy/stats/tests/test_sigma_clipping.py",
    "astropy/time/tests/test_methods.py",
    "astropy/uncertainty/tests/test_distribution.py",
]

def run_tests():
    """Run pytest for all affected modules."""
    print("=" * 80)
    print("Running Tests for np.rollaxis() -> np.moveaxis() Migration")
    print("=" * 80)
    print()

    failed_tests = []
    passed_tests = []

    for module in TEST_MODULES:
        print(f"Testing {module}...")
        print("-" * 80)

        result = subprocess.run(
            [sys.executable, "-m", "pytest", module, "-v", "--tb=short"],
            capture_output=False
        )

        if result.returncode == 0:
            print(f"[PASS] {module}")
            passed_tests.append(module)
        else:
            print(f"[FAIL] {module}")
            failed_tests.append(module)

        print()

    # Summary
    print("=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    print(f"Passed: {len(passed_tests)}/{len(TEST_MODULES)}")
    print(f"Failed: {len(failed_tests)}/{len(TEST_MODULES)}")

    if passed_tests:
        print("\nPassed tests:")
        for test in passed_tests:
            print(f"  [PASS] {test}")

    if failed_tests:
        print("\nFailed tests:")
        for test in failed_tests:
            print(f"  [FAIL] {test}")
        return 1

    print("\n[PASS] All tests passed!")
    return 0

if __name__ == "__main__":
    sys.exit(run_tests())
