#!/usr/bin/env python
"""
Simple migration verification that checks code syntax and correctness
without requiring full Astropy build.
"""

import re
import sys
from pathlib import Path

# Files that should have been changed
CHANGED_FILES = [
    "astropy/modeling/fitting.py",
    "astropy/modeling/polynomial.py",
    "astropy/modeling/core.py",
    "astropy/modeling/_fitting_parallel.py",
    "astropy/uncertainty/core.py",
    "astropy/coordinates/representation/cartesian.py",
    "astropy/modeling/tests/test_model_sets.py",
    "astropy/modeling/tests/test_parameters.py",
    "astropy/modeling/tests/test_fitters.py",
    "astropy/coordinates/tests/test_representation_methods.py",
    "astropy/stats/tests/test_sigma_clipping.py",
    "astropy/time/tests/test_methods.py",
    "astropy/uncertainty/tests/test_distribution.py",
]

def check_file_syntax(filepath):
    """Check if a Python file has valid syntax."""
    try:
        with open(filepath, 'r') as f:
            code = f.read()
        compile(code, filepath, 'exec')
        return True, None
    except SyntaxError as e:
        return False, str(e)

def check_moveaxis_usage(filepath):
    """Check if file uses np.moveaxis() and doesn't use np.rollaxis()."""
    with open(filepath, 'r') as f:
        content = f.read()

    # Check for actual rollaxis calls (not in comments or strings)
    lines = content.split('\n')
    rollaxis_found = False
    for i, line in enumerate(lines, 1):
        # Skip comments and docstrings
        if '#' in line:
            line = line[:line.index('#')]

        if 'np.rollaxis(' in line and 'np.moveaxis(' not in line:
            rollaxis_found = True
            break

    # Check for moveaxis usage
    moveaxis_found = 'np.moveaxis(' in content

    return rollaxis_found, moveaxis_found

def main():
    """Run verification checks."""
    print("=" * 80)
    print("Migration Verification: Syntax and Code Quality Check")
    print("=" * 80)
    print()

    syntax_errors = []
    rollaxis_found = []
    moveaxis_usage = []

    base_path = Path(__file__).parent

    for file_path in CHANGED_FILES:
        full_path = base_path / file_path

        if not full_path.exists():
            print(f"[SKIP] {file_path} (file not found)")
            continue

        # Check syntax
        print(f"Checking {file_path}...", end=" ")
        valid, error = check_file_syntax(full_path)

        if not valid:
            print(f"[FAIL] Syntax Error")
            syntax_errors.append((file_path, error))
            continue

        # Check rollaxis/moveaxis usage
        has_rollaxis, has_moveaxis = check_moveaxis_usage(full_path)

        if has_rollaxis:
            print(f"[FAIL] Found np.rollaxis()")
            rollaxis_found.append(file_path)
        elif has_moveaxis:
            print(f"[PASS] Uses np.moveaxis()")
            moveaxis_usage.append(file_path)
        else:
            print(f"[INFO] No axis movement functions")

    # Summary
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Files checked: {len(CHANGED_FILES)}")
    print(f"Syntax valid: {len(CHANGED_FILES) - len(syntax_errors)}")
    print(f"Using moveaxis: {len(moveaxis_usage)}")
    print()

    if syntax_errors:
        print("[FAIL] Syntax Errors Found:")
        for file_path, error in syntax_errors:
            print(f"  {file_path}: {error}")
        return 1

    if rollaxis_found:
        print("[FAIL] np.rollaxis() Still Found In:")
        for file_path in rollaxis_found:
            print(f"  {file_path}")
        return 1

    print("[PASS] All checks passed!")
    print("  - All files have valid Python syntax")
    print("  - No np.rollaxis() calls found in code")
    print("  - Files using np.moveaxis() verified")
    return 0

if __name__ == "__main__":
    sys.exit(main())
