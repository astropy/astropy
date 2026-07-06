import subprocess
import sys
from pathlib import Path

def run_git_grep(pattern, exclude_patterns=None):
    """Run git grep to find patterns in the codebase."""
    cmd = ['git', 'grep', '-n', pattern, '--', '*.py']
    
    result = subprocess.run(
        cmd,
        cwd='/home/thiago-d-coqueiro/Documents/Projects/Astropy/astropy',
        capture_output=True,
        text=True
    )
    
    lines = result.stdout.strip().split('\n') if result.stdout.strip() else []
    
    # Filter out excluded patterns
    if exclude_patterns:
        filtered = []
        for line in lines:
            if not any(exc in line for exc in exclude_patterns):
                filtered.append(line)
        lines = filtered
    
    return lines

def main():
    print("=" * 80)
    print("VERIFICATION: np.rollaxis() → np.moveaxis() Migration")
    print("=" * 80)
    print()
    
    # Check for actual rollaxis() function calls (not in comments or function lists)
    print("1. Searching for direct np.rollaxis() function calls...")
    print("-" * 80)
    
    rollaxis_calls = run_git_grep(
        r'np\.rollaxis\(',
        exclude_patterns=['#', 'helps={', 'SUBCLASS_SAFE_FUNCTIONS']
    )
    
    if rollaxis_calls:
        print("[FAIL] FOUND DIRECT np.rollaxis() CALLS (need migration):")
        for call in rollaxis_calls:
            print(f"  {call}")
        print()
    else:
        print("[PASS] No direct np.rollaxis() function calls found in code")
        print()
    
    # Check for moveaxis usage
    print("2. Verifying np.moveaxis() is being used...")
    print("-" * 80)
    
    moveaxis_calls = run_git_grep(r'np\.moveaxis\(')
    print(f"[INFO] Found {len(moveaxis_calls)} occurrences of np.moveaxis()")
    print()
    
    # Check for remaining rollaxis in helper functions (expected)
    print("3. Checking for rollaxis in helper/compatibility lists (expected)...")
    print("-" * 80)
    
    helper_rollaxis = run_git_grep('rollaxis')
    expected_files = {
        'astropy/units/quantity_helper/function_helpers.py',
        'astropy/utils/masked/function_helpers.py',
        'astropy/utils/shapes.py'
    }
    
    print("Expected locations (for backward compatibility):")
    for line in helper_rollaxis:
        if any(expected in line for expected in expected_files):
            print(f"{line.split(':')[0].split('/')[-1]}: {line.split(':')[1]}")
    print()
    
    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    if not rollaxis_calls:
        print("[PASS] Migration appears to be COMPLETE")
        print("   - All np.rollaxis() calls have been replaced with np.moveaxis()")
        print("   - Helper/compatibility lists retain rollaxis for backward compatibility")
        print("   - Ready for testing")
        return 0
    else:
        print("[FAIL] Migration INCOMPLETE")
        print(f"   - Found {len(rollaxis_calls)} remaining np.rollaxis() calls")
        return 1

if __name__ == '__main__':
    sys.exit(main())
