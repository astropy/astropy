import sys
import subprocess
import pytest
from pathlib import Path

# Specifies the binary path, which is ../.tmp/
BINARY_NAME = "test_convolve_1d_binary"
if sys.platform.startswith("win"):
    BINARY_NAME +=".exe"
BINARY_PATH = Path(__file__).parent / ".tmp" / BINARY_NAME

@pytest.mark.skipif(not BINARY_PATH.exists(), 
                    reason=f"Compiled C test binary not found at {BINARY_PATH}.")
def test_convolution_c_results():
    """Executes the standalone C test binary via subprocess."""
    
    result = subprocess.run(
        [BINARY_PATH],
        capture_output=True,
        text=True
    )
    
    error_msg = (f"C Test Suite Failed.\nSTDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}")
    assert result.returncode == 0, error_msg