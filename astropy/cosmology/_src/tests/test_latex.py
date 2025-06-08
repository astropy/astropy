import os
import tempfile
import pytest

from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import read_latex, write_latex, read_latex_from_string


def _create_temp_latex(content: str) -> str:

    """
    Helper Function: write the given LaTeX content to a temporary .tex file,
    the return the filename so tests can read from it.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".tex", mode="w")
    tmp.write(content)
    tmp.close()
    return tmp.name


def test_read_latex():
    # Create a minimal LaTeX document with H0 and Om0 parameters
    sample = (
        "\\documentclass{article}\n"
        "\\begin{document}\n"
        "H0 = 65.0\n"
        "Om0 = 0.25\n"
        "\\end{document}\n"
    )
    path = _create_temp_latex(sample)
    cosmo = read_latex(path)
    os.remove(path)

    assert abs(cosmo.H0.value - 65.0) < 1e-8
    assert abs(cosmo.Om0 - 0.25) < 1e-8


def test_write_and_read_roundtrip():
    # Start from a known FlatLambdaCDM instance
    original = FlatLambdaCDM(H0=68.0, Om0=0.32)
    tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix=".tex")
    tmpfile.close()

    # Write it out, then read it back in
    write_latex(original, tmpfile.name)
    loaded = read_latex(tmpfile.name)
    os.remove(tmpfile.name)

    assert abs(loaded.H0.value - 68.0) < 1e-8
    assert abs(loaded.Om0 - 0.32) < 1e-8


def test_read_latex_from_string_function():
    # Test the convenience function that parses a raw string
    latex_str = "H0 = 72.0\nOm0 = 0.28\n"
    cosmo = read_latex_from_string(latex_str)

    assert abs(cosmo.H0.value - 72.0) < 1e-8
    assert abs(cosmo.Om0 - 0.28) < 1e-8


if __name__ == "__main__":
    pytest.main(["-q", __file__])
