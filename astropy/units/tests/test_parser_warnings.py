# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Tests for detailed unit parser warning messages."""

from astropy.units import OGIPInvalidMultiplicationWarning


def test_OGIPInvalidMultiplicationWarning() -> None:
    assert str(OGIPInvalidMultiplicationWarning("m", "(s)")) == (
        "if 'm(s)' was meant to be a multiplication, it should have "
        "been written as 'm (s)'."
    )
