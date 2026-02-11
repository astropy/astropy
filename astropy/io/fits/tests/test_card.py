# Licensed under a 3-clause BSD style license - see PYFITS.rst

import pytest

from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning


def test_fits_strict_string_valid():
    """
    Strict FITS string with doubled quotes must parse cleanly.
    """
    card = fits.Card.fromstring("TEST    = 'a '' b'")
    assert card.value == "a ' b"


def test_fits_malformed_string_warns_and_falls_back():
    """
    Malformed FITS string:
    - strict parsing fails
    - tolerant parsing is used
    - warning is emitted when value is accessed
    """
    card = fits.Card.fromstring("TEST    = 'a ' b ' /c'")
    with pytest.warns(VerifyWarning, match="Non-standard FITS string detected"):
        val = card.value
    assert val == "a ' b"


def test_fits_invalid_string_does_not_raise():
    """
    FITS is tolerant: even badly malformed strings must not raise,
    but must emit a VerifyWarning when parsed.
    """
    card = fits.Card.fromstring("TEST    = 'a ' b '")
    with pytest.warns(VerifyWarning, match="Non-standard FITS string detected"):
        val = card.value
    assert isinstance(val, str)
