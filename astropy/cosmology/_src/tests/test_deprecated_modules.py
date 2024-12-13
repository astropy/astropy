"""Testing deprecated modules in `astropy.cosmology`."""

import re
import warnings

import pytest


def test_funcs():
    """Test `astropy.cosmology.core`."""
    from astropy.cosmology import funcs

    try:
        del funcs.z_at_value
    except Exception:
        pass

    with (
        warnings.catch_warnings(),
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.funcs` is deprecated")
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.funcs import z_at_value  # noqa: F401


def test_parameter():
    """Test `astropy.cosmology.parameter`."""
    from astropy.cosmology import parameter

    try:
        del parameter.Parameter
    except Exception:
        pass

    with (
        warnings.catch_warnings(),
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.parameter` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.parameter import Parameter  # noqa: F401
