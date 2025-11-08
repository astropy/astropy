"""Testing deprecated modules in `astropy.cosmology`."""

import re
import warnings

import pytest


def test_parameter():
    """Test `astropy.cosmology.parameter`."""
    from astropy.cosmology import parameter

    try:
        del parameter.Parameter
    except Exception:
        pass

    with (
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.parameter` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.parameter import Parameter  # noqa: F401
