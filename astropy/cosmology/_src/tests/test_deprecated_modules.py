"""Testing deprecated modules in `astropy.cosmology`."""

import re
import warnings

import pytest


def test_connect():
    """Test `astropy.cosmology.connect`."""
    from astropy.cosmology import connect

    try:
        del connect.CosmologyFromFormat
    except Exception:
        pass

    with (
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.connect` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.connect import CosmologyFromFormat  # noqa: F401


def test_core():
    """Test `astropy.cosmology.core`."""
    from astropy.cosmology import core

    try:
        del core.Cosmology
    except Exception:
        pass

    with (
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.core` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.core import Cosmology  # noqa: F401


def test_flrw():
    """Test `astropy.cosmology.core`."""
    from astropy.cosmology import flrw

    try:
        del flrw.FLRW
    except Exception:
        pass

    with (
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.flrw` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.flrw import FLRW  # noqa: F401


def test_funcs():
    """Test `astropy.cosmology.core`."""
    from astropy.cosmology import funcs

    try:
        del funcs.z_at_value
    except Exception:
        pass

    with (
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
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
        warnings.catch_warnings(),  # Always raise warning so 2x test sees it too.
        pytest.deprecated_call(
            match=re.escape("The module `astropy.cosmology.parameter` is deprecated"),
        ),
    ):
        warnings.simplefilter("always")
        from astropy.cosmology.parameter import Parameter  # noqa: F401
