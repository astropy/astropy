# Licensed under a 3-clause BSD style license - see LICENSE.rst

# STDLIB
from types import MappingProxyType

# THIRD PARTY
import numpy as np
import pytest

# LOCAL
from astropy.cosmology import parameters, realizations


def test_realizations_in_dir():
    """Test the realizations are in ``dir`` of :mod:`astropy.cosmology.parameters`."""
    d = dir(parameters)

    assert set(d) == set(parameters.__all__)
    for n in parameters.available:
        assert n in d


@pytest.mark.parametrize("name", parameters.available)
def test_getting_parameters(name):
    """
    Test getting 'parameters' and that it is derived from the corresponding
    realization.
    """
    params = getattr(parameters, name)

    assert isinstance(params, MappingProxyType)
    assert params["name"] == name

    # Check parameters have the right keys and values
    cosmo = getattr(realizations, name)
    assert params["name"] == cosmo.name
    assert params["cosmology"] == cosmo.__class__.__qualname__
    # All the cosmology parameters are equal
    for n in cosmo.__parameters__:
        assert np.array_equal(params[n], getattr(cosmo, n))
    # All the metadata is included. Parameter values take precedence, so only
    # checking the keys.
    assert set(cosmo.meta.keys()).issubset(params.keys())

    # Lastly, check the generation process.
    m = cosmo.to_format("mapping", cosmology_as_str=True, move_from_meta=True)
    assert params == m
