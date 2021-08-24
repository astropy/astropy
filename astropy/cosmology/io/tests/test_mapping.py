# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.cosmology.core import _COSMOLOGY_CLASSES
from astropy.table import QTable, vstack
from astropy.utils.compat import optional_deps
from astropy.utils.exceptions import AstropyUserWarning


class CosmologyWithKwargs(Cosmology):
    def __init__(self, name="cosmology with kwargs", meta=None, **kwargs):
        super().__init__(name=name, meta=meta, **kwargs)


cosmo_instances = [
    getattr(cosmology.realizations, name) for name in cosmology.parameters.available
]
cosmo_instances.append(CosmologyWithKwargs())


def teardown_module(module):
    # pop CosmologyWithKwargs from registered classes
    # but don't error b/c it fails in parallel
    _COSMOLOGY_CLASSES.pop(CosmologyWithKwargs.__qualname__, None)


###############################################################################

@pytest.mark.parametrize("expected", cosmo_instances)
def test_to_from_mapping_instance(expected):
    # ------------
    # To Mapping
    params = expected.to_format('mapping')

    assert isinstance(params, dict)
    assert params["cosmology"] is expected.__class__
    assert params["name"] == expected.name

    # ------------
    # From Mapping
    params["mismatching"] = "will error"

    # tests are different if the last argument is a **kwarg
    if tuple(expected._init_signature.parameters.values())[-1].kind == 4:
        got = Cosmology.from_format(params, format="mapping")

        assert got.__class__ == expected.__class__
        assert got.name == expected.name
        assert "mismatching" not in got.meta

        return  # don't continue testing

    # read with mismatching parameters errors
    with pytest.raises(TypeError, match="there are unused parameters"):
        Cosmology.from_format(params, format="mapping")

    # unless mismatched are moved to meta
    got = Cosmology.from_format(params, format="mapping", move_to_meta=True)
    assert got.__class__ == expected.__class__
    assert got == expected
    assert got.meta["mismatching"] == "will error"

    # it won't error if everything matches up
    params.pop("mismatching")
    got = Cosmology.from_format(params, format="mapping")
    assert got.__class__ == expected.__class__
    assert got == expected

    # and it will also work if the cosmology is a string
    params["cosmology"] = params["cosmology"].__name__
    got = Cosmology.from_format(params, format="mapping")
    assert got == expected

    # also it auto-identifies 'format'
    got = Cosmology.from_format(params)
    assert got == expected
