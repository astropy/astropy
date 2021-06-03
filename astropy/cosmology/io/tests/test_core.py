# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.cosmology import io
from astropy.table import QTable, vstack
from astropy.utils.compat import optional_deps
from astropy.utils.exceptions import AstropyUserWarning


class CosmologyWithKwargs(Cosmology):
    def __init__(self, name="cosmology with kwargs", **kwargs):
        super().__init__(name=name, **kwargs)


cosmo_instances = [
    getattr(cosmology.realizations, name) for name in cosmology.parameters.available
]
cosmo_instances.append(CosmologyWithKwargs())

# -----------------------------------------------------------------------------


@pytest.mark.parametrize("expected", cosmo_instances)
def test_to_from_mapping_instance(expected):
    # ------------
    # To Mapping
    params = io.to_mapping(expected)

    assert isinstance(params, dict)
    assert params["cosmology"] is expected.__class__
    assert params["name"] == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo

    # ------------
    # From Mapping
    params["mismatching"] = "will error"

    # tests are different if the last argument is a **kwarg
    if tuple(expected._init_signature.parameters.values())[-1].kind == 4:
        got = io.from_mapping(params)

        assert got.__class__ == expected.__class__
        assert got.name == expected.name
        assert "mismatching" not in got.meta

        return  # don't continue teesting

    # read with mismatching parameters errors
    with pytest.raises(TypeError, match="There are unused parameters"):
        io.from_mapping(params)

    # unless mismatched are moved to meta
    got = io.from_mapping(params, move_to_meta=True)
    assert got.__class__ == expected.__class__
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo
    assert got.meta["mismatching"] == "will error"

    # it won't error if everything matches up
    params.pop("mismatching")
    got = io.from_mapping(params)
    assert got.__class__ == expected.__class__
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo

    # and it will also work if the cosmology is a string
    params["cosmology"] = params["cosmology"].__name__
    got = io.from_mapping(params)
    assert got.__class__ == expected.__class__
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo


@pytest.mark.parametrize("expected", cosmo_instances)
def test_to_from_table_instance(expected):
    # ------------
    # To Table
    t = io.to_table(expected)

    assert isinstance(t, QTable)
    assert t.meta["cosmology"] is expected.__class__.__name__
    assert t["name"] == expected.name

    # ------------
    # From Table
    # simple roundtrip
    got = io.from_table(t)  # index : None -> 0

    # Now test index is string
    # the table is not indexed, so this also tests adding an index column
    got = io.from_table(t, index=expected.name)
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo


@pytest.mark.parametrize("expected", cosmo_instances)
def test_ND_table(expected):
    t = vstack([expected.write.to_table(), expected.write.to_table()])
    t["name"][1] = "Other"

    # error for no index
    with pytest.raises(ValueError, match="row index for N-D"):
        io.from_table(t)

    got = io.from_table(t, index=0)
    assert got.name == expected.name
    # assert got == expected  # FIXME! no __eq__ on cosmo

    got = io.from_table(t, index=1)
    assert got.name == "Other"
    # assert got == expected  # FIXME! no __eq__ on cosmo

    # Now test index is string
    # the table is not indexed, so this also tests adding an index column
    got = io.from_table(t, index="Other")
    assert got.name == "Other"
    # assert got == expected  # FIXME! no __eq__ on cosmo
