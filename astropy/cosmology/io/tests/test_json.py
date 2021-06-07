# Licensed under a 3-clause BSD style license - see LICENSE.rst

import json
import os

import pytest

from astropy import cosmology
from astropy.cosmology.io.json import read_json, write_json, json_identify
from astropy.utils.exceptions import AstropyUserWarning

cosmo_instances = [
    getattr(cosmology.realizations, name) for name in cosmology.parameters.available
]


# make a common directory for reading / writing cosmologies
@pytest.fixture(scope="session")
def cosmo_dir(tmpdir_factory):
    drct = tmpdir_factory.mktemp("cosmo")
    return drct


# -----------------------------------------------------------------------------


@pytest.mark.parametrize("cosmo", cosmo_instances)
def test_write_then_read_file(cosmo_dir, cosmo):
    """Read tests happen later."""
    fname = cosmo_dir / f"{cosmo.name}.json"
    write_json(cosmo, fname)

    # Also test kwarg "overwrite"
    assert os.path.exists(fname)  # file exists
    with pytest.raises(IOError):
        write_json(cosmo, fname, overwrite=False)

    assert os.path.exists(fname)  # overwrite file existing file
    write_json(cosmo, fname, overwrite=True)

    # Read back
    got = read_json(fname)
    assert got.__class__ is cosmo.__class__
    assert got.name == cosmo.name
    # assert got == expected  # FIXME! no __eq__ on cosmo


@pytest.mark.parametrize("cosmo", cosmo_instances)
def test_ND(cosmo_dir, cosmo):
    """Read tests happen later."""
    m = cosmo.write.to_mapping()
    m["cosmology"] = m["cosmology"].__name__

    print("mapping = ", m)

    data = {cosmo.name: m, "Other": m}

    fname = cosmo_dir / f"{cosmo.name}_ND.json"
    with open(fname, "w") as write_file:
        json.dump(data, write_file)

    # error for no index
    with pytest.raises(KeyError):
        read_json(fname)

    got = read_json(fname, key=cosmo.name)
    assert got.__class__ == cosmo.__class__
    assert got.name == cosmo.name
    # assert got == cosmo  # FIXME! no __eq__ on cosmo

    got = read_json(fname, key="Other")
    assert got.__class__ == cosmo.__class__
    assert got.name == cosmo.name
    # assert got == cosmo  # FIXME! no __eq__ on cosmo
