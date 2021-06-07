# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import pytest

from astropy import cosmology
from astropy.cosmology.io.ecsv import read_ecsv, write_ecsv, ecsv_identify
from astropy.table import QTable, vstack
from astropy.utils.compat.optional_deps import HAS_YAML


pytestmark = pytest.mark.skipif(not HAS_YAML, reason="Needs PyYAML")

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
    fname = cosmo_dir / f"{cosmo.name}.ecsv"
    write_ecsv(cosmo, fname)

    # Also test kwarg "overwrite"
    assert os.path.exists(fname)  # file exists
    with pytest.raises(IOError):
        write_ecsv(cosmo, fname, overwrite=False)

    assert os.path.exists(fname)  # overwrite file existing file
    write_ecsv(cosmo, fname, overwrite=True)

    # Read back
    got = read_ecsv(fname)
    assert got.__class__ is cosmo.__class__
    assert got == cosmo
    assert got.meta == cosmo.meta  # == doesn't check meta


@pytest.mark.parametrize("cosmo", cosmo_instances)
def test_ND(cosmo_dir, cosmo):
    """Read tests happen later."""
    cosmo1 = cosmo.clone(name="Other")
    t = vstack([cosmo.write.to_table(), cosmo1.write.to_table()])

    fname = cosmo_dir / f"{cosmo.name}_ND.ecsv"
    t.write(fname, overwrite=True)

    # error for no index
    with pytest.raises(ValueError, match="row index for N-D"):
        read_ecsv(fname)

    got = read_ecsv(fname, index=0)
    assert got.name == cosmo.name
    assert got == cosmo
    assert got.meta == cosmo.meta  # == doesn't check meta

    got = read_ecsv(fname, index=1)
    assert got.name == "Other"
    assert got == cosmo1
    assert got.meta == cosmo1.meta  # == doesn't check meta

    # Now test index is string
    # the table is not indexed, so this also tests adding an index column
    got = read_ecsv(fname, index="Other")
    assert got.name == "Other"
    assert got == cosmo1
    assert got.meta == cosmo1.meta  # == doesn't check meta
