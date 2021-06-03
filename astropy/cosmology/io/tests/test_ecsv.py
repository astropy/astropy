# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

import pytest

from astropy import cosmology
from astropy.cosmology import Cosmology
from astropy.cosmology.io.ecsv import read_ecsv, write_ecsv, ecsv_identify
from astropy.table import QTable, vstack
from astropy.utils.compat import optional_deps
from astropy.utils.exceptions import AstropyUserWarning

cosmo_instances = [
    getattr(cosmology.realizations, name) for name in cosmology.parameters.available
]

# TODO! remove in astropy v5.0
if not getattr(optional_deps, "HAS_YAML"):
    save_formats.remove("ascii.ecsv")


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
    assert got.name == cosmo.name
    # assert got == expected  # FIXME! no __eq__ on cosmo


@pytest.mark.parametrize("cosmo", cosmo_instances)
def test_ND(cosmo_dir, cosmo):
    """Read tests happen later."""
    t = vstack([cosmo.write.to_table(), cosmo.write.to_table()])
    t["name"][1] = "Other"

    fname = cosmo_dir / f"{cosmo.name}_ND.ecsv"
    t.write(fname, overwrite=True)

    # error for no index
    with pytest.raises(ValueError, match="row index for N-D"):
        read_ecsv(fname)

    got = read_ecsv(fname, index=0)
    assert got.name == cosmo.name
    # assert got == cosmo  # FIXME! no __eq__ on cosmo

    got = read_ecsv(fname, index=1)
    assert got.name == "Other"
    # assert got == cosmo  # FIXME! no __eq__ on cosmo

    # Now test index is string
    # the table is not indexed, so this also tests adding an index column
    got = read_ecsv(fname, index="Other")
    assert got.name == "Other"
    # assert got == cosmo  # FIXME! no __eq__ on cosmo
