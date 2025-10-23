# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io

from astropy import units as u
from astropy.table import Table


def test_mrt_like_units_parsing():
    # MRT-like snippet with CDS units in header
    content = """
Title: Example MRT with CDS units
Authors: Example
Table: Example
================================================================================
Byte-by-byte Description of file: datafile.txt
--------------------------------------------------------------------------------
 Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
 1-  5 F5.1   10+3J/m/s/kpc2 SBCONT  Surface brightness continuum
 7- 11 F5.1   10-7J/s/kpc2   SBLINE  Surface brightness line
--------------------------------------------------------------------------------
  1.0   2.0
""".strip()

    fh = io.StringIO(content)
    t = Table.read(fh, format="ascii.cds")

    assert t["SBCONT"].unit == (1e3 * u.J / u.m / u.s / (u.kpc**2))
    assert t["SBLINE"].unit == (1e-7 * u.J / u.s / (u.kpc**2))
