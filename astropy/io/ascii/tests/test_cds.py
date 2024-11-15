# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``CDS`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""

from io import StringIO

from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename


def test_roundtrip_cds_table():
    """
    Tests whether or not the CDS writer can roundtrip a table,
    i.e. read a table to ``Table`` object and write it exactly
    as it is back to a file. Since, presently CDS uses a
    MRT format template while writing, only the Byte-By-Byte
    and the data section of the table can be compared between
    original and the newly written table.

    Further, the CDS Reader does not have capability to recognize
    column format from the header of a CDS/MRT table, so this test
    can work for a limited set of simple tables, which don't have
    whitespaces in the column values or mix-in columns. Because of
    this the written table output cannot be directly matched with
    the original file and have to be checked against a list of lines.
    Masked columns are read properly though, and thus are being tested
    during round-tripping.

    The difference between ``cdsFunctional2.dat`` file and ``exp_output``
    is the following:
        * Metadata is different because MRT template is used for writing.
        * Spacing between ``Label`` and ``Explanations`` column in the
            Byte-By-Byte.
        * Units are written as ``[cm.s-2]`` and not ``[cm/s2]``, since both
            are valid according to CDS/MRT standard.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 7  A7       ---    ID       Star ID                              ",
        " 9-12  I4       K      Teff     [4337/4654] Effective temperature    ",
        "14-17  F4.2   [cm.s-2] logg     [0.77/1.28] Surface gravity          ",
        "19-22  F4.2     km.s-1 vturb    [1.23/1.82] Micro-turbulence velocity",
        "24-28  F5.2     [-]    [Fe/H]   [-2.11/-1.5] Metallicity             ",
        "30-33  F4.2     [-]    e_[Fe/H] ? rms uncertainty on [Fe/H]          ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "S05-5   4337 0.77 1.80 -2.07     ",
        "S08-229 4625 1.23 1.23 -1.50     ",
        "S05-10  4342 0.91 1.82 -2.11 0.14",
        "S05-47  4654 1.28 1.74 -1.64 0.16",
    ]
    dat = get_pkg_data_filename(
        "data/cdsFunctional2.dat", package="astropy.io.ascii.tests"
    )
    t = Table.read(dat, format="ascii.cds")
    out = StringIO()
    # TODO: Write with CDS
    t.write(out, format="ascii.cds")
    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output
