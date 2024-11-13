# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``CDS`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""

from io import StringIO

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.table import Column, MaskedColumn, QTable, Table
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyWarning

test_dat = [
    "names e d s i",
    "HD81809 1E-7 22.25608 +2 67",
    "HD103095 -31.6e5 +27.2500 -9E34 -30",
]


def test_roundtrip_mrt_table():
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
    t = Table.read(dat, format="ascii.mrt")
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output


def test_write_byte_by_byte_units():
    t = ascii.read(test_dat)
    col_units = [None, u.C, u.kg, u.m / u.s, u.year]
    t._set_column_attribute("unit", col_units)
    # Add a column with magnitude units.
    # Note that magnitude has to be assigned for each value explicitly.
    t["magnitude"] = [u.Magnitude(25), u.Magnitude(-9)]
    col_units.append(u.mag)
    out = StringIO()
    t.write(out, format="ascii.mrt")
    # Read written table.
    columns = ascii.read(out.getvalue(), format="cds").itercols()
    assert [col.unit for col in columns] == col_units


def test_write_readme_with_default_options():
    exp_output = [
        "Title:",
        "Authors:",
        "Table:",
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names              ",
        "10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e",
        "16-23  F8.5   ---    d       [22.25/27.25] Description of d    ",
        "25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     ",
        "33-35  I3     ---    i       [-30/67] Description of i         ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07  22.25608   2e+00  67",
        "HD103095 -3e+06 27.25000  -9e+34 -30",
    ]
    t = ascii.read(test_dat)
    out = StringIO()
    t.write(out, format="ascii.mrt")
    assert out.getvalue().splitlines() == exp_output


def test_write_empty_table():
    out = StringIO()
    import pytest

    with pytest.raises(NotImplementedError):
        Table().write(out, format="ascii.mrt")


def test_write_null_data_values():
    exp_output = [
        "HD81809  1e-07  22.25608  2.0e+00  67",
        "HD103095 -3e+06 27.25000 -9.0e+34 -30",
        "Sun                       5.3e+27    ",
    ]
    t = ascii.read(test_dat)
    t.add_row(
        ["Sun", "3.25", "0", "5.3e27", "2"], mask=[False, True, True, False, True]
    )
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines) if s.startswith(("------", "======="))]
    lines = lines[i_secs[-1] + 1 :]  # Last section is the data.
    assert lines == exp_output


def test_write_byte_by_byte_for_masked_column():
    """
    This test differs from the ``test_write_null_data_values``
    above in that it tests the column value limits in the Byte-By-Byte
    description section for columns whose values are masked.
    It also checks the description for columns with same values.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names          ",
        "10-14  E5.1   ---    e       [0.0/0.01]? Description of e  ",
        "16-17  F2.0   ---    d       ? Description of d            ",
        "19-25  E7.1   ---    s       [-9e+34/2.0] Description of s ",
        "27-29  I3     ---    i       [-30/67] Description of i     ",
        "31-33  F3.1   ---    sameF   [5.0/5.0] Description of sameF",
        "35-36  I2     ---    sameI   [20] Description of sameI     ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07    2e+00  67 5.0 20",
        "HD103095         -9e+34 -30 5.0 20",
    ]
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name="sameF")
    t.add_column([20, 20], name="sameI")
    t["e"] = MaskedColumn(t["e"], mask=[False, True])
    t["d"] = MaskedColumn(t["d"], mask=[True, True])
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output


exp_coord_cols_output = {
    "generic": [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names              ",
        "10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e",
        "16-23  F8.5   ---    d       [22.25/27.25] Description of d    ",
        "25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     ",
        "33-35  I3     ---    i       [-30/67] Description of i         ",
        "37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    ",
        "41-42  I2     ---    sameI   [20] Description of sameI         ",
        "44-45  I2     h      RAh     Right Ascension (hour)            ",
        "47-48  I2     min    RAm     Right Ascension (minute)          ",
        "50-62  F13.10 s      RAs     Right Ascension (second)          ",
        "   64  A1     ---    DE-     Sign of Declination               ",
        "65-66  I2     deg    DEd     Declination (degree)              ",
        "68-69  I2     arcmin DEm     Declination (arcmin)              ",
        "71-82  F12.9  arcsec DEs     Declination (arcsec)              ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07  22.25608   2e+00  67 5.0 20 22 02 15.4500000000 -61 39 34.599996000",
        "HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12 48 15.2244072000 +17 46 26.496624000",
    ],
    "positive_de": [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names              ",
        "10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e",
        "16-23  F8.5   ---    d       [22.25/27.25] Description of d    ",
        "25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     ",
        "33-35  I3     ---    i       [-30/67] Description of i         ",
        "37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    ",
        "41-42  I2     ---    sameI   [20] Description of sameI         ",
        "44-45  I2     h      RAh     Right Ascension (hour)            ",
        "47-48  I2     min    RAm     Right Ascension (minute)          ",
        "50-62  F13.10 s      RAs     Right Ascension (second)          ",
        "   64  A1     ---    DE-     Sign of Declination               ",
        "65-66  I2     deg    DEd     Declination (degree)              ",
        "68-69  I2     arcmin DEm     Declination (arcmin)              ",
        "71-82  F12.9  arcsec DEs     Declination (arcsec)              ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07  22.25608   2e+00  67 5.0 20 12 48 15.2244072000 +17 46 26.496624000",
        "HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12 48 15.2244072000 +17 46 26.496624000",
    ],
    "galactic": [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names              ",
        "10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e",
        "16-23  F8.5   ---    d       [22.25/27.25] Description of d    ",
        "25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     ",
        "33-35  I3     ---    i       [-30/67] Description of i         ",
        "37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    ",
        "41-42  I2     ---    sameI   [20] Description of sameI         ",
        "44-59  F16.12 deg    GLON    Galactic Longitude                ",
        "61-76  F16.12 deg    GLAT    Galactic Latitude                 ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07  22.25608   2e+00  67 5.0 20 330.071639591690 -45.548080484609",
        "HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 330.071639591690 -45.548080484609",
    ],
    "ecliptic": [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names   Description of names                       ",
        "10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e         ",
        "16-23  F8.5   ---    d       [22.25/27.25] Description of d             ",
        "25-31  E7.1   ---    s       [-9e+34/2.0] Description of s              ",
        "33-35  I3     ---    i       [-30/67] Description of i                  ",
        "37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF             ",
        "41-42  I2     ---    sameI   [20] Description of sameI                  ",
        "44-59  F16.12 deg    ELON    Ecliptic Longitude (geocentrictrueecliptic)",
        "61-76  F16.12 deg    ELAT    Ecliptic Latitude (geocentrictrueecliptic) ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809  1e-07  22.25608   2e+00  67 5.0 20 306.224208650096 -45.621789850825",
        "HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 306.224208650096 -45.621789850825",
    ],
}


def test_write_coord_cols():
    """
    There can only be one such coordinate column in a single table,
    because division of columns into individual component columns requires
    iterating over the table columns, which will have to be done again
    if additional such coordinate columns are present.
    """
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name="sameF")
    t.add_column([20, 20], name="sameI")

    # Coordinates of ASASSN-15lh
    coord = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    # Coordinates of ASASSN-14li
    coordp = SkyCoord(192.06343503, 17.77402684, unit=u.deg)
    cols = [
        Column([coord, coordp]),  # Generic coordinate column
        coordp,  # Coordinate column with positive DEC
        coord.galactic,  # Galactic coordinates
        coord.geocentrictrueecliptic,  # Ecliptic coordinates
    ]

    # Loop through different types of coordinate columns.
    for col, coord_type in zip(cols, exp_coord_cols_output):
        exp_output = exp_coord_cols_output[coord_type]
        t["coord"] = col
        out = StringIO()
        t.write(out, format="ascii.mrt")
        lines = out.getvalue().splitlines()
        i_bbb = lines.index("=" * 80)
        lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
        # Check the written table.
        assert lines == exp_output

        # Check if the original table columns remains unmodified.
        assert t.colnames == ["names", "e", "d", "s", "i", "sameF", "sameI", "coord"]


def test_write_byte_by_byte_bytes_col_format():
    """
    Tests the alignment of Byte counts with respect to hyphen
    in the Bytes column of Byte-By-Byte. The whitespace around the
    hyphen is govered by the number of digits in the total Byte
    count. Single Byte columns should have a single Byte count
    without the hyphen.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 8  A8     ---    names         Description of names              ",
        "10-21  E12.6  ---    e             [-3160000.0/0.01] Description of e",
        "23-30  F8.5   ---    d             [22.25/27.25] Description of d    ",
        "32-38  E7.1   ---    s             [-9e+34/2.0] Description of s     ",
        "40-42  I3     ---    i             [-30/67] Description of i         ",
        "44-46  F3.1   ---    sameF         [5.0/5.0] Description of sameF    ",
        "48-49  I2     ---    sameI         [20] Description of sameI         ",
        "   51  I1     ---    singleByteCol [2] Description of singleByteCol  ",
        "53-54  I2     h      RAh           Right Ascension (hour)            ",
        "56-57  I2     min    RAm           Right Ascension (minute)          ",
        "59-71  F13.10 s      RAs           Right Ascension (second)          ",
        "   73  A1     ---    DE-           Sign of Declination               ",
        "74-75  I2     deg    DEd           Declination (degree)              ",
        "77-78  I2     arcmin DEm           Declination (arcmin)              ",
        "80-91  F12.9  arcsec DEs           Declination (arcsec)              ",
        "--------------------------------------------------------------------------------",
    ]
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name="sameF")
    t.add_column([20, 20], name="sameI")
    t["coord"] = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t["singleByteCol"] = [2, 2]
    t["e"].format = ".5E"
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines) if s.startswith(("------", "======="))]
    # Select only the Byte-By-Byte section.
    lines = lines[i_secs[0] : i_secs[-2]]
    lines.append("-" * 80)  # Append a separator line.
    assert lines == exp_output


def test_write_byte_by_byte_wrapping():
    """
    Test line wrapping in the description column of the
    Byte-By-Byte section of the ReadMe.
    """
    exp_output = """\
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8  A8     ---    thisIsALongColumnLabel This is a tediously long
                                           description. But they do sometimes
                                           have them. Better to put extra
                                           details in the notes. This is a
                                           tediously long description. But they
                                           do sometimes have them. Better to put
                                           extra details in the notes.
10-14  E5.1   ---    e                      [-3160000.0/0.01] Description of e
16-23  F8.5   ---    d                      [22.25/27.25] Description of d
--------------------------------------------------------------------------------
"""
    t = ascii.read(test_dat)
    t.remove_columns(["s", "i"])
    description = (
        "This is a tediously long description."
        " But they do sometimes have them."
        " Better to put extra details in the notes. "
    )
    t["names"].description = description * 2
    t["names"].name = "thisIsALongColumnLabel"
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines) if s.startswith(("------", "======="))]
    # Select only the Byte-By-Byte section.
    lines = lines[i_secs[0] : i_secs[-2]]
    lines.append("-" * 80)  # Append a separator line.
    assert lines == exp_output.splitlines()


def test_write_mixin_and_broken_cols():
    """
    Tests conversion to string values for ``mix-in`` columns other than
    ``SkyCoord`` and for columns with only partial ``SkyCoord`` values.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        "  1-  7  A7     ---    name    Description of name       ",
        "  9- 74  A66    ---    Unknown Description of Unknown    ",
        " 76-114  A39    ---    cart    Description of cart       ",
        "116-138  A23    ---    time    Description of time       ",
        "140-142  F3.1   m      q       [1.0/1.0] Description of q",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD81809 <SkyCoord (ICRS): (ra, dec) in deg",
        "    (330.564375, -61.65961111)> (0.41342785, -0.23329341, -0.88014294)  2019-01-01 00:00:00.000 1.0",
        "random  12                                                                 (0.41342785, -0.23329341, -0.88014294)  2019-01-01 00:00:00.000 1.0",
    ]
    t = Table()
    t["name"] = ["HD81809"]
    coord = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t["coord"] = Column(coord)
    t.add_row(["random", 12])
    t["cart"] = coord.cartesian
    t["time"] = Time("2019-1-1")
    t["q"] = u.Quantity(1.0, u.m)
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    # Check the written table.
    assert lines == exp_output


def test_write_extra_skycoord_cols():
    """
    Tests output for cases when table contains multiple ``SkyCoord`` columns.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 7  A7     ---    name    Description of name     ",
        " 9-10  I2     h      RAh     Right Ascension (hour)  ",
        "12-13  I2     min    RAm     Right Ascension (minute)",
        "15-27  F13.10 s      RAs     Right Ascension (second)",
        "   29  A1     ---    DE-     Sign of Declination     ",
        "30-31  I2     deg    DEd     Declination (degree)    ",
        "33-34  I2     arcmin DEm     Declination (arcmin)    ",
        "36-47  F12.9  arcsec DEs     Declination (arcsec)    ",
        "49-62  A14    ---    coord2  Description of coord2   ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD4760   0 49 39.9000000000 +06 24 07.999200000 12.4163 6.407 ",
        "HD81809 22 02 15.4500000000 -61 39 34.599996000 330.564 -61.66",
    ]
    t = Table()
    t["name"] = ["HD4760", "HD81809"]
    t["coord1"] = SkyCoord([12.41625, 330.564375], [6.402222, -61.65961111], unit=u.deg)
    t["coord2"] = SkyCoord([12.41630, 330.564400], [6.407, -61.66], unit=u.deg)
    out = StringIO()
    with pytest.warns(
        UserWarning,
        match=r"column 2 is being skipped with designation of a "
        r"string valued column `coord2`",
    ):
        t.write(out, format="ascii.mrt")

    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and following lines.
    # Check the written table.
    assert lines[:-2] == exp_output[:-2]

    for a, b in zip(lines[-2:], exp_output[-2:]):
        assert a[:18] == b[:18]
        assert a[30:42] == b[30:42]
        assert_allclose(np.fromstring(a[2:], sep=" "), np.fromstring(b[2:], sep=" "))


def test_write_skycoord_with_format():
    """
    Tests output with custom setting for ``SkyCoord`` (second) columns.
    """
    exp_output = [
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        " 1- 7  A7     ---    name    Description of name     ",
        " 9-10  I2     h      RAh     Right Ascension (hour)  ",
        "12-13  I2     min    RAm     Right Ascension (minute)",
        "15-19  F5.2   s      RAs     Right Ascension (second)",
        "   21  A1     ---    DE-     Sign of Declination     ",
        "22-23  I2     deg    DEd     Declination (degree)    ",
        "25-26  I2     arcmin DEm     Declination (arcmin)    ",
        "28-31  F4.1   arcsec DEs     Declination (arcsec)    ",
        "--------------------------------------------------------------------------------",
        "Notes:",
        "--------------------------------------------------------------------------------",
        "HD4760   0 49 39.90 +06 24 08.0",
        "HD81809 22 02 15.45 -61 39 34.6",
    ]
    t = Table()
    t["name"] = ["HD4760", "HD81809"]
    t["coord"] = SkyCoord([12.41625, 330.564375], [6.402222, -61.65961111], unit=u.deg)

    out = StringIO()
    # This will raise a warning because `formats` is checked before the writer creating the
    # final list of columns is called.
    with pytest.warns(
        AstropyWarning,
        match=r"The key.s. {'[RD][AE]s', '[RD][AE]s'} specified in "
        r"the formats argument do not match a column name.",
    ):
        t.write(out, format="ascii.mrt", formats={"RAs": "05.2f", "DEs": "04.1f"})

    lines = out.getvalue().splitlines()
    i_bbb = lines.index("=" * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and following lines.
    # Check the written table.
    assert lines == exp_output


def test_write_qtable():
    # Regression test for gh-12804
    qt = QTable([np.arange(4) * u.m, ["a", "b", "c", "ddd"]], names=["a", "b"])
    out = StringIO()
    qt.write(out, format="mrt")
    result = out.getvalue()
    assert "Description of a" in result
    assert "Description of b" in result
