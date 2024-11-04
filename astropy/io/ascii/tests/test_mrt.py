# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
from io import StringIO

import numpy as np
import pytest
from numpy.testing import assert_allclose

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io.ascii.mrt import MAX_SIZE_README_LINE, MrtHeader
from astropy.table import Column, MaskedColumn, QTable, Table
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyWarning

test_dat = [
    "names e d s i",
    "HD81809 1E-7 22.25608 +2 67",
    "HD103095 -31.6e5 +27.2500 -9E34 -30",
]


def test_find_map_note():
    good_twolines = (
        "Note (1)",
        ["", "    M = some description of M", "        with a continuation lines"],
    )
    text_with_equal = (
        "Note (1)",
        ["", "    M = 10.0 gives a magnitude and has", "    a continuation lines"],
    )
    # TODO: How should this case be treated? Only one key like this unlikely,
    # but a single-line note on a new line is also unlikely
    # I guess it boils down to which one is least likely
    single_compatible_line = ("Note (1)", ["", "    M = some flag"])

    hdr = MrtHeader()
    assert isinstance(hdr._find_map_note(*good_twolines), dict)
    assert (
        hdr._find_map_note(*good_twolines)["M"]
        == "some description of M\nwith a continuation lines"
    )
    assert hdr._find_map_note(*text_with_equal) is None
    assert isinstance(hdr._find_map_note(*single_compatible_line), dict)


def test_extract_metadata_lines():
    pattern_not_enough_groups = "A = B"
    good_pattern = "([A-Z][a-z]+):(.*)$"
    lines_with_gaps = [
        "Title: A title",
        "Authors: Someone",
        "Table: A table",
        "" "Notes (1): A note",
    ]

    hdr = MrtHeader()
    with pytest.raises(ValueError):
        hdr._extract_metadata_lines([], pattern_not_enough_groups)
    hdr._extract_metadata_lines(lines_with_gaps, good_pattern, allow_gaps=True)
    hdr._extract_metadata_lines(lines_with_gaps[:3], good_pattern, allow_gaps=False)
    with pytest.raises(ValueError):
        hdr._extract_metadata_lines(lines_with_gaps, good_pattern, allow_gaps=False)


def test_roundtrip_mrt_meta():
    """
    Tests whether or not the MRT writer can roundtrip a table,
    i.e. read a table to ``Table`` object and write it exactly
    as it is back to a file. Since, presently CDS uses a
    MRT format template while writing, only the Byte-By-Byte
    and the data section of the table can be compared between
    original and the newly written table.

    Further, the MRT Reader does not have capability to recognize
    column format from the header of a CDS/MRT table, so this test
    can work for a limited set of simple tables, which don't have
    whitespaces in the column values or mix-in columns. Because of
    this the written table output cannot be directly matched with
    the original file and have to be checked against a list of lines.
    Masked columns are read properly though, and thus are being tested
    during round-tripping.
    """
    exp_output = [
        "Title: The Taurus Spitzer Survey: New Candidate Taurus Members Selected Using",
        "    Sensitive Mid-Infrared Photometry",
        "Authors: Rebull L.M., Padgett D.L., McCabe C.-E., Hillenbrand L.A.,",
        "    Stapelfeldt K.R., Noriega-Crespo A., Carey S.J., Brooke T., Huard T.,",
        "    Terebey S., Audard M., Monin J.-L., Fukagawa M., Gudel M., Knapp G.R.,",
        "    Menard F., Allen L.E., Angione J.R., Baldovin-Saavedra C., Bouvier J.,",
        "    Briggs K., Dougados C., Evans N.J., Flagey N., Guieu S., Grosso N.,",
        "    Glauser A.M., Harvey P., Hines D., Latter W.B., Skinner S.L., Strom S.,",
        "    Tromp J., Wolf S.",
        "Table: Spitzer measurements for sample of previously identified Taurus members",
        "================================================================================",
        "Byte-by-byte Description of file: table.dat",
        "--------------------------------------------------------------------------------",
        " Bytes Format Units  Label     Explanations",
        "--------------------------------------------------------------------------------",
        "  1- 15  A15    ---    SST      Spitzer Tau name",
        " 17- 39  A23    ---    CName    Common name",
        "     41  A1     ---    l_3.6mag ? Limit flag on 3.6mag",
        " 43- 47  F5.2   mag    3.6mag   [6.62/15.33] Spitzer/IRAC 3.6 micron band",
        "                                 magnitude (1)",
        " 49- 52  F4.2   mag    e_3.6mag [0.05/0.15]? Uncertainty in 3.6mag",
        "     54  A1     ---    l_4.5mag ? Limit flag on 4.5mag",
        " 56- 60  F5.2   mag    4.5mag   [6.1/14.25]? Spitzer/IRAC 4.5 micron band",
        "                                 magnitude (1)",
        " 62- 65  F4.2   mag    e_4.5mag [0.05/0.22]? Uncertainty in 4.5mag",
        "     67  A1     ---    l_5.8mag ? Limit flag on 5.8mag",
        " 69- 73  F5.2   mag    5.8mag   [3.49/13.7] Spitzer/IRAC 5.8 micron band",
        "                                 magnitude (1)",
        " 75- 78  F4.2   mag    e_5.8mag [0.05/0.12]? Uncertainty in 5.8mag",
        "     80  A1     ---    l_8mag   ? Limit flag on 8.0mag",
        " 82- 86  F5.2   mag    8mag     [3.52/13.41]? Spitzer/IRAC 8.0 micron band",
        "                                 magnitude (1)",
        " 88- 91  F4.2   mag    e_8mag   [0.04/0.17]? Uncertainty in 8mag",
        "     93  A1     ---    l_24mag  ? Limit flag on 24mag",
        " 95- 99  F5.2   mag    24mag    [0.45/11.05]? Spitzer/MIPS 24 micron band",
        "                                 magnitude (1)",
        "101-104  F4.2   mag    e_24mag  [0.01/0.33]? Uncertainty in 24mag",
        "    106  A1     ---    l_70mag  ? Limit flag on 70mag",
        "108-112  F5.2   mag    70mag    [-2.54/3.14]? Spitzer/MIPS 70 micron band",
        "                                 magnitude (1)",
        "114-117  F4.2   mag    e_70mag  ? Uncertainty in 70mag",
        "    119  A1     ---    l_160mag ? Limit flag on 160mag",
        "121-125  F5.2   mag    160mag   ? Spitzer/MIPS 160 micron band magnitude (1)",
        "127-130  F4.2   mag    e_160mag ? Uncertainty in 160mag",
        "132-134  A3     ---    ID24/70  ? Identification in 24/70 micron color-magnitude",
        "                                 diagram",
        "136-144  A9     ---    IDKS/24  ? Identification in Ks/70 micron color-magnitude",
        "                                 diagram",
        "146-154  A9     ---    ID8/24   ? Identification in 8/24 micron color-magnitude",
        "                                 diagram",
        "156-164  A9     ---    ID4.5/8  ? Identification in 4.5/8 micron color-magnitude",
        "                                 diagram",
        "166-168  A3     ---    IDIRAC   ? Identification in IRAC color-color diagram",
        "170-172  A3     ---    Note     ? Additional note (2)",
        "--------------------------------------------------------------------------------",
        "Note (1): To convert between magnitudes and flux densities, we use",
        "    M= 2.5 log(F_zeropt_/F) where the zero-point flux densities for the",
        "    seven Spitzer bands are 280.9, 179.7, 115.0, and 64.13 Jy for IRAC",
        "    and 7.14, 0.775, and 0.159 Jy for MIPS.  IRAC effective wavelengths",
        "    are 3.6, 4.5, 5.8, and 8.0 microns; MIPS effective wavelengths are",
        "    24, 70, and 160 microns.",
        "Note (2):",
        "    b = MIPS-160 flux density for this object is subject to confusion with a",
        "        nearby source or sources.",
        "    c = MIPS-160 flux density for this object is compromised by missing and/or",
        "        saturated data.",
        "    d = MIPS-160 flux density for this object is hard saturated.",
        "    e = IRAC flux densities for 043835.4+261041=HV Tau C do not appear in our",
        "        automatically-extracted catalog. Flux densities here are those from",
        "        Hartmann et al. (2005); since their observations have more redundancy",
        "        at IRAC bands, they are able to obtain reliable flux densities for",
        "        this object at IRAC bands.  MIPS flux densities are determined from",
        "        our data.",
        "    f = The image morphology around 041426.2+280603 is complex; careful PSF",
        "        subtraction and modeling will be required to apportion flux densities",
        "        among the three local maxima seen in close proximity in the IRAC",
        "        images, which may or may not be three physically distinct sources.",
    ]
    dat = get_pkg_data_filename("data/cds2.dat", package="astropy.io.ascii.tests")
    t = Table.read(dat, format="ascii.mrt")
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    reversed_index = lines[::-1].index("-" * 80)
    meta_lines = lines[: len(lines) - reversed_index - 1]
    assert all(len(line) <= MAX_SIZE_README_LINE for line in meta_lines)
    assert meta_lines == exp_output


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
