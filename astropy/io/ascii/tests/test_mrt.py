# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
from io import StringIO

import pytest

from astropy.io.ascii.cds import MAX_SIZE_README_LINE
from astropy.io.ascii.mrt import MrtHeader
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename


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
    # I guess it boils down to which one is least likely?
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
    dat = get_pkg_data_filename("data/mrt2.dat", package="astropy.io.ascii.tests")
    t = Table.read(dat, format="ascii.mrt")
    out = StringIO()
    t.write(out, format="ascii.mrt")
    lines = out.getvalue().splitlines()
    reversed_index = lines[::-1].index("-" * 80)
    meta_lines = lines[: len(lines) - reversed_index - 1]
    assert all(len(line) <= MAX_SIZE_README_LINE for line in meta_lines)
    assert meta_lines == exp_output
