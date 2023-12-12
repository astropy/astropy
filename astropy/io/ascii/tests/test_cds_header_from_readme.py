# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest

from astropy import units as u
from astropy.io import ascii

from .common import (
    assert_almost_equal,
    assert_equal,
    setup_function,  # noqa: F401
    teardown_function,  # noqa: F401
)


def read_table1(readme, data):
    reader = ascii.Cds(readme)
    return reader.read(data)


def read_table2(readme, data):
    reader = ascii.get_reader(reader_cls=ascii.Cds, readme=readme)
    reader.outputter = ascii.TableOutputter()
    return reader.read(data)


def read_table3(readme, data):
    return ascii.read(data, readme=readme)


def test_description():
    readme = "data/cds/description/ReadMe"
    data = "data/cds/description/table.dat"
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 2)
        assert_equal(table["Cluster"].description, "Cluster name")
        assert_equal(table["Star"].description, "")
        assert_equal(table["Wave"].description, "wave ? Wavelength in Angstroms")
        assert_equal(table["El"].description, "a")
        assert_equal(
            table["ion"].description, "- Ionization stage (1 for neutral element)"
        )
        assert_equal(
            table["loggf"].description,
            "log10 of the gf value - logarithm base 10 of stat. weight times "
            "oscillator strength",
        )
        assert_equal(table["EW"].description, "Equivalent width (in mA)")
        assert_equal(
            table["Q"].description, "DAOSPEC quality parameter Q (large values are bad)"
        )


def test_multi_header():
    readme = "data/cds/multi/ReadMe"
    data = "data/cds/multi/lhs2065.dat"
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 18)
        assert_almost_equal(table["Lambda"][-1], 6479.32)
        assert_equal(table["Fnu"][-1], "0.285937")
    data = "data/cds/multi/lp944-20.dat"
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 18)
        assert_almost_equal(table["Lambda"][0], 6476.09)
        assert_equal(table["Fnu"][-1], "0.489005")


def test_glob_header():
    readme = "data/cds/glob/ReadMe"
    data = "data/cds/glob/lmxbrefs.dat"
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 291)
        assert_equal(table["Name"][-1], "J1914+0953")
        assert_equal(table["BibCode"][-2], "2005A&A...432..235R")


def test_header_from_readme():
    r = ascii.Cds("data/vizier/ReadMe")
    table = r.read("data/vizier/table1.dat")
    assert len(r.data.data_lines) == 15
    assert len(table) == 15
    assert len(table.keys()) == 18
    Bmag = [
        14.79,
        15.00,
        14.80,
        12.38,
        12.36,
        12.24,
        13.75,
        13.65,
        13.41,
        11.59,
        11.68,
        11.53,
        13.92,
        14.03,
        14.18,
    ]
    for i, val in enumerate(table.field("Bmag")):
        assert val == Bmag[i]

    table = r.read("data/vizier/table5.dat")
    assert len(r.data.data_lines) == 49
    assert len(table) == 49
    assert len(table.keys()) == 10
    Q = [
        0.289,
        0.325,
        0.510,
        0.577,
        0.539,
        0.390,
        0.957,
        0.736,
        1.435,
        1.117,
        1.473,
        0.808,
        1.416,
        2.209,
        0.617,
        1.046,
        1.604,
        1.419,
        1.431,
        1.183,
        1.210,
        1.005,
        0.706,
        0.665,
        0.340,
        0.323,
        0.391,
        0.280,
        0.343,
        0.369,
        0.495,
        0.828,
        1.113,
        0.499,
        1.038,
        0.260,
        0.863,
        1.638,
        0.479,
        0.232,
        0.627,
        0.671,
        0.371,
        0.851,
        0.607,
        -9.999,
        1.958,
        1.416,
        0.949,
    ]
    for i, val in enumerate(table.field("Q")):
        if val is np.ma.masked:
            # text value for a missing value in that table
            assert Q[i] == -9.999
        else:
            assert val == Q[i]


@pytest.mark.parametrize("reader_cls", (ascii.Cds, ascii.Mrt))
def test_cds_units(reader_cls):
    from astropy import units

    data_and_readme = "data/cds.dat"
    reader = ascii.get_reader(reader_cls)
    table = reader.read(data_and_readme)
    # column unit is GMsun (giga solar masses)
    # make sure this is parsed correctly, not as a "string" unit
    assert table["Fit"].to(units.solMass).unit == units.solMass


@pytest.mark.parametrize("reader_cls", (ascii.Cds, ascii.Mrt))
def test_cds_function_units(reader_cls):
    data_and_readme = "data/cdsFunctional.dat"
    reader = ascii.get_reader(reader_cls)
    table = reader.read(data_and_readme)
    assert table["logg"].unit == u.dex(u.cm / u.s**2)
    assert table["logTe"].unit == u.dex(u.K)
    assert table["Mass"].unit == u.Msun
    assert table["e_Mass"].unit == u.Msun
    assert table["Age"].unit == u.Myr
    assert table["e_Age"].unit == u.Myr


@pytest.mark.parametrize("reader_cls", (ascii.Cds, ascii.Mrt))
def test_cds_function_units2(reader_cls):
    # This one includes some dimensionless dex.
    data_and_readme = "data/cdsFunctional2.dat"
    reader = ascii.get_reader(reader_cls)
    table = reader.read(data_and_readme)
    assert table["Teff"].unit == u.K
    assert table["logg"].unit == u.dex(u.cm / u.s**2)
    assert table["vturb"].unit == u.km / u.s
    assert table["[Fe/H]"].unit == u.dex(u.one)
    assert table["e_[Fe/H]"].unit == u.dex(u.one)
    assert_almost_equal(
        table["[Fe/H]"].to(u.one), 10.0 ** (np.array([-2.07, -1.50, -2.11, -1.64]))
    )


def test_cds_ignore_nullable():
    # Make sure CDS reader_cls does not ignore nullabilty for columns
    # with a limit specifier
    readme = "data/cds/null/ReadMe"
    data = "data/cds/null/table.dat"
    r = ascii.Cds(readme)
    r.read(data)
    assert_equal(r.header.cols[6].description, "Temperature class codified (10)")
    assert_equal(r.header.cols[8].description, "Luminosity class codified (11)")
    assert_equal(r.header.cols[5].description, "Pericenter position angle (18)")


def test_cds_no_whitespace():
    # Make sure CDS reader_cls only checks null values when an '=' symbol is present,
    # and read description text even if there is no whitespace after '?'.
    readme = "data/cds/null/ReadMe1"
    data = "data/cds/null/table.dat"
    r = ascii.Cds(readme)
    r.read(data)
    assert_equal(r.header.cols[6].description, "Temperature class codified (10)")
    assert_equal(r.header.cols[6].null, "")
    assert_equal(r.header.cols[7].description, "Equivalent width (in mA)")
    assert_equal(r.header.cols[7].null, "-9.9")
    assert_equal(
        r.header.cols[10].description,
        "DAOSPEC quality parameter Q (large values are bad)",
    )
    assert_equal(r.header.cols[10].null, "-9.999")


def test_cds_order():
    # Make sure CDS reader_cls does not ignore order specifier that maybe present after
    # the null specifier '?'
    readme = "data/cds/null/ReadMe1"
    data = "data/cds/null/table.dat"
    r = ascii.Cds(readme)
    r.read(data)
    assert_equal(r.header.cols[5].description, "Catalogue Identification Number")
    assert_equal(r.header.cols[8].description, "Equivalent width (in mA)")
    assert_equal(r.header.cols[9].description, "Luminosity class codified (11)")


if __name__ == "__main__":  # run from main directory; not from test/
    test_header_from_readme()
    test_multi_header()
    test_glob_header()
    test_description()
    test_cds_units()
    test_cds_ignore_nullable()
    test_cds_no_whitespace()
    test_cds_order()
