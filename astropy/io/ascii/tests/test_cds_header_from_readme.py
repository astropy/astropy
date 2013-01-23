# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ....tests.helper import pytest
from ... import ascii as asciitable

from .common import (raises, numpy_lt_1p5,
                     assert_equal, assert_almost_equal, assert_true,
                     setup_function, teardown_function, has_isnan)


def read_table1(readme, data):
    reader = asciitable.Cds(readme)
    return reader.read(data)


def read_table2(readme, data):
    reader = asciitable.get_reader(Reader=asciitable.Cds, readme=readme)
    reader.outputter = asciitable.TableOutputter()
    return reader.read(data)


def read_table3(readme, data):
    return asciitable.read(data, readme=readme)


def test_multi_header():
    readme = 't/cds/multi/ReadMe'
    data = 't/cds/multi/lhs2065.dat'
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 18)
        assert_almost_equal(table['Lambda'][-1], 6479.32)
        assert_equal(table['Fnu'][-1], '0.285937')


def test_glob_header():
    readme = 't/cds/glob/ReadMe'
    data = 't/cds/glob/lmxbrefs.dat'
    for read_table in (read_table1, read_table2, read_table3):
        table = read_table(readme, data)
        assert_equal(len(table), 291)
        assert_equal(table['Name'][-1], 'J1914+0953')
        assert_equal(table['BibCode'][-2], '2005A&A...432..235R')


@pytest.mark.xfail('numpy_lt_1p5')
def test_header_from_readme():
    r = asciitable.Cds("t/vizier/ReadMe")
    table = r.read("t/vizier/table1.dat")
    assert len(r.data.data_lines) == 15
    assert len(table) == 15
    assert len(table.keys()) == 18
    Bmag = [14.79,
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
            14.18]
    for i, val in enumerate(table.field('Bmag')):
        assert val == Bmag[i]

    table = r.read("t/vizier/table5.dat")
    assert len(r.data.data_lines) == 49
    assert len(table) == 49
    assert len(table.keys()) == 10
    Q = [ 0.289,
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
          0.949]
    if has_isnan:
        from .common import isnan
        for i, val in enumerate(table.field('Q')):
            if isnan(val):
                assert Q[i] == -9.999  # text value for a missing value in that
                                       # table
            else:
                assert val == Q[i]

if __name__ == "__main__":  # run from main directory; not from test/
    test_header_from_readme()
    test_multi_header()
    test_glob_header()
