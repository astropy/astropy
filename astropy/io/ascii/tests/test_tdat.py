# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``tdat`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""

import copy
from io import StringIO

import numpy as np
import pytest

from astropy.io import ascii
from astropy.io.ascii.tdat import TdatFormatWarning
from astropy.table import Table
from astropy.table.table_helpers import simple_table
from astropy.units import allclose as quantity_allclose

test_dat = [
    "<HEADER>",
    "#",
    "# TABLE: heasarc_simple",
    "# TOTAL ROWS: 7",
    "#",
    "table_name = heasarc_simple",
    'table_description = "Test table"',
    "table_security = public",
    "#",
    "# Table Parameters",
    "#",
    "field[record_number] = int4  [meta.id] (key) // Unique Identifier for Entry",
    "field[id] = int4  [meta.id] (index) // Source ID Number",
    "field[name] = char12  [meta.id;meta.main] (index) // String Name",
    "field[ra] = float8:.4f_degree [pos.eq.ra] (index) // Right Ascension",
    "field[dec] = float8:.4f_degree [pos.eq.dec] (index) // Declination",
    "field[empty] = float8:.4f // Empty // Comment",
    "#",
    "parameter_defaults = name ra dec",
    "#",
    "# Virtual Parameters",
    "#",
    "frequency_regime = Gamma-ray",
    "observatory_name = GAMMA-RAY BURSTS",
    "row_type = GRB",
    "table_author = Example et al.",
    "table_priority = 3.01",
    "table_type = Observation",
    "unique_key = record_number",
    "#",
    "# Data Format Specification",
    "#",
    "line[1] = record_number id name ra dec empty",
    "#",
    "<DATA>",
    "1|10|aaa|1.0|1.0||",
    "2|20|b|2.0|||",
    "3|30|c||3.0||",
    "4|20|||||",
    "5||||||",
    "|60|f|6.0|6.0||",
    "7| 70| g | 7.0 |7.0||",
    "<END>",
]
test_table = Table.read(test_dat, format="ascii.tdat")
# Corresponds to simple_table()
SIMPLE_LINES = [
    '<HEADER>',
    'table_name = astropy_table',
    'table_description = "A table created via astropy"',
    '#',
    '# Table Parameters',
    '#',
    'field[a] = integer',
    'field[b] = float',
    'field[c] = char1',
    '#',
    '# Data Format Specification',
    '#',
    'line[1] = a b c',
    '#',
    '<DATA>',
    '1|1.0|c|',
    '2|2.0|d|',
    '3|3.0|e|',
    '<END>'
]


def test_write_simple():
    """
    Write a simple table with common types.  This shows the compact version
    of serialization with one line per column.
    """
    t = simple_table()

    out = StringIO()
    with pytest.raises(TdatFormatWarning) as err:
        t.write(out, format="ascii.tdat")
        assert out.getvalue().splitlines() == SIMPLE_LINES


def test_write_full():
    """
    Write a full-featured table with common types and explicitly checkout output
    """
    t = test_table
    lines = ['<HEADER>',
             '# TABLE: heasarc_simple',
             '# TOTAL ROWS: 7',
             'table_name = heasarc_simple',
             'table_description = "Test table"',
             'table_security = public',
             '#',
             '# Table Parameters',
             '#',
             'field[record_number] = int4  [meta.id] (key) // Unique Identifier for Entry',
             'field[id] = int4  [meta.id] (index) // Source ID Number',
             'field[name] = char12  [meta.id;meta.main] (index) // String Name',
             'field[ra] = float8:.4f_degree [pos.eq.ra] (index) // Right Ascension',
             'field[dec] = float8:.4f_degree [pos.eq.dec] (index) // Declination',
             'field[empty] = float8:.4f // Empty // Comment',
             '#',
             'parameter_defaults = name ra dec',
             '#',
             '# Virtual Parameters',
             '#',
             'frequency_regime = Gamma-ray',
             'observatory_name = GAMMA-RAY BURSTS',
             'row_type = GRB',
             'table_author = Example et al.',
             'table_priority = 3.01',
             'table_type = Observation',
             'unique_key = record_number',
             '#',
             '# Data Format Specification',
             '#',
             'line[1] = record_number id name ra dec empty',
             '#',
             '<DATA>',
             '1|10|aaa|1.0000|1.0000||',
             '2|20|b|2.0000|||',
             '3|30|c||3.0000||',
             '4|20|||||',
             '5||||||',
             '|60|f|6.0000|6.0000||',
             '7|70|g|7.0000|7.0000||',
             '<END>']

    out = StringIO()
    t.write(out, format="ascii.tdat")
    assert out.getvalue().splitlines() == lines


def test_write_read_roundtrip():
    """
    Write a full-featured table with all types and see that it round-trips on
    readback.
    """
    t = test_table
    out = StringIO()
    t.write(out, format="ascii.tdat")

    t2s = [
        Table.read(out.getvalue(), format="ascii.tdat"),
        ascii.read(out.getvalue(), format="tdat"),
    ]
    for t2 in t2s:
        assert t.meta == t2.meta
        for name in t.colnames:
            assert t[name].attrs_equal(t2[name])
            assert (np.all(t[name] == t2[name])
                    or np.all(t[name].mask == t2[name].mask))


def test_write_read_roundtrip_empty_table(tmp_path):
    # see https://github.com/astropy/astropy/issues/13191
    with pytest.raises(TdatFormatWarning) as err:
        sfile = tmp_path / "x.tdat"
        Table().write(sfile)
        t = Table.read(sfile)
        assert len(t) == 0
        assert len(t.colnames) == 0


def test_bad_delimiter():
    """
    Passing a delimiter other than | (pipe) gives an exception
    """
    out = StringIO()
    with pytest.raises(ValueError) as err:
        test_table.write(out, format="ascii.tdat", delimiter=",")
        assert "only pipe and space delimitter is allowed in tdat format" in str(err.value)


def test_bad_header_start():
    """
    Bad header without initial <HEADER>
    """
    lines = copy.copy(SIMPLE_LINES)
    lines[0] = "<DATA>"
    with pytest.raises(ascii.tdat.TdatFormatError):
        Table.read("\n".join(lines), format="ascii.tdat")


def assert_objects_equal(obj1, obj2, attrs, compare_class=True):
    if compare_class:
        assert obj1.__class__ is obj2.__class__

    assert obj1.shape == obj2.shape

    info_attrs = [
        "info.name",
        "info.format",
        "info.unit",
        "info.description",
        "info.dtype",
    ]
    for attr in attrs + info_attrs:
        a1 = obj1
        a2 = obj2
        for subattr in attr.split("."):
            try:
                a1 = getattr(a1, subattr)
                a2 = getattr(a2, subattr)
            except AttributeError:
                a1 = a1[subattr]
                a2 = a2[subattr]

        if isinstance(a1, np.ndarray) and a1.dtype.kind == "f":
            assert quantity_allclose(a1, a2, rtol=1e-10)
        else:
            assert np.all(a1 == a2)

    # For no attrs that means we just compare directly.
    if not attrs:
        if isinstance(obj1, np.ndarray) and obj1.dtype.kind == "f":
            assert quantity_allclose(obj1, obj2, rtol=1e-15)
        else:
            assert np.all(obj1 == obj2)


def test_round_trip_masked_table_default(tmp_path):
    """Test (mostly) round-trip of MaskedColumn through tdat using default serialization
    that uses an empty string "" to mark NULL values.  Note:

    >>> simple_table(masked=True)
    <Table masked=True length=3>
      a      b     c
    int64 float64 str1
    ----- ------- ----
       --     1.0    c
        2     2.0   --
        3      --    e
    """
    filename = tmp_path / "x.tdat"

    t = simple_table(masked=True)  # int, float, and str cols with one masked element
    with pytest.raises(TdatFormatWarning) as err:
        t.write(filename)

        t2 = Table.read(filename)
        assert t2.masked is False
        assert t2.colnames == t.colnames
        for name in t2.colnames:
            # From formal perspective the round-trip columns are the "same"
            assert np.all(t2[name].mask == t[name].mask)
            assert np.all(t2[name] == t[name])
    
            # But peeking under the mask shows that the underlying data are changed
            # because by default ECSV uses "" to represent masked elements.
            t[name].mask = False
            t2[name].mask = False
            assert not np.all(t2[name] == t[name])  # Expected diff
