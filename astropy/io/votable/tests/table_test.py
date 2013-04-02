# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test the conversion to/from astropy.table
"""

import io
import os
import shutil
import tempfile
from distutils import version

import numpy as np
from ....tests.helper import pytest

from ....utils.data import get_pkg_data_filename, get_pkg_data_fileobj
from ..table import parse, writeto
from .. import tree

numpy_lt_1p5 = version.LooseVersion(np.__version__) < version.LooseVersion('1.5')

TMP_DIR = None
def setup_module():
    global TMP_DIR
    TMP_DIR = tempfile.mkdtemp()


def teardown_module():
    shutil.rmtree(TMP_DIR)


@pytest.mark.xfail('numpy_lt_1p5')
def test_table():
    # Read the VOTABLE
    votable = parse(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False)
    table = votable.get_first_table()
    astropy_table = table.to_table()

    for name in table.array.dtype.names:
        assert np.all(astropy_table.mask[name] == table.array.mask[name])

    votable2 = tree.VOTableFile.from_table(astropy_table)
    t = votable2.get_first_table()

    field_types = [
        ('string_test', {'datatype': 'char', 'arraysize': '*'}),
        ('string_test_2', {'datatype': 'char', 'arraysize': '10'}),
        ('unicode_test', {'datatype': 'unicodeChar', 'arraysize': '*'}),
        ('fixed_unicode_test', {'datatype': 'unicodeChar', 'arraysize': '10'}),
        ('string_array_test', {'datatype': 'char', 'arraysize': '4'}),
        ('unsignedByte', {'datatype': 'unsignedByte'}),
        ('short', {'datatype': 'short'}),
        ('int', {'datatype': 'int'}),
        ('long', {'datatype': 'long'}),
        ('double', {'datatype': 'double'}),
        ('float', {'datatype': 'float'}),
        ('array', {'datatype': 'long', 'arraysize': '2*'}),
        ('bit', {'datatype': 'bit'}),
        ('bitarray', {'datatype': 'bit', 'arraysize': '3x2'}),
        ('bitvararray', {'datatype': 'bit', 'arraysize': '*'}),
        ('bitvararray2', {'datatype': 'bit', 'arraysize': '3x2*'}),
        ('floatComplex', {'datatype': 'floatComplex'}),
        ('doubleComplex', {'datatype': 'doubleComplex'}),
        ('doubleComplexArray', {'datatype': 'doubleComplex', 'arraysize': '*'}),
        ('doubleComplexArrayFixed', {'datatype': 'doubleComplex', 'arraysize': '2'}),
        ('boolean', {'datatype': 'bit'}),
        ('booleanArray', {'datatype': 'bit', 'arraysize': '4'}),
        ('nulls', {'datatype': 'int'}),
        ('nulls_array', {'datatype': 'int', 'arraysize': '2x2'}),
        ('precision1', {'datatype': 'double'}),
        ('precision2', {'datatype': 'double'}),
        ('doublearray', {'datatype': 'double', 'arraysize': '*'}),
        ('bitarray2', {'datatype': 'bit', 'arraysize': '16'})]

    for field, type in zip(t.fields, field_types):
        name, d = type
        assert field.ID == name
        assert field.datatype == d['datatype']
        if 'arraysize' in d:
            assert field.arraysize == d['arraysize']

    writeto(votable2, os.path.join(TMP_DIR, "through_table.xml"))


@pytest.mark.xfail('numpy_lt_1p5')
def test_read_through_table_interface():
    from ....table import Table

    with get_pkg_data_fileobj('data/regression.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable', table_id='main_table')

    assert len(t) == 5

    fn = os.path.join(TMP_DIR, "table_interface.xml")
    t.write(fn, table_id='FOO', format='votable')

    with open(fn, 'rb') as fd:
        t2 = Table.read(fd, format='votable', table_id='FOO')

    assert len(t2) == 5


@pytest.mark.xfail('numpy_lt_1p5')
def test_read_through_table_interface2():
    from ....table import Table

    with get_pkg_data_fileobj('data/regression.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable', table_id='last_table')

    assert len(t) == 0


@pytest.mark.xfail('numpy_lt_1p5')
def test_from_table_without_mask():
    from ....table import Table, Column
    t = Table()
    c = Column(data=[1,2,3], name='a')
    t.add_column(c)
    output = io.BytesIO()
    t.write(output, format='votable')


@pytest.mark.xfail('numpy_lt_1p5')
def test_table_read_with_unnamed_tables():
    """
    Issue #927
    """
    from ....table import Table

    with get_pkg_data_fileobj('data/names.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable')

    assert len(t) == 1
