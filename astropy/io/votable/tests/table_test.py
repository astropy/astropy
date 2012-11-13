"""
Test the conversion to/from astropy.table
"""

import os
import shutil
import tempfile

import numpy as np

from ....config import get_data_filename
from ..table import parse, writeto
from .. import tree


TMP_DIR = None
def setup_module():
    global TMP_DIR
    TMP_DIR = tempfile.mkdtemp()


def teardown_module():
    shutil.rmtree(TMP_DIR)


def test_table():
    # Read the VOTABLE
    votable = parse(
        get_data_filename('data/regression.xml'),
        pedantic=False)
    table = votable.get_first_table()
    astropy_table = table.to_table()
    mask = np.array([ (False, False, False, False, False, False, False, False, False,
                       False, False, False, False,
                       [[False, False], [False, False], [False, False]], False, False,
                       True, True, False, [False, False], False, [False, False, False, False],
                       False, [[True, True], [True, True]], False, False, False,
                       [False, False, False, False, False, False, False, False, False, False,
                        False, False, False, False, False, False]),
                      (False, False, False, False, False, False, False, False, True, False,
                       False, False, False, [[False, False], [False, False], [False, False]],
                       False, False, False, False, False, [False, False], False,
                       [False, False, False, False], True, [[False, False], [False, False]],
                       False, False, False, [True, True, True, True, True, True, True, True,
                                             True, True, True, True, True, True, True, True]),
                      (False, False, False, False, False, False, False, False, False, False,
                       False, False, False, [[False, False], [False, False], [False, False]],
                       False, False, False, False, False, [False, False], False,
                       [False, False, True, False], False, [[True, False], [True, False]],
                       False, False, False, [True, True, True, True, True, True, True, True,
                                             True, True, True, True, True, True, True, True]),
                      (False, False, False, False, False, False, False, False, False, True,
                       False, False, False, [[True, True], [True, True], [True, True]],
                       False, False, True, True, False, [False, False], False,
                       [True, True, True, True], True, [[False, True], [False, True]],
                       True, True, False, [True, True, True, True, True, True, True, True, True,
                                           True, True, True, True, True, True, True]),
                      (False, False, False, False, False, False, False, True, True, False, True,
                       False, True, [[True, True], [True, True], [True, True]],
                       False, False, True, True, False, [False, False], True,
                       [True, True, True, True], True, [[True, True], [True, True]],
                       True, True, False, [True, True, True, True, True, True, True, True, True,
                                           True, True, True, True, True, True, True])],
                    dtype=[('string_test', '|b1'), ('string_test_2', '|b1'),
                           ('unicode_test', '|b1'), ('fixed_unicode_test', '|b1'),
                           ('string_array_test', '|b1'), ('unsignedByte', '|b1'),
                           ('short', '|b1'), ('int', '|b1'), ('long', '|b1'), ('double', '|b1'),
                           ('float', '|b1'), ('array', '|b1'), ('bit', '|b1'),
                           ('bitarray', '|b1', (3, 2)), ('bitvararray', '|b1'),
                           ('bitvararray2', '|b1'), ('floatComplex', '|b1'),
                           ('doubleComplex', '|b1'), ('doubleComplexArray', '|b1'),
                           ('doubleComplexArrayFixed', '|b1', (2,)), ('boolean', '|b1'),
                           ('booleanArray', '|b1', (4,)), ('nulls', '|b1'),
                           ('nulls_array', '|b1', (2, 2)), ('precision1', '|b1'),
                           ('precision2', '|b1'), ('doublearray', '|b1'),
                           ('bitarray2', '|b1', (16,))])

    assert np.all(astropy_table.mask == mask)

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
