# -*- coding: utf-8 -*-

# TEST_UNICODE_LITERALS

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is a set of regression tests for vo.
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from ....extern import six
from ....extern.six.moves import xrange

# STDLIB
import difflib
import io
import os
import shutil
import sys
import tempfile

# THIRD-PARTY
from numpy.testing import assert_array_equal
import numpy as np

# LOCAL
from ..table import parse, parse_single_table, validate
from .. import tree
from ..exceptions import VOTableSpecError, VOWarning
from ..xmlutil import validate_schema
from ....utils.data import get_pkg_data_filename, get_pkg_data_filenames
from ....tests.helper import pytest, raises, catch_warnings
from ....utils.compat import gzip

# Determine the kind of float formatting in this build of Python
if hasattr(sys, 'float_repr_style'):
    legacy_float_repr = (sys.float_repr_style == 'legacy')
else:
    legacy_float_repr = sys.platform.startswith('win')


def assert_validate_schema(filename, version):
    if sys.platform.startswith('win'):
        return

    try:
        rc, stdout, stderr = validate_schema(filename, version)
    except OSError:
        # If xmllint is not installed, we want the test to pass anyway
        return
    assert rc == 0, 'File did not validate against VOTable schema'


def test_parse_single_table():
    table = parse_single_table(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False)
    assert isinstance(table, tree.Table)
    assert len(table.array) == 5


def test_parse_single_table2():
    table2 = parse_single_table(
        get_pkg_data_filename('data/regression.xml'),
        table_number=1,
        pedantic=False)
    assert isinstance(table2, tree.Table)
    assert len(table2.array) == 1
    assert len(table2.array.dtype.names) == 28


@raises(IndexError)
def test_parse_single_table3():
    table2 = parse_single_table(
        get_pkg_data_filename('data/regression.xml'),
        table_number=3, pedantic=False)


def _test_regression(tmpdir, _python_based=False, binary_mode=1):
    # Read the VOTABLE
    votable = parse(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False,
        _debug_python_based_parser=_python_based)
    table = votable.get_first_table()

    dtypes = [
        ((str('string test'), str('string_test')), str('|O8')),
        ((str('fixed string test'), str('string_test_2')), str('|S10')),
        (str('unicode_test'), str('|O8')),
        ((str('unicode test'), str('fixed_unicode_test')), str('<U10')),
        ((str('string array test'), str('string_array_test')), str('|S4')),
        (str('unsignedByte'), str('|u1')),
        (str('short'), str('<i2')),
        (str('int'), str('<i4')),
        (str('long'), str('<i8')),
        (str('double'), str('<f8')),
        (str('float'), str('<f4')),
        (str('array'), str('|O8')),
        (str('bit'), str('|b1')),
        (str('bitarray'), str('|b1'), (3, 2)),
        (str('bitvararray'), str('|O8')),
        (str('bitvararray2'), str('|O8')),
        (str('floatComplex'), str('<c8')),
        (str('doubleComplex'), str('<c16')),
        (str('doubleComplexArray'), str('|O8')),
        (str('doubleComplexArrayFixed'), str('<c16'), (2,)),
        (str('boolean'), str('|b1')),
        (str('booleanArray'), str('|b1'), (4,)),
        (str('nulls'), str('<i4')),
        (str('nulls_array'), str('<i4'), (2, 2)),
        (str('precision1'), str('<f8')),
        (str('precision2'), str('<f8')),
        (str('doublearray'), str('|O8')),
        (str('bitarray2'), str('|b1'), (16,))
        ]
    if sys.byteorder == 'big':
        new_dtypes = []
        for dtype in dtypes:
            dtype = list(dtype)
            dtype[1] = dtype[1].replace(str('<'), str('>'))
            new_dtypes.append(tuple(dtype))
        dtypes = new_dtypes
    assert table.array.dtype == dtypes

    votable.to_xml(str(tmpdir.join("regression.tabledata.xml")),
                   _debug_python_based_parser=_python_based)
    assert_validate_schema(str(tmpdir.join("regression.tabledata.xml")),
                           votable.version)

    if binary_mode == 1:
        votable.get_first_table().format = 'binary'
        votable.version = '1.1'
    elif binary_mode == 2:
        votable.get_first_table()._config['version_1_3_or_later'] = True
        votable.get_first_table().format = 'binary2'
        votable.version = '1.3'

    # Also try passing a file handle
    with open(str(tmpdir.join("regression.binary.xml")), "wb") as fd:
        votable.to_xml(fd, _debug_python_based_parser=_python_based)
    assert_validate_schema(str(tmpdir.join("regression.binary.xml")),
                           votable.version)
    # Also try passing a file handle
    with open(str(tmpdir.join("regression.binary.xml")), "rb") as fd:
        votable2 = parse(fd, pedantic=False,
                         _debug_python_based_parser=_python_based)
    votable2.get_first_table().format = 'tabledata'
    votable2.to_xml(str(tmpdir.join("regression.bin.tabledata.xml")),
                    _astropy_version="testing",
                    _debug_python_based_parser=_python_based)
    assert_validate_schema(str(tmpdir.join("regression.bin.tabledata.xml")),
                           votable.version)

    with io.open(
        get_pkg_data_filename(
            'data/regression.bin.tabledata.truth.{0}.xml'.format(
                votable.version)),
        'rt', encoding='utf-8') as fd:
        truth = fd.readlines()
    with io.open(str(tmpdir.join("regression.bin.tabledata.xml")),
                 'rt', encoding='utf-8') as fd:
        output = fd.readlines()

    # If the lines happen to be different, print a diff
    # This is convenient for debugging
    for line in difflib.unified_diff(truth, output):
        sys.stdout.write(
            line.
            encode('unicode_escape').
            replace('\\n', '\n'))

    assert truth == output

    # Test implicit gzip saving
    votable2.to_xml(
        str(tmpdir.join("regression.bin.tabledata.xml.gz")),
        _astropy_version="testing",
        _debug_python_based_parser=_python_based)
    with gzip.GzipFile(
        str(tmpdir.join("regression.bin.tabledata.xml.gz")), 'rb') as gzfd:
        output = gzfd.readlines()
    output = [x.decode('utf-8').rstrip() for x in output]
    truth = [x.rstrip() for x in truth]

    assert truth == output


@pytest.mark.xfail(str('legacy_float_repr'))
def test_regression(tmpdir):
    _test_regression(tmpdir, False)


@pytest.mark.xfail(str('legacy_float_repr'))
def test_regression_python_based_parser(tmpdir):
    _test_regression(tmpdir, True)


@pytest.mark.xfail(str('legacy_float_repr'))
def test_regression_binary2(tmpdir):
    _test_regression(tmpdir, False, 2)


class TestFixups:
    def setup_class(self):
        self.table = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False).get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    def test_implicit_id(self):
        assert_array_equal(self.array['string_test_2'],
                           self.array['fixed string test'])


class TestReferences:
    def setup_class(self):
        self.votable = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False)
        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    def test_fieldref(self):
        fieldref = self.table.groups[1].entries[0]
        assert isinstance(fieldref, tree.FieldRef)
        assert fieldref.get_ref().name == 'boolean'
        assert fieldref.get_ref().datatype == 'boolean'

    def test_paramref(self):
        paramref = self.table.groups[0].entries[0]
        assert isinstance(paramref, tree.ParamRef)
        assert paramref.get_ref().name == 'INPUT'
        assert paramref.get_ref().datatype == 'float'

    def test_iter_fields_and_params_on_a_group(self):
        assert len(list(self.table.groups[1].iter_fields_and_params())) == 2

    def test_iter_groups_on_a_group(self):
        assert len(list(self.table.groups[1].iter_groups())) == 1

    def test_iter_groups(self):
        # Because of the ref'd table, there are more logical groups
        # than actually exist in the file
        assert len(list(self.votable.iter_groups())) == 9

    def test_ref_table(self):
        tables = list(self.votable.iter_tables())
        for x, y in zip(tables[0].array.data[0], tables[1].array.data[0]):
            assert_array_equal(x, y)

    def test_iter_coosys(self):
        assert len(list(self.votable.iter_coosys())) == 1


def test_select_columns_by_index():
    columns = [0, 5, 13]
    table = parse(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False, columns=columns).get_first_table()
    array = table.array
    mask = table.array.mask
    assert array['string_test'][0] == b"String & test"
    columns = ['string_test', 'unsignedByte', 'bitarray']
    for c in columns:
        assert not np.all(mask[c])
    assert np.all(mask['unicode_test'])


def test_select_columns_by_name():
    columns = ['string_test', 'unsignedByte', 'bitarray']
    table = parse(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False, columns=columns).get_first_table()
    array = table.array
    mask = table.array.mask
    assert array['string_test'][0] == b"String & test"
    for c in columns:
        assert not np.all(mask[c])
    assert np.all(mask['unicode_test'])


class TestParse:
    def setup_class(self):
        self.votable = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False)
        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    def test_string_test(self):
        assert issubclass(self.array['string_test'].dtype.type,
                          np.object_)
        assert_array_equal(
            self.array['string_test'],
            [b'String & test', b'String &amp; test', b'XXXX',
             b'', b''])

    def test_fixed_string_test(self):
        assert issubclass(self.array['string_test_2'].dtype.type,
                          np.string_)
        assert_array_equal(
            self.array['string_test_2'],
            [b'Fixed stri', b'0123456789', b'XXXX', b'', b''])

    def test_unicode_test(self):
        assert issubclass(self.array['unicode_test'].dtype.type,
                          np.object_)
        assert_array_equal(self.array['unicode_test'],
                           ["Ceçi n'est pas un pipe",
                            'வணக்கம்',
                            'XXXX', '', ''])

    def test_fixed_unicode_test(self):
        assert issubclass(self.array['fixed_unicode_test'].dtype.type,
                          np.unicode_)
        assert_array_equal(self.array['fixed_unicode_test'],
                           ["Ceçi n'est",
                            'வணக்கம்',
                            '0123456789', '', ''])

    def test_unsignedByte(self):
        assert issubclass(self.array['unsignedByte'].dtype.type,
                          np.uint8)
        assert_array_equal(self.array['unsignedByte'],
                           [128, 255, 0, 255, 255])
        assert not np.any(self.mask['unsignedByte'])

    def test_short(self):
        assert issubclass(self.array['short'].dtype.type,
                          np.int16)
        assert_array_equal(self.array['short'],
                           [4096, 32767, -4096, 32767, 32767])
        assert not np.any(self.mask['short'])

    def test_int(self):
        assert issubclass(self.array['int'].dtype.type,
                          np.int32)
        assert_array_equal(
            self.array['int'],
            [268435456, 2147483647, -268435456, 268435455, 123456789])
        assert_array_equal(self.mask['int'],
                           [False, False, False, False, True])

    def test_long(self):
        assert issubclass(self.array['long'].dtype.type,
                          np.int64)
        assert_array_equal(
            self.array['long'],
            [922337203685477, 123456789, -1152921504606846976,
             1152921504606846975, 123456789])
        assert_array_equal(self.mask['long'],
                           [False, True, False, False, True])

    def test_double(self):
        assert issubclass(self.array['double'].dtype.type,
                          np.float64)
        assert_array_equal(self.array['double'],
                           [8.999999, 0.0, np.inf, np.nan, -np.inf])
        assert_array_equal(self.mask['double'],
                           [False, False, False, True, False])

    def test_float(self):
        assert issubclass(self.array['float'].dtype.type,
                          np.float32)
        assert_array_equal(self.array['float'],
                           [1.0, 0.0, np.inf, np.inf, np.nan])
        assert_array_equal(self.mask['float'],
                           [False, False, False, False, True])

    def test_array(self):
        assert issubclass(self.array['array'].dtype.type,
                          np.object_)
        match = [[],
                 [[42, 32], [12, 32]],
                 [[12, 34], [56, 78], [87, 65], [43, 21]],
                 [[-1, 23]],
                 [[31, -1]]]
        for a, b in zip(self.array['array'], match):
            # assert issubclass(a.dtype.type, np.int64)
            # assert a.shape[1] == 2
            for a0, b0 in zip(a, b):
                assert issubclass(a0.dtype.type, np.int64)
                assert_array_equal(a0, b0)
        assert self.array.data['array'][3].mask[0][0]
        assert self.array.data['array'][4].mask[0][1]

    def test_bit(self):
        assert issubclass(self.array['bit'].dtype.type,
                          np.bool_)
        assert_array_equal(self.array['bit'],
                           [True, False, True, False, False])

    def test_bit_mask(self):
        assert_array_equal(self.mask['bit'],
                           [False, False, False, False, True])

    def test_bitarray(self):
        assert issubclass(self.array['bitarray'].dtype.type,
                          np.bool_)
        assert self.array['bitarray'].shape == (5, 3, 2)
        assert_array_equal(self.array['bitarray'],
                           [[[ True, False],
                             [ True,  True],
                             [False,  True]],

                            [[False,  True],
                             [False, False],
                             [ True,  True]],

                            [[ True,  True],
                             [ True, False],
                             [False, False]],

                            [[False, False],
                             [False, False],
                             [False, False]],

                            [[False, False],
                             [False, False],
                             [False, False]]])

    def test_bitarray_mask(self):
        assert_array_equal(self.mask['bitarray'],
                           [[[False, False],
                             [False, False],
                             [False, False]],

                            [[False, False],
                             [False, False],
                             [False, False]],

                            [[False, False],
                             [False, False],
                             [False, False]],

                            [[ True,  True],
                             [ True,  True],
                             [ True,  True]],

                            [[ True,  True],
                             [ True,  True],
                             [ True,  True]]])

    def test_bitvararray(self):
        assert issubclass(self.array['bitvararray'].dtype.type,
                          np.object_)
        match = [[ True,  True,  True],
                 [False, False, False, False, False],
                 [ True, False,  True, False,  True],
                 [], []]
        for a, b in zip(self.array['bitvararray'], match):
            assert_array_equal(a, b)
        match_mask = [[False, False, False],
                      [False, False, False, False, False],
                      [False, False, False, False, False],
                      False, False]
        for a, b in zip(self.array['bitvararray'], match_mask):
            assert_array_equal(a.mask, b)

    def test_bitvararray2(self):
        assert issubclass(self.array['bitvararray2'].dtype.type,
                          np.object_)
        match = [[],

                 [[[False,  True],
                   [False, False],
                   [ True, False]],
                  [[ True, False],
                   [ True, False],
                   [ True, False]]],

                 [[[ True,  True],
                   [ True,  True],
                   [ True,  True]]],

                 [],

                 []]
        for a, b in zip(self.array['bitvararray2'], match):
            for a0, b0 in zip(a, b):
                assert a0.shape == (3, 2)
                assert issubclass(a0.dtype.type, np.bool_)
                assert_array_equal(a0, b0)

    def test_floatComplex(self):
        assert issubclass(self.array['floatComplex'].dtype.type,
                          np.complex64)
        assert_array_equal(self.array['floatComplex'],
                           [np.nan+0j, 0+0j, 0+-1j, np.nan+0j, np.nan+0j])
        assert_array_equal(self.mask['floatComplex'],
                           [True, False, False, True, True])

    def test_doubleComplex(self):
        assert issubclass(self.array['doubleComplex'].dtype.type,
                          np.complex128)
        assert_array_equal(
            self.array['doubleComplex'],
            [np.nan+0j, 0+0j, 0+-1j, np.nan+(np.inf*1j), np.nan+0j])
        assert_array_equal(self.mask['doubleComplex'],
                           [True, False, False, True, True])

    def test_doubleComplexArray(self):
        assert issubclass(self.array['doubleComplexArray'].dtype.type,
                          np.object_)
        assert ([len(x) for x in self.array['doubleComplexArray']] ==
                [0, 2, 2, 0, 0])

    def test_boolean(self):
        assert issubclass(self.array['boolean'].dtype.type,
                          np.bool_)
        assert_array_equal(self.array['boolean'],
                           [True, False, True, False, False])

    def test_boolean_mask(self):
        assert_array_equal(self.mask['boolean'],
                           [False, False, False, False, True])

    def test_boolean_array(self):
        assert issubclass(self.array['booleanArray'].dtype.type,
                          np.bool_)
        assert_array_equal(self.array['booleanArray'],
                           [[ True,  True,  True,  True],
                            [ True,  True, False,  True],
                            [ True,  True, False,  True],
                            [False, False, False, False],
                            [False, False, False, False]])

    def test_boolean_array_mask(self):
        assert_array_equal(self.mask['booleanArray'],
                           [[False, False, False, False],
                            [False, False, False, False],
                            [False, False,  True, False],
                            [ True,  True,  True,  True],
                            [ True,  True,  True,  True]])

    def test_nulls(self):
        assert_array_equal(self.array['nulls'],
                           [0, -9, 2, -9, -9])
        assert_array_equal(self.mask['nulls'],
                           [False, True, False, True, True])

    def test_nulls_array(self):
        assert_array_equal(self.array['nulls_array'],
                           [[[-9, -9], [-9, -9]],
                            [[0, 1], [2, 3]],
                            [[-9, 0], [-9, 1]],
                            [[0, -9], [1, -9]],
                            [[-9, -9], [-9, -9]]])
        assert_array_equal(self.mask['nulls_array'],
                           [[[ True,  True],
                             [ True,  True]],

                            [[False, False],
                             [False, False]],

                            [[ True, False],
                             [ True, False]],

                            [[False,  True],
                             [False,  True]],

                            [[ True,  True],
                             [ True,  True]]])

    def test_double_array(self):
        assert issubclass(self.array['doublearray'].dtype.type,
                          np.object_)
        assert len(self.array['doublearray'][0]) == 0
        assert_array_equal(self.array['doublearray'][1],
                           [0, 1, np.inf, -np.inf, np.nan, 0, -1])
        assert_array_equal(self.array.data['doublearray'][1].mask,
                           [False, False, False, False, False, False, True])

    def test_bit_array2(self):
        assert_array_equal(self.array['bitarray2'][0],
                           [True, True, True, True,
                            False, False, False, False,
                            True, True, True, True,
                            False, False, False, False])

    def test_bit_array2_mask(self):
        assert not np.any(self.mask['bitarray2'][0])
        assert np.all(self.mask['bitarray2'][1:])

    def test_get_coosys_by_id(self):
        coosys = self.votable.get_coosys_by_id('J2000')
        assert coosys.system == 'eq_FK5'

    def test_get_field_by_utype(self):
        fields = list(self.votable.get_fields_by_utype("myint"))
        assert fields[0].name == "int"
        assert fields[0].values.min == -1000

    def test_get_info_by_id(self):
        info = self.votable.get_info_by_id('QUERY_STATUS')
        assert info.value == 'OK'

        if self.votable.version != '1.1':
            info = self.votable.get_info_by_id("ErrorInfo")
            assert info.value == "One might expect to find some INFO here, too..."

    def test_repr(self):
        assert '3 tables' in repr(self.votable)
        assert repr(list(self.votable.iter_fields_and_params())[0]) == \
            '<PARAM ID="awesome" arraysize="*" datatype="float" name="INPUT" unit="deg" value="[0.0 0.0]"/>'
        # Smoke test
        repr(list(self.votable.iter_groups()))


class TestThroughTableData(TestParse):
    def setup_class(self):
        votable = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False)

        self.xmlout = bio = io.BytesIO()
        votable.to_xml(bio)
        bio.seek(0)
        self.votable = parse(bio, pedantic=False)
        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    def test_bit_mask(self):
        assert_array_equal(self.mask['bit'],
                           [False, False, False, False, False])

    def test_bitarray_mask(self):
        assert not np.any(self.mask['bitarray'])

    def test_bit_array2_mask(self):
        assert not np.any(self.mask['bitarray2'])

    def test_schema(self, tmpdir):
        # have to use an actual file because assert_validate_schema only works
        # on filenames, not file-like objects
        fn = str(tmpdir.join("test_through_tabledata.xml"))
        with open(fn, 'wb') as f:
            f.write(self.xmlout.getvalue())
        assert_validate_schema(fn, '1.1')


class TestThroughBinary(TestParse):
    def setup_class(self):
        votable = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False)
        votable.get_first_table().format = 'binary'

        self.xmlout = bio = io.BytesIO()
        votable.to_xml(bio)
        bio.seek(0)
        self.votable = parse(bio, pedantic=False)

        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    # Masked values in bit fields don't roundtrip through the binary
    # representation -- that's not a bug, just a limitation, so
    # override the mask array checks here.
    def test_bit_mask(self):
        assert not np.any(self.mask['bit'])

    def test_bitarray_mask(self):
        assert not np.any(self.mask['bitarray'])

    def test_bit_array2_mask(self):
        assert not np.any(self.mask['bitarray2'])


class TestThroughBinary2(TestParse):
    def setup_class(self):
        votable = parse(
            get_pkg_data_filename('data/regression.xml'),
            pedantic=False)
        votable.version = '1.3'
        votable.get_first_table()._config['version_1_3_or_later'] = True
        votable.get_first_table().format = 'binary2'

        self.xmlout = bio = io.BytesIO()
        votable.to_xml(bio)
        bio.seek(0)
        self.votable = parse(bio, pedantic=False)

        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.array.mask

    def test_get_coosys_by_id(self):
        # No COOSYS in VOTable 1.2 or later
        pass


def table_from_scratch():
    from ..tree import VOTableFile, Resource, Table, Field

    # Create a new VOTable file...
    votable = VOTableFile()

    # ...with one resource...
    resource = Resource()
    votable.resources.append(resource)

    # ... with one table
    table = Table(votable)
    resource.tables.append(table)

    # Define some fields
    table.fields.extend([
            Field(votable, ID="filename", datatype="char"),
            Field(votable, ID="matrix", datatype="double", arraysize="2x2")])

    # Now, use those field definitions to create the numpy record arrays, with
    # the given number of rows
    table.create_arrays(2)

    # Now table.array can be filled with data
    table.array[0] = ('test1.xml', [[1, 0], [0, 1]])
    table.array[1] = ('test2.xml', [[0.5, 0.3], [0.2, 0.1]])

    # Now write the whole thing to a file.
    # Note, we have to use the top-level votable file object
    out = io.StringIO()
    votable.to_xml(out)


def test_open_files():
    def test_file(filename):
        parse(filename, pedantic=False)

    for filename in get_pkg_data_filenames('data', pattern='*.xml'):
        if filename.endswith('custom_datatype.xml'):
            continue
        yield test_file, filename


@raises(VOTableSpecError)
def test_too_many_columns():
    votable = parse(
        get_pkg_data_filename('data/too_many_columns.xml.gz'),
        pedantic=False)


def test_build_from_scratch(tmpdir):
    # Create a new VOTable file...
    votable = tree.VOTableFile()

    # ...with one resource...
    resource = tree.Resource()
    votable.resources.append(resource)

    # ... with one table
    table = tree.Table(votable)
    resource.tables.append(table)

    # Define some fields
    table.fields.extend([
        tree.Field(votable, ID="filename", datatype="char"),
        tree.Field(votable, ID="matrix", datatype="double", arraysize="2x2")])

    # Now, use those field definitions to create the numpy record arrays, with
    # the given number of rows
    table.create_arrays(2)

    # Now table.array can be filled with data
    table.array[0] = ('test1.xml', [[1, 0], [0, 1]])
    table.array[1] = ('test2.xml', [[0.5, 0.3], [0.2, 0.1]])

    # Now write the whole thing to a file.
    # Note, we have to use the top-level votable file object
    votable.to_xml(str(tmpdir.join("new_votable.xml")))

    votable = parse(str(tmpdir.join("new_votable.xml")))

    table = votable.get_first_table()
    assert_array_equal(
        table.array.mask, np.array([(False, [[False, False], [False, False]]),
                                    (False, [[False, False], [False, False]])],
                                    dtype=[(str('filename'), str('?')),
                                           (str('matrix'), str('?'), (2, 2))]))


def test_validate():
    output = io.StringIO()

    # We can't test xmllint, because we can't rely on it being on the
    # user's machine.
    with catch_warnings():
        result = validate(get_pkg_data_filename('data/regression.xml'),
                          output, xmllint=False)

    assert result == False

    output.seek(0)
    output = output.readlines()

    # Uncomment to generate new groundtruth
    # with io.open('validation.txt', 'wt', encoding='utf-8') as fd:
    #     fd.write(u''.join(output))

    with io.open(
        get_pkg_data_filename('data/validation.txt'),
        'rt', encoding='utf-8') as fd:
        truth = fd.readlines()

    truth = truth[1:]
    output = output[1:-1]

    for line in difflib.unified_diff(truth, output):
        if six.PY3:
            sys.stdout.write(
                line.replace('\\n', '\n'))
        else:
            sys.stdout.write(
                line.encode('unicode_escape').
                replace('\\n', '\n'))

    assert truth == output


def test_gzip_filehandles(tmpdir):
    votable = parse(
        get_pkg_data_filename('data/regression.xml'),
        pedantic=False)


    with open(str(tmpdir.join("regression.compressed.xml")), 'wb') as fd:
        votable.to_xml(
            fd,
            compressed=True,
            _astropy_version="testing")

    with open(str(tmpdir.join("regression.compressed.xml")), 'rb') as fd:
        votable = parse(
            fd,
            pedantic=False)


def test_from_scratch_example():
    with catch_warnings(VOWarning) as warning_lines:
        try:
            _run_test_from_scratch_example()
        except ValueError as e:
            warning_lines.append(str(e))

    assert len(warning_lines) == 0


def _run_test_from_scratch_example():
    from ..tree import VOTableFile, Resource, Table, Field

    # Create a new VOTable file...
    votable = VOTableFile()

    # ...with one resource...
    resource = Resource()
    votable.resources.append(resource)

    # ... with one table
    table = Table(votable)
    resource.tables.append(table)

    # Define some fields
    table.fields.extend([
        Field(votable, name="filename", datatype="char", arraysize="*"),
        Field(votable, name="matrix", datatype="double", arraysize="2x2")])

    # Now, use those field definitions to create the numpy record arrays, with
    # the given number of rows
    table.create_arrays(2)

    # Now table.array can be filled with data
    table.array[0] = ('test1.xml', [[1, 0], [0, 1]])
    table.array[1] = ('test2.xml', [[0.5, 0.3], [0.2, 0.1]])

    assert table.array[0][0] == 'test1.xml'


def test_fileobj():
    # Assert that what we get back is a raw C file pointer
    # so it will be super fast in the C extension.
    from ....utils.xml import iterparser
    filename = get_pkg_data_filename('data/regression.xml')
    with iterparser._convert_to_fd_or_read_function(filename) as fd:
        if sys.platform == 'win32':
            fd()
        else:
            if six.PY3:
                assert isinstance(fd, io.FileIO)
            elif six.PY2:
                assert isinstance(fd, file)


def test_nonstandard_units():
    from .... import units as u

    votable = parse(
        get_pkg_data_filename('data/nonstandard_units.xml'),
        pedantic=False)

    assert isinstance(
        votable.get_first_table().fields[0].unit, u.UnrecognizedUnit)

    votable = parse(
        get_pkg_data_filename('data/nonstandard_units.xml'),
        pedantic=False,
        unit_format='generic')

    assert not isinstance(
        votable.get_first_table().fields[0].unit, u.UnrecognizedUnit)


def test_resource_structure():
    # Based on issue #1223, as reported by @astro-friedel and @RayPlante
    from astropy.io.votable import tree as vot

    vtf = vot.VOTableFile()

    r1 = vot.Resource()
    vtf.resources.append(r1)
    t1 = vot.Table(vtf)
    t1.name = "t1"
    t2 = vot.Table(vtf)
    t2.name = 't2'
    r1.tables.append(t1)
    r1.tables.append(t2)

    r2 = vot.Resource()
    vtf.resources.append(r2)
    t3 = vot.Table(vtf)
    t3.name = "t3"
    t4 = vot.Table(vtf)
    t4.name = "t4"
    r2.tables.append(t3)
    r2.tables.append(t4)

    r3 = vot.Resource()
    vtf.resources.append(r3)
    t5 = vot.Table(vtf)
    t5.name = "t5"
    t6 = vot.Table(vtf)
    t6.name = "t6"
    r3.tables.append(t5)
    r3.tables.append(t6)

    buff = io.BytesIO()
    vtf.to_xml(buff)

    buff.seek(0)
    vtf2 = parse(buff)

    assert len(vtf2.resources) == 3

    for r in xrange(len(vtf2.resources)):
        res = vtf2.resources[r]
        assert len(res.tables) == 2
        assert len(res.resources) == 0


def test_no_resource_check():
    output = io.StringIO()

    with catch_warnings():
        # We can't test xmllint, because we can't rely on it being on the
        # user's machine.
        result = validate(get_pkg_data_filename('data/no_resource.xml'),
                          output, xmllint=False)

    assert result == False

    output.seek(0)
    output = output.readlines()

    # Uncomment to generate new groundtruth
    # with io.open('no_resource.txt', 'wt', encoding='utf-8') as fd:
    #     fd.write(u''.join(output))

    with io.open(
        get_pkg_data_filename('data/no_resource.txt'),
        'rt', encoding='utf-8') as fd:
        truth = fd.readlines()

    truth = truth[1:]
    output = output[1:-1]

    for line in difflib.unified_diff(truth, output):
        if six.PY3:
            sys.stdout.write(
                line.replace('\\n', '\n'))
        else:
            sys.stdout.write(
                line.encode('unicode_escape').
                replace('\\n', '\n'))

    assert truth == output


def test_instantiate_vowarning():
    # This used to raise a deprecation exception on Python 2.6.
    # See https://github.com/astropy/astroquery/pull/276
    VOWarning(())


def test_custom_datatype():
    votable = parse(
        get_pkg_data_filename('data/custom_datatype.xml'),
        pedantic=False,
        datatype_mapping={'bar': 'int'}
    )

    table = votable.get_first_table()
    assert table.array.dtype['foo'] == np.int32
