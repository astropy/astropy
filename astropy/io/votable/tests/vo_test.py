"""
This is a set of regression tests for vo.
"""

from __future__ import absolute_import, print_function

# STDLIB
import difflib
from distutils import version
import glob
import io
import os
import shutil
import subprocess
import sys
import tempfile
import warnings

# THIRD-PARTY
from numpy.testing import assert_array_equal
import numpy as np

# LOCAL
from ..table import parse, parse_single_table, validate
from .. import tree
from ..util import IS_PY3K
from ..exceptions import VOTableSpecError, VOWarning
from ..xmlutil import validate_schema
from ....config import get_data_filename, get_data_fileobj, get_data_filenames
from ....tests.helper import pytest, raises
from ....utils.compat import gzip

numpy_has_complex_bug = (
    version.LooseVersion(np.__version__) < version.LooseVersion('1.5')
    )

join = os.path.join


# Update this to use py.test's tmpdir functionality
def setup_module():
    global TMP_DIR
    TMP_DIR = tempfile.mkdtemp()


def teardown_module():
    shutil.rmtree(TMP_DIR)


def assert_validate_schema(filename):
    if sys.platform.startswith('win'):
        return

    try:
        rc, stdout, stderr = validate_schema(filename, '1.1')
    except OSError:
        # If xmllint is not installed, we want the test to pass anyway
        return
    assert rc == 0, 'File did not validate against VOTable schema'


def test_parse_single_table():
    table = parse_single_table(
        get_data_filename('data/regression.xml'),
        pedantic=False)
    assert isinstance(table, tree.Table)
    assert len(table.array) == 5


def test_parse_single_table2():
    table2 = parse_single_table(
        get_data_filename('data/regression.xml'),
        table_number=1,
        pedantic=False)
    assert isinstance(table2, tree.Table)
    assert len(table2.array) == 1
    assert len(table2.array.dtype.names) == 28


@raises(IndexError)
def test_parse_single_table3():
    table2 = parse_single_table(
        get_data_filename('data/regression.xml'),
        table_number=2, pedantic=False)


def _test_regression(_python_based=False):
    # Read the VOTABLE
    votable = parse(
        get_data_filename('data/regression.xml'),
        pedantic=False,
        _debug_python_based_parser=_python_based)
    table = votable.get_first_table()

    assert table.array.dtype == [
        (('string test', 'string_test'), '|O8'),
        (('fixed string test', 'string_test_2'), '|S10'),
        ('unicode_test', '|O8'),
        (('unicode test', 'fixed_unicode_test'), '<U10'),
        (('string array test', 'string_array_test'), '|S4'),
        ('unsignedByte', '|u1'),
        ('short', '<i2'),
        ('int', '<i4'),
        ('long', '<i8'),
        ('double', '<f8'),
        ('float', '<f4'),
        ('array', '|O8'),
        ('bit', '|b1'),
        ('bitarray', '|b1', (3, 2)),
        ('bitvararray', '|O8'),
        ('bitvararray2', '|O8'),
        ('floatComplex', '<c8'),
        ('doubleComplex', '<c16'),
        ('doubleComplexArray', '|O8'),
        ('doubleComplexArrayFixed', '<c16', (2,)),
        ('boolean', '|b1'),
        ('booleanArray', '|b1', (4,)),
        ('nulls', '<i4'),
        ('nulls_array', '<i4', (2, 2)),
        ('precision1', '<f8'),
        ('precision2', '<f8'),
        ('doublearray', '|O8'),
        ('bitarray2', '|b1', (16,))
        ]

    votable.to_xml(join(TMP_DIR, "regression.tabledata.xml"),
                   _debug_python_based_parser=_python_based)
    assert_validate_schema(join(TMP_DIR, "regression.tabledata.xml"))
    votable.get_first_table().format = 'binary'
    # Also try passing a file handle
    with open(join(TMP_DIR, "regression.binary.xml"), "wb") as fd:
        votable.to_xml(fd, _debug_python_based_parser=_python_based)
    assert_validate_schema(join(TMP_DIR, "regression.binary.xml"))
    # Also try passing a file handle
    with open(join(TMP_DIR, "regression.binary.xml"), "rb") as fd:
        votable2 = parse(fd, pedantic=False,
                         _debug_python_based_parser=_python_based)
    votable2.get_first_table().format = 'tabledata'
    votable2.to_xml(join(TMP_DIR, "regression.bin.tabledata.xml"),
                    _astropy_version="testing",
                    _debug_python_based_parser=_python_based)
    assert_validate_schema(join(TMP_DIR, "regression.bin.tabledata.xml"))

    with io.open(
        get_data_filename('data/regression.bin.tabledata.truth.xml'),
        'rt', encoding='utf-8') as fd:
        truth = fd.readlines()
    with io.open(
        join(TMP_DIR, "regression.bin.tabledata.xml"),
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
        join(TMP_DIR, "regression.bin.tabledata.xml.gz"),
        _astropy_version="testing",
        _debug_python_based_parser=_python_based)
    with gzip.GzipFile(
        join(TMP_DIR, "regression.bin.tabledata.xml.gz"), 'rb') as gzfd:
        output = gzfd.readlines()
    output = [x.decode('utf-8').rstrip() for x in output]
    truth = [x.rstrip() for x in truth]

    assert truth == output


def test_regression():
    _test_regression(False)


def test_regression_python_based_parser():
    _test_regression(True)


class TestFixups:
    def setup_class(self):
        self.table = parse(
            get_data_filename('data/regression.xml'),
            pedantic=False).get_first_table()
        self.array = self.table.array
        self.mask = self.table.mask

    def test_implicit_id(self):
        assert_array_equal(self.array['string_test_2'],
                           self.array['fixed string test'])


class TestReferences:
    def setup_class(self):
        self.votable = parse(
            get_data_filename('data/regression.xml'),
            pedantic=False)
        self.table = self.votable.get_first_table()
        self.array = self.table.array
        self.mask = self.table.mask

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
        assert len(list(self.votable.iter_groups())) == 6

    def test_ref_table(self):
        tables = list(self.votable.iter_tables())
        for x, y in zip(tables[0].array[0], tables[1].array[0]):
            assert_array_equal(x, y)

    def test_iter_coosys(self):
        assert len(list(self.votable.iter_coosys())) == 1


def test_select_columns_by_index():
    columns = [0, 5, 13]
    table = parse(
        get_data_filename('data/regression.xml'),
        pedantic=False, columns=columns).get_first_table()
    array = table.array
    mask = table.mask
    assert array['string_test'][0] == b"String & test"
    columns = ['string_test', 'unsignedByte', 'bitarray']
    for c in columns:
        assert not np.all(mask[c])
    assert np.all(mask['unicode_test'])


def test_select_columns_by_name():
    columns = ['string_test', 'unsignedByte', 'bitarray']
    table = parse(
        get_data_filename('data/regression.xml'),
        pedantic=False, columns=columns).get_first_table()
    array = table.array
    mask = table.mask
    assert array['string_test'][0] == b"String & test"
    for c in columns:
        assert not np.all(mask[c])
    assert np.all(mask['unicode_test'])


class TestParse:
    def setup_class(self):
        self.table = parse(
            get_data_filename('data/regression.xml'),
            pedantic=False).get_first_table()
        self.array = self.table.array
        self.mask = self.table.mask

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
                           [u"Ce\xe7i n'est pas un pipe",
                            u'\u0bb5\u0ba3\u0b95\u0bcd\u0b95\u0bae\u0bcd',
                            u'XXXX', u'', u''])

    def test_fixed_unicode_test(self):
        assert issubclass(self.array['fixed_unicode_test'].dtype.type,
                          np.unicode_)
        assert_array_equal(self.array['fixed_unicode_test'],
                           [u"Ce\xe7i n'est",
                            u'\u0bb5\u0ba3\u0b95\u0bcd\u0b95\u0bae\u0bcd',
                            u'0123456789', u'', u''])

    def test_unsignedByte(self):
        assert issubclass(self.array['unsignedByte'].dtype.type,
                          np.uint8)
        assert_array_equal(self.array['unsignedByte'],
                           [128, 0, 233, 255, 0])
        assert not np.any(self.mask['unsignedByte'])

    def test_short(self):
        assert issubclass(self.array['short'].dtype.type,
                          np.int16)
        assert_array_equal(self.array['short'],
                           [4096, 0, -4096, -1, 0])
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
                           [1.0, 0.0, np.inf, np.nan, -np.inf])
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
        assert self.mask['array'][3][0][0]
        assert self.mask['array'][4][0][1]

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
        self.array['bitarray'].shape == (5, 3, 2)
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
                      True, True]
        for a, b in zip(self.mask['bitvararray'], match_mask):
            assert_array_equal(a, b)

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

    @pytest.mark.xfail('numpy_has_complex_bug')
    def test_floatComplex(self):
        assert issubclass(self.array['floatComplex'].dtype.type,
                          np.complex64)
        assert_array_equal(self.array['floatComplex'],
                           [np.nan+0j, 0+0j, 0+-1j, np.nan+0j, np.nan+0j])
        assert_array_equal(self.mask['floatComplex'],
                           [True, False, False, True, True])

    @pytest.mark.xfail('numpy_has_complex_bug')
    def test_doubleComplex(self):
        assert issubclass(self.array['doubleComplex'].dtype.type,
                          np.complex128)
        assert_array_equal(
            self.array['doubleComplex'],
            [np.nan+0j, 0+0j, 0+-1j, np.nan+(np.inf*1j), np.nan+0j])
        assert_array_equal(self.mask['doubleComplex'],
                           [True, False, False, True, True])

    @pytest.mark.xfail('numpy_has_complex_bug')
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
        assert_array_equal(self.mask['doublearray'][1],
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


class TestThroughTableData(TestParse):
    def setup_class(self):
        votable = parse(
            get_data_filename('data/regression.xml'),
            pedantic=False)
        votable.to_xml(join(TMP_DIR, "test_through_tabledata.xml"))
        self.table = parse(join(TMP_DIR, "test_through_tabledata.xml"),
                           pedantic=False).get_first_table()
        self.array = self.table.array
        self.mask = self.table.mask

    def test_schema(self):
        assert_validate_schema(join(TMP_DIR, "test_through_tabledata.xml"))


class TestThroughBinary(TestParse):
    def setup_class(self):
        votable = parse(
            get_data_filename('data/regression.xml'),
            pedantic=False)
        votable.get_first_table().format = 'binary'
        votable.to_xml(join(TMP_DIR, "test_through_binary.xml"))
        self.table = parse(join(TMP_DIR, "test_through_binary.xml"),
                           pedantic=False).get_first_table()
        self.array = self.table.array
        self.mask = self.table.mask

    # Masked values in bit fields don't roundtrip through the binary
    # representation -- that's not a bug, just a limitation, so
    # override the mask array checks here.
    def test_bit_mask(self):
        assert not np.any(self.mask['bit'])

    def test_bitarray_mask(self):
        assert not np.any(self.mask['bitarray'])

    def test_bit_array2_mask(self):
        assert not np.any(self.mask['bitarray2'])


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

    for filename in get_data_filenames('data', '*.xml'):
        yield test_file, filename


@raises(VOTableSpecError)
def test_too_many_columns():
    votable = parse(
        get_data_filename('data/too_many_columns.xml.gz'),
        pedantic=False)


def test_build_from_scratch():
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
    votable.to_xml(os.path.join(TMP_DIR, "new_votable.xml"))

    votable = parse(os.path.join(TMP_DIR, "new_votable.xml"))

    table = votable.get_first_table()
    assert_array_equal(
        table.mask, np.array([(False, [[False, False], [False, False]]),
                              (False, [[False, False], [False, False]])],
                             dtype=[('filename', '?'),
                                    ('matrix', '?', (2, 2))]))


def test_validate():
    output = io.StringIO()

    # We can't test xmllint, because we can't rely on it being on the
    # user's machine.
    result = validate(get_data_filename('data/regression.xml'),
                      output, xmllint=False)

    assert result == False

    output.seek(0)
    output = output.readlines()

    # Uncomment to generate new groundtruth
    # with io.open('validation.txt', 'wt', encoding='utf-8') as fd:
    #     fd.write(u''.join(output))

    with io.open(
        get_data_filename('data/validation.txt'),
        'rt', encoding='utf-8') as fd:
        truth = fd.readlines()

    truth = truth[1:]
    output = output[1:-1]

    for line in difflib.unified_diff(truth, output):
        if IS_PY3K:
            sys.stdout.write(
                line.replace('\\n', '\n'))
        else:
            sys.stdout.write(
                line.encode('unicode_escape').
                replace('\\n', '\n'))

    assert truth == output


def test_gzip_filehandles():
    votable = parse(
        get_data_filename('data/regression.xml'),
        pedantic=False)

    with open(join(TMP_DIR, "regression.compressed.xml"), 'wb') as fd:
        votable.to_xml(
            fd,
            compressed=True,
            _astropy_version="testing")

    with open(join(TMP_DIR, "regression.compressed.xml"), 'rb') as fd:
        votable = parse(
            fd,
            pedantic=False)


def test_from_scratch_example():
    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        warnings.simplefilter("always", VOWarning, append=True)
        try:
            _run_test_from_scratch_example()
        except ValueError as e:
            lines.append(str(e))

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
