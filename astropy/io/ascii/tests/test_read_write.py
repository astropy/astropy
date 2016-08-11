# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy

from ....extern.six.moves import cStringIO as StringIO
from ....tests.helper import pytest
from ... import ascii
from ....table import Table

# importing .common sets the current directory
#   so that test files are read from the right place
from .common import setup_function, teardown_function

from .test_write import test_defs

debug = False


def create_reader_kwargs_from_writer_kwargs(test_write_def):
    """Create a test definition with appropriate kwargs for a Reader."""

    remove_col_for_read = ('formats', 'include_names', 'exclude_names', 'strip_whitespace')

    test_read_def = copy.deepcopy(test_write_def)
    if 'Writer' in test_read_def['kwargs']:
        test_read_def['kwargs']['Reader'] = test_write_def['kwargs']['Writer']
        del test_read_def['kwargs']['Writer']
    # Delimiter is a special case.  In general, we want to keep the delimiter as set
    # But in the special case of 'delimiter=None', we want to remove 'delimiter' from 'kwargs'.
    if 'delimiter' in test_read_def['kwargs']:
        if test_read_def['kwargs']['delimiter'] is None:
            del test_read_def['kwargs']['delimiter']
    for name in remove_col_for_read:
        if name in test_read_def['kwargs']:
            del test_read_def['kwargs'][name]

    return test_read_def


def check_read_write_table(test_def, table, fast_writer):
    """Check reading and writing using the astropy.io.ascii interface."""
    in_out = StringIO()

    test_read_def = create_reader_kwargs_from_writer_kwargs(test_def)
    try:
        if debug:
            print(test_def['kwargs'])
            print(test_read_def['kwargs'])
        ascii.write(table, in_out, fast_writer=fast_writer,
                    **test_def['kwargs'])
        foo = ascii.read(in_out, **test_read_def['kwargs'])
    except ValueError as e:  # if format doesn't have a fast writer, ignore
        if 'not in the list of formats with fast writers' not in str(e):
            raise e
        return
    if debug:
        print('Expected:\n%s' % test_def['out'])
        print('Actual:\n%s' % in_out.getvalue())
        print('Read:\n%s' % str(foo))
    assert [x.strip() for x in in_out.getvalue().strip().splitlines()] == \
           [x.strip() for x in test_def['out'].strip().splitlines()]


def check_read_write_table_via_table(test_def, table, fast_writer):
    """Check reading and writing using the astropy.table.Table interface."""
    in_out = StringIO()

    test_write_def = copy.deepcopy(test_def)
    # Table would like 'format' text string instead of 'Writer' and 'Reader'.
    if 'Writer' in test_def['kwargs']:
        format_name = 'ascii.{0}'.format(test_write_def['kwargs']['Writer']._format_name)
        test_write_def['kwargs']['format'] = format_name
        del test_write_def['kwargs']['Writer']

    test_read_def = create_reader_kwargs_from_writer_kwargs(test_write_def)

    try:
        table.write(in_out, fast_writer=fast_writer,
                    **test_write_def['kwargs'])
        foo = Table.read(in_out, **test_read_def['kwargs'])
    except ValueError as e:  # if format doesn't have a fast writer, ignore
        if 'not in the list of formats with fast writers' not in str(e):
            raise e
        return
    if debug:
        print('Expected:\n%s' % test_def['out'])
        print('Actual:\n%s' % in_out.getvalue())
        print('Read:\n%s' % str(foo))
    assert [x.strip() for x in in_out.getvalue().strip().splitlines()] == \
           [x.strip() for x in test_def['out'].strip().splitlines()]


@pytest.mark.parametrize("fast_writer", [True, False])
def test_write_table(fast_writer):
    table = ascii.get_reader(Reader=ascii.Daophot)
    data = table.read('t/daophot.dat')

    for test_def in test_defs:
        check_read_write_table(test_def, data, fast_writer)
        check_read_write_table_via_table(test_def, data, fast_writer)
