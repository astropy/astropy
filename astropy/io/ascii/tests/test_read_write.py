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

test_write_defs = copy.deepcopy(test_defs)
test_read_defs = copy.deepcopy(test_defs)
read_kwarg_names = ('format')

for td in test_read_defs:
    if 'Writer' in td['kwargs']:
        td = {'Reader': td['kwargs']['Writer']}


def check_read_write_table(test_def, table, fast_writer):
    in_out = StringIO()

    test_read_def = {}
    if 'Writer' in test_def['kwargs']:
        test_read_def['kwargs'] = {'Reader': test_def['kwargs']['Writer']}
    else:
        test_read_def['kwargs'] = {}

#    if 'kwargs' in test_def and 'format' in test_def['kwargs']:
#        # astropy.io.ascii expects 'basic',
#        # while astropy.table expects 'ascii.basic'
#        test_def['kwargs']['format'] = test_def['kwargs']['format'].replace('ascii.', '')

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
    in_out = StringIO()

    test_write_def = copy.deepcopy(test_def)
    if 'Writer' in test_def['kwargs']:
        format_name = 'ascii.{0}'.format(test_write_def['kwargs']['Writer']._format_name)
        test_write_def['kwargs']['format'] = format_name
        test_read_def = {'kwargs': {'format': format_name}}

        del test_write_def['kwargs']['Writer']

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
