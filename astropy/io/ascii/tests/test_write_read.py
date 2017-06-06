# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy

from ....extern import six
from ....extern.six.moves import cStringIO as StringIO
from ....tests.helper import pytest
from ... import ascii
from ....table import Table

# importing .common sets the current directory
#   so that test files are read from the right place
from .common import setup_function, teardown_function

# We will test that test tables that can be written (test_write.test_defs)
#  can also be read in.
from .test_write import test_defs

debug = False

# BeautifulSoup is required to *read* HTML Files
# Check to see if the BeautifulSoup dependency is present.
try:
    from bs4 import BeautifulSoup, FeatureNotFound
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False


def create_reader_kwargs_from_writer_kwargs(test_write_def):
    """Create a test definition with appropriate kwargs for a Reader."""
    remove_col_for_read = ('formats', 'include_names', 'exclude_names',
                           'strip_whitespace')

    test_read_def = copy.deepcopy(test_write_def)
    if 'Writer' in test_read_def['kwargs']:
        test_read_def['kwargs']['Reader'] = test_write_def['kwargs']['Writer']
        del test_read_def['kwargs']['Writer']
    # Delimiter is a special case.
    # In general, we want to keep the 'delimiter' option unchanged.
    # But in the special case of 'delimiter=None',
    #   we want to remove 'delimiter' from 'kwargs'
    #   because it can confuse the Reader.
    if 'delimiter' in test_read_def['kwargs']:
        if test_read_def['kwargs']['delimiter'] is None:
            del test_read_def['kwargs']['delimiter']
    for name in remove_col_for_read:
        if name in test_read_def['kwargs']:
            del test_read_def['kwargs'][name]

    return test_read_def


def check_read_write_table_via_ascii(test_def, table, fast_writer):
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
        format_name = 'ascii.{0}'.format(
            test_write_def['kwargs']['Writer']._format_name)
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


def writer_or_format_is_html(def_kwargs):
    """Identify a test_def that will write/read HTML.

    If 'Writer': ascii.HTML or 'format': 'ascii.html'
    then return False
    Otherwise return True
    >>> print(writer_or_format_is_html({'Writer': ascii.HTML}))
    True
    >>> print(writer_or_format_is_html({'format': 'ascii.html'}))
    True
    >>> print(writer_or_format_is_html({'Writer': ascii.Csv}))
    False
    """
    if 'Writer' in def_kwargs and def_kwargs['Writer'] == ascii.HTML:
        return True
    if 'format' in def_kwargs and def_kwargs['format'] == 'ascii.html':
        return True

    return False


check_functions = \
    (check_read_write_table_via_ascii, check_read_write_table_via_table)


@pytest.mark.parametrize("test_def", test_defs)
@pytest.mark.parametrize("check_function", check_functions)
@pytest.mark.parametrize("fast_writer", [True, False])
def test_read_write_table(check_function, test_def, fast_writer):
    table = ascii.get_reader(Reader=ascii.Daophot)
    data = table.read('t/daophot.dat')  # A small sample of test data.

    # If BeautifulSoup is not imported.
    # Then skip this test_def if this is an HTML write/read
    if not HAS_BEAUTIFUL_SOUP and writer_or_format_is_html(test_def['kwargs']):
        msg = "BeautifulSoup not installed.  Skipping write+read HTML test."
        pytest.skip(msg)

    check_function(test_def, data, fast_writer)
