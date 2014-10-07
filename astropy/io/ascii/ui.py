# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

ui.py:
  Provides the main user functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import re
import os
import sys

from . import core
from . import basic
from . import cds
from . import daophot
from . import sextractor
from . import ipac
from . import latex
from . import html
from . import fastbasic
from . import cparser

from ...table import Table

# Default setting for guess parameter in read()
_GUESS = True


def set_guess(guess):
    """Set the default value of the ``guess`` parameter for read()

    :param guess: New default ``guess`` value (True|False)
    """
    global _GUESS
    _GUESS = guess


def get_reader(Reader=None, Inputter=None, Outputter=None, **kwargs):
    """Initialize a table reader allowing for common customizations.  Most of the
    default behavior for various parameters is determined by the Reader class.

    :param Reader: Reader class (DEPRECATED) (default= :class:`Basic`)
    :param Inputter: Inputter class
    :param Outputter: Outputter class
    :param delimiter: column delimiter string
    :param comment: regular expression defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param header_start: line index for the header line not counting comment lines
    :param data_start: line index for the start of data not counting comment lines
    :param data_end: line index for the end of data (can be negative to count from end)
    :param converters: dict of converters
    :param data_Splitter: Splitter class to split data columns
    :param header_Splitter: Splitter class to split header columns
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    :param fill_values: specification of fill values for bad or missing table values
    :param fill_include_names: list of names to include in fill_values (default=None selects all names)
    :param fill_exclude_names: list of names to exlude from fill_values (applied after ``fill_include_names``)
    """
    # This function is a light wrapper around core._get_reader to provide a public interface
    # with a default Reader.
    if Reader is None:
        Reader = basic.Basic
    reader = core._get_reader(Reader, Inputter=Inputter, Outputter=Outputter, **kwargs)
    return reader


def _get_format_class(format, ReaderWriter, label):
    if format is not None and ReaderWriter is not None:
        raise ValueError('Cannot supply both format and {0} keywords'.format(label))

    if format is not None:
        if format in core.FORMAT_CLASSES:
            ReaderWriter = core.FORMAT_CLASSES[format]
        else:
            raise ValueError('ASCII format {0!r} not in allowed list {1}'
                             .format(format, sorted(core.FORMAT_CLASSES)))
    return ReaderWriter


def read(table, guess=None, **kwargs):
    """Read the input ``table`` and return the table.  Most of
    the default behavior for various parameters is determined by the Reader
    class.

    :param table: input table (file name, file-like object, list of strings, or single newline-separated string)
    :param guess: try to guess the table format (default=True)
    :param format: input table format
    :param Inputter: Inputter class
    :param Outputter: Outputter class (default=TableOutputter)
    :param delimiter: column delimiter string
    :param comment: regular expression defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param header_start: line index for the header line not counting comment lines
    :param data_start: line index for the start of data not counting comment lines
    :param data_end: line index for the end of data (can be negative to count from end)
    :param converters: dict of converters
    :param data_Splitter: Splitter class to split data columns
    :param header_Splitter: Splitter class to split header columns
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    :param fill_values: specification of fill values for bad or missing table values (default=('', '0'))
    :param fill_include_names: list of names to include in fill_values (default=None selects all names)
    :param fill_exclude_names: list of names to exlude from fill_values (applied after ``fill_include_names``)
    :param fast_reader: whether to use the C engine, can also be a dict with options which default to False (default=True)
    :param Reader: Reader class (DEPRECATED) (default=``ascii.Basic``)
    """

    if 'fill_values' not in kwargs:
        kwargs['fill_values'] = [('', '0')]

    # If an Outputter is supplied in kwargs that will take precedence.
    new_kwargs = {}
    new_kwargs['Outputter'] = core.TableOutputter
    fast_reader_param = kwargs.get('fast_reader', True)
    if 'Outputter' in kwargs: # user specified Outputter, not supported for fast reading
        fast_reader_param = False
    format = kwargs.get('format')
    new_kwargs.update(kwargs)

    # Get the Reader class based on possible format and Reader kwarg inputs.
    Reader = _get_format_class(format, kwargs.get('Reader'), 'Reader')
    if Reader is not None:
        new_kwargs['Reader'] = Reader
        format = Reader._format_name

    # Remove format keyword if there, this is only allowed in read() not get_reader()
    if 'format' in new_kwargs:
        del new_kwargs['format']

    if guess is None:
        guess = _GUESS
    if guess:
        dat = _guess(table, new_kwargs, format, fast_reader_param)
    else:
        reader = get_reader(**new_kwargs)
        # Try the fast reader first if applicable
        if fast_reader_param and format is not None and 'fast_{0}'.format(format) \
                                                        in core.FAST_CLASSES:
            new_kwargs['Reader'] = core.FAST_CLASSES['fast_{0}'.format(format)]
            fast_reader = get_reader(**new_kwargs)
            try:
                return fast_reader.read(table)
            except (core.ParameterError, cparser.CParserError) as e:
                # special testing value to avoid falling back on the slow reader
                if fast_reader_param == 'force':
                    raise e
                # If the fast reader doesn't work, try the slow version
                dat = reader.read(table)
        else:
            dat = reader.read(table)

    return dat


def _guess(table, read_kwargs, format, fast_reader):
    """Try to read the table using various sets of keyword args. First try the
    original args supplied in the read() call. Then try the standard guess
    keyword args. For each key/val pair specified explicitly in the read()
    call make sure that if there is a corresponding definition in the guess
    then it must have the same val.  If not then skip this guess."""

    # Keep a trace of all failed guesses kwarg
    failed_kwargs = []
    fast_kwargs = []

    first_kwargs = [read_kwargs.copy()]
    if fast_reader and format is not None and 'fast_{0}'.format(format) in \
                                                         core.FAST_CLASSES:
        # If a fast version of the reader is available, try that before the slow version
        fast_kwargs = read_kwargs.copy()
        fast_kwargs['Reader'] = core.FAST_CLASSES['fast_{0}'.format(format)]
        first_kwargs = [fast_kwargs] + first_kwargs

    # First try guessing
    for guess_kwargs in first_kwargs + _get_guess_kwargs_list():
        guess_kwargs_ok = True  # guess_kwargs are consistent with user_kwargs?
        for key, val in read_kwargs.items():
            # Do guess_kwargs.update(read_kwargs) except that if guess_args has
            # a conflicting key/val pair then skip this guess entirely.
            if key not in guess_kwargs:
                guess_kwargs[key] = val
            elif val != guess_kwargs[key] and guess_kwargs != fast_kwargs:
                guess_kwargs_ok = False
                break

        if not guess_kwargs_ok:
            # User-supplied kwarg is inconsistent with the guess-supplied kwarg, e.g.
            # user supplies delimiter="|" but the guess wants to try delimiter=" ",
            # so skip the guess entirely.
            continue

        try:
            # If guessing will try all Readers then use strict req'ts on column names
            if 'Reader' not in read_kwargs:
                guess_kwargs['strict_names'] = True

            reader = get_reader(**guess_kwargs)
            dat = reader.read(table)

            # When guessing require at least two columns
            if len(dat.colnames) <= 1:
                del dat
                raise ValueError

            return dat

        except (core.InconsistentTableError, ValueError, TypeError,
                core.OptionalTableImportError, core.ParameterError, cparser.CParserError) as e:
            failed_kwargs.append(guess_kwargs)
    else:
        # failed all guesses, try the original read_kwargs without column requirements
        try:
            reader = get_reader(**read_kwargs)
            return reader.read(table)
        except (core.InconsistentTableError, ValueError, ImportError,
                core.OptionalTableImportError, core.ParameterError, cparser.CParserError):
            failed_kwargs.append(read_kwargs)
            lines = ['\nERROR: Unable to guess table format with the guesses listed below:']
            for kwargs in failed_kwargs:
                sorted_keys = sorted([x for x in sorted(kwargs)
                                      if x not in ('Reader', 'Outputter')])
                reader_repr = repr(kwargs.get('Reader', basic.Basic))
                keys_vals = ['Reader:' + re.search(r"\.(\w+)'>", reader_repr).group(1)]
                kwargs_sorted = ((key, kwargs[key]) for key in sorted_keys)
                keys_vals.extend(['%s: %s' % (key, repr(val)) for key, val in kwargs_sorted])
                lines.append(' '.join(keys_vals))
            lines.append('ERROR: Unable to guess table format with the guesses listed above.')
            lines.append('Check the table and try with guess=False '
                         'and appropriate arguments to read()')
            raise core.InconsistentTableError('\n'.join(lines))


def _get_guess_kwargs_list():
    guess_kwargs_list = [dict(Reader=basic.Rdb),
                         dict(Reader=fastbasic.FastTab),
                         dict(Reader=basic.Tab),
                         dict(Reader=cds.Cds),
                         dict(Reader=daophot.Daophot),
                         dict(Reader=sextractor.SExtractor),
                         dict(Reader=ipac.Ipac),
                         dict(Reader=latex.Latex),
                         dict(Reader=latex.AASTex),
                         dict(Reader=html.HTML)
                         ]
    for Reader in (basic.CommentedHeader, fastbasic.FastBasic, basic.Basic,
                   fastbasic.FastNoHeader, basic.NoHeader):
        for delimiter in ("|", ",", " ", "\s"):
            for quotechar in ('"', "'"):
                guess_kwargs_list.append(dict(
                    Reader=Reader, delimiter=delimiter, quotechar=quotechar))
    return guess_kwargs_list

extra_writer_pars = ('delimiter', 'comment', 'quotechar', 'formats',
                     'names', 'include_names', 'exclude_names', 'strip_whitespace')


def get_writer(Writer=None, fast_writer=True, **kwargs):
    """Initialize a table writer allowing for common customizations.  Most of the
    default behavior for various parameters is determined by the Writer class.

    :param Writer: Writer class (DEPRECATED) (default=``ascii.Basic``)
    :param delimiter: column delimiter string
    :param write_comment: string defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param formats: dict of format specifiers or formatting functions
    :param strip_whitespace: strip surrounding whitespace from column values (default=True)
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    :param fast_writer: whether to use the fast Cython writer (default=True)
    """
    if Writer is None:
        Writer = basic.Basic
    if 'strip_whitespace' not in kwargs:
        kwargs['strip_whitespace'] = True
    writer = core._get_writer(Writer, fast_writer, **kwargs)
    return writer


def write(table, output=None,  format=None, Writer=None, fast_writer=True, **kwargs):
    """Write the input ``table`` to ``filename``.  Most of the default behavior
    for various parameters is determined by the Writer class.

    :param table: input table (Reader object, NumPy struct array, list of lists, etc)
    :param output: output [filename, file-like object] (default = sys.stdout)
    :param format: output format (default=``basic``)
    :param delimiter: column delimiter string
    :param write_comment: string defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param formats: dict of format specifiers or formatting functions
    :param strip_whitespace: strip surrounding whitespace from column values (default=True)
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    :param fast_writer: whether to use the fast Cython writer (default=True)
    :param Writer: Writer class (DEPRECATED) (default=``ascii.Basic``)
    """
    if output is None:
        output = sys.stdout

    table = Table(table, names=kwargs.get('names'))

    Writer = _get_format_class(format, Writer, 'Writer')
    writer = get_writer(Writer=Writer, fast_writer=fast_writer, **kwargs)
    if writer._format_name in core.FAST_CLASSES:
        writer.write(table, output)
        return

    lines = writer.write(table)

    # Write the lines to output
    outstr = os.linesep.join(lines)
    if not hasattr(output, 'write'):
        output = open(output, 'w')
        output.write(outstr)
        output.write(os.linesep)
        output.close()
    else:
        output.write(outstr)
        output.write(os.linesep)
