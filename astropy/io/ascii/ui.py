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
from . import ecsv
from . import sextractor
from . import ipac
from . import latex
from . import html
from . import fastbasic
from . import cparser
from . import fixedwidth

from ...table import Table
from ...utils.data import get_readable_fileobj

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# Default setting for guess parameter in read()
_GUESS = True


def set_guess(guess):
    """
    Set the default value of the ``guess`` parameter for read()

    Parameters
    ----------
    guess : bool
        New default ``guess`` value (e.g., True or False)

    """
    global _GUESS
    _GUESS = guess


def get_reader(Reader=None, Inputter=None, Outputter=None, **kwargs):
    """
    Initialize a table reader allowing for common customizations.  Most of the
    default behavior for various parameters is determined by the Reader class.

    Parameters
    ----------
    Reader : `~astropy.io.ascii.BaseReader`
        Reader class (DEPRECATED) (default= :class:`Basic`)
    Inputter : `~astropy.io.ascii.BaseInputter`
        Inputter class
    Outputter : `~astropy.io.ascii.BaseOutputter`
        Outputter class
    delimiter : str
        Column delimiter string
    comment : str
        Regular expression defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    header_start : int
        Line index for the header line not counting comment lines
    data_start : int
        Line index for the start of data not counting comment lines
    data_end : int
        Line index for the end of data (can be negative to count from end)
    converters : dict
        Dictionary of converters
    data_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split data columns
    header_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split header columns
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output (default= ``None`` selects all names)
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fill_values : dict
        specification of fill values for bad or missing table values
    fill_include_names : list
        List of names to include in fill_values (default= ``None`` selects all names)
    fill_exclude_names : list
        List of names to exclude from fill_values (applied after ``fill_include_names``)

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
    """
    Read the input ``table`` and return the table.  Most of
    the default behavior for various parameters is determined by the Reader
    class.

    Parameters
    ----------
    table : str, file-like, list
        Input table as a file name, file-like object, list of strings, or
        single newline-separated string.
    guess : bool
        Try to guess the table format (default= ``True``)
    format : str, `~astropy.io.ascii.BaseReader`
        Input table format
    Inputter : `~astropy.io.ascii.BaseInputter`
        Inputter class
    Outputter : `~astropy.io.ascii.BaseOutputter`
        Outputter class
    delimiter : str
        Column delimiter string
    comment : str
        Regular expression defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    header_start : int
        Line index for the header line not counting comment lines
    data_start : int
        Line index for the start of data not counting comment lines
    data_end : int
        Line index for the end of data (can be negative to count from end)
    converters : dict
        Dictionary of converters
    data_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split data columns
    header_Splitter : `~astropy.io.ascii.BaseSplitter`
        Splitter class to split header columns
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output (default= ``None`` selects all names)
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fill_values : dict
        specification of fill values for bad or missing table values
    fill_include_names : list
        List of names to include in fill_values (default= ``None`` selects all names)
    fill_exclude_names : list
        List of names to exclude from fill_values (applied after ``fill_include_names``)
    fast_reader : bool
        Whether to use the C engine, can also be a dict with options which default to ``False``
        (default= ``True``)
    Reader : `~astropy.io.ascii.BaseReader`
        Reader class (DEPRECATED) (default= :class:`Basic`).
    """

    if 'fill_values' not in kwargs:
        kwargs['fill_values'] = [('', '0')]

    # If an Outputter is supplied in kwargs that will take precedence.
    new_kwargs = {}
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

    # If `table` is a readable file object then read in the file now.  This
    # prevents problems in Python 3 with the file object getting closed or
    # left at the file end.  See #3132, #3013, #3109, #2001.  If a `readme`
    # arg was passed that implies CDS format, in which case the original
    # `table` as the data filename must be left intact.
    if 'readme' not in read_kwargs:
        try:
            with get_readable_fileobj(table) as fileobj:
                table = fileobj.read()
        except:
            pass

    # Keep a trace of all failed guesses kwarg
    failed_kwargs = []
    full_list_guess = _get_guess_kwargs_list(read_kwargs)

    if fast_reader and format is not None and 'fast_{0}'.format(format) in \
                                                         core.FAST_CLASSES:
        # If a fast version of the reader is available, try that before the slow version
        fast_kwargs = read_kwargs.copy()
        fast_kwargs['Reader'] = core.FAST_CLASSES['fast_{0}'.format(format)]
        full_list_guess = [fast_kwargs] + full_list_guess
    else:
        fast_kwargs = None

    # First try guessing
    for guess_kwargs in full_list_guess:
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
            reader.guessing = True
            return reader.read(table)

        except (core.InconsistentTableError, ValueError, TypeError, AttributeError,
                core.OptionalTableImportError, core.ParameterError, cparser.CParserError):
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

            msg = ['',
                   '************************************************************************',
                   '** ERROR: Unable to guess table format with the guesses listed above. **',
                   '**                                                                    **',
                   '** To figure out why the table did not read, use guess=False and      **',
                   '** appropriate arguments to read().  In particular specify the format **',
                   '** and any known attributes like the delimiter.                       **',
                   '************************************************************************']
            lines.extend(msg)
            raise core.InconsistentTableError('\n'.join(lines))

def _get_guess_kwargs_list(read_kwargs):
    guess_kwargs_list = []

    # Start with ECSV because an ECSV file will be read by Basic.  This format
    # has very specific header requirements and fails out quickly.
    if HAS_YAML:
        guess_kwargs_list.append(dict(Reader=ecsv.Ecsv))

    # Now try readers that accept the common arguments with the input arguments
    # (Unless there are not arguments - we try that in the next step anyway.)
    # FixedWidthTwoLine would also be read by Basic, so it needs to come first.
    if len(read_kwargs) > 0:
        for reader in [fixedwidth.FixedWidthTwoLine,
                       basic.Basic]:
            first_kwargs = read_kwargs.copy()
            first_kwargs.update(dict(Reader=reader))
            guess_kwargs_list.append(first_kwargs)

    # Then try a list of readers with default arguments
    guess_kwargs_list.extend([dict(Reader=fixedwidth.FixedWidthTwoLine),
                              dict(Reader=fastbasic.FastBasic),
                              dict(Reader=basic.Basic),
                              dict(Reader=basic.Rdb),
                              dict(Reader=fastbasic.FastTab),
                              dict(Reader=basic.Tab),
                              dict(Reader=cds.Cds),
                              dict(Reader=daophot.Daophot),
                              dict(Reader=sextractor.SExtractor),
                              dict(Reader=ipac.Ipac),
                              dict(Reader=latex.Latex),
                              dict(Reader=latex.AASTex),
                              dict(Reader=html.HTML)
                              ])

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
    """
    Initialize a table writer allowing for common customizations.  Most of the
    default behavior for various parameters is determined by the Writer class.

    Parameters
    ----------
    Writer : ``Writer``
        Writer class (DEPRECATED) (default= :class:`Basic`)
    delimiter : str
        Column delimiter string
    write_comment : str
        String defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    formats : dict
        Dictionary of format specifiers or formatting functions
    strip_whitespace : bool
        Strip surrounding whitespace from column values (default= ``True``)
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output (default= ``None`` selects all names)
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fast_writer : bool
        Whether to use the fast Cython writer (default= ``True``)

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

    Parameters
    ----------
    table : `~astropy.io.ascii.BaseReader`, array_like, str, file_like, list
        Input table as a Reader object, Numpy struct array, file name,
        file-like object, list of strings, or single newline-separated string.
    output : str, file_like
        Output [filename, file-like object] (default = ``sys.stdout``)
    format : str
        Output table format (default= ``basic``)
    delimiter : str
        Column delimiter string
    write_comment : str
        String defining a comment line in table
    quotechar : str
        One-character string to quote fields containing special characters
    formats : dict
        Dictionary of format specifiers or formatting functions
    strip_whitespace : bool
        Strip surrounding whitespace from column values (default= ``True``)
    names : list
        List of names corresponding to each data column
    include_names : list
        List of names to include in output (default= ``None`` selects all names)
    exclude_names : list
        List of names to exclude from output (applied after ``include_names``)
    fast_writer : bool
        Whether to use the fast Cython writer (default= ``True``)
    Writer : ``Writer``
        Writer class (DEPRECATED) (default= :class:`Basic`)
    """
    if output is None:
        output = sys.stdout

    table = Table(table, names=kwargs.get('names'))

    if table.has_mixin_columns:
        fast_writer = False

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
