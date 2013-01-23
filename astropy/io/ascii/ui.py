# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

ui.py:
  Provides the main user functions for reading and writing tables.

:Copyright: Smithsonian Astrophysical Observatory (2010)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

## Copyright (c) 2010, Smithsonian Astrophysical Observatory
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import re
import os
import sys

from . import core
from . import basic
from . import cds
from . import daophot
from . import sextractor
from . import ipac
from .core import next, izip, any
from . import latex

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

    :param Reader: Reader class (default= :class:`Basic`)
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

def read(table, guess=None, **kwargs):
    """Read the input ``table`` and return the table.  Most of
    the default behavior for various parameters is determined by the Reader
    class.

    :param table: input table (file name, list of strings, or single newline-separated string)
    :param guess: try to guess the table format (default=True)
    :param Reader: Reader class (default=``ascii.Basic``)
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

    # Provide a simple way to choose between the two common outputters.  If an
    # Outputter is supplied in kwargs that will take precedence.
    new_kwargs = {}
    new_kwargs['Outputter'] = core.TableOutputter
    new_kwargs.update(kwargs)

    if guess is None:
        guess = _GUESS
    if guess:
        dat = _guess(table, new_kwargs)
    else:
        reader = get_reader(**new_kwargs)
        dat = reader.read(table)
    return dat


def _is_number(x):
    try:
        x = float(x)
        return True
    except ValueError:
        pass
    return False


def _guess(table, read_kwargs):
    """Try to read the table using various sets of keyword args. First try the
    original args supplied in the read() call. Then try the standard guess
    keyword args. For each key/val pair specified explicitly in the read()
    call make sure that if there is a corresponding definition in the guess
    then it must have the same val.  If not then skip this guess."""

    # Keep a trace of all failed guesses kwarg
    failed_kwargs = []

    # First try guessing
    for guess_kwargs in [read_kwargs.copy()] + _get_guess_kwargs_list():
        guess_kwargs_ok = True  # guess_kwargs are consistent with user_kwargs?
        for key, val in read_kwargs.items():
            # Do guess_kwargs.update(read_kwargs) except that if guess_args has
            # a conflicting key/val pair then skip this guess entirely.
            if key not in guess_kwargs:
                guess_kwargs[key] = val
            elif val != guess_kwargs[key]:
                guess_kwargs_ok = False
                break

        if not guess_kwargs_ok:
            # User-supplied kwarg is inconsistent with the guess-supplied kwarg, e.g.
            # user supplies delimiter="|" but the guess wants to try delimiter=" ", 
            # so skip the guess entirely.
            continue

        try:
            reader = get_reader(**guess_kwargs)
            dat = reader.read(table)
            # When guessing impose additional requirements on column names and number of cols
            bads = [" ", ",", "|", "\t", "'", '"']
            if (len(reader.cols) <= 1 or
                any(_is_number(col.name) or 
                     len(col.name) == 0 or 
                     col.name[0] in bads or 
                     col.name[-1] in bads for col in reader.cols)):
                raise ValueError
            return dat
        except (core.InconsistentTableError, ValueError, TypeError):
            failed_kwargs.append(guess_kwargs)
            pass
    else:
        # failed all guesses, try the original read_kwargs without column requirements
        try:
            reader = get_reader(**read_kwargs)
            return reader.read(table)
        except (core.InconsistentTableError, ValueError):
            failed_kwargs.append(read_kwargs)
            lines = ['\nERROR: Unable to guess table for with the guesses listed below:']
            for kwargs in failed_kwargs:
                sorted_keys = sorted([x for x in sorted(kwargs) if x not in ('Reader', 'Outputter')])
                reader_repr = repr(kwargs.get('Reader', basic.Basic))
                keys_vals = ['Reader:' + re.search(r"\.(\w+)'>", reader_repr).group(1)]
                kwargs_sorted = ((key, kwargs[key]) for key in sorted_keys)
                keys_vals.extend(['%s: %s' % (key, repr(val)) for key, val in kwargs_sorted])
                lines.append(' '.join(keys_vals))
            lines.append('ERROR: Unable to guess table for with the guesses listed above.')
            lines.append('Check the table and try with guess=False and appropriate arguments to read()')
            raise core.InconsistentTableError('\n'.join(lines))
    
def _get_guess_kwargs_list():
    guess_kwargs_list = [dict(Reader=basic.Rdb),
                         dict(Reader=basic.Tab),
                         dict(Reader=cds.Cds),
                         dict(Reader=daophot.Daophot),
                         dict(Reader=sextractor.SExtractor),
                         dict(Reader=ipac.Ipac),
                         dict(Reader=latex.Latex),
                         dict(Reader=latex.AASTex)
                         ]
    for Reader in (basic.CommentedHeader, basic.Basic, basic.NoHeader):
        for delimiter in ("|", ",", " ", "\s"):
            for quotechar in ('"', "'"):
                guess_kwargs_list.append(dict(
                    Reader=Reader, delimiter=delimiter, quotechar=quotechar))
    return guess_kwargs_list

extra_writer_pars = ('delimiter', 'comment', 'quotechar', 'formats',
                     'names', 'include_names', 'exclude_names', 'strip_whitespace')

def get_writer(Writer=None, **kwargs):
    """Initialize a table writer allowing for common customizations.  Most of the
    default behavior for various parameters is determined by the Writer class.

    :param Writer: Writer class (default=``ascii.Basic``)
    :param delimiter: column delimiter string
    :param write_comment: string defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param formats: dict of format specifiers or formatting functions
    :param strip_whitespace: strip surrounding whitespace from column values (default=True)
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    """
    if Writer is None:
        Writer = basic.Basic
    if 'strip_whitespace' not in kwargs:
        kwargs['strip_whitespace'] = True
    writer = core._get_writer(Writer, **kwargs)
    return writer

def write(table, output=sys.stdout,  Writer=None, **kwargs):
    """Write the input ``table`` to ``filename``.  Most of the default behavior
    for various parameters is determined by the Writer class.

    :param table: input table (Reader object, NumPy struct array, list of lists, etc)
    :param output: output [filename, file-like object] (default = sys.stdout)
    :param Writer: Writer class (default=``ascii.Basic``)
    :param delimiter: column delimiter string
    :param write_comment: string defining a comment line in table
    :param quotechar: one-character string to quote fields containing special characters
    :param formats: dict of format specifiers or formatting functions
    :param strip_whitespace: strip surrounding whitespace from column values (default=True)
    :param names: list of names corresponding to each data column
    :param include_names: list of names to include in output (default=None selects all names)
    :param exclude_names: list of names to exlude from output (applied after ``include_names``)
    """

    table = Table(table, names=kwargs.get('names'))

    names = set(table.colnames)
    if 'include_names' in kwargs:
        names.intersection_update(kwargs['include_names'])
    if 'exclude_names' in kwargs:
        names.difference_update(kwargs['exclude_names'])
    if names != set(table.colnames):
        remove_names = set(table.colnames) - set(names)
        table.remove_columns(remove_names)

    table.cols = table.columns.values()

    writer = get_writer(Writer=Writer, **kwargs)
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

