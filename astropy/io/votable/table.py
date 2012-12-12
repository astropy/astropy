"""
This file contains a contains the high-level functions to read a
VOTable file.
"""

from __future__ import division, absolute_import

# STDLIB
import io
import os
import sys
import warnings

# LOCAL
from . import exceptions
from . import tree
from ...utils.xml import iterparser
from ...utils import data
from ...config import ConfigurationItem


__all__ = ['parse', 'parse_single_table', 'from_table', 'writeto', 'validate']


PEDANTIC = ConfigurationItem(
    'pedantic',
    False,
    'When True, treat fixable violations of the VOTable spec as exceptions.')


def parse(source, columns=None, invalid='exception', pedantic=None,
          chunk_size=tree.DEFAULT_CHUNK_SIZE, table_number=None,
          filename=None,
          _debug_python_based_parser=False):
    """
    Parses a VOTABLE_ xml file (or file-like object), and returns a
    `~astropy.io.votable.tree.VOTable` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a VOTABLE_ xml file.

    columns : sequence of str, optional
        List of field names to include in the output.  The default is
        to include all fields.

    invalid : str, optional
        One of the following values:

            - 'exception': throw an exception when an invalid value is
              encountered (default)

            - 'mask': mask out invalid values

    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        When not provided, uses the configuration setting
        `astropy.io.votable.pedantic`, which defaults to True.

    chunk_size : int, optional
        The number of rows to read before converting to an array.
        Higher numbers are likely to be faster, but will consume more
        memory.

    table_number : int, optional
        The number of table in the file to read in.  If `None`, all
        tables will be read.  If a number, 0 refers to the first table
        in the file, and only that numbered table will be parsed and
        read in.

    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    votable : `astropy.io.votable.tree.VOTableFile` object

    See also
    --------
    astropy.io.votable.exceptions : The exceptions this function may raise.
    """
    invalid = invalid.lower()
    assert invalid in ('exception', 'mask')

    if pedantic is None:
        pedantic = PEDANTIC()

    config = {
        'columns'      :      columns,
        'invalid'      :      invalid,
        'pedantic'     :     pedantic,
        'chunk_size'   :   chunk_size,
        'table_number' : table_number,
        'filename'     :     filename}

    if filename is None and isinstance(source, basestring):
        config['filename'] = source

    with iterparser.get_xml_iterator(
        source,
        _debug_python_based_parser=_debug_python_based_parser) as iterator:
        return tree.VOTableFile(
            config=config, pos=(1, 1)).parse(iterator, config)


def parse_single_table(source, **kwargs):
    """
    Parses a VOTABLE_ xml file (or file-like object), reading and
    returning only the first `~astropy.io.votable.tree.Table`
    instance.

    See `parse` for a description of the keyword arguments.

    Returns
    -------
    votable : `astropy.io.votable.tree.Table` object
    """
    if kwargs.get('table_number') is None:
        kwargs['table_number'] = 0

    votable = parse(source, **kwargs)

    return votable.get_first_table()


def writeto(table, file):
    """
    Writes a `astropy.io.vo.VOTableFile` to a VOTABLE_ xml file.

    Parameters
    ----------
    table : `astropy.io.vo.VOTableFile` or `astropy.table.Table` instance.

    file : str or writable file-like object
        Path or file object to write to
    """
    from ...table import Table
    if isinstance(table, Table):
        table = tree.VOTableFile.from_table(table)
    elif not isinstance(table, tree.VOTableFile):
        raise TypeError(
            "first argument must be astropy.io.vo.VOTableFile or "
            "astropy.table.Table instance")
    table.to_xml(file, _debug_python_based_parser=True)


def validate(source, output=sys.stdout, xmllint=False, filename=None):
    """
    Prints a validation report for the given file.

    Parameters
    ----------
    source : str or readable file-like object
        Path to a VOTABLE_ xml file.

    output : writable file-like object, optional
        Where to output the report.  Defaults to `sys.stdout`.
        If `None`, the output will be returned as a string.

    xmllint : bool, optional
        When `True`, also send the file to `xmllint` for schema and
        DTD validation.  Requires that `xmllint` is installed.  The
        default is `False`.  `source` must be a file on the local
        filesystem in order for `xmllint` to work.

    filename : str, optional
        A filename to use in the error messages.  If not provided, one
        will be automatically determined from ``source``.

    Returns
    -------
    is_valid : bool or str
        Returns `True` if no warnings were found.  If `output` is
        `None`, the return value will be a string.
    """
    import textwrap
    from . import converters, xmlutil
    from ...utils.console import print_code_line, color_print

    return_as_str = False
    if output is None:
        output = io.StringIO()

    lines = []
    votable = None

    # This is a special variable used by the Python warnings
    # infrastructure to keep track of warnings that have already been
    # seen.  Since we want to get every single warning out of this, we
    # have to delete all of them first.
    for module in (exceptions, converters, tree, xmlutil):
        if hasattr(module, '__warningregistry__'):
            del module.__warningregistry__

    with data.get_readable_fileobj(source, encoding='binary') as fd:
        content = fd.read()
    content_buffer = io.BytesIO(content)
    content_buffer.seek(0)

    if filename is None:
        if isinstance(source, basestring):
            filename = source
        elif hasattr(source, 'name'):
            filename = source.name
        elif hasattr(source, 'url'):
            filename = source.url
        else:
            filename = "<unknown>"

    with warnings.catch_warnings(record=True) as warning_lines:
        warnings.resetwarnings()
        warnings.simplefilter("always", exceptions.VOWarning, append=True)
        try:
            votable = parse(content_buffer, pedantic=False, filename=filename)
        except ValueError as e:
            lines.append(str(e))
    lines = [str(x.message) for x in warning_lines] + lines

    content_buffer.seek(0)
    output.write(u"Validation report for {0}\n\n".format(filename))

    if len(lines):
        xml_lines = iterparser.xml_readlines(content_buffer)

        for warning in lines:
            w = exceptions.parse_vowarning(warning)

            if not w['is_something']:
                output.write(warning)
                output.write(u'\n\n')
            else:
                line = xml_lines[w['nline'] - 1]
                warning = w['warning']
                if w['is_warning']:
                    color = 'yellow'
                else:
                    color = 'red'
                color_print(
                    u'{0:d}: '.format(w['nline']), '',
                    warning or 'EXC', color,
                    u': ', '',
                    textwrap.fill(
                        w['message'],
                        initial_indent=u'          ',
                        subsequent_indent=u'  ').lstrip(),
                    file=output)
                print_code_line(line, w['nchar'], file=output)
            output.write(u'\n')
    else:
        output.write(u'astropy.io.votable found no violations.\n\n')

    success = 0
    if xmllint and os.path.exists(filename):
        from ...utils.xml import validate

        if votable is None:
            version = "1.1"
        else:
            version = votable.version
        success, stdout, stderr = validate.validate_schema(
            filename, version)

        if success != 0:
            output.write(
                u'xmllint schema violations:\n\n')
            output.write(stderr)
        else:
            output.write(u'xmllint passed\n')

    if return_as_str:
        return output.getvalue()
    return len(lines) == 0 and success == 0


def from_table(table):
    """
    Given an `astropy.table.Table` object, return a
    `~astropy.io.votable.tree.VOTableFile` file structure containing
    just that single table.

    Parameters
    ----------
    table : `astropy.table.Table` instance

    Returns
    -------
    votable : `astropy.io.votable.tree.VOTableFile` instance
    """
    return tree.VOTableFile.from_table(table)
