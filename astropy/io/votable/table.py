# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This file contains a contains the high-level functions to read a
VOTable file.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...extern import six

# STDLIB
import io
import os
import sys
import textwrap
import warnings

# LOCAL
from . import exceptions
from . import tree
from ...utils.xml import iterparser
from ...utils import data
from ...config import ConfigAlias


__all__ = ['parse', 'parse_single_table', 'from_table', 'writeto', 'validate',
           'reset_vo_warnings']


PEDANTIC = ConfigAlias(
    '0.4', 'PEDANTIC', 'pedantic',
    'astropy.io.votable.table', 'astropy.io.votable')


def parse(source, columns=None, invalid='exception', pedantic=None,
          chunk_size=tree.DEFAULT_CHUNK_SIZE, table_number=None,
          table_id=None, filename=None, unit_format=None,
          datatype_mapping=None, _debug_python_based_parser=False):
    """
    Parses a VOTABLE_ xml file (or file-like object), and returns a
    `~astropy.io.votable.tree.VOTableFile` object.

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
        ``astropy.io.votable.pedantic``, which defaults to False.

    chunk_size : int, optional
        The number of rows to read before converting to an array.
        Higher numbers are likely to be faster, but will consume more
        memory.

    table_number : int, optional
        The number of table in the file to read in.  If `None`, all
        tables will be read.  If a number, 0 refers to the first table
        in the file, and only that numbered table will be parsed and
        read in.  Should not be used with ``table_id``.

    table_id : str, optional
        The ID of the table in the file to read in.  Should not be
        used with ``table_number``.

    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    unit_format : str, astropy.units.format.Base instance or None, optional
        The unit format to use when parsing unit attributes.  If a
        string, must be the name of a unit formatter. The built-in
        formats include ``generic``, ``fits``, ``cds``, and
        ``vounit``.  A custom formatter may be provided by passing a
        `~astropy.units.UnitBase` instance.  If `None` (default),
        the unit format to use will be the one specified by the
        VOTable specification (which is ``cds`` up to version 1.2 of
        VOTable, and (probably) ``vounit`` in future versions of the
        spec).

    datatype_mapping : dict of str to str, optional
        A mapping of datatype names to valid VOTable datatype names.
        For example, if the file being read contains the datatype
        "unsignedInt" (an invalid datatype in VOTable), include the
        mapping ``{"unsignedInt": "long"}``.

    Returns
    -------
    votable : `~astropy.io.votable.tree.VOTableFile` object

    See also
    --------
    astropy.io.votable.exceptions : The exceptions this function may raise.
    """
    from . import conf

    invalid = invalid.lower()
    assert invalid in ('exception', 'mask')

    if pedantic is None:
        pedantic = conf.pedantic

    if datatype_mapping is None:
        datatype_mapping = {}

    config = {
        'columns'          : columns,
        'invalid'          : invalid,
        'pedantic'         : pedantic,
        'chunk_size'       : chunk_size,
        'table_number'     : table_number,
        'filename'         : filename,
        'unit_format'      : unit_format,
        'datatype_mapping' : datatype_mapping
    }

    if filename is None and isinstance(source, six.string_types):
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
    votable : `~astropy.io.votable.tree.Table` object
    """
    if kwargs.get('table_number') is None:
        kwargs['table_number'] = 0

    votable = parse(source, **kwargs)

    return votable.get_first_table()


def writeto(table, file, tabledata_format=None):
    """
    Writes a `~astropy.io.votable.tree.VOTableFile` to a VOTABLE_ xml file.

    Parameters
    ----------
    table : `~astropy.io.votable.tree.VOTableFile` or `~astropy.table.Table` instance.

    file : str or writable file-like object
        Path or file object to write to

    tabledata_format : str, optional
        Override the format of the table(s) data to write.  Must be
        one of ``tabledata`` (text representation), ``binary`` or
        ``binary2``.  By default, use the format that was specified in
        each ``table`` object as it was created or read in.  See
        :ref:`votable-serialization`.
    """
    from ...table import Table
    if isinstance(table, Table):
        table = tree.VOTableFile.from_table(table)
    elif not isinstance(table, tree.VOTableFile):
        raise TypeError(
            "first argument must be astropy.io.vo.VOTableFile or "
            "astropy.table.Table instance")
    table.to_xml(file, tabledata_format=tabledata_format,
                 _debug_python_based_parser=True)


def validate(source, output=None, xmllint=False, filename=None):
    """
    Prints a validation report for the given file.

    Parameters
    ----------
    source : str or readable file-like object
        Path to a VOTABLE_ xml file.

    output : writable file-like object, optional
        Where to output the report.  Defaults to ``sys.stdout``.
        If `None`, the output will be returned as a string.

    xmllint : bool, optional
        When `True`, also send the file to ``xmllint`` for schema and
        DTD validation.  Requires that ``xmllint`` is installed.  The
        default is `False`.  ``source`` must be a file on the local
        filesystem in order for ``xmllint`` to work.

    filename : str, optional
        A filename to use in the error messages.  If not provided, one
        will be automatically determined from ``source``.

    Returns
    -------
    is_valid : bool or str
        Returns `True` if no warnings were found.  If ``output`` is
        `None`, the return value will be a string.
    """

    from ...utils.console import print_code_line, color_print

    if output is None:
        output = sys.stdout

    return_as_str = False
    if output is None:
        output = io.StringIO()

    lines = []
    votable = None

    reset_vo_warnings()

    with data.get_readable_fileobj(source, encoding='binary') as fd:
        content = fd.read()
    content_buffer = io.BytesIO(content)
    content_buffer.seek(0)

    if filename is None:
        if isinstance(source, six.string_types):
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

    lines = [str(x.message) for x in warning_lines if
             issubclass(x.category, exceptions.VOWarning)] + lines

    content_buffer.seek(0)
    output.write("Validation report for {0}\n\n".format(filename))

    if len(lines):
        xml_lines = iterparser.xml_readlines(content_buffer)

        for warning in lines:
            w = exceptions.parse_vowarning(warning)

            if not w['is_something']:
                output.write(w['message'])
                output.write('\n\n')
            else:
                line = xml_lines[w['nline'] - 1]
                warning = w['warning']
                if w['is_warning']:
                    color = 'yellow'
                else:
                    color = 'red'
                color_print(
                    '{0:d}: '.format(w['nline']), '',
                    warning or 'EXC', color,
                    ': ', '',
                    textwrap.fill(
                        w['message'],
                        initial_indent='          ',
                        subsequent_indent='  ').lstrip(),
                    file=output)
                print_code_line(line, w['nchar'], file=output)
            output.write('\n')
    else:
        output.write('astropy.io.votable found no violations.\n\n')

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
                'xmllint schema violations:\n\n')
            output.write(stderr)
        else:
            output.write('xmllint passed\n')

    if return_as_str:
        return output.getvalue()
    return len(lines) == 0 and success == 0


def from_table(table, table_id=None):
    """
    Given an `~astropy.table.Table` object, return a
    `~astropy.io.votable.tree.VOTableFile` file structure containing
    just that single table.

    Parameters
    ----------
    table : `~astropy.table.Table` instance

    table_id : str, optional
        If not `None`, set the given id on the returned
        `~astropy.io.votable.tree.Table` instance.

    Returns
    -------
    votable : `~astropy.io.votable.tree.VOTableFile` instance
    """
    return tree.VOTableFile.from_table(table, table_id=table_id)


def is_votable(source):
    """
    Reads the header of a file to determine if it is a VOTable file.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a VOTABLE_ xml file.

    Returns
    -------
    is_votable : bool
        Returns `True` if the given file is a VOTable file.
    """
    try:
        with iterparser.get_xml_iterator(source) as iterator:
            for start, tag, data, pos in iterator:
                if tag != 'xml':
                    return False
                break

            for start, tag, data, pos in iterator:
                if tag != 'VOTABLE':
                    return False
                break

            return True
    except ValueError:
        return False


def reset_vo_warnings():
    """
    Resets all of the vo warning state so that warnings that
    have already been emitted will be emitted again. This is
    used, for example, by `validate` which must emit all
    warnings each time it is called.

    """
    from . import converters, xmlutil

    #-----------------------------------------------------------#
    # This is a special variable used by the Python warnings    #
    # infrastructure to keep track of warnings that have        #
    # already been seen.  Since we want to get every single     #
    # warning out of this, we have to delete all of them first. #
    #-----------------------------------------------------------#
    for module in (converters, exceptions, tree, xmlutil):
        if hasattr(module, '__warningregistry__'):
            del module.__warningregistry__
