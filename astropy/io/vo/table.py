"""
This file contains a contains the high-level functions to read a
VOTable file.
"""

from __future__ import division, absolute_import

# STDLIB
import io
import sys
import warnings

# LOCAL
from . import exceptions
from . import tree
from . import util
from ...utils.xml import iterparser


__all__ = ['parse', 'parse_single_table', 'validate']


def parse(source, columns=None, invalid='exception', pedantic=True,
          chunk_size=tree.DEFAULT_CHUNK_SIZE, table_number=None,
          filename=None,
          _debug_python_based_parser=False):
    """
    Parses a VOTABLE_ xml file (or file-like object), and returns a
    `~astropy.io.vo.tree.VOTable` object, with nested
    `~astropy.io.vo.tree.Resource` instances and
    `~astropy.io.vo.tree.Table` instances.

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
    votable : `astropy.io.vo.tree.VOTableFile` object

    See also
    --------
    astropy.io.vo.exceptions : The exceptions this function may raise.
    """
    invalid = invalid.lower()
    assert invalid in ('exception', 'mask')

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
    Parses a VOTABLE_ xml file (or file-like object), reads only the
    first TABLE_ element, and returns a `~astropy.io.vo.tree.Table`
    instance.

    See `parse` for a description of the keyword arguments.

    Returns
    -------
    votable : `astropy.io.vo.tree.Table` object
    """
    if kwargs.get('table_number') is None:
        kwargs['table_number'] = 0

    votable = parse(source, **kwargs)

    return votable.get_first_table()


def validate(filename, output=sys.stdout, xmllint=False):
    """
    Prints a validation report for the given file.

    Parameters
    ----------
    filename : str path
        Path to a VOTABLE_ xml file.

    output : writable file-like object, optional
        Where to output the report.  Defaults to `sys.stdout`.
        If `None`, the output will be returned as a string.

    xmllint : bool, optional
        When `True`, also send the file to `xmllint` for schema and
        DTD validation.  Requires that `xmllint` is installed.  The
        default is `False`.

    Returns
    -------
    is_valid : bool or str
        Returns `True` if no warnings were found.  If `output` is
        `None`, the return value will be a string.
    """
    import textwrap
    from ...utils.console import print_code_line, color_print

    return_as_str = False
    if output is None:
        output = io.StringIO()

    lines = []
    votable = None

    # This is a special variable used by the Python warnings
    # infrastructure to keep track of warnings that have already been
    # seen.  Since we want to get every single warning out of this, we
    # have to delete it first.
    if hasattr(exceptions, '__warningregistry__'):
        del exceptions.__warningregistry__

    with io.open(filename, 'rb') as input:
        with warnings.catch_warnings(record=True) as warning_lines:
            warnings.resetwarnings()
            warnings.simplefilter("always", append=True)
            try:
                votable = parse(input, pedantic=False, filename=filename)
            except ValueError as e:
                lines.append(str(e))
    lines = [str(x.message) for x in warning_lines] + lines

    output.write(u"Validation report for {0}\n\n".format(filename))

    if len(lines):
        xml_lines = iterparser.xml_readlines(filename)

        for warning in lines:
            w = exceptions.parse_vowarning(warning)

            if not w['is_something']:
                output.write(warning)
                output.write(u'\n\n')
            else:
                line = xml_lines[w['nline'] - 1]
                warning = w['warning']
                if warning.startswith('W'):
                    color = 'yellow'
                else:
                    color = 'red'
                color_print(
                    u'{0:d}: '.format(w['nline']), '',
                    warning, color,
                    u': ', '',
                    textwrap.fill(
                        w['message'],
                        initial_indent=u'          ',
                        subsequent_indent=u'  ').lstrip(),
                    file=output)
                print_code_line(line, w['nchar'], file=output)
            output.write(u'\n')
    else:
        output.write(u'astropy.io.vo found no violations.\n\n')

    success = 0
    if xmllint:
        if votable is None:
            version = "1.1"
        else:
            version = votable.version
        success, stdout, stderr = xmlutil.validate_schema(
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

