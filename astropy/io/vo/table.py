"""
This file contains a contains the high-level functions to read a
VOTable file.
"""

from __future__ import division, absolute_import

# LOCAL
from . import tree
from . import util
from . import xmlutil


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
    votable : `astropy.io.vo.tree.VOTable` object

    Raises
    ------
    `astropy.io.vo.exceptions.VOTableSpecError` : The table violates the specification.
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

    source = util.convert_to_fd_or_read_function(source)

    iterator = xmlutil.get_xml_iterator(
        source, _debug_python_based_parser=_debug_python_based_parser)

    return tree.VOTableFile(config=config, pos=(1, 1)).parse(iterator, config)


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
