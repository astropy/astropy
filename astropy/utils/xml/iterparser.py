# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes a fast iterator-based XML parser.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six

# STDLIB
import contextlib
import io
import sys

# ASTROPY
from .. import data


__all__ = ['get_xml_iterator', 'get_xml_encoding', 'xml_readlines']


@contextlib.contextmanager
def _convert_to_fd_or_read_function(fd):
    """
    Returns a function suitable for streaming input, or a file object.

    This function is only useful if passing off to C code where:

       - If it's a real file object, we want to use it as a real
         C file object to avoid the Python overhead.

       - If it's not a real file object, it's much handier to just
         have a Python function to call.

    This is somewhat quirky behavior, of course, which is why it is
    private.  For a more useful version of similar behavior, see
    `astropy.utils.misc.get_readable_fileobj`.

    Parameters
    ----------
    fd : object
        May be:

            - a file object.  If the file is uncompressed, this raw
              file object is returned verbatim.  Otherwise, the read
              method is returned.

            - a function that reads from a stream, in which case it is
              returned verbatim.

            - a file path, in which case it is opened.  Again, like a
              file object, if it's uncompressed, a raw file object is
              returned, otherwise its read method.

            - an object with a :meth:`read` method, in which case that
              method is returned.

    Returns
    -------
    fd : context-dependent
        See above.
    """
    if six.callable(fd):
        yield fd
        return

    with data.get_readable_fileobj(fd, encoding='binary') as new_fd:
        if sys.platform.startswith('win'):
            yield new_fd.read
        else:
            if six.PY3:
                if isinstance(new_fd, io.FileIO):
                    yield new_fd
                else:
                    yield new_fd.read
            elif six.PY2:
                if isinstance(new_fd, file):
                    yield new_fd
                else:
                    yield new_fd.read


def _fast_iterparse(fd, buffersize=2 ** 10):
    from xml.parsers import expat

    if not six.callable(fd):
        read = fd.read
    else:
        read = fd

    queue = []
    text = []

    def start(name, attr):
        queue.append((True, name, attr,
                      (parser.CurrentLineNumber, parser.CurrentColumnNumber)))
        del text[:]

    if sys.version_info[:3] < (2, 6, 5):  # pragma py2
        # Due to Python issue #4978, convert all keys to byte strings
        _start = start
        def start(name, attr):
            attr = dict((k.encode('utf-8'), v) for (k, v) in six.iteritems(attr))
            return _start(name, attr)

    def end(name):
        queue.append((False, name, ''.join(text).strip(),
                      (parser.CurrentLineNumber, parser.CurrentColumnNumber)))

    parser = expat.ParserCreate()
    if six.PY2:
        parser.returns_unicode = True
    parser.specified_attributes = True
    parser.StartElementHandler = start
    parser.EndElementHandler = end
    parser.CharacterDataHandler = text.append
    Parse = parser.Parse

    data = read(buffersize)
    while data:
        Parse(data, False)
        for elem in queue:
            yield elem
        del queue[:]
        data = read(buffersize)

    Parse('', True)
    for elem in queue:
        yield elem


# Try to import the C version of the iterparser, otherwise fall back
# to the Python implementation above.
_slow_iterparse = _fast_iterparse
try:
    from . import _iterparser
    _fast_iterparse = _iterparser.IterParser
except ImportError:
    pass


@contextlib.contextmanager
def get_xml_iterator(source, _debug_python_based_parser=False):
    """
    Returns an iterator over the elements of an XML file.

    The iterator doesn't ever build a tree, so it is much more memory
    and time efficient than the alternative in ``cElementTree``.

    Parameters
    ----------
    fd : readable file-like object or read function

    Returns
    -------
    parts : iterator

        The iterator returns 4-tuples (*start*, *tag*, *data*, *pos*):

            - *start*: when `True` is a start element event, otherwise
              an end element event.

            - *tag*: The name of the element

            - *data*: Depends on the value of *event*:

                - if *start* == `True`, data is a dictionary of
                  attributes

                - if *start* == `False`, data is a string containing
                  the text content of the element

            - *pos*: Tuple (*line*, *col*) indicating the source of the
              event.
    """
    with _convert_to_fd_or_read_function(source) as fd:
        if _debug_python_based_parser:
            context = _slow_iterparse(fd)
        else:
            context = _fast_iterparse(fd)
        yield iter(context)


def get_xml_encoding(source):
    """
    Determine the encoding of an XML file by reading its header.

    Parameters
    ----------
    source : readable file-like object, read function or str path

    Returns
    -------
    encoding : str
    """
    with get_xml_iterator(source) as iterator:
        start, tag, data, pos = six.next(iterator)
        if not start or tag != 'xml':
            raise IOError('Invalid XML file')

    # The XML spec says that no encoding === utf-8
    return data.get('encoding') or 'utf-8'


def xml_readlines(source):
    """
    Get the lines from a given XML file.  Correctly determines the
    encoding and always returns unicode.

    Parameters
    ----------
    source : readable file-like object, read function or str path

    Returns
    -------
    lines : list of unicode
    """
    encoding = get_xml_encoding(source)

    with data.get_readable_fileobj(source, encoding=encoding) as input:
        input.seek(0)
        xml_lines = input.readlines()

    return xml_lines
