# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes a fast iterator-based XML parser.
"""

# STDLIB
import collections
import contextlib
import io
import sys


__all__ = ['get_xml_iterator', 'get_xml_encoding', 'xml_readlines']


############################################################
# TODO: Refactor this into a py3k compatibility module
IS_PY3K = (sys.version_info[0] >= 3)

if IS_PY3K:
    def is_callable(o):
        """
        Abstracts away the different ways to test for a callable object in
        Python 2.x and 3.x.
        """
        return isinstance(o, collections.Callable)
else:
    def is_callable(o):
        """
        Abstracts away the different ways to test for a callable object in
        Python 2.x and 3.x.
        """
        return callable(o)
############################################################


@contextlib.contextmanager
def _convert_to_fd_or_read_function(fd):
    """
    Returns a function suitable for streaming input, or a file object.

    This function is only useful if passing off to C code where:

       - If it's a real file object, we want to use it as a real
         C file object to avoid the Python overhead.

       - If it's not a real file object, it's much handier to just
         have a Python function to call.

    Parameters
    ----------
    fd : object
        May be:

            - a file object, in which case it is returned verbatim.

            - a function that reads from a stream, in which case it is
              returned verbatim.

            - a file path, in which case it is opened.  If it ends in
              `.gz`, it is assumed to be a gzipped file, and the
              :meth:`read` method on the file object is returned.
              Otherwise, the raw file object is returned.

           - an object with a :meth:`read` method, in which case that
             method is returned.

    Returns
    -------
    fd : context-dependent
        See above.
    """
    if not sys.platform.startswith('win'):
        if IS_PY3K:
            if isinstance(fd, io.IOBase):
                yield fd
                return
        else:
            if isinstance(fd, file):
                yield fd
                return
    if is_callable(fd):
        yield fd
        return
    elif isinstance(fd, basestring):
        if fd.endswith('.gz'):
            from ...utils.compat import gzip
            with gzip.GzipFile(fd, 'rb') as real_fd:
                yield real_fd.read
                return
        else:
            with open(fd, 'rb') as real_fd:
                if sys.platform.startswith('win'):
                    # On Windows, we can't pass a real file descriptor
                    # to the C level, so we pass the read method
                    yield real_fd.read
                    return
                yield real_fd
                return
    elif hasattr(fd, 'read'):
        assert is_callable(fd.read)
        yield fd.read
        return
    else:
        raise TypeError("Can not be coerced to read function")


def _fast_iterparse(fd, buffersize=2 ** 10):
    from xml.parsers import expat

    if not is_callable(fd):
        read = fd.read
    else:
        read = fd

    queue = []
    text = []

    def start(name, attr):
        queue.append((True, name, attr,
                      (parser.CurrentLineNumber, parser.CurrentColumnNumber)))
        del text[:]

    def end(name):
        queue.append((False, name, ''.join(text).strip(),
                      (parser.CurrentLineNumber, parser.CurrentColumnNumber)))

    parser = expat.ParserCreate()
    if not IS_PY3K:
        parser.returns_unicode = False
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
    and time efficient than the alternative in `cElementTree`.

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
        start, tag, data, pos = iterator.next()
        if not start or tag != 'xml':
            raise IOError('Invalid XML file')

    return data['encoding']


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

    with io.open(source, 'rt', encoding=encoding) as input:
        xml_lines = input.readlines()

    return xml_lines
