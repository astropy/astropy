"""
Various XML-related utilities
"""

from __future__ import division, absolute_import

# STDLIB
import os

# ASTROPY
from ...utils.xml import check as xml_check
from ...utils.xml import validate

# LOCAL
from . import util
from .util import IS_PY3K
from .exceptions import (warn_or_raise, vo_warn, vo_raise,
     VOTableChangeWarning, W02, W03, W04, W05, W34)


def fast_iterparse(fd, buffersize=2 ** 10):
    """
    Based on :class:`cElementTree.iterparse`, but doesn't ever build a
    tree at all.  This makes things much faster and more memory
    efficient.

    Parameters
    ----------
    fd : readable file-like object, read function or str path

    Returns
    -------
    parts : iterator

        The iterator returns 3-tuples (*event*, *tag*, *data*):

            - *event*: may be either 'start' (when an element begins) or 'end'
              (when an element is completed).

            - *tag*: The name of the element

            - *data*: Depends on the value of *event*:

                - if *event* == 'start', data is a dictionary of
                  attributes

                - if *event* == 'end', data is a list of strings, that
                  when joined is the text content of the element
    """
    from xml.parsers import expat

    close_at_end = False
    if not util.is_callable(fd):
        if not hasattr(fd, 'read'):
            close_at_end = True
            fd = open(fd, 'rb')
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

    if close_at_end:
        fd.close()
    Parse('', True)
    for elem in queue:
        yield elem

# Try to import the C version of the iterparser, otherwise fall back
# to the Python implementation above.
slow_iterparse = fast_iterparse
try:
    from . import iterparser
    fast_iterparse = iterparser.IterParser
except ImportError:
    pass


def get_xml_iterator(source, _debug_python_based_parser=False):
    """
    Returns an iterator over the elements of an XML file.  See
    :class:`fast_iterparse` for more information.
    """
    if _debug_python_based_parser:
        context = slow_iterparse(source)
    else:
        context = fast_iterparse(source)
    iterator = iter(context)
    return iterator


def check_id(ID, name='ID', config={}, pos=None):
    """
    Raises a `astropy.io.vo.exceptions.VOTableSpecError` if *ID* is not a
    valid XML ID_.  *name* is the name of the attribute being checked
    (used only for error messages).
    """
    if (ID is not None and not xml_check.check_id(ID)):
        warn_or_raise(W02, W02, (name, ID), config, pos)
        return False
    return True


def fix_id(ID, config={}, pos=None):
    """
    Given an arbitrary string, create one that can be used as an xml
    id.  This is rather simplistic at the moment, since it just
    replaces non-valid characters with underscores.
    """
    if ID is None:
        return None
    corrected = xml_check.fix_id(ID)
    if corrected != ID:
        vo_warn(W03, (ID, corrected), config, pos)
    return corrected


_token_regex = r"(?![\r\l\t ])[^\r\l\t]*(?![\r\l\t ])"


def check_token(token, attr_name, config={}, pos=None):
    """
    Raises a `ValueError` if *token* is not a valid XML token, as
    defined by XML Schema Part 2.
    """
    if (token is not None and not xml_check.check_token(token)):
        return False
    return True


def check_mime_content_type(content_type, config={}, pos=None):
    """
    Raises a `ValueError` if *content_type* is not a valid MIME
    content type (syntactically at least), as defined by RFC 2045.
    """
    if (content_type is not None and
        not xml_check.check_mime_content_type(content_type)):
        warn_or_raise(W04, W04, content_type, config, pos)
        return False
    return True


def check_anyuri(uri, config={}, pos=None):
    """
    Raises a `ValueError` if *uri* is not a valid URI as defined in RFC
    2396.
    """
    if (uri is not None and not xml_check.check_anyuri(uri)):
        warn_or_raise(W05, W05, uri, config, pos)
        return False
    return True


def validate_schema(filename, version='1.2'):
    """
    Validates the given file against the appropriate VOTable schema.

    Parameters
    ----------
    filename : str
        The path to the XML file to validate

    version : str
        The VOTABLE version to check, which must be a string \"1.0\",
        \"1.1\", or \"1.2\".

        For version \"1.0\", it is checked against a DTD, since that
        version did not have an XML Schema.

    Returns
    -------
    returncode, stdout, stderr : int, str, str
        Returns the returncode from xmllint and the stdout and stderr
        as strings
    """
    import subprocess

    assert version in ('1.0', '1.1', '1.2')

    if version in ('1.1', '1.2'):
        schema_path = os.path.join(
            os.path.dirname(__file__),
            "data", "VOTable.v%s.xsd" % version)
    else:
        schema_path = os.path.join(
            os.path.dirname(__file__),
            "data", "VOTable.dtd")

    return validate.validate_schema(filename, schema_path)
