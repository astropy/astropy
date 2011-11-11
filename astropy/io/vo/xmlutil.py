"""
Various XML-related utilities
"""

from __future__ import division, absolute_import

# STDLIB
import contextlib
import os
import re
import sys
import textwrap
import urlparse

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
    if (ID is not None and
        re.match(r"^[A-Za-z_][A-Za-z0-9_\.\-]*$", ID) is None):
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
    if re.match(r"^[A-Za-z_][A-Za-z0-9_\.\-]*$", ID):
        return ID
    if len(ID):
        corrected = ID
        if not len(corrected) or re.match('^[^A-Za-z_]$', corrected[0]):
            corrected = '_' + corrected
        corrected = (re.sub(r"[^A-Za-z_]", '_', corrected[0]) +
                     re.sub(r"[^A-Za-z0-9_\.\-]", "_", corrected[1:]))
        vo_warn(W03, (ID, corrected), config, pos)
        return corrected
    return ''


_token_regex = r"(?![\r\l\t ])[^\r\l\t]*(?![\r\l\t ])"


def check_token(token, attr_name, config={}, pos=None):
    """
    Raises a `ValueError` if *token* is not a valid XML token, as
    defined by XML Schema Part 2.
    """
    if (token is not None and
        not (token == '' or
             re.match(
                "[^\r\n\t ]?([^\r\n\t ]| [^\r\n\t ])*[^\r\n\t ]?$", token))):
        warn_or_raise(W34, W34, (token, attr_name), config, pos)
    return True


def check_mime_content_type(content_type, config={}, pos=None):
    """
    Raises a `ValueError` if *content_type* is not a valid MIME
    content type (syntactically at least), as defined by RFC 2045.
    """
    ctrls = ''.join(chr(x) for x in xrange(0, 0x20))
    token_regex = '[^()<>@,;:\\\"/[\]?= %s\x7f]+' % ctrls
    if (content_type is not None and
        re.match(
            r'(?P<type>%s)/(?P<subtype>%s)$' % (token_regex, token_regex),
            content_type) is None):
        warn_or_raise(W04, W04, content_type, config, pos)
        return False
    return True


def check_anyuri(uri, config={}, pos=None):
    """
    Raises a `ValueError` if *uri* is not a valid URI as defined in RFC
    2396.
    """
    if uri is not None:
        if (re.match(
            r"(([a-zA-Z][0-9a-zA-Z+\-\.]*:)?/{0,2}" +
            r"[0-9a-zA-Z;/?:@&=+$\.\-_!~*'()%]+)?" +
            r"(#[0-9a-zA-Z;/?:@&=+$\.\-_!~*'()%]+)?",
            uri) is None):
            warn_or_raise(W05, W05, uri, config, pos)
            return False
        try:
            urlparse.urlparse(uri)
        except:
            warn_or_raise(W05, W05, uri, config, pos)
            return False
    return True


def validate_schema(filename, version='1.2'):
    """
    Validates the given file against the appropriate VOTable schema
    corresponding to the given *version*, which must be a string
    "1.0", "1.1", or "1.2".

    For version "1.0", it is checked against a DTD, since that version
    did not have an XML Schema.
    """
    import subprocess

    assert version in ('1.0', '1.1', '1.2')

    if version in ('1.1', '1.2'):
        schema_path = os.path.join(
            os.path.dirname(__file__),
            "data", "VOTable.v%s.xsd" % version)
        schema_part = '--schema %s' % schema_path
    else:
        dtd_path = os.path.join(
            os.path.dirname(__file__),
            "data", "VOTable.dtd")
        schema_part = '--dtdvalid %s' % dtd_path

    p = subprocess.Popen(
        "xmllint --noout --nonet %s %s" %
        (schema_part, filename),
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if p.returncode == 127:
        raise OSError(
            "xmllint not found, so can not validate schema")

    return p.returncode, stdout, stderr

