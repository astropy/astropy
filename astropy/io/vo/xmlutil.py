# Copyright (C) 2008-2010 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

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
from .voexceptions import warn_or_raise, vo_warn, vo_raise, \
     VOTableChangeWarning, W02, W03, W04, W05, W34


def xml_escape_cdata(s):
    s = s.replace(u"&", u"&amp;")
    s = s.replace(u"<", u"&lt;")
    s = s.replace(u">", u"&gt;")
    return s


def xml_escape(s):
    s = s.replace(u"&", u"&amp;")
    s = s.replace(u"'", u"&apos;")
    s = s.replace(u"\"", u"&quot;")
    s = s.replace(u"<", u"&lt;")
    s = s.replace(u">", u"&gt;")
    return s


def fast_iterparse(fd, buffersize = 2 ** 10):
    """
    Based on :class:`cElementTree.iterparse`, but doesn't ever build a
    tree at all.  This makes things much faster and more memory
    efficient.

    The iterator returns 3-tuples (*event*, *tag*, *data*):

    - *event*: may be either 'start' (when an element begins) or 'end'
      (when an element is completed).

    - *tag*: The name of the element

    - *data*: Depends on the value of *event*:

      - if *event* == 'start', data is a dictionary of attributes

      - if *event* == 'end', data is a list of strings, that when
        joined is the text content of the element
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
    xml_escape = iterparser.escape_xml
    xml_escape_cdata = iterparser.escape_xml_cdata
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


def object_attrs(obj, attrs):
    """
    Converts an object with a bunch of attributes on an object into a
    dictionary for use by the XMLWriter.

    *obj* is any Python object

    *attrs* is a sequence of attribute names to pull from the object

    If any of the attributes is ``None``, it will not appear in the
    output dictionary.
    """
    d = {}
    for attr in attrs:
        if getattr(obj, attr) is not None:
            d[attr.replace('_', '-')] = str(getattr(obj, attr))
    return d


def check_id(ID, name='ID', config={}, pos=None):
    """
    Raises a :exc:`~vo.voexceptions.VOTableSpecError` if *ID* is not a
    valid XML ID_.  *name* is the name of the attribute being checked
    (used only for error messages).
    """
    if ID is not None and re.match(r"^[A-Za-z_][A-Za-z0-9_\.\-]*$", ID) is None:
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
    Raises a :exc:`ValueError` if *token* is not a valid XML token, as
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
    Raises a :exc:`ValueError` if *content_type* is not a valid MIME
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
    Raises a :exc:`ValueError` if *uri* is not a valid URI as defined in RFC
    2396.
    """
    if uri is not None:
        if (re.match(
            r"(([a-zA-Z][0-9a-zA-Z+\-\.]*:)?/{0,2}[0-9a-zA-Z;/?:@&=+$\.\-_!~*'()%]+)?(#[0-9a-zA-Z;/?:@&=+$\.\-_!~*'()%]+)?",
            uri) is None):
            warn_or_raise(W05, W05, uri, config, pos)
            return False
        try:
            urlparse.urlparse(uri)
        except:
            warn_or_raise(W05, W05, uri, config, pos)
            return False
    return True


class XMLWriter:
    def __init__(self, file):
        """
        *file* is a writable file-like object.
        """
        self.write = file.write
        if hasattr(file, "flush"):
            self.flush = file.flush
        self._open = 0 # true if start tag is open
        self._tags = []
        self._data = []
        self._indentation = u" " * 64

    def _flush(self, indent=True, wrap=False):
        """
        Flush internal buffers.
        """
        if self._open:
            if indent:
                self.write(u">\n")
            else:
                self.write(u">")
            self._open = 0
        if self._data:
            data = u''.join(self._data)
            if wrap:
                indent = self.get_indentation_spaces(1)
                data = textwrap.fill(
                    data,
                    initial_indent=indent,
                    subsequent_indent=indent)
                self.write('\n')
                self.write(xml_escape_cdata(data))
                self.write('\n')
                self.write(self.get_indentation_spaces())
            else:
                self.write(xml_escape_cdata(data))
            self._data = []

    def start(self, tag, attrib={}, **extra):
        """
        Opens a new element.  Attributes can be given as keyword
        arguments, or as a string/string dictionary. The method
        returns an opaque identifier that can be passed to the
        :meth:`close` method, to close all open elements up to and
        including this one.

        *tag*: Element tag.

        *attrib*: Attribute dictionary.  Alternatively, attributes can
        be given as keyword arguments.

        Returns an element identifier.
        """
        self._flush()
        # This is just busy work -- we know our tag names are clean
        # tag = xml_escape_cdata(tag)
        self._data = []
        self._tags.append(tag)
        self.write(self.get_indentation_spaces(-1))
        self.write(u"<%s" % tag)
        if attrib or extra:
            attrib = attrib.copy()
            attrib.update(extra)
            attrib = attrib.items()
            attrib.sort()
            for k, v in attrib:
                if not v == '' and v is not None:
                    # This is just busy work -- we know our keys are clean
                    # k = xml_escape_cdata(k)
                    v = xml_escape(v)
                    self.write(u" %s=\"%s\"" % (k, v))
        self._open = 1

        return len(self._tags)

    @contextlib.contextmanager
    def tag(self, tag, attrib={}, **extra):
        """
        A convenience method for use with the `with` statement::

            with writer.tag('foo'):
                writer.element('bar')
            # </foo> is implicitly closed here

        Parameters are the same as to `start`.
        """
        self.start(tag, attrib, **extra)
        yield
        self.end(tag)

    def comment(self, comment):
        """
        Adds a comment to the output stream.

        *comment*: Comment text, as a Unicode string.
        """
        self._flush()
        self.write(self.get_indentation_spaces())
        self.write(u"<!-- %s -->\n" % escape_cdata(comment))

    def data(self, text):
        """
        Adds character data to the output stream.

        *text*: Character data, as a Unicode string.
        """
        self._data.append(text)

    def end(self, tag=None, indent=True, wrap=False):
        """
        Closes the current element (opened by the most recent call to
        `start`).

        *tag*: Element tag.  If given, the tag must match the start
        tag.  If omitted, the current element is closed.
        """
        if tag:
            assert self._tags, "unbalanced end(%s)" % tag
            assert tag == self._tags[-1],\
                   "expected end(%s), got %s" % (self._tags[-1], tag)
        else:
            assert self._tags, "unbalanced end()"
        tag = self._tags.pop()
        if self._data:
            self._flush(indent, wrap)
        elif self._open:
            self._open = 0
            self.write(u"/>\n")
            return
        if indent:
            self.write(self.get_indentation_spaces())
        self.write(u"</%s>\n" % tag)

    def close(self, id):
        """
        Closes open elements, up to (and including) the element identified
        by the given identifier.

        *id*: Element identifier, as returned by the `start` method.
        """
        while len(self._tags) > id:
            self.end()

    def element(self, tag, text=None, wrap=False, attrib={}, **extra):
        """
        Adds an entire element.  This is the same as calling `start`,
        `data`, and `end` in sequence. The `text` argument
        can be omitted.
        """
        self.start(tag, attrib, **extra)
        if text:
            self.data(text)
        self.end(indent=False, wrap=wrap)

    def flush(self):
        pass # replaced by the constructor

    def get_indentation(self):
        """
        Returns the number of indentation levels the file is currently
        in.
        """
        return len(self._tags)

    def get_indentation_spaces(self, offset=0):
        """
        Returns a string of spaces that matches the current
        indentation level.
        """
        return self._indentation[:len(self._tags)+offset]


def validate_schema(filename, version='1.2'):
    """
    Validates the given file against the appropriate VOTable schema
    corresponding to the given *version*, which must be a string "1.0",
    "1.1", or "1.2".

    For version "1.0", it is checked against a DTD, since that version did
    not have an XML Schema.
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
    return p.returncode, stdout, stderr

