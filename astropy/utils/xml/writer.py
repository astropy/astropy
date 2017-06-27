# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Contains a class that makes it simple to stream out well-formed and
nicely-indented XML.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from ...extern import six

# STDLIB
import contextlib
import textwrap

try:
    import bleach
    HAS_BLEACH = True
except ImportError:
    HAS_BLEACH = False

try:
    from . import _iterparser
except ImportError:
    def xml_escape_cdata(s):
        """
        Escapes &, < and > in an XML CDATA string.
        """
        s = s.replace("&", "&amp;")
        s = s.replace("<", "&lt;")
        s = s.replace(">", "&gt;")
        return s

    def xml_escape(s):
        """
        Escapes &, ', ", < and > in an XML attribute value.
        """
        s = s.replace("&", "&amp;")
        s = s.replace("'", "&apos;")
        s = s.replace("\"", "&quot;")
        s = s.replace("<", "&lt;")
        s = s.replace(">", "&gt;")
        return s
else:
    xml_escape_cdata = _iterparser.escape_xml_cdata
    xml_escape = _iterparser.escape_xml


class XMLWriter:
    """
    A class to write well-formed and nicely indented XML.

    Use like this::

        w = XMLWriter(fh)
        with w.tag('html'):
            with w.tag('body'):
                w.data('This is the content')

    Which produces::

        <html>
         <body>
          This is the content
         </body>
        </html>
    """

    def __init__(self, file):
        """
        Parameters
        ----------
        file : writable file-like object.
        """
        self.write = file.write
        if hasattr(file, "flush"):
            self.flush = file.flush
        self._open = 0  # true if start tag is open
        self._tags = []
        self._data = []
        self._indentation = " " * 64

        self.xml_escape_cdata = xml_escape_cdata
        self.xml_escape = xml_escape

    def _flush(self, indent=True, wrap=False):
        """
        Flush internal buffers.
        """
        if self._open:
            if indent:
                self.write(">\n")
            else:
                self.write(">")
            self._open = 0
        if self._data:
            data = ''.join(self._data)
            if wrap:
                indent = self.get_indentation_spaces(1)
                data = textwrap.fill(
                    data,
                    initial_indent=indent,
                    subsequent_indent=indent)
                self.write('\n')
                self.write(self.xml_escape_cdata(data))
                self.write('\n')
                self.write(self.get_indentation_spaces())
            else:
                self.write(self.xml_escape_cdata(data))
            self._data = []

    def start(self, tag, attrib={}, **extra):
        """
        Opens a new element.  Attributes can be given as keyword
        arguments, or as a string/string dictionary.  The method
        returns an opaque identifier that can be passed to the
        :meth:`close` method, to close all open elements up to and
        including this one.

        Parameters
        ----------
        tag : str
            The element name

        attrib : dict of str -> str
            Attribute dictionary.  Alternatively, attributes can
            be given as keyword arguments.

        Returns
        -------
        id : int
            Returns an element identifier.
        """
        self._flush()
        # This is just busy work -- we know our tag names are clean
        # tag = xml_escape_cdata(tag)
        self._data = []
        self._tags.append(tag)
        self.write(self.get_indentation_spaces(-1))
        self.write("<{}".format(tag))
        if attrib or extra:
            attrib = attrib.copy()
            attrib.update(extra)
            attrib = list(six.iteritems(attrib))
            attrib.sort()
            for k, v in attrib:
                if v is not None:
                    # This is just busy work -- we know our keys are clean
                    # k = xml_escape_cdata(k)
                    v = self.xml_escape(v)
                    self.write(" {}=\"{}\"".format(k, v))
        self._open = 1

        return len(self._tags)

    @contextlib.contextmanager
    def xml_cleaning_method(self, method='escape_xml', **clean_kwargs):
        """Context manager to control how XML data tags are cleaned (escaped) to
        remove potentially unsafe characters or constructs.

        The default (``method='escape_xml'``) applies brute-force escaping of
        certain key XML characters like ``<``, ``>``, and ``&`` to ensure that
        the output is not valid XML.

        In order to explicitly allow certain XML tags (e.g. link reference or
        emphasis tags), use ``method='bleach_clean'``.  This sanitizes the data
        string using the ``clean`` function of the
        `http://bleach.readthedocs.io/en/latest/clean.html <bleach>`_ package.
        Any additional keyword arguments will be passed directly to the
        ``clean`` function.

        Finally, use ``method='none'`` to disable any sanitization. This should
        be used sparingly.

        Example::

          w = writer.XMLWriter(ListWriter(lines))
          with w.xml_cleaning_method('bleach_clean'):
              w.start('td')
              w.data('<a href="http://google.com">google.com</a>')
              w.end()

        Parameters
        ----------
        method : str
            Cleaning method.  Allowed values are "escape_xml",
            "bleach_clean", and "none".

        **clean_kwargs : keyword args
            Additional keyword args that are passed to the
            bleach.clean() function.
        """
        current_xml_escape_cdata = self.xml_escape_cdata

        if method == 'bleach_clean':
            if HAS_BLEACH:
                if clean_kwargs is None:
                    clean_kwargs = {}
                self.xml_escape_cdata = lambda x: bleach.clean(x, **clean_kwargs)
            else:
                raise ValueError('bleach package is required when HTML escaping is disabled.\n'
                                 'Use "pip install bleach".')
        elif method == "none":
            self.xml_escape_cdata = lambda x: x
        elif method != 'escape_xml':
            raise ValueError('allowed values of method are "escape_xml", "bleach_clean", and "none"')

        yield

        self.xml_escape_cdata = current_xml_escape_cdata

    @contextlib.contextmanager
    def tag(self, tag, attrib={}, **extra):
        """
        A convenience method for creating wrapper elements using the
        ``with`` statement.

        Examples
        --------

        >>> with writer.tag('foo'):  # doctest: +SKIP
        ...     writer.element('bar')
        ... # </foo> is implicitly closed here
        ...

        Parameters are the same as to `start`.
        """
        self.start(tag, attrib, **extra)
        yield
        self.end(tag)

    def comment(self, comment):
        """
        Adds a comment to the output stream.

        Parameters
        ----------
        comment : str
            Comment text, as a Unicode string.
        """
        self._flush()
        self.write(self.get_indentation_spaces())
        self.write("<!-- {} -->\n".format(self.xml_escape_cdata(comment)))

    def data(self, text):
        """
        Adds character data to the output stream.

        Parameters
        ----------
        text : str
            Character data, as a Unicode string.
        """
        self._data.append(text)

    def end(self, tag=None, indent=True, wrap=False):
        """
        Closes the current element (opened by the most recent call to
        `start`).

        Parameters
        ----------
        tag : str
            Element name.  If given, the tag must match the start tag.
            If omitted, the current element is closed.
        """
        if tag:
            if not self._tags:
                raise ValueError("unbalanced end({})".format(tag))
            if tag != self._tags[-1]:
                raise ValueError("expected end({}), got {}".format(
                        self._tags[-1], tag))
        else:
            if not self._tags:
                raise ValueError("unbalanced end()")
        tag = self._tags.pop()
        if self._data:
            self._flush(indent, wrap)
        elif self._open:
            self._open = 0
            self.write("/>\n")
            return
        if indent:
            self.write(self.get_indentation_spaces())
        self.write("</{}>\n".format(tag))

    def close(self, id):
        """
        Closes open elements, up to (and including) the element identified
        by the given identifier.

        Parameters
        ----------
        id : int
            Element identifier, as returned by the `start` method.
        """
        while len(self._tags) > id:
            self.end()

    def element(self, tag, text=None, wrap=False, attrib={}, **extra):
        """
        Adds an entire element.  This is the same as calling `start`,
        `data`, and `end` in sequence. The ``text`` argument
        can be omitted.
        """
        self.start(tag, attrib, **extra)
        if text:
            self.data(text)
        self.end(indent=False, wrap=wrap)

    def flush(self):
        pass  # replaced by the constructor

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
        return self._indentation[:len(self._tags) + offset]

    @staticmethod
    def object_attrs(obj, attrs):
        """
        Converts an object with a bunch of attributes on an object
        into a dictionary for use by the `XMLWriter`.

        Parameters
        ----------
        obj : object
            Any Python object

        attrs : sequence of str
            Attribute names to pull from the object

        Returns
        -------
        attrs : dict
            Maps attribute names to the values retrieved from
            ``obj.attr``.  If any of the attributes is `None`, it will
            not appear in the output dictionary.
        """
        d = {}
        for attr in attrs:
            if getattr(obj, attr) is not None:
                d[attr.replace('_', '-')] = six.text_type(getattr(obj, attr))
        return d
