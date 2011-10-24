# -*- coding: utf-8 -*-

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
Exceptions and warnings used by the vo package.
"""

from __future__ import division, absolute_import

# STDLIB
import re
import sys
from warnings import warn

# LOCAL
from .util import IS_PY3K

MAX_WARNINGS = 10

def _format_message(message, name, config={}, pos=None):
    if pos is None:
        pos = ('?', '?')
    filename = config.get('filename', '?')
    return '%s:%s:%s: %s: %s' % (filename, pos[0], pos[1], name, message)

def _suppressed_warning(warning, config):
    warning_class = type(warning)
    config.setdefault('_warning_counts', {}).setdefault(warning_class, 0)
    config['_warning_counts'][warning_class] += 1
    message_count = config['_warning_counts'][warning_class]
    if message_count <= MAX_WARNINGS:
        if message_count == MAX_WARNINGS:
            warning.formatted_message += ' (suppressing further warnings of this type...)'
        warn(str(warning))

def warn_or_raise(warning_class, exception_class=None, args=(), config={}, pos=None):
    """
    Warn or raise an exception, depending on the pedantic setting.
    """
    if config.get('pedantic'):
        if exception_class is None:
            exception_class = warning_class
        vo_raise(exception_class, args, config, pos)
    else:
        vo_warn(warning_class, args, config, pos)

def vo_raise(exception_class, args=(), config={}, pos=None):
    """
    Raise an exception, with proper position information if available.
    """
    raise exception_class(args, config, pos)

def vo_reraise(exc, config={}, pos=None, additional=''):
    """
    Raise an exception, with proper position information if available.

    Restores the original traceback of the exception, and should only
    be called within an "except:" block of code.
    """
    message = _format_message(str(exc), exc.__class__.__name__, config, pos)
    if message.split()[0] == str(exc).split()[0]:
        message = str(exc)
    if len(additional):
        message += ' ' + additional
    exc.args = (message,)
    raise exc, None, sys.exc_info()[2]

def vo_warn(warning_class, args=(), config={}, pos=None):
    """
    Warn, with proper position information if available.
    """
    warning = warning_class(args, config, pos)
    _suppressed_warning(warning, config)

def warn_unknown_attrs(element, attrs, config, pos, good_attr=[]):
    for attr in attrs:
        if attr not in good_attr:
            vo_warn(W48, (attr, element), config, pos)

warning_pat = re.compile(
        r":?(?P<nline>[0-9?]+):(?P<nchar>[0-9?]+): ((?P<warning>[WE]\d+): )?(?P<rest>.*)$")
def parse_vowarning(line):
    """
    Parses the vo warning string back into its parts.
    """
    result = {}
    match = warning_pat.search(line)
    if match:
        result['warning'] = warning = match.group('warning')
        if warning is not None:
            result['is_warning'] = (warning[0].upper() == 'W')
            result['is_exception'] = (warning[0].upper() == 'E')
            result['number'] = int(match.group('warning')[1:])
            if warning.startswith('W'):
                result['doc_url'] = "warnings.html#%s" % warning.lower()
            else:
                result['doc_url'] = "exceptions.html#%s" % warning.lower()
        else:
            result['is_warning'] = False
            result['is_exception'] = False
            result['is_other'] = True
            result['number'] = None
            result['doc_url'] = None
        result['nline'] = int(match.group('nline'))
        result['nchar'] = int(match.group('nchar'))
        result['message'] = match.group('rest')
        result['is_something'] = True
    else:
        result['warning'] = None
        result['is_warning'] = False
        result['is_exception'] = False
        result['is_other'] = False
        result['is_something'] = False

    return result

class VOWarning(object):
    """
    The base class of all VO warnings and exceptions.  Handles the
    formatting of the message with a warning or exception code,
    filename, line and column number.
    """
    default_args = ()

    def __init__(self, args, config={}, pos=None):
        msg = self.message % args
        self.formatted_message = _format_message(
            msg, self.__class__.__name__, config, pos)

    def __str__(self):
        return self.formatted_message

    @classmethod
    def get_short_name(cls):
        if len(cls.default_args):
            return cls.message % cls.default_args
        return cls.message

class VOTableChangeWarning(VOWarning, SyntaxWarning):
    """
    Warning emitted when a change has been made to the input XML file.
    """
    pass

class VOTableSpecWarning(VOWarning, SyntaxWarning):
    """
    Warning emitted when the input XML file violates the spec, but
    there is an obvious workaround.
    """
    pass

class UnimplementedWarning(VOWarning, SyntaxWarning):
    """
    Warning emitted when the a feature of VOTABLE is not implemented.
    """
    pass

class IOWarning(VOWarning, RuntimeWarning):
    """
    Warning emitted when a network or IO error occurred, but the
    system is able to recover using a cached copy of the data or some
    other fallback.
    """
    pass

class VOTableSpecError(VOWarning, ValueError):
    """
    Error raised when the input XML file violates the spec and there
    is no good workaround.
    """
    pass

class W01(VOTableSpecWarning):
    """
    The VOTable spec states:

        If a cell contains an array or complex number, it should be
        encoded as multiple numbers separated by whitespace.

    Many VOTable files in the wild use commas as a separator instead,
    and ``vo.table`` supports this convention when not in
    :ref:`pedantic-mode`.

    `vo.table` always outputs files using only spaces, regardless of
    how they were input.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#toc-header-35>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:TABLEDATA>`_
    """

    message = "Array uses commas rather than whitespace"

class W02(VOTableSpecWarning):
    """
    XML ids must match the following regular expression::

        ^[A-Za-z_][A-Za-z0-9_\.\-]*$

    The VOTable 1.1 says the following:

        According to the XML standard, the attribute ``ID`` is a
        string beginning with a letter or underscore (``_``), followed
        by a sequence of letters, digits, or any of the punctuation
        characters ``.`` (dot), ``-`` (dash), ``_`` (underscore), or
        ``:`` (colon).

    However, this is in conflict with the XML standard, which says
    colons may not be used.  VOTable 1.1's own schema does not allow a
    colon here.  Therefore, ``vo.table`` disallows the colon.

    VOTable 1.2 corrects this error in the specification.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `XML Names <http://www.w3.org/TR/REC-xml/#NT-Name>`_
    """

    message = "%s attribute '%s' is invalid.  Must be a standard XML id"
    default_args = ('x', 'y')

class W03(VOTableChangeWarning):
    """
    The VOTable 1.1 spec says the following about ``name`` vs. ``ID``
    on ``FIELD`` and ``VALUE`` elements:

        ``ID`` and ``name`` attributes have a different role in
        VOTable: the ``ID`` is meant as a *unique identifier* of an
        element seen as a VOTable component, while the ``name`` is
        meant for presentation purposes, and need not to be unique
        throughout the VOTable document. The ``ID`` attribute is
        therefore required in the elements which have to be
        referenced, but in principle any element may have an ``ID``
        attribute. ... In summary, the ``ID`` is different from the
        ``name`` attribute in that (a) the ``ID`` attribute is made
        from a restricted character set, and must be unique throughout
        a VOTable document whereas names are standard XML attributes
        and need not be unique; and (b) there should be support in the
        parsing software to look up references and extract the
        relevant element with matching ``ID``.

    It is further recommended in the VOTable 1.2 spec:

        While the ``ID`` attribute has to be unique in a VOTable
        document, the ``name`` attribute need not. It is however
        recommended, as a good practice, to assign unique names within
        a ``TABLE`` element. This recommendation means that, between a
        ``TABLE`` and its corresponding closing ``TABLE`` tag,
        ``name`` attributes of ``FIELD``, ``PARAM`` and optional
        ``GROUP`` elements should be all different.

    Since ``vo.table`` requires a unique identifier for each of its
    columns, ``ID`` is used for the column name when present.
    However, when ``ID`` is not present, (since it is not required by
    the specification) ``name`` is used instead.  However, ``name``
    must be cleansed by replacing invalid characters (such as
    whitespace) with underscores.

    .. note::
        This warning does not indicate that the input file is invalid
        with respect to the VOTable specification, only that the
        column names in the record array may not match exactly the
        ``name`` attributes specified in the file.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "Implictly generating an ID from a name '%s' -> '%s'"
    default_args = ('x', 'y')

class W04(VOTableSpecWarning):
    """
    The ``content-type`` attribute must use MIME content-type syntax as
    defined in `RFC 2046 <http://tools.ietf.org/html/rfc2046>`_.

    The current check for validity is somewhat over-permissive.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:link>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:link>`_
    """

    message = "content-type '%s' must be a valid MIME content type"
    default_args = ('x',)

class W05(VOTableSpecWarning):
    """
    The attribute must be a valid URI as defined in `RFC 2396
    <http://www.ietf.org/rfc/rfc2396.txt>`_.
    """

    message = "'%s' is not a valid URI"
    default_args = ('x',)

class W06(VOTableSpecWarning):
    """
    This warning is emitted when a ``ucd`` attribute does not match
    the syntax of a `unified content descriptor
    <http://vizier.u-strasbg.fr/doc/UCD.htx>`_.

    If the VOTable version is 1.2 or later, the UCD will also be
    checked to ensure it conforms to the controlled vocabulary defined
    by UCD1+.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:ucd>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:ucd>`_
    """

    message = "Invalid UCD '%s': %s"
    default_args = ('x', 'explanation')

class W07(VOTableSpecWarning):
    """
    As astro year field is a Besselian or Julian year matching the
    regular expression::

        ^[JB]?[0-9]+([.][0-9]*)?$

    Defined in this XML Schema snippet::

        <xs:simpleType  name="astroYear">
          <xs:restriction base="xs:token">
            <xs:pattern  value="[JB]?[0-9]+([.][0-9]*)?"/>
          </xs:restriction>
        </xs:simpleType>
    """

    message = "Invalid astroYear in %s: '%s'"
    default_args = ('x', 'y')

class W08(VOTableSpecWarning):
    """
    To avoid local-dependent number parsing differences, ``vo.table``
    may require a string or unicode string where a numeric type may
    make more sense.
    """

    if IS_PY3K:
        message = "'%s' must be a str or bytes object"
    else:
        message = "'%s' must be a str or unicode object"
    default_args = ('x',)

class W09(VOTableSpecWarning):
    """
    The VOTable specification uses the attribute name ``ID`` (with
    uppercase letters) to specify unique identifiers.  Some
    VOTable-producing tools use the more standard lowercase ``id``
    instead.  ``vo.table`` accepts ``id`` and emits this warning when
    not in ``pedantic`` mode.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "ID attribute not capitalized"

class W10(VOTableSpecWarning):
    """
    The parser has encountered an element that does not exist in the
    specification, or appears in an invalid context.  Check the file
    against the VOTable schema (with a tool such as `xmllint
    <http://xmlsoft.org/xmllint.html>`_.  If the file validates
    against the schema, and you still receive this warning, this may
    indicate a bug in ``vo.table``.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_
    """

    message = "Unknown tag '%s'.  Ignoring"
    default_args = ('x',)

class W11(VOTableSpecWarning):
    """
    Earlier versions of the VOTable specification used a ``gref``
    attribute on the ``LINK`` element to specify a `GLU reference
    <http://simbad3.u-strasbg.fr/glu/glu.htx>`_.  New files should
    specify a ``glu:`` protocol using the ``href`` attribute.

    Since ``vo.table`` does not currently support GLU references, it
    likewise does not automatically convert the ``gref`` attribute to
    the new form.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:link>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:link>`_
    """

    message = "The gref attribute on LINK is deprecated in VOTable 1.1"

class W12(VOTableChangeWarning):
    """
    In order to name the columns of the Numpy record array, each
    ``FIELD`` element must have either an ``ID`` or ``name`` attribute
    to derive a name from.  Strictly speaking, according to the
    VOTable schema, the ``name`` attribute is required.  However, if
    ``name`` is not present by ``ID`` is, and *pedantic mode* is off,
    ``vo.table`` will continue without a ``name`` defined.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "'%s' element must have at least one of 'ID' or 'name' attributes"
    default_args = ('x',)

class W13(VOTableSpecWarning):
    """
    Some VOTable files in the wild use non-standard datatype names.  These
    are mapped to standard ones using the following mapping::

       string        -> char
       unicodeString -> unicodeChar
       int16         -> short
       int32         -> int
       int64         -> long
       float32       -> float
       float64       -> double

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "'%s' is not a valid VOTable datatype, should be '%s'"
    default_args = ('x', 'y')

# W14: Deprecated

class W15(VOTableSpecWarning):
    """
    The ``name`` attribute is required on every ``FIELD`` element.
    However, many VOTable files in the wild omit it and provide only
    an ``ID`` instead.  In this case, when *pedantic mode* is off,
    ``vo.table`` will copy the ``name`` attribute to a new ``ID``
    attribute.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "%s element missing required 'name' attribute"
    default_args = ('x',)

# W16: Deprecated

class W17(VOTableSpecWarning):
    """
    A ``DESCRIPTION`` element can only appear once within its parent
    element.

    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_)

    However, it is a `proposed extension
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:addesc>`_
    to VOTable 1.2.
    """

    message = "%s element contains more than one DESCRIPTION element"
    default_args = ('x',)

class W18(VOTableSpecWarning):
    """
    The number of rows explicitly specified in the ``nrows`` attribute
    does not match the actual number of rows (``TR`` elements) present
    in the ``TABLE``.  This may indicate truncation of the file, or an
    internal error in the tool that produced it.  If *pedantic mode*
    is off, parsing will proceed, with the loss of some performance.

    **References:** `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC10>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC10>`_
    """

    message = 'TABLE specified nrows=%s, but table contains %s rows'
    default_args = ('x', 'y')

class W19(VOTableSpecWarning):
    """
    The column fields as defined using ``FIELD`` elements do not match
    those in the headers of the embedded FITS file.  If *pedantic
    mode* is off, the embedded FITS file will take precedence.
    """

    message = 'The fields defined in the VOTable do not match those in the embedded FITS file'

class W20(VOTableSpecWarning):
    """
    If no version number is explicitly given in the VOTable file, the
    parser assumes it is written to the VOTable 1.1 specification.
    """

    message = 'No version number specified in file.  Assuming %s'
    default_args = ('1.1',)

class W21(UnimplementedWarning):
    """
    Unknown issues may arise using ``vo.table`` with VOTable files
    from a version other than 1.1 or 1.2.
    """

    message = 'vo.table is designed for VOTable version 1.1 and 1.2, but this file is %s'
    default_args = ('x',)

class W22(VOTableSpecWarning):
    """
    Version 1.0 of the VOTable specification used the ``DEFINITIONS``
    element to define coordinate systems.  Version 1.1 now uses
    ``COOSYS`` elements throughout the document.

    **References:** `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:definitions>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:definitions>`_
    """

    message = 'The DEFINITIONS element is deprecated in VOTable 1.1.  Ignoring'

class W23(IOWarning):
    """
    Raised when the VO service database can not be updated (possibly
    due to a network outage).  This is only a warning, since an older
    and possible out-of-date VO service database was available
    locally.
    """

    message = "Unable to update service information for '%s'"
    default_args = ('x',)

class W24(VOWarning, FutureWarning):
    """
    The VO catalog database retrieved from the www is designed for a
    newer version of vo.table.  This may cause problems or limited
    features performing service queries.  Consider upgrading vo.table
    to the latest version.
    """

    message = "The VO catalog database is for a later version of vo.table"

class W25(IOWarning):
    """
    A VO service query failed due to a network error or malformed
    arguments.  Another alternative service may be attempted.  If all
    services fail, an exception will be raised.
    """

    message = "'%s' failed with: %s"
    default_args = ('service', '...')

class W26(VOTableSpecWarning):
    """
    The given element was not supported inside of the given element
    until the specified VOTable version, however the version declared
    in the file is for an earlier version.  These attributes may not
    be written out to the file.
    """

    message = "'%s' inside '%s' added in VOTable %s"
    default_args = ('child', 'parent', 'X.X')

class W27(VOTableSpecWarning):
    """
    The ``COOSYS`` element was deprecated in VOTABLE version 1.2 in
    favor of a reference to the Space-Time Coordinate (STC) data
    model (see `utype
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:utype>`_
    and the IVOA note `referencing STC in VOTable
    <http://ivoa.net/Documents/latest/VOTableSTC.html>`_.
    """

    message = "COOSYS deprecated in VOTable 1.2"

class W28(VOTableSpecWarning):
    """
    The given attribute was not supported on the given element until the
    specified VOTable version, however the version declared in the file is
    for an earlier version.  These attributes may not be written out to
    the file.
    """

    message = "'%s' on '%s' added in VOTable %s"
    default_args = ('attribute', 'element', 'X.X')

class W29(VOTableSpecWarning):
    """
    Some VOTable files specify their version number in the form "v1.0",
    when the only supported forms in the spec are "1.0".

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_
    """

    message = "Version specified in non-standard form '%s'"
    default_args = ('v1.0',)

class W30(VOTableSpecWarning):
    """
    Some VOTable files write missing floating-point values in non-standard
    ways, such as "null" and "-".  In non-pedantic mode, any non-standard
    floating-point literals are treated as missing values.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "Invalid literal for float '%s'.  Treating as empty."
    default_args = ('x',)

class W31(VOTableSpecWarning):
    """
    Since NaN's can not be represented in integer fields directly, a null
    value must be specified in the FIELD descriptor to support reading
    NaN's from the tabledata.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "NaN given in an integral field without a specified null value"

class W32(VOTableSpecWarning):
    """
    Each field in a table must have a unique ID.  If two or more fields
    have the same ID, some will be renamed to ensure that all IDs are
    unique.

    From the VOTable 1.2 spec:

        The ``ID`` and ``ref`` attributes are defined as XML types
        ``ID`` and ``IDREF`` respectively. This means that the
        contents of ``ID`` is an identifier which must be unique
        throughout a VOTable document, and that the contents of the
        ``ref`` attribute represents a reference to an identifier
        which must exist in the VOTable document.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "Duplicate ID '%s' renamed to '%s' to ensure uniqueness"
    default_args = ('x', 'x_2')

class W33(VOTableChangeWarning):
    """
    Each field in a table must have a unique name.  If two or more
    fields have the same name, some will be renamed to ensure that all
    names are unique.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_
    """

    message = "Column name '%s' renamed to '%s' to ensure uniqueness"
    default_args = ('x', 'x_2')

class W34(VOTableSpecWarning):
    """
    The attribute requires the value to be a valid XML token, as
    defined by `XML 1.0
    <http://www.w3.org/TR/2000/WD-xml-2e-20000814#NT-Nmtoken>`_.
    """

    message = "'%s' is an invalid token for attribute '%s'"
    default_args = ('x', 'y')

class W35(VOTableSpecWarning):
    """
    The ``name`` and ``value`` attributes are required on all ``INFO``
    elements.

    **References:** `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC32>`
    """

    message = "'%s' attribute required for INFO elements"
    default_args = ('x',)

class W36(VOTableSpecWarning):
    """
    If the field specifies a ``null`` value, that value must conform
    to the given ``datatype``.

    **References:** `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`
    """

    message = "null value '%s' does not match field datatype, setting to 0"
    default_args = ('x',)

class W37(UnimplementedWarning):
    """
    The 3 datatypes defined in the VOTable specification and supported by
    vo.table are ``TABLEDATA``, ``BINARY`` and ``FITS``.

    **References:** `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:data>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:data>`
    """

    message = "Unsupported data format '%s'"
    default_args = ('x',)

class W38(VOTableSpecWarning):
    """
    The only encoding for local binary data supported by the VOTable
    specification is base64.
    """

    message = "Inline binary data must be base64 encoded, got '%s'"
    default_args = ('x',)

class W39(VOTableSpecWarning):
    """
    Bit values do not support masking.  This warning is raised upon
    setting masked data in a bit column.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "Bit values can not be masked"

class W40(VOTableSpecWarning):
    """
    This is a terrible hack to support Simple Image Access Protocol
    results from `<archive.noao.edu>`_.  It creates a field for the
    coordinate projection type of type "double", which
    actually contains character data.  We have to hack the field
    to store character data, or we can't read it in.  A warning
    will be raised when this happens.
    """

    message = "'cprojection' datatype repaired"

class W41(VOTableSpecWarning):
    """
    An XML namespace was specified on the ``VOTABLE`` element, but the
    namespace does not match what is expected for a ``VOTABLE`` file.

    The ``VOTABLE`` namespace is::

      http://www.ivoa.net/xml/VOTable/vX.X

    where "X.X" is the version number.

    Some files in the wild set the namespace to the location of the
    VOTable schema, which is not correct and will not pass some
    validating parsers.
    """

    message = "An XML namespace is specified, but is incorrect.  Expected '%s', got '%s'"
    default_args = ('x', 'y')

class W42(VOTableSpecWarning):
    """
    The root element should specify a namespace.

    The ``VOTABLE`` namespace is::

        http://www.ivoa.net/xml/VOTable/vX.X

    where "X.X" is the version number.
    """

    message = "No XML namespace specified"

class W43(VOTableSpecWarning):
    """
    Referenced elements should be defined before referees.  From the
    VOTable 1.2 spec:

       In VOTable1.2, it is further recommended to place the ID
       attribute prior to referencing it whenever possible.
    """

    message = "%s ref='%s' which has not already been defined"
    default_args = ('element', 'x',)

class W44(VOTableSpecWarning):
    """
    ``VALUES`` elements that reference another element should not have
    their own content.

    From the VOTable 1.2 spec:

        The ``ref`` attribute of a ``VALUES`` element can be used to
        avoid a repetition of the domain definition, by referring to a
        previously defined ``VALUES`` element having the referenced
        ``ID`` attribute. When specified, the ``ref`` attribute
        defines completely the domain without any other element or
        attribute, as e.g. ``<VALUES ref="RAdomain"/>``
    """

    message = "VALUES element with ref attribute has content ('%s')"
    default_args = ('element',)

class W45(VOWarning, ValueError):
    """
    The ``content-role`` attribute on the ``LINK`` element must be one of
    the following::

        query, hints, doc, location

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_
    """

    message = "content-role attribute '%s' invalid"
    default_args = ('x',)

class W46(VOTableSpecWarning):
    """
    The given char or unicode string is too long for the specified
    field length.
    """

    message = "%s value is too long for specified length of %s"
    default_args = ('char or unicode', 'x')

class W47(VOTableSpecWarning):
    """
    If no arraysize is specified on a char field, the default of '1'
    is implied, but this is rarely what is intended.
    """

    message = "Missing arraysize indicates length 1"

class W48(VOTableSpecWarning):
    """
    The attribute is not defined in the specification.
    """

    message = "Unknown attribute '%s' on %s"
    default_args = ('attribute', 'element')

class W49(VOTableSpecWarning):
    """
    Empty cell illegal for integer fields.

    If a "null" value was specified for the cell, it will be used for the value,
    otherwise, 0 will be used.
    """

    message = "Empty cell illegal for integer fields."
    
class E01(VOWarning, ValueError):
    """
    The size specifier for a ``char`` or ``unicode`` field must be
    only a number followed, optionally, by an asterisk.
    Multi-dimensional size specifiers are not supported for these
    datatypes.

    Strings, which are defined as a set of characters, can be
    represented in VOTable as a fixed- or variable-length array of
    characters::

        <FIELD name="unboundedString" datatype="char" arraysize="*"/>

    A 1D array of strings can be represented as a 2D array of
    characters, but given the logic above, it is possible to define a
    variable-length array of fixed-length strings, but not a
    fixed-length array of variable-length strings.
    """

    message = "Invalid size specifier '%s' for a %s field (in field '%s')"
    default_args = ('x', 'char/unicode', 'y')

class E02(VOWarning, ValueError):
    """
    The number of array elements in the data does not match that specified
    in the FIELD specifier.
    """

    message = "Incorrect number of elements in array.  Expected multiple of %s, got %s"
    default_args = ('x', 'y')

class E03(VOWarning, ValueError):
    """
    Complex numbers should be two values separated by whitespace.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "'%s' does not parse as a complex number"
    default_args = ('x',)

class E04(VOWarning, ValueError):
    """
    A ``bit`` array should be a string of '0's and '1's.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "Invalid bit value '%s'"
    default_args = ('x',)

class E05(VOWarning, ValueError):
    """
    A ``boolean`` value should be one of the following strings (case
    insensitive) in the ``TABLEDATA`` format::

        'TRUE', 'FALSE', '1', '0', 'T', 'F', '\\0', ' ', '?'

    and in ``BINARY`` format::

        'T', 'F', '1', '0', '\\0', ' ', '?'

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "Invalid boolean value '%s'"
    default_args = ('x',)

class E06(VOWarning, ValueError):
    """
    The supported datatypes are::

        double, float, bit, boolean, unsignedByte, short, int, long,
        floatComplex, doubleComplex, char, unicodeChar

    The following non-standard aliases are also supported, but in
    these case :ref:`W13 <W13>` will be raised::

        string        -> char
        unicodeString -> unicodeChar
        int16         -> short
        int32         -> int
        int64         -> long
        float32       -> float
        float64       -> double

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_
    """

    message = "Unknown datatype '%s' on field '%s'"
    default_args = ('x', 'y')

# E07: Deprecated

class E08(VOWarning, ValueError):
    """
    The ``type`` attribute on the ``VALUES`` element must be either
    ``legal`` or ``actual``.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`_
    """

    message = "type must be 'legal' or 'actual', but is '%s'"
    default_args = ('x',)

class E09(VOWarning, ValueError):
    """
    The ``MIN``, ``MAX`` and ``OPTION`` elements must always have a
    ``value`` attribute.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`_
    """

    message = "'%s' must have a value attribute"
    default_args = ('x',)

class E10(VOWarning, ValueError):
    """
    From VOTable 1.1 and later, ``FIELD`` and ``PARAM`` elements must have
    a ``datatype`` field.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`_
    """

    message = "'datatype' attribute required on all '%s' elements"
    default_args = ('FIELD',)

class E11(VOWarning, ValueError):
    """
    The precision attribute is meant to express the number of significant
    digits, either as a number of decimal places (e.g. ``precision="F2"`` or
    equivalently ``precision="2"`` to express 2 significant figures
    after the decimal point), or as a number of significant figures
    (e.g. ``precision="E5"`` indicates a relative precision of 10-5).

    It is validated using the following regular expression::

        [EF]?[1-9][0-9]*

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:form>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:form>`_
    """

    message = "precision '%s' is invalid"
    default_args = ('x',)

class E12(VOWarning, ValueError):
    """
    The width attribute is meant to indicate to the application the
    number of characters to be used for input or output of the
    quantity.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:form>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:form>`_
    """

    message = "width must be a positive integer, got '%s'"
    default_args = ('x',)

class E13(VOWarning, ValueError):
    u"""
    From the VOTable 1.2 spec:

        A table cell can contain an array of a given primitive type,
        with a fixed or variable number of elements; the array may
        even be multidimensional. For instance, the position of a
        point in a 3D space can be defined by the following::

            <FIELD ID="point_3D" datatype="double" arraysize="3"/>

        and each cell corresponding to that definition must contain
        exactly 3 numbers. An asterisk (\*) may be appended to
        indicate a variable number of elements in the array, as in::

            <FIELD ID="values" datatype="int" arraysize="100*"/>

        where it is specified that each cell corresponding to that
        definition contains 0 to 100 integer numbers. The number may
        be omitted to specify an unbounded array (in practice up to
        =~2×10⁹ elements).

        A table cell can also contain a multidimensional array of a
        given primitive type. This is specified by a sequence of
        dimensions separated by the ``x`` character, with the first
        dimension changing fastest; as in the case of a simple array,
        the last dimension may be variable in length. As an example,
        the following definition declares a table cell which may
        contain a set of up to 10 images, each of 64×64 bytes::

            <FIELD ID="thumbs" datatype="unsignedByte" arraysize="64×64×10*"/>

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:dim>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:dim>`_
    """

    message = "Invalid arraysize attribute '%s'"
    default_args = ('x',)

class E14(VOWarning, ValueError):
    """
    All ``PARAM`` elements must have a ``value`` attribute.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`_
    """

    message = "value attribute is required for all PARAM elements"

class E15(VOWarning, ValueError):
    """
    All ``COOSYS`` elements must have an ``ID`` attribute.

    Note that the VOTable 1.1 specification says this attribute is
    optional, but its corresponding schema indicates it is required.

    In VOTable 1.2, the ``COOSYS`` element is deprecated.
    """

    message = "ID attribute is required for all COOSYS elements"

class E16(VOTableSpecWarning):
    """
    The ``system`` attribute on the ``COOSYS`` element must be one of the
    following::

      'eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5', 'galactic',
      'supergalactic', 'xy', 'barycentric', 'geo_app'

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:COOSYS>`_
    """

    message = "Invalid system attribute '%s'"
    default_args = ('x',)

class E17(VOWarning, ValueError):
    """
    ``extnum`` attribute must be a positive integer.

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_
    """

    message = "extnum must be a positive integer"

class E18(VOWarning, ValueError):
    """
    The ``type`` attribute of the ``RESOURCE`` element must be one of
    "results" or "meta".

    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_
    """

    message = "type must be 'results' or 'meta', not '%s'"
    default_args = ('x',)

class E19(VOWarning, ValueError):
    """
    Raised either when the file doesn't appear to be XML, or the root
    element is not VOTABLE.
    """

    message = "File does not appear to be a VOTABLE"

class E20(VOTableSpecError):
    """
    The table had only *x* fields defined, but the data itself has more
    columns than that.
    """

    message = "Data has more columns than are defined in the header (%s)"
    default_args = ('x',)

class E21(VOWarning, ValueError):
    """
    The table had *x* fields defined, but the data itself has only *y*
    columns.
    """

    message = "Data has fewer columns (%s) than are defined in the header (%s)"
    default_args = ('x', 'y')
