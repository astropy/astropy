# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
.. _warnings:

Warnings
--------

.. note::
    Most of the following warnings indicate violations of the VOTable
    specification.  They should be reported to the authors of the
    tools that produced the VOTable file.

    To control the warnings emitted, use the standard Python
    :mod:`warnings` module.  Most of these are of the type
    `VOTableSpecWarning`.

{warnings}

.. _exceptions:

Exceptions
----------

.. note::

    This is a list of many of the fatal exceptions emitted by ``astropy.io.votable``
    when the file does not conform to spec.  Other exceptions may be
    raised due to unforeseen cases or bugs in ``astropy.io.votable`` itself.

{exceptions}
"""


# STDLIB
import io
import re

from textwrap import dedent
from warnings import warn

from astropy.utils.exceptions import AstropyWarning


__all__ = [
    'warn_or_raise', 'vo_raise', 'vo_reraise', 'vo_warn',
    'warn_unknown_attrs', 'parse_vowarning', 'VOWarning',
    'VOTableChangeWarning', 'VOTableSpecWarning',
    'UnimplementedWarning', 'IOWarning', 'VOTableSpecError']


MAX_WARNINGS = 10


def _format_message(message, name, config=None, pos=None):
    if config is None:
        config = {}
    if pos is None:
        pos = ('?', '?')
    filename = config.get('filename', '?')
    return '{}:{}:{}: {}: {}'.format(filename, pos[0], pos[1], name, message)


def _suppressed_warning(warning, config, stacklevel=2):
    warning_class = type(warning)
    config.setdefault('_warning_counts', dict()).setdefault(warning_class, 0)
    config['_warning_counts'][warning_class] += 1
    message_count = config['_warning_counts'][warning_class]
    if message_count <= MAX_WARNINGS:
        if message_count == MAX_WARNINGS:
            warning.formatted_message += \
                ' (suppressing further warnings of this type...)'
        warn(warning, stacklevel=stacklevel+1)


def warn_or_raise(warning_class, exception_class=None, args=(), config=None,
                  pos=None, stacklevel=1):
    """
    Warn or raise an exception, depending on the verify setting.
    """
    if config is None:
        config = {}
    # NOTE: the default here is deliberately warn rather than ignore, since
    # one would expect that calling warn_or_raise without config should not
    # silence the warnings.
    config_value = config.get('verify', 'warn')
    if config_value == 'exception':
        if exception_class is None:
            exception_class = warning_class
        vo_raise(exception_class, args, config, pos)
    elif config_value == 'warn':
        vo_warn(warning_class, args, config, pos, stacklevel=stacklevel+1)


def vo_raise(exception_class, args=(), config=None, pos=None):
    """
    Raise an exception, with proper position information if available.
    """
    if config is None:
        config = {}
    raise exception_class(args, config, pos)


def vo_reraise(exc, config=None, pos=None, additional=''):
    """
    Raise an exception, with proper position information if available.

    Restores the original traceback of the exception, and should only
    be called within an "except:" block of code.
    """
    if config is None:
        config = {}
    message = _format_message(str(exc), exc.__class__.__name__, config, pos)
    if message.split()[0] == str(exc).split()[0]:
        message = str(exc)
    if len(additional):
        message += ' ' + additional
    exc.args = (message,)
    raise exc


def vo_warn(warning_class, args=(), config=None, pos=None, stacklevel=1):
    """
    Warn, with proper position information if available.
    """
    if config is None:
        config = {}
    # NOTE: the default here is deliberately warn rather than ignore, since
    # one would expect that calling warn_or_raise without config should not
    # silence the warnings.
    if config.get('verify', 'warn') != 'ignore':
        warning = warning_class(args, config, pos)
        _suppressed_warning(warning, config, stacklevel=stacklevel+1)


def warn_unknown_attrs(element, attrs, config, pos, good_attr=[], stacklevel=1):
    for attr in attrs:
        if attr not in good_attr:
            vo_warn(W48, (attr, element), config, pos, stacklevel=stacklevel+1)


_warning_pat = re.compile(
    r":?(?P<nline>[0-9?]+):(?P<nchar>[0-9?]+): " +
     r"((?P<warning>[WE]\d+): )?(?P<rest>.*)$")


def parse_vowarning(line):
    """
    Parses the vo warning string back into its parts.
    """
    result = {}
    match = _warning_pat.search(line)
    if match:
        result['warning'] = warning = match.group('warning')
        if warning is not None:
            result['is_warning'] = (warning[0].upper() == 'W')
            result['is_exception'] = not result['is_warning']
            result['number'] = int(match.group('warning')[1:])
            result['doc_url'] = "io/votable/api_exceptions.html#{}".format(
                warning.lower())
        else:
            result['is_warning'] = False
            result['is_exception'] = False
            result['is_other'] = True
            result['number'] = None
            result['doc_url'] = None
        try:
            result['nline'] = int(match.group('nline'))
        except ValueError:
            result['nline'] = 0
        try:
            result['nchar'] = int(match.group('nchar'))
        except ValueError:
            result['nchar'] = 0
        result['message'] = match.group('rest')
        result['is_something'] = True
    else:
        result['warning'] = None
        result['is_warning'] = False
        result['is_exception'] = False
        result['is_other'] = False
        result['is_something'] = False
        if not isinstance(line, str):
            line = line.decode('utf-8')
        result['message'] = line

    return result


class VOWarning(AstropyWarning):
    """
    The base class of all VO warnings and exceptions.

    Handles the formatting of the message with a warning or exception
    code, filename, line and column number.
    """
    default_args = ()
    message_template = ''

    def __init__(self, args, config=None, pos=None):
        if config is None:
            config = {}
        if not isinstance(args, tuple):
            args = (args, )
        msg = self.message_template.format(*args)

        self.formatted_message = _format_message(
            msg, self.__class__.__name__, config, pos)
        Warning.__init__(self, self.formatted_message)

    def __str__(self):
        return self.formatted_message

    @classmethod
    def get_short_name(cls):
        if len(cls.default_args):
            return cls.message_template.format(*cls.default_args)
        return cls.message_template


class VOTableChangeWarning(VOWarning, SyntaxWarning):
    """
    A change has been made to the input XML file.
    """


class VOTableSpecWarning(VOWarning, SyntaxWarning):
    """
    The input XML file violates the spec, but there is an obvious workaround.
    """


class UnimplementedWarning(VOWarning, SyntaxWarning):
    """
    A feature of the VOTABLE_ spec is not implemented.
    """


class IOWarning(VOWarning, RuntimeWarning):
    """
    A network or IO error occurred, but was recovered using the cache.
    """


class VOTableSpecError(VOWarning, ValueError):
    """
    The input XML file violates the spec and there is no good workaround.
    """


class W01(VOTableSpecWarning):
    """
    The VOTable spec states:

        If a cell contains an array or complex number, it should be
        encoded as multiple numbers separated by whitespace.

    Many VOTable files in the wild use commas as a separator instead,
    and ``astropy.io.votable`` supports this convention when not in
    :ref:`pedantic-mode`.

    ``astropy.io.votable`` always outputs files using only spaces, regardless of
    how they were input.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#toc-header-35>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:TABLEDATA>`__
    """

    message_template = "Array uses commas rather than whitespace"


class W02(VOTableSpecWarning):
    r"""
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
    colon here.  Therefore, ``astropy.io.votable`` disallows the colon.

    VOTable 1.2 corrects this error in the specification.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `XML Names <http://www.w3.org/TR/REC-xml/#NT-Name>`__
    """

    message_template = "{} attribute '{}' is invalid.  Must be a standard XML id"
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

    Since ``astropy.io.votable`` requires a unique identifier for each of its
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
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = "Implicitly generating an ID from a name '{}' -> '{}'"
    default_args = ('x', 'y')


class W04(VOTableSpecWarning):
    """
    The ``content-type`` attribute must use MIME content-type syntax as
    defined in `RFC 2046 <https://tools.ietf.org/html/rfc2046>`__.

    The current check for validity is somewhat over-permissive.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:link>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:link>`__
    """

    message_template = "content-type '{}' must be a valid MIME content type"
    default_args = ('x',)


class W05(VOTableSpecWarning):
    """
    The attribute must be a valid URI as defined in `RFC 2396
    <https://www.ietf.org/rfc/rfc2396.txt>`_.
    """

    message_template = "'{}' is not a valid URI"
    default_args = ('x',)


class W06(VOTableSpecWarning):
    """
    This warning is emitted when a ``ucd`` attribute does not match
    the syntax of a `unified content descriptor
    <http://vizier.u-strasbg.fr/doc/UCD.htx>`__.

    If the VOTable version is 1.2 or later, the UCD will also be
    checked to ensure it conforms to the controlled vocabulary defined
    by UCD1+.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:ucd>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:ucd>`__
    """

    message_template = "Invalid UCD '{}': {}"
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

    message_template = "Invalid astroYear in {}: '{}'"
    default_args = ('x', 'y')


class W08(VOTableSpecWarning):
    """
    To avoid local-dependent number parsing differences, ``astropy.io.votable``
    may require a string or unicode string where a numeric type may
    make more sense.
    """

    message_template = "'{}' must be a str or bytes object"

    default_args = ('x',)


class W09(VOTableSpecWarning):
    """
    The VOTable specification uses the attribute name ``ID`` (with
    uppercase letters) to specify unique identifiers.  Some
    VOTable-producing tools use the more standard lowercase ``id``
    instead. ``astropy.io.votable`` accepts ``id`` and emits this warning if
    ``verify`` is ``'warn'``.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = "ID attribute not capitalized"


class W10(VOTableSpecWarning):
    """
    The parser has encountered an element that does not exist in the
    specification, or appears in an invalid context.  Check the file
    against the VOTable schema (with a tool such as `xmllint
    <http://xmlsoft.org/xmllint.html>`__.  If the file validates
    against the schema, and you still receive this warning, this may
    indicate a bug in ``astropy.io.votable``.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__
    """

    message_template = "Unknown tag '{}'.  Ignoring"
    default_args = ('x',)


class W11(VOTableSpecWarning):
    """
    Earlier versions of the VOTable specification used a ``gref``
    attribute on the ``LINK`` element to specify a `GLU reference
    <http://aladin.u-strasbg.fr/glu/>`__.  New files should
    specify a ``glu:`` protocol using the ``href`` attribute.

    Since ``astropy.io.votable`` does not currently support GLU references, it
    likewise does not automatically convert the ``gref`` attribute to
    the new form.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:link>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:link>`__
    """

    message_template = "The gref attribute on LINK is deprecated in VOTable 1.1"


class W12(VOTableChangeWarning):
    """
    In order to name the columns of the Numpy record array, each
    ``FIELD`` element must have either an ``ID`` or ``name`` attribute
    to derive a name from.  Strictly speaking, according to the
    VOTable schema, the ``name`` attribute is required.  However, if
    ``name`` is not present by ``ID`` is, and ``verify`` is not ``'exception'``,
    ``astropy.io.votable`` will continue without a ``name`` defined.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = (
        "'{}' element must have at least one of 'ID' or 'name' attributes")
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
       unsignedInt   -> long
       unsignedShort -> int

    To add more datatype mappings during parsing, use the
    ``datatype_mapping`` keyword to `astropy.io.votable.parse`.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "'{}' is not a valid VOTable datatype, should be '{}'"
    default_args = ('x', 'y')


# W14: Deprecated


class W15(VOTableSpecWarning):
    """
    The ``name`` attribute is required on every ``FIELD`` element.
    However, many VOTable files in the wild omit it and provide only
    an ``ID`` instead.  In this case, when ``verify`` is not ``'exception'``
    ``astropy.io.votable`` will copy the ``name`` attribute to a new ``ID``
    attribute.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = "{} element missing required 'name' attribute"
    default_args = ('x',)

# W16: Deprecated


class W17(VOTableSpecWarning):
    """
    A ``DESCRIPTION`` element can only appear once within its parent
    element.

    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__)

    However, it is a `proposed extension
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:addesc>`__
    to VOTable 1.2.
    """

    message_template = "{} element contains more than one DESCRIPTION element"
    default_args = ('x',)


class W18(VOTableSpecWarning):
    """
    The number of rows explicitly specified in the ``nrows`` attribute
    does not match the actual number of rows (``TR`` elements) present
    in the ``TABLE``.  This may indicate truncation of the file, or an
    internal error in the tool that produced it.  If ``verify`` is not
    ``'exception'``, parsing will proceed, with the loss of some performance.

    **References:** `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC10>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC10>`__
    """

    message_template = 'TABLE specified nrows={}, but table contains {} rows'
    default_args = ('x', 'y')


class W19(VOTableSpecWarning):
    """
    The column fields as defined using ``FIELD`` elements do not match
    those in the headers of the embedded FITS file.  If ``verify`` is not
    ``'exception'``, the embedded FITS file will take precedence.
    """

    message_template = (
        'The fields defined in the VOTable do not match those in the ' +
        'embedded FITS file')


class W20(VOTableSpecWarning):
    """
    If no version number is explicitly given in the VOTable file, the
    parser assumes it is written to the VOTable 1.1 specification.
    """

    message_template = 'No version number specified in file.  Assuming {}'
    default_args = ('1.1',)


class W21(UnimplementedWarning):
    """
    Unknown issues may arise using ``astropy.io.votable`` with VOTable files
    from a version other than 1.1, 1.2, 1.3, or 1.4.
    """

    message_template = (
        'astropy.io.votable is designed for VOTable version 1.1, 1.2, 1.3,'
        ' and 1.4, but this file is {}')
    default_args = ('x',)


class W22(VOTableSpecWarning):
    """
    Version 1.0 of the VOTable specification used the ``DEFINITIONS``
    element to define coordinate systems.  Version 1.1 now uses
    ``COOSYS`` elements throughout the document.

    **References:** `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:definitions>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:definitions>`__
    """

    message_template = 'The DEFINITIONS element is deprecated in VOTable 1.1.  Ignoring'


class W23(IOWarning):
    """
    Raised when the VO service database can not be updated (possibly
    due to a network outage).  This is only a warning, since an older
    and possible out-of-date VO service database was available
    locally.
    """

    message_template = "Unable to update service information for '{}'"
    default_args = ('x',)


class W24(VOWarning, FutureWarning):
    """
    The VO catalog database retrieved from the www is designed for a
    newer version of ``astropy.io.votable``.  This may cause problems or limited
    features performing service queries.  Consider upgrading ``astropy.io.votable``
    to the latest version.
    """

    message_template = "The VO catalog database is for a later version of astropy.io.votable"


class W25(IOWarning):
    """
    A VO service query failed due to a network error or malformed
    arguments.  Another alternative service may be attempted.  If all
    services fail, an exception will be raised.
    """

    message_template = "'{}' failed with: {}"
    default_args = ('service', '...')


class W26(VOTableSpecWarning):
    """
    The given element was not supported inside of the given element
    until the specified VOTable version, however the version declared
    in the file is for an earlier version.  These attributes may not
    be written out to the file.
    """

    message_template = "'{}' inside '{}' added in VOTable {}"
    default_args = ('child', 'parent', 'X.X')


class W27(VOTableSpecWarning):
    """
    The ``COOSYS`` element was deprecated in VOTABLE version 1.2 in
    favor of a reference to the Space-Time Coordinate (STC) data
    model (see `utype
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:utype>`__
    and the IVOA note `referencing STC in VOTable
    <http://ivoa.net/Documents/latest/VOTableSTC.html>`__.
    """

    message_template = "COOSYS deprecated in VOTable 1.2"


class W28(VOTableSpecWarning):
    """
    The given attribute was not supported on the given element until the
    specified VOTable version, however the version declared in the file is
    for an earlier version.  These attributes may not be written out to
    the file.
    """

    message_template = "'{}' on '{}' added in VOTable {}"
    default_args = ('attribute', 'element', 'X.X')


class W29(VOTableSpecWarning):
    """
    Some VOTable files specify their version number in the form "v1.0",
    when the only supported forms in the spec are "1.0".

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__
    """

    message_template = "Version specified in non-standard form '{}'"
    default_args = ('v1.0',)


class W30(VOTableSpecWarning):
    """
    Some VOTable files write missing floating-point values in non-standard ways,
    such as "null" and "-".  If ``verify`` is not ``'exception'``, any
    non-standard floating-point literals are treated as missing values.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "Invalid literal for float '{}'.  Treating as empty."
    default_args = ('x',)


class W31(VOTableSpecWarning):
    """
    Since NaN's can not be represented in integer fields directly, a null
    value must be specified in the FIELD descriptor to support reading
    NaN's from the tabledata.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "NaN given in an integral field without a specified null value"


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
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = "Duplicate ID '{}' renamed to '{}' to ensure uniqueness"
    default_args = ('x', 'x_2')


class W33(VOTableChangeWarning):
    """
    Each field in a table must have a unique name.  If two or more
    fields have the same name, some will be renamed to ensure that all
    names are unique.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`__
    """

    message_template = "Column name '{}' renamed to '{}' to ensure uniqueness"
    default_args = ('x', 'x_2')


class W34(VOTableSpecWarning):
    """
    The attribute requires the value to be a valid XML token, as
    defined by `XML 1.0
    <http://www.w3.org/TR/2000/WD-xml-2e-20000814#NT-Nmtoken>`__.
    """

    message_template = "'{}' is an invalid token for attribute '{}'"
    default_args = ('x', 'y')


class W35(VOTableSpecWarning):
    """
    The ``name`` and ``value`` attributes are required on all ``INFO``
    elements.

    **References:** `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC32>`__
    """

    message_template = "'{}' attribute required for INFO elements"
    default_args = ('x',)


class W36(VOTableSpecWarning):
    """
    If the field specifies a ``null`` value, that value must conform
    to the given ``datatype``.

    **References:** `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`__
    """

    message_template = "null value '{}' does not match field datatype, setting to 0"
    default_args = ('x',)


class W37(UnimplementedWarning):
    """
    The 3 datatypes defined in the VOTable specification and supported by
    ``astropy.io.votable`` are ``TABLEDATA``, ``BINARY`` and ``FITS``.

    **References:** `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:data>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:data>`__
    """

    message_template = "Unsupported data format '{}'"
    default_args = ('x',)


class W38(VOTableSpecWarning):
    """
    The only encoding for local binary data supported by the VOTable
    specification is base64.
    """

    message_template = "Inline binary data must be base64 encoded, got '{}'"
    default_args = ('x',)


class W39(VOTableSpecWarning):
    """
    Bit values do not support masking.  This warning is raised upon
    setting masked data in a bit column.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "Bit values can not be masked"


class W40(VOTableSpecWarning):
    """
    This is a terrible hack to support Simple Image Access Protocol
    results from `archive.noao.edu <http://archive.noao.edu>`__.  It
    creates a field for the coordinate projection type of type "double",
    which actually contains character data.  We have to hack the field
    to store character data, or we can't read it in.  A warning will be
    raised when this happens.
    """

    message_template = "'cprojection' datatype repaired"


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

    message_template = (
        "An XML namespace is specified, but is incorrect.  Expected " +
        "'{}', got '{}'")
    default_args = ('x', 'y')


class W42(VOTableSpecWarning):
    """
    The root element should specify a namespace.

    The ``VOTABLE`` namespace is::

        http://www.ivoa.net/xml/VOTable/vX.X

    where "X.X" is the version number.
    """

    message_template = "No XML namespace specified"


class W43(VOTableSpecWarning):
    """
    Referenced elements should be defined before referees.  From the
    VOTable 1.2 spec:

       In VOTable1.2, it is further recommended to place the ID
       attribute prior to referencing it whenever possible.
    """

    message_template = "{} ref='{}' which has not already been defined"
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

    message_template = "VALUES element with ref attribute has content ('{}')"
    default_args = ('element',)


class W45(VOWarning, ValueError):
    """
    The ``content-role`` attribute on the ``LINK`` element must be one of
    the following::

        query, hints, doc, location

    And in VOTable 1.3, additionally::

        type

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__
    `1.3
    <http://www.ivoa.net/documents/VOTable/20130315/PR-VOTable-1.3-20130315.html#sec:link>`__
    """

    message_template = "content-role attribute '{}' invalid"
    default_args = ('x',)


class W46(VOTableSpecWarning):
    """
    The given char or unicode string is too long for the specified
    field length.
    """

    message_template = "{} value is too long for specified length of {}"
    default_args = ('char or unicode', 'x')


class W47(VOTableSpecWarning):
    """
    If no arraysize is specified on a char field, the default of '1'
    is implied, but this is rarely what is intended.
    """

    message_template = "Missing arraysize indicates length 1"


class W48(VOTableSpecWarning):
    """
    The attribute is not defined in the specification.
    """

    message_template = "Unknown attribute '{}' on {}"
    default_args = ('attribute', 'element')


class W49(VOTableSpecWarning):
    """
    Prior to VOTable 1.3, the empty cell was illegal for integer
    fields.

    If a \"null\" value was specified for the cell, it will be used
    for the value, otherwise, 0 will be used.
    """

    message_template = "Empty cell illegal for integer fields."


class W50(VOTableSpecWarning):
    """
    Invalid unit string as defined in the `Standards for Astronomical
    Catalogues, Version 2.0
    <http://cdsarc.u-strasbg.fr/doc/catstd-3.2.htx>`_.

    Consider passing an explicit ``unit_format`` parameter if the units
    in this file conform to another specification.
    """

    message_template = "Invalid unit string '{}'"
    default_args = ('x',)


class W51(VOTableSpecWarning):
    """
    The integer value is out of range for the size of the field.
    """

    message_template = "Value '{}' is out of range for a {} integer field"
    default_args = ('x', 'n-bit')


class W52(VOTableSpecWarning):
    """
    The BINARY2 format was introduced in VOTable 1.3.  It should
    not be present in files marked as an earlier version.
    """

    message_template = ("The BINARY2 format was introduced in VOTable 1.3, but "
               "this file is declared as version '{}'")
    default_args = ('1.2',)


class W53(VOTableSpecWarning):
    """
    The VOTABLE element must contain at least one RESOURCE element.
    """

    message_template = ("VOTABLE element must contain at least one RESOURCE element.")
    default_args = ()


class W54(VOTableSpecWarning):
    """
    The TIMESYS element was introduced in VOTable 1.4.  It should
    not be present in files marked as an earlier version.
    """

    message_template = (
        "The TIMESYS element was introduced in VOTable 1.4, but "
        "this file is declared as version '{}'")
    default_args = ('1.3',)


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

    message_template = "Invalid size specifier '{}' for a {} field (in field '{}')"
    default_args = ('x', 'char/unicode', 'y')


class E02(VOWarning, ValueError):
    """
    The number of array elements in the data does not match that specified
    in the FIELD specifier.
    """

    message_template = (
        "Incorrect number of elements in array. " +
        "Expected multiple of {}, got {}")
    default_args = ('x', 'y')


class E03(VOWarning, ValueError):
    """
    Complex numbers should be two values separated by whitespace.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "'{}' does not parse as a complex number"
    default_args = ('x',)


class E04(VOWarning, ValueError):
    """
    A ``bit`` array should be a string of '0's and '1's.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "Invalid bit value '{}'"
    default_args = ('x',)


class E05(VOWarning, ValueError):
    r"""
    A ``boolean`` value should be one of the following strings (case
    insensitive) in the ``TABLEDATA`` format::

        'TRUE', 'FALSE', '1', '0', 'T', 'F', '\0', ' ', '?'

    and in ``BINARY`` format::

        'T', 'F', '1', '0', '\0', ' ', '?'

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "Invalid boolean value '{}'"
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
        unsignedInt   -> long
        unsignedShort -> int

    To add more datatype mappings during parsing, use the
    ``datatype_mapping`` keyword to `astropy.io.votable.parse`.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`__
    """

    message_template = "Unknown datatype '{}' on field '{}'"
    default_args = ('x', 'y')

# E07: Deprecated


class E08(VOWarning, ValueError):
    """
    The ``type`` attribute on the ``VALUES`` element must be either
    ``legal`` or ``actual``.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`__
    """

    message_template = "type must be 'legal' or 'actual', but is '{}'"
    default_args = ('x',)


class E09(VOWarning, ValueError):
    """
    The ``MIN``, ``MAX`` and ``OPTION`` elements must always have a
    ``value`` attribute.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`__
    """

    message_template = "'{}' must have a value attribute"
    default_args = ('x',)


class E10(VOWarning, ValueError):
    """
    From VOTable 1.1 and later, ``FIELD`` and ``PARAM`` elements must have
    a ``datatype`` field.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`__
    """

    message_template = "'datatype' attribute required on all '{}' elements"
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
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:form>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:form>`__
    """

    message_template = "precision '{}' is invalid"
    default_args = ('x',)


class E12(VOWarning, ValueError):
    """
    The width attribute is meant to indicate to the application the
    number of characters to be used for input or output of the
    quantity.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:form>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:form>`__
    """

    message_template = "width must be a positive integer, got '{}'"
    default_args = ('x',)


class E13(VOWarning, ValueError):
    r"""
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
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:dim>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#sec:dim>`__
    """

    message_template = "Invalid arraysize attribute '{}'"
    default_args = ('x',)


class E14(VOWarning, ValueError):
    """
    All ``PARAM`` elements must have a ``value`` attribute.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`__
    """

    message_template = "value attribute is required for all PARAM elements"


class E15(VOWarning, ValueError):
    """
    All ``COOSYS`` elements must have an ``ID`` attribute.

    Note that the VOTable 1.1 specification says this attribute is
    optional, but its corresponding schema indicates it is required.

    In VOTable 1.2, the ``COOSYS`` element is deprecated.
    """

    message_template = "ID attribute is required for all COOSYS elements"


class E16(VOTableSpecWarning):
    """
    The ``system`` attribute on the ``COOSYS`` element must be one of the
    following::

      'eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5', 'galactic',
      'supergalactic', 'xy', 'barycentric', 'geo_app'

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:COOSYS>`__
    """

    message_template = "Invalid system attribute '{}'"
    default_args = ('x',)


class E17(VOWarning, ValueError):
    """
    ``extnum`` attribute must be a positive integer.

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__
    """

    message_template = "extnum must be a positive integer"


class E18(VOWarning, ValueError):
    """
    The ``type`` attribute of the ``RESOURCE`` element must be one of
    "results" or "meta".

    **References**: `1.1
    <http://www.ivoa.net/documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`__,
    `1.2
    <http://www.ivoa.net/documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`__
    """

    message_template = "type must be 'results' or 'meta', not '{}'"
    default_args = ('x',)


class E19(VOWarning, ValueError):
    """
    Raised either when the file doesn't appear to be XML, or the root
    element is not VOTABLE.
    """

    message_template = "File does not appear to be a VOTABLE"


class E20(VOTableSpecError):
    """
    The table had only *x* fields defined, but the data itself has more
    columns than that.
    """

    message_template = "Data has more columns than are defined in the header ({})"
    default_args = ('x',)


class E21(VOWarning, ValueError):
    """
    The table had *x* fields defined, but the data itself has only *y*
    columns.
    """

    message_template = "Data has fewer columns ({}) than are defined in the header ({})"
    default_args = ('x', 'y')


class E22(VOWarning, ValueError):
    """
    All ``TIMESYS`` elements must have an ``ID`` attribute.
    """

    message_template = "ID attribute is required for all TIMESYS elements"


class E23(VOTableSpecWarning):
    """
    The ``timeorigin`` attribute on the ``TIMESYS`` element must be
    either a floating point literal specifiying a valid Julian Date,
    or, for convenience, the string "MJD-origin" (standing for 2400000.5)
    or the string "JD-origin" (standing for 0).

    **References**: `1.4
    <http://www.ivoa.net/documents/VOTable/20191021/REC-VOTable-1.4-20191021.html#ToC21>`__
    """

    message_template = "Invalid timeorigin attribute '{}'"
    default_args = ('x',)


class E25(VOTableSpecWarning):
    """
    A VOTable cannot have a DATA section without any defined FIELD; DATA will be ignored.
    """

    message_template = "No FIELDs are defined; DATA section will be ignored."


def _get_warning_and_exception_classes(prefix):
    classes = []
    for key, val in globals().items():
        if re.match(prefix + "[0-9]{2}", key):
            classes.append((key, val))
    classes.sort()
    return classes


def _build_doc_string():
    def generate_set(prefix):
        classes = _get_warning_and_exception_classes(prefix)

        out = io.StringIO()

        for name, cls in classes:
            out.write(f".. _{name}:\n\n")
            msg = "{}: {}".format(cls.__name__, cls.get_short_name())
            if not isinstance(msg, str):
                msg = msg.decode('utf-8')
            out.write(msg)
            out.write('\n')
            out.write('~' * len(msg))
            out.write('\n\n')
            doc = cls.__doc__
            if not isinstance(doc, str):
                doc = doc.decode('utf-8')
            out.write(dedent(doc))
            out.write('\n\n')

        return out.getvalue()

    warnings = generate_set('W')
    exceptions = generate_set('E')

    return {'warnings': warnings,
            'exceptions': exceptions}


if __doc__ is not None:
    __doc__ = __doc__.format(**_build_doc_string())

__all__.extend([x[0] for x in _get_warning_and_exception_classes('W')])
__all__.extend([x[0] for x in _get_warning_and_exception_classes('E')])
