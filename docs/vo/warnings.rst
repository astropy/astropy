
.. _warnings:

Warnings
========

.. note::
    Most of the following warnings indicate violations of the VOTable
    specification.  They should be reported to the authors of the
    tools that produced the VOTable file.

    To control the warnings emitted, use the standard Python
    :mod:`warnings` module.  Most of these are of the type
    `VOTableSpecWarning`.


.. _W01:

W01: Array uses commas rather than whitespace
---------------------------------------------

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

.. _W02:

W02: x attribute 'y' is invalid.  Must be a standard XML id
-----------------------------------------------------------

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

.. _W03:

W03: Implictly generating an ID from a name 'x' -> 'y'
------------------------------------------------------

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

.. _W04:

W04: content-type 'x' must be a valid MIME content type
-------------------------------------------------------

The ``content-type`` attribute must use MIME content-type syntax as
defined in `RFC 2046 <http://tools.ietf.org/html/rfc2046>`_.

The current check for validity is somewhat over-permissive.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:link>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:link>`_

.. _W05:

W05: 'x' is not a valid URI
---------------------------

The attribute must be a valid URI as defined in `RFC 2396
<http://www.ietf.org/rfc/rfc2396.txt>`_.

.. _W06:

W06: Invalid UCD 'x': explanation
---------------------------------

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

.. _W07:

W07: Invalid astroYear in x: 'y'
--------------------------------

As astro year field is a Besselian or Julian year matching the
regular expression::

    ^[JB]?[0-9]+([.][0-9]*)?$

Defined in this XML Schema snippet::

    <xs:simpleType  name="astroYear">
      <xs:restriction base="xs:token">
        <xs:pattern  value="[JB]?[0-9]+([.][0-9]*)?"/>
      </xs:restriction>
    </xs:simpleType>

.. _W08:

W08: 'x' must be a str or unicode object
----------------------------------------

To avoid local-dependent number parsing differences, ``vo.table``
may require a string or unicode string where a numeric type may
make more sense.

.. _W09:

W09: ID attribute not capitalized
---------------------------------

The VOTable specification uses the attribute name ``ID`` (with
uppercase letters) to specify unique identifiers.  Some
VOTable-producing tools use the more standard lowercase ``id``
instead.  ``vo.table`` accepts ``id`` and emits this warning when
not in ``pedantic`` mode.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_

.. _W10:

W10: Unknown tag 'x'.  Ignoring
-------------------------------

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

.. _W11:

W11: The gref attribute on LINK is deprecated in VOTable 1.1
------------------------------------------------------------

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

.. _W12:

W12: 'x' element must have at least one of 'ID' or 'name' attributes
--------------------------------------------------------------------

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

.. _W13:

W13: 'x' is not a valid VOTable datatype, should be 'y'
-------------------------------------------------------

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

.. _W15:

W15: x element missing required 'name' attribute
------------------------------------------------

The ``name`` attribute is required on every ``FIELD`` element.
However, many VOTable files in the wild omit it and provide only
an ``ID`` instead.  In this case, when *pedantic mode* is off,
``vo.table`` will copy the ``name`` attribute to a new ``ID``
attribute.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_

.. _W17:

W17: x element contains more than one DESCRIPTION element
---------------------------------------------------------

A ``DESCRIPTION`` element can only appear once within its parent
element.

According to the schema, it may only occur once (`1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_)

However, it is a `proposed extension
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:addesc>`_
to VOTable 1.2.

.. _W18:

W18: TABLE specified nrows=x, but table contains y rows
-------------------------------------------------------

The number of rows explicitly specified in the ``nrows`` attribute
does not match the actual number of rows (``TR`` elements) present
in the ``TABLE``.  This may indicate truncation of the file, or an
internal error in the tool that produced it.  If *pedantic mode*
is off, parsing will proceed, with the loss of some performance.

**References:** `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC10>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC10>`_

.. _W19:

W19: The fields defined in the VOTable do not match those in the embedded FITS file
-----------------------------------------------------------------------------------

The column fields as defined using ``FIELD`` elements do not match
those in the headers of the embedded FITS file.  If *pedantic
mode* is off, the embedded FITS file will take precedence.

.. _W20:

W20: No version number specified in file.  Assuming 1.1
-------------------------------------------------------

If no version number is explicitly given in the VOTable file, the
parser assumes it is written to the VOTable 1.1 specification.

.. _W21:

W21: vo.table is designed for VOTable version 1.1 and 1.2, but this file is x
-----------------------------------------------------------------------------

Unknown issues may arise using ``vo.table`` with VOTable files
from a version other than 1.1 or 1.2.

.. _W22:

W22: The DEFINITIONS element is deprecated in VOTable 1.1.  Ignoring
--------------------------------------------------------------------

Version 1.0 of the VOTable specification used the ``DEFINITIONS``
element to define coordinate systems.  Version 1.1 now uses
``COOSYS`` elements throughout the document.

**References:** `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:definitions>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:definitions>`_

.. _W23:

W23: Unable to update service information for 'x'
-------------------------------------------------

Raised when the VO service database can not be updated (possibly
due to a network outage).  This is only a warning, since an older
and possible out-of-date VO service database was available
locally.

.. _W24:

W24: The VO catalog database is for a later version of vo.table
---------------------------------------------------------------

The VO catalog database retrieved from the www is designed for a
newer version of vo.table.  This may cause problems or limited
features performing service queries.  Consider upgrading vo.table
to the latest version.

.. _W25:

W25: 'service' failed with: ...
-------------------------------

A VO service query failed due to a network error or malformed
arguments.  Another alternative service may be attempted.  If all
services fail, an exception will be raised.

.. _W26:

W26: 'child' inside 'parent' added in VOTable X.X
-------------------------------------------------

The given element was not supported inside of the given element
until the specified VOTable version, however the version declared
in the file is for an earlier version.  These attributes may not
be written out to the file.

.. _W27:

W27: COOSYS deprecated in VOTable 1.2
-------------------------------------

The ``COOSYS`` element was deprecated in VOTABLE version 1.2 in
favor of a reference to the Space-Time Coordinate (STC) data
model (see `utype
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:utype>`_
and the IVOA note `referencing STC in VOTable
<http://ivoa.net/Documents/latest/VOTableSTC.html>`_.

.. _W28:

W28: 'attribute' on 'element' added in VOTable X.X
--------------------------------------------------

The given attribute was not supported on the given element until the
specified VOTable version, however the version declared in the file is
for an earlier version.  These attributes may not be written out to
the file.

.. _W29:

W29: Version specified in non-standard form 'v1.0'
--------------------------------------------------

Some VOTable files specify their version number in the form "v1.0",
when the only supported forms in the spec are "1.0".

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_

.. _W30:

W30: Invalid literal for float 'x'.  Treating as empty.
-------------------------------------------------------

Some VOTable files write missing floating-point values in non-standard
ways, such as "null" and "-".  In non-pedantic mode, any non-standard
floating-point literals are treated as missing values.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _W31:

W31: NaN given in an integral field without a specified null value
------------------------------------------------------------------

Since NaN's can not be represented in integer fields directly, a null
value must be specified in the FIELD descriptor to support reading
NaN's from the tabledata.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _W32:

W32: Duplicate ID 'x' renamed to 'x_2' to ensure uniqueness
-----------------------------------------------------------

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

.. _W33:

W33: Column name 'x' renamed to 'x_2' to ensure uniqueness
----------------------------------------------------------

Each field in a table must have a unique name.  If two or more
fields have the same name, some will be renamed to ensure that all
names are unique.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:name>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:name>`_

.. _W34:

W34: 'x' is an invalid token for attribute 'y'
----------------------------------------------

The attribute requires the value to be a valid XML token, as
defined by `XML 1.0
<http://www.w3.org/TR/2000/WD-xml-2e-20000814#NT-Nmtoken>`_.

.. _W35:

W35: 'x' attribute required for INFO elements
---------------------------------------------

The ``name`` and ``value`` attributes are required on all ``INFO``
elements.

**References:** `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC32>`

.. _W36:

W36: null value 'x' does not match field datatype, setting to 0
---------------------------------------------------------------

If the field specifies a ``null`` value, that value must conform
to the given ``datatype``.

**References:** `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`

.. _W37:

W37: Unsupported data format 'x'
--------------------------------

The 3 datatypes defined in the VOTable specification and supported by
vo.table are ``TABLEDATA``, ``BINARY`` and ``FITS``.

**References:** `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:data>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:data>`

.. _W38:

W38: Inline binary data must be base64 encoded, got 'x'
-------------------------------------------------------

The only encoding for local binary data supported by the VOTable
specification is base64.

.. _W39:

W39: Bit values can not be masked
---------------------------------

Bit values do not support masking.  This warning is raised upon
setting masked data in a bit column.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _W40:

W40: 'cprojection' datatype repaired
------------------------------------

This is a terrible hack to support Simple Image Access Protocol
results from `<archive.noao.edu>`_.  It creates a field for the
coordinate projection type of type "double", which
actually contains character data.  We have to hack the field
to store character data, or we can't read it in.  A warning
will be raised when this happens.

.. _W41:

W41: An XML namespace is specified, but is incorrect.  Expected 'x', got 'y'
----------------------------------------------------------------------------

An XML namespace was specified on the ``VOTABLE`` element, but the
namespace does not match what is expected for a ``VOTABLE`` file.

The ``VOTABLE`` namespace is::

  http://www.ivoa.net/xml/VOTable/vX.X

where "X.X" is the version number.

Some files in the wild set the namespace to the location of the
VOTable schema, which is not correct and will not pass some
validating parsers.

.. _W42:

W42: No XML namespace specified
-------------------------------

The root element should specify a namespace.

The ``VOTABLE`` namespace is::

    http://www.ivoa.net/xml/VOTable/vX.X

where "X.X" is the version number.

.. _W43:

W43: element ref='x' which has not already been defined
-------------------------------------------------------

Referenced elements should be defined before referees.  From the
VOTable 1.2 spec:

   In VOTable1.2, it is further recommended to place the ID
   attribute prior to referencing it whenever possible.

.. _W44:

W44: VALUES element with ref attribute has content ('element')
--------------------------------------------------------------

``VALUES`` elements that reference another element should not have
their own content.

From the VOTable 1.2 spec:

    The ``ref`` attribute of a ``VALUES`` element can be used to
    avoid a repetition of the domain definition, by referring to a
    previously defined ``VALUES`` element having the referenced
    ``ID`` attribute. When specified, the ``ref`` attribute
    defines completely the domain without any other element or
    attribute, as e.g. ``<VALUES ref="RAdomain"/>``

.. _W45:

W45: content-role attribute 'x' invalid
---------------------------------------

The ``content-role`` attribute on the ``LINK`` element must be one of
the following::

    query, hints, doc, location

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_

