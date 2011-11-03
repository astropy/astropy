
.. _exceptions:

Exceptions
==========

.. note::

    This is a list of many of the fatal exceptions emitted by vo.table
    when the file does not conform to spec.  Other exceptions may be
    raised due to unforeseen cases or bugs in vo.table itself.


.. _E01:

E01: Invalid size specifier 'x' for a char/unicode field (in field 'y')
-----------------------------------------------------------------------

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

.. _E02:

E02: Incorrect number of elements in array.  Expected multiple of x, got y
--------------------------------------------------------------------------

The number of array elements in the data does not match that specified
in the FIELD specifier.

.. _E03:

E03: 'x' does not parse as a complex number
-------------------------------------------

Complex numbers should be two values separated by whitespace.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _E04:

E04: Invalid bit value 'x'
--------------------------

A ``bit`` array should be a string of '0's and '1's.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _E05:

E05: Invalid boolean value 'x'
------------------------------

A ``boolean`` value should be one of the following strings (case
insensitive) in the ``TABLEDATA`` format::

    'TRUE', 'FALSE', '1', '0', 'T', 'F', '\0', ' ', '?'

and in ``BINARY`` format::

    'T', 'F', '1', '0', '\0', ' ', '?'

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:datatypes>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:datatypes>`_

.. _E06:

E06: Unknown datatype 'x' on field 'y'
--------------------------------------

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

.. _E08:

E08: type must be 'legal' or 'actual', but is 'x'
-------------------------------------------------

The ``type`` attribute on the ``VALUES`` element must be either
``legal`` or ``actual``.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`_

.. _E09:

E09: 'x' must have a value attribute
------------------------------------

The ``MIN``, ``MAX`` and ``OPTION`` elements must always have a
``value`` attribute.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:values>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:values>`_

.. _E10:

E10: 'datatype' attribute required on all 'FIELD' elements
----------------------------------------------------------

From VOTable 1.1 and later, ``FIELD`` and ``PARAM`` elements must have
a ``datatype`` field.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`_

.. _E11:

E11: precision 'x' is invalid
-----------------------------

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

.. _E12:

E12: width must be a positive integer, got 'x'
----------------------------------------------

The width attribute is meant to indicate to the application the
number of characters to be used for input or output of the
quantity.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:form>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:form>`_

.. _E13:

E13: Invalid arraysize attribute 'x'
------------------------------------

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

.. _E14:

E14: value attribute is required for all PARAM elements
-------------------------------------------------------

All ``PARAM`` elements must have a ``value`` attribute.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:FIELD>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#elem:FIELD>`_

.. _E15:

E15: ID attribute is required for all COOSYS elements
-----------------------------------------------------

All ``COOSYS`` elements must have an ``ID`` attribute.

Note that the VOTable 1.1 specification says this attribute is
optional, but its corresponding schema indicates it is required.

In VOTable 1.2, the ``COOSYS`` element is deprecated.

.. _E16:

E16: Invalid system attribute 'x'
---------------------------------

The ``system`` attribute on the ``COOSYS`` element must be one of the
following::

  'eq_FK4', 'eq_FK5', 'ICRS', 'ecl_FK4', 'ecl_FK5', 'galactic',
  'supergalactic', 'xy', 'barycentric', 'geo_app'

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#elem:COOSYS>`_

.. _E17:

E17: extnum must be a positive integer
--------------------------------------

``extnum`` attribute must be a positive integer.

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_

.. _E18:

E18: type must be 'results' or 'meta', not 'x'
----------------------------------------------

The ``type`` attribute of the ``RESOURCE`` element must be one of
"results" or "meta".

**References**: `1.1
<http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC54>`_,
`1.2
<http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#ToC58>`_

.. _E19:

E19: File does not appear to be a VOTABLE
-----------------------------------------

Raised either when the file doesn't appear to be XML, or the root
element is not VOTABLE.

.. _E20:

E20: Data has more columns than are defined in the header (x)
-------------------------------------------------------------

The table had only *x* fields defined, but the data itself has more
columns than that.

.. _E21:

E21: Data has fewer columns (x) than are defined in the header (y)
------------------------------------------------------------------

The table had *x* fields defined, but the data itself has only *y*
columns.

