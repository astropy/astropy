# Licensed under a 3-clause BSD style license - see PYFITS.rst

import copy
import numbers
import operator
import re
import sys
import warnings
import weakref
from collections import OrderedDict
from contextlib import suppress
from functools import reduce

import numpy as np
from numpy import char as chararray

from astropy.utils import indent, isiterable, lazyproperty
from astropy.utils.exceptions import AstropyUserWarning

from .card import CARD_LENGTH, Card
from .util import NotifierMixin, _convert_array, _is_int, cmp, encode_ascii, pairwise
from .verify import VerifyError, VerifyWarning

__all__ = ["Column", "ColDefs", "Delayed"]


# mapping from TFORM data type to numpy data type (code)
# L: Logical (Boolean)
# B: Unsigned Byte
# I: 16-bit Integer
# J: 32-bit Integer
# K: 64-bit Integer
# E: Single-precision Floating Point
# D: Double-precision Floating Point
# C: Single-precision Complex
# M: Double-precision Complex
# A: Character
FITS2NUMPY = {
    "L": "i1",
    "B": "u1",
    "I": "i2",
    "J": "i4",
    "K": "i8",
    "E": "f4",
    "D": "f8",
    "C": "c8",
    "M": "c16",
    "A": "S",
}

# the inverse dictionary of the above
NUMPY2FITS = {val: key for key, val in FITS2NUMPY.items()}
# Normally booleans are represented as ints in Astropy, but if passed in a numpy
# boolean array, that should be supported
NUMPY2FITS["b1"] = "L"
# Add unsigned types, which will be stored as signed ints with a TZERO card.
NUMPY2FITS["u2"] = "I"
NUMPY2FITS["u4"] = "J"
NUMPY2FITS["u8"] = "K"
# Add half precision floating point numbers which will be up-converted to
# single precision.
NUMPY2FITS["f2"] = "E"

# This is the order in which values are converted to FITS types
# Note that only double precision floating point/complex are supported
FORMATORDER = ["L", "B", "I", "J", "K", "D", "M", "A"]

# Convert single precision floating point/complex to double precision.
FITSUPCONVERTERS = {"E": "D", "C": "M"}

# mapping from ASCII table TFORM data type to numpy data type
# A: Character
# I: Integer (32-bit)
# J: Integer (64-bit; non-standard)
# F: Float (64-bit; fixed decimal notation)
# E: Float (64-bit; exponential notation)
# D: Float (64-bit; exponential notation, always 64-bit by convention)
ASCII2NUMPY = {"A": "S", "I": "i4", "J": "i8", "F": "f8", "E": "f8", "D": "f8"}

# Maps FITS ASCII column format codes to the appropriate Python string
# formatting codes for that type.
ASCII2STR = {"A": "", "I": "d", "J": "d", "F": "f", "E": "E", "D": "E"}

# For each ASCII table format code, provides a default width (and decimal
# precision) for when one isn't given explicitly in the column format
ASCII_DEFAULT_WIDTHS = {
    "A": (1, 0),
    "I": (10, 0),
    "J": (15, 0),
    "E": (15, 7),
    "F": (16, 7),
    "D": (25, 17),
}

# TDISPn for both ASCII and Binary tables
TDISP_RE_DICT = {}
TDISP_RE_DICT["F"] = re.compile(
    r"(?:(?P<formatc>[F])(?:(?P<width>[0-9]+)\.{1}(?P<precision>[0-9]+))+)|"
)
TDISP_RE_DICT["A"] = TDISP_RE_DICT["L"] = re.compile(
    r"(?:(?P<formatc>[AL])(?P<width>[0-9]+)+)|"
)
TDISP_RE_DICT["I"] = TDISP_RE_DICT["B"] = TDISP_RE_DICT["O"] = TDISP_RE_DICT[
    "Z"
] = re.compile(
    r"(?:(?P<formatc>[IBOZ])(?:(?P<width>[0-9]+)"
    r"(?:\.{0,1}(?P<precision>[0-9]+))?))|"
)
TDISP_RE_DICT["E"] = TDISP_RE_DICT["G"] = TDISP_RE_DICT["D"] = re.compile(
    r"(?:(?P<formatc>[EGD])(?:(?P<width>[0-9]+)\."
    r"(?P<precision>[0-9]+))+)"
    r"(?:E{0,1}(?P<exponential>[0-9]+)?)|"
)
TDISP_RE_DICT["EN"] = TDISP_RE_DICT["ES"] = re.compile(
    r"(?:(?P<formatc>E[NS])(?:(?P<width>[0-9]+)\.{1}(?P<precision>[0-9]+))+)"
)

# mapping from TDISP format to python format
# A: Character
# L: Logical (Boolean)
# I: 16-bit Integer
#    Can't predefine zero padding and space padding before hand without
#    knowing the value being formatted, so grabbing precision and using that
#    to zero pad, ignoring width. Same with B, O, and Z
# B: Binary Integer
# O: Octal Integer
# Z: Hexadecimal Integer
# F: Float (64-bit; fixed decimal notation)
# EN: Float (engineering fortran format, exponential multiple of thee
# ES: Float (scientific, same as EN but non-zero leading digit
# E: Float, exponential notation
#    Can't get exponential restriction to work without knowing value
#    before hand, so just using width and precision, same with D, G, EN, and
#    ES formats
# D: Double-precision Floating Point with exponential
#    (E but for double precision)
# G: Double-precision Floating Point, may or may not show exponent
TDISP_FMT_DICT = {
    "I": "{{:{width}d}}",
    "B": "{{:{width}b}}",
    "O": "{{:{width}o}}",
    "Z": "{{:{width}x}}",
    "F": "{{:{width}.{precision}f}}",
    "G": "{{:{width}.{precision}g}}",
}
TDISP_FMT_DICT["A"] = TDISP_FMT_DICT["L"] = "{{:>{width}}}"
TDISP_FMT_DICT["E"] = TDISP_FMT_DICT["D"] = TDISP_FMT_DICT["EN"] = TDISP_FMT_DICT[
    "ES"
] = "{{:{width}.{precision}e}}"

# tuple of column/field definition common names and keyword names, make
# sure to preserve the one-to-one correspondence when updating the list(s).
# Use lists, instead of dictionaries so the names can be displayed in a
# preferred order.
KEYWORD_NAMES = (
    "TTYPE",
    "TFORM",
    "TUNIT",
    "TNULL",
    "TSCAL",
    "TZERO",
    "TDISP",
    "TBCOL",
    "TDIM",
    "TCTYP",
    "TCUNI",
    "TCRPX",
    "TCRVL",
    "TCDLT",
    "TRPOS",
)
KEYWORD_ATTRIBUTES = (
    "name",
    "format",
    "unit",
    "null",
    "bscale",
    "bzero",
    "disp",
    "start",
    "dim",
    "coord_type",
    "coord_unit",
    "coord_ref_point",
    "coord_ref_value",
    "coord_inc",
    "time_ref_pos",
)
"""This is a list of the attributes that can be set on `Column` objects."""


KEYWORD_TO_ATTRIBUTE = OrderedDict(zip(KEYWORD_NAMES, KEYWORD_ATTRIBUTES))

ATTRIBUTE_TO_KEYWORD = OrderedDict(zip(KEYWORD_ATTRIBUTES, KEYWORD_NAMES))


# TODO: Define a list of default comments to associate with each table keyword

# TFORMn regular expression
TFORMAT_RE = re.compile(
    r"(?P<repeat>^[0-9]*)(?P<format>[LXBIJKAEDCMPQ])(?P<option>[!-~]*)", re.I
)

# TFORMn for ASCII tables; two different versions depending on whether
# the format is floating-point or not; allows empty values for width
# in which case defaults are used
TFORMAT_ASCII_RE = re.compile(
    r"(?:(?P<format>[AIJ])(?P<width>[0-9]+)?)|"
    r"(?:(?P<formatf>[FED])"
    r"(?:(?P<widthf>[0-9]+)(?:\."
    r"(?P<precision>[0-9]+))?)?)"
)

TTYPE_RE = re.compile(r"[0-9a-zA-Z_]+")
"""
Regular expression for valid table column names.  See FITS Standard v3.0 section 7.2.2.
"""

# table definition keyword regular expression
TDEF_RE = re.compile(r"(?P<label>^T[A-Z]*)(?P<num>[1-9][0-9 ]*$)")

# table dimension keyword regular expression (fairly flexible with whitespace)
TDIM_RE = re.compile(r"\(\s*(?P<dims>(?:\d+\s*)(?:,\s*\d+\s*)*\s*)\)\s*")

# value for ASCII table cell with value = TNULL
# this can be reset by user.
ASCIITNULL = 0

# The default placeholder to use for NULL values in ASCII tables when
# converting from binary to ASCII tables
DEFAULT_ASCII_TNULL = "---"


class Delayed:
    """Delayed file-reading data."""

    def __init__(self, hdu=None, field=None):
        self.hdu = weakref.proxy(hdu)
        self.field = field

    def __getitem__(self, key):
        # This forces the data for the HDU to be read, which will replace
        # the corresponding Delayed objects in the Tables Columns to be
        # transformed into ndarrays.  It will also return the value of the
        # requested data element.
        return self.hdu.data[key][self.field]


class _BaseColumnFormat(str):
    """
    Base class for binary table column formats (just called _ColumnFormat)
    and ASCII table column formats (_AsciiColumnFormat).
    """

    def __eq__(self, other):
        if not other:
            return False

        if isinstance(other, str):
            if not isinstance(other, self.__class__):
                try:
                    other = self.__class__(other)
                except ValueError:
                    return False
        else:
            return False

        return self.canonical == other.canonical

    def __hash__(self):
        return hash(self.canonical)

    @lazyproperty
    def dtype(self):
        """
        The Numpy dtype object created from the format's associated recformat.
        """
        return np.dtype(self.recformat)

    @classmethod
    def from_column_format(cls, format):
        """Creates a column format object from another column format object
        regardless of their type.

        That is, this can convert a _ColumnFormat to an _AsciiColumnFormat
        or vice versa at least in cases where a direct translation is possible.
        """
        return cls.from_recformat(format.recformat)


class _ColumnFormat(_BaseColumnFormat):
    """
    Represents a FITS binary table column format.

    This is an enhancement over using a normal string for the format, since the
    repeat count, format code, and option are available as separate attributes,
    and smart comparison is used.  For example 1J == J.
    """

    def __new__(cls, format):
        self = super().__new__(cls, format)
        self.repeat, self.format, self.option = _parse_tformat(format)
        self.format = self.format.upper()
        if self.format in ("P", "Q"):
            # TODO: There should be a generic factory that returns either
            # _FormatP or _FormatQ as appropriate for a given TFORMn
            if self.format == "P":
                recformat = _FormatP.from_tform(format)
            else:
                recformat = _FormatQ.from_tform(format)
            # Format of variable length arrays
            self.p_format = recformat.format
        else:
            self.p_format = None
        return self

    @classmethod
    def from_recformat(cls, recformat):
        """Creates a column format from a Numpy record dtype format."""
        return cls(_convert_format(recformat, reverse=True))

    @lazyproperty
    def recformat(self):
        """Returns the equivalent Numpy record format string."""
        return _convert_format(self)

    @lazyproperty
    def canonical(self):
        """
        Returns a 'canonical' string representation of this format.

        This is in the proper form of rTa where T is the single character data
        type code, a is the optional part, and r is the repeat.  If repeat == 1
        (the default) it is left out of this representation.
        """
        if self.repeat == 1:
            repeat = ""
        else:
            repeat = str(self.repeat)

        return f"{repeat}{self.format}{self.option}"


class _AsciiColumnFormat(_BaseColumnFormat):
    """Similar to _ColumnFormat but specifically for columns in ASCII tables.

    The formats of ASCII table columns and binary table columns are inherently
    incompatible in FITS.  They don't support the same ranges and types of
    values, and even reuse format codes in subtly different ways.  For example
    the format code 'Iw' in ASCII columns refers to any integer whose string
    representation is at most w characters wide, so 'I' can represent
    effectively any integer that will fit in a FITS columns.  Whereas for
    binary tables 'I' very explicitly refers to a 16-bit signed integer.

    Conversions between the two column formats can be performed using the
    ``to/from_binary`` methods on this class, or the ``to/from_ascii``
    methods on the `_ColumnFormat` class.  But again, not all conversions are
    possible and may result in a `ValueError`.
    """

    def __new__(cls, format, strict=False):
        self = super().__new__(cls, format)
        self.format, self.width, self.precision = _parse_ascii_tformat(format, strict)

        # If no width has been specified, set the dtype here to default as well
        if format == self.format:
            self.recformat = ASCII2NUMPY[format]

        # This is to support handling logical (boolean) data from binary tables
        # in an ASCII table
        self._pseudo_logical = False
        return self

    @classmethod
    def from_column_format(cls, format):
        inst = cls.from_recformat(format.recformat)
        # Hack
        if format.format == "L":
            inst._pseudo_logical = True
        return inst

    @classmethod
    def from_recformat(cls, recformat):
        """Creates a column format from a Numpy record dtype format."""
        return cls(_convert_ascii_format(recformat, reverse=True))

    @lazyproperty
    def recformat(self):
        """Returns the equivalent Numpy record format string."""
        return _convert_ascii_format(self)

    @lazyproperty
    def canonical(self):
        """
        Returns a 'canonical' string representation of this format.

        This is in the proper form of Tw.d where T is the single character data
        type code, w is the width in characters for this field, and d is the
        number of digits after the decimal place (for format codes 'E', 'F',
        and 'D' only).
        """
        if self.format in ("E", "F", "D"):
            return f"{self.format}{self.width}.{self.precision}"

        return f"{self.format}{self.width}"


class _FormatX(str):
    """For X format in binary tables."""

    def __new__(cls, repeat=1):
        nbytes = ((repeat - 1) // 8) + 1
        # use an array, even if it is only ONE u1 (i.e. use tuple always)
        obj = super().__new__(cls, repr((nbytes,)) + "u1")
        obj.repeat = repeat
        return obj

    def __getnewargs__(self):
        return (self.repeat,)

    @property
    def tform(self):
        return f"{self.repeat}X"


# TODO: Table column formats need to be verified upon first reading the file;
# as it is, an invalid P format will raise a VerifyError from some deep,
# unexpected place
class _FormatP(str):
    """For P format in variable length table."""

    # As far as I can tell from my reading of the FITS standard, a type code is
    # *required* for P and Q formats; there is no default
    _format_re_template = (
        r"(?P<repeat>\d+)?{}(?P<dtype>[LXBIJKAEDCM])(?:\((?P<max>\d*)\))?"
    )
    _format_code = "P"
    _format_re = re.compile(_format_re_template.format(_format_code))
    _descriptor_format = "2i4"

    def __new__(cls, dtype, repeat=None, max=None):
        obj = super().__new__(cls, cls._descriptor_format)
        obj.format = NUMPY2FITS[dtype]
        obj.dtype = dtype
        obj.repeat = repeat
        obj.max = max
        return obj

    def __getnewargs__(self):
        return (self.dtype, self.repeat, self.max)

    @classmethod
    def from_tform(cls, format):
        m = cls._format_re.match(format)
        if not m or m.group("dtype") not in FITS2NUMPY:
            raise VerifyError(f"Invalid column format: {format}")
        repeat = m.group("repeat")
        array_dtype = m.group("dtype")
        max = m.group("max")
        if not max:
            max = None
        return cls(FITS2NUMPY[array_dtype], repeat=repeat, max=max)

    @property
    def tform(self):
        repeat = "" if self.repeat is None else self.repeat
        max = "" if self.max is None else self.max
        return f"{repeat}{self._format_code}{self.format}({max})"


class _FormatQ(_FormatP):
    """Carries type description of the Q format for variable length arrays.

    The Q format is like the P format but uses 64-bit integers in the array
    descriptors, allowing for heaps stored beyond 2GB into a file.
    """

    _format_code = "Q"
    _format_re = re.compile(_FormatP._format_re_template.format(_format_code))
    _descriptor_format = "2i8"


class ColumnAttribute:
    """
    Descriptor for attributes of `Column` that are associated with keywords
    in the FITS header and describe properties of the column as specified in
    the FITS standard.

    Each `ColumnAttribute` may have a ``validator`` method defined on it.
    This validates values set on this attribute to ensure that they meet the
    FITS standard.  Invalid values will raise a warning and will not be used in
    formatting the column.  The validator should take two arguments--the
    `Column` it is being assigned to, and the new value for the attribute, and
    it must raise an `AssertionError` if the value is invalid.

    The `ColumnAttribute` itself is a decorator that can be used to define the
    ``validator`` for each column attribute.  For example::

        @ColumnAttribute('TTYPE')
        def name(col, name):
            if not isinstance(name, str):
                raise AssertionError

    The actual object returned by this decorator is the `ColumnAttribute`
    instance though, not the ``name`` function.  As such ``name`` is not a
    method of the class it is defined in.

    The setter for `ColumnAttribute` also updates the header of any table
    HDU this column is attached to in order to reflect the change.  The
    ``validator`` should ensure that the value is valid for inclusion in a FITS
    header.
    """

    def __init__(self, keyword):
        self._keyword = keyword
        self._validator = None

        # The name of the attribute associated with this keyword is currently
        # determined from the KEYWORD_NAMES/ATTRIBUTES lists.  This could be
        # make more flexible in the future, for example, to support custom
        # column attributes.
        self._attr = "_" + KEYWORD_TO_ATTRIBUTE[self._keyword]

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        else:
            return getattr(obj, self._attr)

    def __set__(self, obj, value):
        if self._validator is not None:
            self._validator(obj, value)

        old_value = getattr(obj, self._attr, None)
        setattr(obj, self._attr, value)
        obj._notify("column_attribute_changed", obj, self._attr[1:], old_value, value)

    def __call__(self, func):
        """
        Set the validator for this column attribute.

        Returns ``self`` so that this can be used as a decorator, as described
        in the docs for this class.
        """
        self._validator = func

        return self

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._keyword}')"


class Column(NotifierMixin):
    """
    Class which contains the definition of one column, e.g.  ``ttype``,
    ``tform``, etc. and the array containing values for the column.
    """

    def __init__(
        self,
        name=None,
        format=None,
        unit=None,
        null=None,
        bscale=None,
        bzero=None,
        disp=None,
        start=None,
        dim=None,
        array=None,
        ascii=None,
        coord_type=None,
        coord_unit=None,
        coord_ref_point=None,
        coord_ref_value=None,
        coord_inc=None,
        time_ref_pos=None,
    ):
        """
        Construct a `Column` by specifying attributes.  All attributes
        except ``format`` can be optional; see :ref:`astropy:column_creation`
        and :ref:`astropy:creating_ascii_table` for more information regarding
        ``TFORM`` keyword.

        Parameters
        ----------
        name : str, optional
            column name, corresponding to ``TTYPE`` keyword

        format : str
            column format, corresponding to ``TFORM`` keyword

        unit : str, optional
            column unit, corresponding to ``TUNIT`` keyword

        null : str, optional
            null value, corresponding to ``TNULL`` keyword

        bscale : int-like, optional
            bscale value, corresponding to ``TSCAL`` keyword

        bzero : int-like, optional
            bzero value, corresponding to ``TZERO`` keyword

        disp : str, optional
            display format, corresponding to ``TDISP`` keyword

        start : int, optional
            column starting position (ASCII table only), corresponding
            to ``TBCOL`` keyword

        dim : str, optional
            column dimension corresponding to ``TDIM`` keyword

        array : iterable, optional
            a `list`, `numpy.ndarray` (or other iterable that can be used to
            initialize an ndarray) providing initial data for this column.
            The array will be automatically converted, if possible, to the data
            format of the column.  In the case were non-trivial ``bscale``
            and/or ``bzero`` arguments are given, the values in the array must
            be the *physical* values--that is, the values of column as if the
            scaling has already been applied (the array stored on the column
            object will then be converted back to its storage values).

        ascii : bool, optional
            set `True` if this describes a column for an ASCII table; this
            may be required to disambiguate the column format

        coord_type : str, optional
            coordinate/axis type corresponding to ``TCTYP`` keyword

        coord_unit : str, optional
            coordinate/axis unit corresponding to ``TCUNI`` keyword

        coord_ref_point : int-like, optional
            pixel coordinate of the reference point corresponding to ``TCRPX``
            keyword

        coord_ref_value : int-like, optional
            coordinate value at reference point corresponding to ``TCRVL``
            keyword

        coord_inc : int-like, optional
            coordinate increment at reference point corresponding to ``TCDLT``
            keyword

        time_ref_pos : str, optional
            reference position for a time coordinate column corresponding to
            ``TRPOS`` keyword
        """
        if format is None:
            raise ValueError("Must specify format to construct Column.")

        # any of the input argument (except array) can be a Card or just
        # a number/string
        kwargs = {"ascii": ascii}
        for attr in KEYWORD_ATTRIBUTES:
            value = locals()[attr]  # get the argument's value

            if isinstance(value, Card):
                value = value.value

            kwargs[attr] = value

        valid_kwargs, invalid_kwargs = self._verify_keywords(**kwargs)

        if invalid_kwargs:
            msg = ["The following keyword arguments to Column were invalid:"]

            for val in invalid_kwargs.values():
                msg.append(indent(val[1]))

            raise VerifyError("\n".join(msg))

        for attr in KEYWORD_ATTRIBUTES:
            setattr(self, attr, valid_kwargs.get(attr))

        # TODO: Try to eliminate the following two special cases
        # for recformat and dim:
        # This is not actually stored as an attribute on columns for some
        # reason
        recformat = valid_kwargs["recformat"]

        # The 'dim' keyword's original value is stored in self.dim, while
        # *only* the tuple form is stored in self._dims.
        self._dims = self.dim
        self.dim = dim

        # Awful hack to use for now to keep track of whether the column holds
        # pseudo-unsigned int data
        self._pseudo_unsigned_ints = False

        # if the column data is not ndarray, make it to be one, i.e.
        # input arrays can be just list or tuple, not required to be ndarray
        # does not include Object array because there is no guarantee
        # the elements in the object array are consistent.
        if not isinstance(array, (np.ndarray, chararray.chararray, Delayed)):
            try:  # try to convert to a ndarray first
                if array is not None:
                    array = np.array(array)
            except Exception:
                try:  # then try to convert it to a strings array
                    itemsize = int(recformat[1:])
                    array = chararray.array(array, itemsize=itemsize)
                except ValueError:
                    # then try variable length array
                    # Note: This includes _FormatQ by inheritance
                    if isinstance(recformat, _FormatP):
                        array = _VLF(array, dtype=recformat.dtype)
                    else:
                        raise ValueError(
                            f"Data is inconsistent with the format `{format}`."
                        )

        array = self._convert_to_valid_data_type(array)

        # We have required (through documentation) that arrays passed in to
        # this constructor are already in their physical values, so we make
        # note of that here
        if isinstance(array, np.ndarray):
            self._physical_values = True
        else:
            self._physical_values = False

        self._parent_fits_rec = None
        self.array = array

    def __repr__(self):
        text = ""
        for attr in KEYWORD_ATTRIBUTES:
            value = getattr(self, attr)
            if value is not None:
                text += attr + " = " + repr(value) + "; "
        return text[:-2]

    def __eq__(self, other):
        """
        Two columns are equal if their name and format are the same.  Other
        attributes aren't taken into account at this time.
        """
        # According to the FITS standard column names must be case-insensitive
        a = (self.name.lower(), self.format)
        b = (other.name.lower(), other.format)
        return a == b

    def __hash__(self):
        """
        Like __eq__, the hash of a column should be based on the unique column
        name and format, and be case-insensitive with respect to the column
        name.
        """
        return hash((self.name.lower(), self.format))

    @property
    def array(self):
        """
        The Numpy `~numpy.ndarray` associated with this `Column`.

        If the column was instantiated with an array passed to the ``array``
        argument, this will return that array.  However, if the column is
        later added to a table, such as via `BinTableHDU.from_columns` as
        is typically the case, this attribute will be updated to reference
        the associated field in the table, which may no longer be the same
        array.
        """
        # Ideally the .array attribute never would have existed in the first
        # place, or would have been internal-only.  This is a legacy of the
        # older design from Astropy that needs to have continued support, for
        # now.

        # One of the main problems with this design was that it created a
        # reference cycle.  When the .array attribute was updated after
        # creating a FITS_rec from the column (as explained in the docstring) a
        # reference cycle was created.  This is because the code in BinTableHDU
        # (and a few other places) does essentially the following:
        #
        # data._coldefs = columns  # The ColDefs object holding this Column
        # for col in columns:
        #     col.array = data.field(col.name)
        #
        # This way each columns .array attribute now points to the field in the
        # table data.  It's actually a pretty confusing interface (since it
        # replaces the array originally pointed to by .array), but it's the way
        # things have been for a long, long time.
        #
        # However, this results, in *many* cases, in a reference cycle.
        # Because the array returned by data.field(col.name), while sometimes
        # an array that owns its own data, is usually like a slice of the
        # original data.  It has the original FITS_rec as the array .base.
        # This results in the following reference cycle (for the n-th column):
        #
        #    data -> data._coldefs -> data._coldefs[n] ->
        #     data._coldefs[n].array -> data._coldefs[n].array.base -> data
        #
        # Because ndarray objects do not handled by Python's garbage collector
        # the reference cycle cannot be broken.  Therefore the FITS_rec's
        # refcount never goes to zero, its __del__ is never called, and its
        # memory is never freed.  This didn't occur in *all* cases, but it did
        # occur in many cases.
        #
        # To get around this, Column.array is no longer a simple attribute
        # like it was previously.  Now each Column has a ._parent_fits_rec
        # attribute which is a weakref to a FITS_rec object.  Code that
        # previously assigned each col.array to field in a FITS_rec (as in
        # the example a few paragraphs above) is still used, however now
        # array.setter checks if a reference cycle will be created.  And if
        # so, instead of saving directly to the Column's __dict__, it creates
        # the ._prent_fits_rec weakref, and all lookups of the column's .array
        # go through that instead.
        #
        # This alone does not fully solve the problem.  Because
        # _parent_fits_rec is a weakref, if the user ever holds a reference to
        # the Column, but deletes all references to the underlying FITS_rec,
        # the .array attribute would suddenly start returning None instead of
        # the array data.  This problem is resolved on FITS_rec's end.  See the
        # note in the FITS_rec._coldefs property for the rest of the story.

        # If the Columns's array is not a reference to an existing FITS_rec,
        # then it is just stored in self.__dict__; otherwise check the
        # _parent_fits_rec reference if it 's still available.
        if "array" in self.__dict__:
            return self.__dict__["array"]
        elif self._parent_fits_rec is not None:
            parent = self._parent_fits_rec()
            if parent is not None:
                return parent[self.name]
        else:
            return None

    @array.setter
    def array(self, array):
        # The following looks over the bases of the given array to check if it
        # has a ._coldefs attribute (i.e. is a FITS_rec) and that that _coldefs
        # contains this Column itself, and would create a reference cycle if we
        # stored the array directly in self.__dict__.
        # In this case it instead sets up the _parent_fits_rec weakref to the
        # underlying FITS_rec, so that array.getter can return arrays through
        # self._parent_fits_rec().field(self.name), rather than storing a
        # hard reference to the field like it used to.
        base = array
        while True:
            if hasattr(base, "_coldefs") and isinstance(base._coldefs, ColDefs):
                for col in base._coldefs:
                    if col is self and self._parent_fits_rec is None:
                        self._parent_fits_rec = weakref.ref(base)

                        # Just in case the user already set .array to their own
                        # array.
                        if "array" in self.__dict__:
                            del self.__dict__["array"]
                        return

            if getattr(base, "base", None) is not None:
                base = base.base
            else:
                break

        self.__dict__["array"] = array

    @array.deleter
    def array(self):
        try:
            del self.__dict__["array"]
        except KeyError:
            pass

        self._parent_fits_rec = None

    @ColumnAttribute("TTYPE")
    def name(col, name):
        if name is None:
            # Allow None to indicate deleting the name, or to just indicate an
            # unspecified name (when creating a new Column).
            return

        # Check that the name meets the recommended standard--other column
        # names are *allowed*, but will be discouraged
        if isinstance(name, str) and not TTYPE_RE.match(name):
            warnings.warn(
                "It is strongly recommended that column names contain only "
                "upper and lower-case ASCII letters, digits, or underscores "
                "for maximum compatibility with other software "
                f"(got {name!r}).",
                VerifyWarning,
            )

        # This ensures that the new name can fit into a single FITS card
        # without any special extension like CONTINUE cards or the like.
        if not isinstance(name, str) or len(str(Card("TTYPE", name))) != CARD_LENGTH:
            raise AssertionError(
                "Column name must be a string able to fit in a single "
                "FITS card--typically this means a maximum of 68 "
                "characters, though it may be fewer if the string "
                "contains special characters like quotes."
            )

    @ColumnAttribute("TCTYP")
    def coord_type(col, coord_type):
        if coord_type is None:
            return

        if not isinstance(coord_type, str) or len(coord_type) > 8:
            raise AssertionError(
                "Coordinate/axis type must be a string of atmost 8 characters."
            )

    @ColumnAttribute("TCUNI")
    def coord_unit(col, coord_unit):
        if coord_unit is not None and not isinstance(coord_unit, str):
            raise AssertionError("Coordinate/axis unit must be a string.")

    @ColumnAttribute("TCRPX")
    def coord_ref_point(col, coord_ref_point):
        if coord_ref_point is not None and not isinstance(
            coord_ref_point, numbers.Real
        ):
            raise AssertionError(
                "Pixel coordinate of the reference point must be real floating type."
            )

    @ColumnAttribute("TCRVL")
    def coord_ref_value(col, coord_ref_value):
        if coord_ref_value is not None and not isinstance(
            coord_ref_value, numbers.Real
        ):
            raise AssertionError(
                "Coordinate value at reference point must be real floating type."
            )

    @ColumnAttribute("TCDLT")
    def coord_inc(col, coord_inc):
        if coord_inc is not None and not isinstance(coord_inc, numbers.Real):
            raise AssertionError("Coordinate increment must be real floating type.")

    @ColumnAttribute("TRPOS")
    def time_ref_pos(col, time_ref_pos):
        if time_ref_pos is not None and not isinstance(time_ref_pos, str):
            raise AssertionError("Time reference position must be a string.")

    format = ColumnAttribute("TFORM")
    unit = ColumnAttribute("TUNIT")
    null = ColumnAttribute("TNULL")
    bscale = ColumnAttribute("TSCAL")
    bzero = ColumnAttribute("TZERO")
    disp = ColumnAttribute("TDISP")
    start = ColumnAttribute("TBCOL")
    dim = ColumnAttribute("TDIM")

    @lazyproperty
    def ascii(self):
        """Whether this `Column` represents a column in an ASCII table."""
        return isinstance(self.format, _AsciiColumnFormat)

    @lazyproperty
    def dtype(self):
        return self.format.dtype

    def copy(self):
        """
        Return a copy of this `Column`.
        """
        tmp = Column(format="I")  # just use a throw-away format
        tmp.__dict__ = self.__dict__.copy()
        return tmp

    @staticmethod
    def _convert_format(format, cls):
        """The format argument to this class's initializer may come in many
        forms.  This uses the given column format class ``cls`` to convert
        to a format of that type.

        TODO: There should be an abc base class for column format classes
        """
        # Short circuit in case we're already a _BaseColumnFormat--there is at
        # least one case in which this can happen
        if isinstance(format, _BaseColumnFormat):
            return format, format.recformat

        if format in NUMPY2FITS:
            with suppress(VerifyError):
                # legit recarray format?
                recformat = format
                format = cls.from_recformat(format)

        try:
            # legit FITS format?
            format = cls(format)
            recformat = format.recformat
        except VerifyError:
            raise VerifyError(f"Illegal format `{format}`.")

        return format, recformat

    @classmethod
    def _verify_keywords(
        cls,
        name=None,
        format=None,
        unit=None,
        null=None,
        bscale=None,
        bzero=None,
        disp=None,
        start=None,
        dim=None,
        ascii=None,
        coord_type=None,
        coord_unit=None,
        coord_ref_point=None,
        coord_ref_value=None,
        coord_inc=None,
        time_ref_pos=None,
    ):
        """
        Given the keyword arguments used to initialize a Column, specifically
        those that typically read from a FITS header (so excluding array),
        verify that each keyword has a valid value.

        Returns a 2-tuple of dicts.  The first maps valid keywords to their
        values.  The second maps invalid keywords to a 2-tuple of their value,
        and a message explaining why they were found invalid.
        """
        valid = {}
        invalid = {}

        try:
            format, recformat = cls._determine_formats(format, start, dim, ascii)
            valid.update(format=format, recformat=recformat)
        except (ValueError, VerifyError) as err:
            msg = (
                f"Column format option (TFORMn) failed verification: {err!s} "
                "The invalid value will be ignored for the purpose of "
                "formatting the data in this column."
            )
            invalid["format"] = (format, msg)
        except AttributeError as err:
            msg = (
                "Column format option (TFORMn) must be a string with a valid "
                f"FITS table format (got {format!s}: {err!s}). "
                "The invalid value will be ignored for the purpose of "
                "formatting the data in this column."
            )
            invalid["format"] = (format, msg)

        # Currently we don't have any validation for name, unit, bscale, or
        # bzero so include those by default
        # TODO: Add validation for these keywords, obviously
        for k, v in [
            ("name", name),
            ("unit", unit),
            ("bscale", bscale),
            ("bzero", bzero),
        ]:
            if v is not None and v != "":
                valid[k] = v

        # Validate null option
        # Note: Enough code exists that thinks empty strings are sensible
        # inputs for these options that we need to treat '' as None
        if null is not None and null != "":
            msg = None
            if isinstance(format, _AsciiColumnFormat):
                null = str(null)
                if len(null) > format.width:
                    msg = (
                        "ASCII table null option (TNULLn) is longer than "
                        "the column's character width and will be truncated "
                        f"(got {null!r})."
                    )
            else:
                tnull_formats = ("B", "I", "J", "K")

                if not _is_int(null):
                    # Make this an exception instead of a warning, since any
                    # non-int value is meaningless
                    msg = (
                        "Column null option (TNULLn) must be an integer for "
                        f"binary table columns (got {null!r}).  The invalid value "
                        "will be ignored for the purpose of formatting "
                        "the data in this column."
                    )

                elif not (
                    format.format in tnull_formats
                    or (
                        format.format in ("P", "Q") and format.p_format in tnull_formats
                    )
                ):
                    # TODO: We should also check that TNULLn's integer value
                    # is in the range allowed by the column's format
                    msg = (
                        "Column null option (TNULLn) is invalid for binary "
                        "table columns of type {!r} (got {!r}).  The invalid "
                        "value will be ignored for the purpose of formatting "
                        "the data in this column.".format(format, null)
                    )

            if msg is None:
                valid["null"] = null
            else:
                invalid["null"] = (null, msg)

        # Validate the disp option
        # TODO: Add full parsing and validation of TDISPn keywords
        if disp is not None and disp != "":
            msg = None
            if not isinstance(disp, str):
                msg = (
                    "Column disp option (TDISPn) must be a string (got "
                    f"{disp!r}). The invalid value will be ignored for the "
                    "purpose of formatting the data in this column."
                )

            elif isinstance(format, _AsciiColumnFormat) and disp[0].upper() == "L":
                # disp is at least one character long and has the 'L' format
                # which is not recognized for ASCII tables
                msg = (
                    "Column disp option (TDISPn) may not use the 'L' format "
                    "with ASCII table columns.  The invalid value will be "
                    "ignored for the purpose of formatting the data in this "
                    "column."
                )

            if msg is None:
                try:
                    _parse_tdisp_format(disp)
                    valid["disp"] = disp
                except VerifyError as err:
                    msg = (
                        "Column disp option (TDISPn) failed verification: "
                        f"{err!s} The invalid value will be ignored for the "
                        "purpose of formatting the data in this column."
                    )
                    invalid["disp"] = (disp, msg)
            else:
                invalid["disp"] = (disp, msg)

        # Validate the start option
        if start is not None and start != "":
            msg = None
            if not isinstance(format, _AsciiColumnFormat):
                # The 'start' option only applies to ASCII columns
                msg = (
                    "Column start option (TBCOLn) is not allowed for binary "
                    f"table columns (got {start!r}).  The invalid keyword will be "
                    "ignored for the purpose of formatting the data in this "
                    "column."
                )
            else:
                try:
                    start = int(start)
                except (TypeError, ValueError):
                    pass

                if not _is_int(start) or start < 1:
                    msg = (
                        "Column start option (TBCOLn) must be a positive integer "
                        f"(got {start!r}).  The invalid value will be ignored for the "
                        "purpose of formatting the data in this column."
                    )

            if msg is None:
                valid["start"] = start
            else:
                invalid["start"] = (start, msg)

        # Process TDIMn options
        # ASCII table columns can't have a TDIMn keyword associated with it;
        # for now we just issue a warning and ignore it.
        # TODO: This should be checked by the FITS verification code
        if dim is not None and dim != "":
            msg = None
            dims_tuple = ()
            # NOTE: If valid, the dim keyword's value in the valid dict is
            # a tuple, not the original string; if invalid just the original
            # string is returned
            if isinstance(format, _AsciiColumnFormat):
                msg = (
                    "Column dim option (TDIMn) is not allowed for ASCII table "
                    f"columns (got {dim!r}).  The invalid keyword will be ignored "
                    "for the purpose of formatting this column."
                )

            elif isinstance(dim, str):
                dims_tuple = _parse_tdim(dim)
            elif isinstance(dim, tuple):
                dims_tuple = dim
            else:
                msg = (
                    "`dim` argument must be a string containing a valid value "
                    "for the TDIMn header keyword associated with this column, "
                    "or a tuple containing the C-order dimensions for the "
                    "column.  The invalid value will be ignored for the purpose "
                    "of formatting this column."
                )

            if dims_tuple:
                if isinstance(recformat, _FormatP):
                    # TDIMs have different meaning for VLA format,
                    # no warning should be thrown
                    msg = None
                elif reduce(operator.mul, dims_tuple) > format.repeat:
                    msg = (
                        "The repeat count of the column format {!r} for column {!r} "
                        "is fewer than the number of elements per the TDIM "
                        "argument {!r}.  The invalid TDIMn value will be ignored "
                        "for the purpose of formatting this column.".format(
                            name, format, dim
                        )
                    )

            if msg is None:
                valid["dim"] = dims_tuple
            else:
                invalid["dim"] = (dim, msg)

        if coord_type is not None and coord_type != "":
            msg = None
            if not isinstance(coord_type, str):
                msg = (
                    "Coordinate/axis type option (TCTYPn) must be a string "
                    "(got {!r}). The invalid keyword will be ignored for the "
                    "purpose of formatting this column.".format(coord_type)
                )
            elif len(coord_type) > 8:
                msg = (
                    "Coordinate/axis type option (TCTYPn) must be a string "
                    f"of atmost 8 characters (got {coord_type!r}). The invalid keyword "
                    "will be ignored for the purpose of formatting this "
                    "column."
                )

            if msg is None:
                valid["coord_type"] = coord_type
            else:
                invalid["coord_type"] = (coord_type, msg)

        if coord_unit is not None and coord_unit != "":
            msg = None
            if not isinstance(coord_unit, str):
                msg = (
                    "Coordinate/axis unit option (TCUNIn) must be a string "
                    "(got {!r}). The invalid keyword will be ignored for the "
                    "purpose of formatting this column.".format(coord_unit)
                )

            if msg is None:
                valid["coord_unit"] = coord_unit
            else:
                invalid["coord_unit"] = (coord_unit, msg)

        for k, v in [
            ("coord_ref_point", coord_ref_point),
            ("coord_ref_value", coord_ref_value),
            ("coord_inc", coord_inc),
        ]:
            if v is not None and v != "":
                msg = None
                if not isinstance(v, numbers.Real):
                    msg = (
                        "Column {} option ({}n) must be a real floating type (got"
                        " {!r}). The invalid value will be ignored for the purpose of"
                        " formatting the data in this column.".format(
                            k, ATTRIBUTE_TO_KEYWORD[k], v
                        )
                    )

                if msg is None:
                    valid[k] = v
                else:
                    invalid[k] = (v, msg)

        if time_ref_pos is not None and time_ref_pos != "":
            msg = None
            if not isinstance(time_ref_pos, str):
                msg = (
                    "Time coordinate reference position option (TRPOSn) must be "
                    "a string (got {!r}). The invalid keyword will be ignored for "
                    "the purpose of formatting this column.".format(time_ref_pos)
                )

            if msg is None:
                valid["time_ref_pos"] = time_ref_pos
            else:
                invalid["time_ref_pos"] = (time_ref_pos, msg)

        return valid, invalid

    @classmethod
    def _determine_formats(cls, format, start, dim, ascii):
        """
        Given a format string and whether or not the Column is for an
        ASCII table (ascii=None means unspecified, but lean toward binary table
        where ambiguous) create an appropriate _BaseColumnFormat instance for
        the column's format, and determine the appropriate recarray format.

        The values of the start and dim keyword arguments are also useful, as
        the former is only valid for ASCII tables and the latter only for
        BINARY tables.
        """
        # If the given format string is unambiguously a Numpy dtype or one of
        # the Numpy record format type specifiers supported by Astropy then that
        # should take priority--otherwise assume it is a FITS format
        if isinstance(format, np.dtype):
            format, _, _ = _dtype_to_recformat(format)

        # check format
        if ascii is None and not isinstance(format, _BaseColumnFormat):
            # We're just give a string which could be either a Numpy format
            # code, or a format for a binary column array *or* a format for an
            # ASCII column array--there may be many ambiguities here.  Try our
            # best to guess what the user intended.
            format, recformat = cls._guess_format(format, start, dim)
        elif not ascii and not isinstance(format, _BaseColumnFormat):
            format, recformat = cls._convert_format(format, _ColumnFormat)
        elif ascii and not isinstance(format, _AsciiColumnFormat):
            format, recformat = cls._convert_format(format, _AsciiColumnFormat)
        else:
            # The format is already acceptable and unambiguous
            recformat = format.recformat

        return format, recformat

    @classmethod
    def _guess_format(cls, format, start, dim):
        if start and dim:
            # This is impossible; this can't be a valid FITS column
            raise ValueError(
                "Columns cannot have both a start (TCOLn) and dim "
                "(TDIMn) option, since the former is only applies to "
                "ASCII tables, and the latter is only valid for binary tables."
            )
        elif start:
            # Only ASCII table columns can have a 'start' option
            guess_format = _AsciiColumnFormat
        elif dim:
            # Only binary tables can have a dim option
            guess_format = _ColumnFormat
        else:
            # If the format is *technically* a valid binary column format
            # (i.e. it has a valid format code followed by arbitrary
            # "optional" codes), but it is also strictly a valid ASCII
            # table format, then assume an ASCII table column was being
            # requested (the more likely case, after all).
            with suppress(VerifyError):
                format = _AsciiColumnFormat(format, strict=True)

            # A safe guess which reflects the existing behavior of previous
            # Astropy versions
            guess_format = _ColumnFormat

        try:
            format, recformat = cls._convert_format(format, guess_format)
        except VerifyError:
            # For whatever reason our guess was wrong (for example if we got
            # just 'F' that's not a valid binary format, but it an ASCII format
            # code albeit with the width/precision omitted
            guess_format = (
                _AsciiColumnFormat if guess_format is _ColumnFormat else _ColumnFormat
            )
            # If this fails too we're out of options--it is truly an invalid
            # format, or at least not supported
            format, recformat = cls._convert_format(format, guess_format)

        return format, recformat

    def _convert_to_valid_data_type(self, array):
        # Convert the format to a type we understand
        if isinstance(array, Delayed):
            return array
        elif array is None:
            return array
        else:
            format = self.format
            dims = self._dims
            if dims and format.format not in "PQ":
                shape = dims[:-1] if "A" in format else dims
                shape = (len(array),) + shape
                array = array.reshape(shape)

            if "P" in format or "Q" in format:
                return array
            elif "A" in format:
                if array.dtype.char in "SU":
                    if dims:
                        # The 'last' dimension (first in the order given
                        # in the TDIMn keyword itself) is the number of
                        # characters in each string
                        fsize = dims[-1]
                    else:
                        fsize = np.dtype(format.recformat).itemsize
                    return chararray.array(array, itemsize=fsize, copy=False)
                else:
                    return _convert_array(array, np.dtype(format.recformat))
            elif "L" in format:
                # boolean needs to be scaled back to storage values ('T', 'F')
                if array.dtype == np.dtype("bool"):
                    return np.where(array == np.False_, ord("F"), ord("T"))
                else:
                    return np.where(array == 0, ord("F"), ord("T"))
            elif "X" in format:
                return _convert_array(array, np.dtype("uint8"))
            else:
                # Preserve byte order of the original array for now; see #77
                numpy_format = array.dtype.byteorder + format.recformat

                # Handle arrays passed in as unsigned ints as pseudo-unsigned
                # int arrays; blatantly tacked in here for now--we need columns
                # to have explicit knowledge of whether they treated as
                # pseudo-unsigned
                bzeros = {
                    2: np.uint16(2**15),
                    4: np.uint32(2**31),
                    8: np.uint64(2**63),
                }
                if (
                    array.dtype.kind == "u"
                    and array.dtype.itemsize in bzeros
                    and self.bscale in (1, None, "")
                    and self.bzero == bzeros[array.dtype.itemsize]
                ):
                    # Basically the array is uint, has scale == 1.0, and the
                    # bzero is the appropriate value for a pseudo-unsigned
                    # integer of the input dtype, then go ahead and assume that
                    # uint is assumed
                    numpy_format = numpy_format.replace("i", "u")
                    self._pseudo_unsigned_ints = True

                # The .base here means we're dropping the shape information,
                # which is only used to format recarray fields, and is not
                # useful for converting input arrays to the correct data type
                dtype = np.dtype(numpy_format).base

                return _convert_array(array, dtype)


class ColDefs(NotifierMixin):
    """
    Column definitions class.

    It has attributes corresponding to the `Column` attributes
    (e.g. `ColDefs` has the attribute ``names`` while `Column`
    has ``name``). Each attribute in `ColDefs` is a list of
    corresponding attribute values from all `Column` objects.
    """

    _padding_byte = "\x00"
    _col_format_cls = _ColumnFormat

    def __new__(cls, input, ascii=False):
        klass = cls

        if hasattr(input, "_columns_type") and issubclass(input._columns_type, ColDefs):
            klass = input._columns_type
        elif hasattr(input, "_col_format_cls") and issubclass(
            input._col_format_cls, _AsciiColumnFormat
        ):
            klass = _AsciiColDefs

        if ascii:  # force ASCII if this has been explicitly requested
            klass = _AsciiColDefs

        return object.__new__(klass)

    def __getnewargs__(self):
        return (self._arrays,)

    def __init__(self, input, ascii=False):
        """
        Parameters
        ----------
        input : sequence of `Column` or `ColDefs` or ndarray or `~numpy.recarray`
            An existing table HDU, an existing `ColDefs`, or any multi-field
            Numpy array or `numpy.recarray`.

        ascii : bool
            Use True to ensure that ASCII table columns are used.

        """
        from .fitsrec import FITS_rec
        from .hdu.table import _TableBaseHDU

        if isinstance(input, ColDefs):
            self._init_from_coldefs(input)
        elif (
            isinstance(input, FITS_rec)
            and hasattr(input, "_coldefs")
            and input._coldefs
        ):
            # If given a FITS_rec object we can directly copy its columns, but
            # only if its columns have already been defined, otherwise this
            # will loop back in on itself and blow up
            self._init_from_coldefs(input._coldefs)
        elif isinstance(input, np.ndarray) and input.dtype.fields is not None:
            # Construct columns from the fields of a record array
            self._init_from_array(input)
        elif isiterable(input):
            # if the input is a list of Columns
            self._init_from_sequence(input)
        elif isinstance(input, _TableBaseHDU):
            # Construct columns from fields in an HDU header
            self._init_from_table(input)
        else:
            raise TypeError(
                "Input to ColDefs must be a table HDU, a list "
                "of Columns, or a record/field array."
            )

        # Listen for changes on all columns
        for col in self.columns:
            col._add_listener(self)

    def _init_from_coldefs(self, coldefs):
        """Initialize from an existing ColDefs object (just copy the
        columns and convert their formats if necessary).
        """
        self.columns = [self._copy_column(col) for col in coldefs]

    def _init_from_sequence(self, columns):
        for idx, col in enumerate(columns):
            if not isinstance(col, Column):
                raise TypeError(f"Element {idx} in the ColDefs input is not a Column.")

        self._init_from_coldefs(columns)

    def _init_from_array(self, array):
        self.columns = []
        for idx in range(len(array.dtype)):
            cname = array.dtype.names[idx]
            ftype = array.dtype.fields[cname][0]

            if ftype.kind == "O":
                dtypes = {np.array(array[cname][i]).dtype for i in range(len(array))}
                if (len(dtypes) > 1) or (np.dtype("O") in dtypes):
                    raise TypeError(
                        f"Column '{cname}' contains unsupported object types or "
                        f"mixed types: {dtypes}"
                    )
                ftype = dtypes.pop()
                format = self._col_format_cls.from_recformat(ftype)
                format = f"P{format}()"
            else:
                format = self._col_format_cls.from_recformat(ftype)

            # Determine the appropriate dimensions for items in the column
            dim = array.dtype[idx].shape[::-1]
            if dim and (len(dim) > 0 or "A" in format):
                if "A" in format:
                    # should take into account multidimensional items in the column
                    dimel = int(re.findall("[0-9]+", str(ftype.subdtype[0]))[0])
                    # n x m string arrays must include the max string
                    # length in their dimensions (e.g. l x n x m)
                    dim = (dimel,) + dim
                dim = "(" + ",".join(str(d) for d in dim) + ")"
            else:
                dim = None

            # Check for unsigned ints.
            bzero = None
            if ftype.base.kind == "u":
                if "I" in format:
                    bzero = np.uint16(2**15)
                elif "J" in format:
                    bzero = np.uint32(2**31)
                elif "K" in format:
                    bzero = np.uint64(2**63)

            c = Column(
                name=cname,
                format=format,
                array=array.view(np.ndarray)[cname],
                bzero=bzero,
                dim=dim,
            )
            self.columns.append(c)

    def _init_from_table(self, table):
        hdr = table._header
        nfields = hdr["TFIELDS"]

        # go through header keywords to pick out column definition keywords
        # definition dictionaries for each field
        col_keywords = [{} for i in range(nfields)]
        for keyword in hdr:
            key = TDEF_RE.match(keyword)
            try:
                label = key.group("label")
            except Exception:
                continue  # skip if there is no match
            if label in KEYWORD_NAMES:
                col = int(key.group("num"))
                if 0 < col <= nfields:
                    attr = KEYWORD_TO_ATTRIBUTE[label]
                    value = hdr[keyword]
                    if attr == "format":
                        # Go ahead and convert the format value to the
                        # appropriate ColumnFormat container now
                        value = self._col_format_cls(value)
                    col_keywords[col - 1][attr] = value

        # Verify the column keywords and display any warnings if necessary;
        # we only want to pass on the valid keywords
        for idx, kwargs in enumerate(col_keywords):
            valid_kwargs, invalid_kwargs = Column._verify_keywords(**kwargs)
            for val in invalid_kwargs.values():
                warnings.warn(
                    f"Invalid keyword for column {idx + 1}: {val[1]}", VerifyWarning
                )
            # Special cases for recformat and dim
            # TODO: Try to eliminate the need for these special cases
            del valid_kwargs["recformat"]
            if "dim" in valid_kwargs:
                valid_kwargs["dim"] = kwargs["dim"]
            col_keywords[idx] = valid_kwargs

        # data reading will be delayed
        for col in range(nfields):
            col_keywords[col]["array"] = Delayed(table, col)

        # now build the columns
        self.columns = [Column(**attrs) for attrs in col_keywords]

        # Add the table HDU is a listener to changes to the columns
        # (either changes to individual columns, or changes to the set of
        # columns (add/remove/etc.))
        self._add_listener(table)

    def __copy__(self):
        return self.__class__(self)

    def __deepcopy__(self, memo):
        return self.__class__([copy.deepcopy(c, memo) for c in self.columns])

    def _copy_column(self, column):
        """Utility function used currently only by _init_from_coldefs
        to help convert columns from binary format to ASCII format or vice
        versa if necessary (otherwise performs a straight copy).
        """
        if isinstance(column.format, self._col_format_cls):
            # This column has a FITS format compatible with this column
            # definitions class (that is ascii or binary)
            return column.copy()

        new_column = column.copy()

        # Try to use the Numpy recformat as the equivalency between the
        # two formats; if that conversion can't be made then these
        # columns can't be transferred
        # TODO: Catch exceptions here and raise an explicit error about
        # column format conversion
        new_column.format = self._col_format_cls.from_column_format(column.format)

        # Handle a few special cases of column format options that are not
        # compatible between ASCII an binary tables
        # TODO: This is sort of hacked in right now; we really need
        # separate classes for ASCII and Binary table Columns, and they
        # should handle formatting issues like these
        if not isinstance(new_column.format, _AsciiColumnFormat):
            # the column is a binary table column...
            new_column.start = None
            if new_column.null is not None:
                # We can't just "guess" a value to represent null
                # values in the new column, so just disable this for
                # now; users may modify it later
                new_column.null = None
        else:
            # the column is an ASCII table column...
            if new_column.null is not None:
                new_column.null = DEFAULT_ASCII_TNULL
            if new_column.disp is not None and new_column.disp.upper().startswith("L"):
                # ASCII columns may not use the logical data display format;
                # for now just drop the TDISPn option for this column as we
                # don't have a systematic conversion of boolean data to ASCII
                # tables yet
                new_column.disp = None

        return new_column

    def __getattr__(self, name):
        """
        Automatically returns the values for the given keyword attribute for
        all `Column`s in this list.

        Implements for example self.units, self.formats, etc.
        """
        cname = name[:-1]
        if cname in KEYWORD_ATTRIBUTES and name[-1] == "s":
            attr = []
            for col in self.columns:
                val = getattr(col, cname)
                attr.append(val if val is not None else "")
            return attr
        raise AttributeError(name)

    @lazyproperty
    def dtype(self):
        # Note: This previously returned a dtype that just used the raw field
        # widths based on the format's repeat count, and did not incorporate
        # field *shapes* as provided by TDIMn keywords.
        # Now this incorporates TDIMn from the start, which makes *this* method
        # a little more complicated, but simplifies code elsewhere (for example
        # fields will have the correct shapes even in the raw recarray).
        formats = []
        offsets = [0]

        for format_, dim in zip(self.formats, self._dims):
            dt = format_.dtype

            if len(offsets) < len(self.formats):
                # Note: the size of the *original* format_ may be greater than
                # one would expect from the number of elements determined by
                # dim.  The FITS format allows this--the rest of the field is
                # filled with undefined values.
                offsets.append(offsets[-1] + dt.itemsize)

            if dim and format_.format not in "PQ":
                # Note: VLA array descriptors should not be reshaped
                # as they are always of shape (2,)
                if format_.format == "A":
                    dt = np.dtype((dt.char + str(dim[-1]), dim[:-1]))
                else:
                    dt = np.dtype((dt.base, dim))

            formats.append(dt)

        return np.dtype({"names": self.names, "formats": formats, "offsets": offsets})

    @lazyproperty
    def names(self):
        return [col.name for col in self.columns]

    @lazyproperty
    def formats(self):
        return [col.format for col in self.columns]

    @lazyproperty
    def _arrays(self):
        return [col.array for col in self.columns]

    @lazyproperty
    def _recformats(self):
        return [fmt.recformat for fmt in self.formats]

    @lazyproperty
    def _dims(self):
        """Returns the values of the TDIMn keywords parsed into tuples."""
        return [col._dims for col in self.columns]

    def __getitem__(self, key):
        if isinstance(key, str):
            key = _get_index(self.names, key)

        x = self.columns[key]
        if _is_int(key):
            return x
        else:
            return ColDefs(x)

    def __len__(self):
        return len(self.columns)

    def __repr__(self):
        rep = "ColDefs("
        if hasattr(self, "columns") and self.columns:
            # The hasattr check is mostly just useful in debugging sessions
            # where self.columns may not be defined yet
            rep += "\n    "
            rep += "\n    ".join([repr(c) for c in self.columns])
            rep += "\n"
        rep += ")"
        return rep

    def __add__(self, other, option="left"):
        if isinstance(other, Column):
            b = [other]
        elif isinstance(other, ColDefs):
            b = list(other.columns)
        else:
            raise TypeError("Wrong type of input.")
        if option == "left":
            tmp = list(self.columns) + b
        else:
            tmp = b + list(self.columns)
        return ColDefs(tmp)

    def __radd__(self, other):
        return self.__add__(other, "right")

    def __sub__(self, other):
        if not isinstance(other, (list, tuple)):
            other = [other]
        _other = [_get_index(self.names, key) for key in other]
        indx = list(range(len(self)))
        for x in _other:
            indx.remove(x)
        tmp = [self[i] for i in indx]
        return ColDefs(tmp)

    def _update_column_attribute_changed(self, column, attr, old_value, new_value):
        """
        Handle column attribute changed notifications from columns that are
        members of this `ColDefs`.

        `ColDefs` itself does not currently do anything with this, and just
        bubbles the notification up to any listening table HDUs that may need
        to update their headers, etc.  However, this also informs the table of
        the numerical index of the column that changed.
        """
        idx = 0
        for idx, col in enumerate(self.columns):
            if col is column:
                break

        if attr == "name":
            del self.names
        elif attr == "format":
            del self.formats

        self._notify(
            "column_attribute_changed", column, idx, attr, old_value, new_value
        )

    def add_col(self, column):
        """
        Append one `Column` to the column definition.
        """
        if not isinstance(column, Column):
            raise AssertionError

        # Ask the HDU object to load the data before we modify our columns
        self._notify("load_data")

        self._arrays.append(column.array)
        # Obliterate caches of certain things
        del self.dtype
        del self._recformats
        del self._dims
        del self.names
        del self.formats

        self.columns.append(column)

        # Listen for changes on the new column
        column._add_listener(self)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._notify("column_added", self, column)
        return self

    def del_col(self, col_name):
        """
        Delete (the definition of) one `Column`.

        col_name : str or int
            The column's name or index
        """
        # Ask the HDU object to load the data before we modify our columns
        self._notify("load_data")

        indx = _get_index(self.names, col_name)
        col = self.columns[indx]

        del self._arrays[indx]
        # Obliterate caches of certain things
        del self.dtype
        del self._recformats
        del self._dims
        del self.names
        del self.formats

        del self.columns[indx]

        col._remove_listener(self)

        # If this ColDefs is being tracked by a table HDU, inform the HDU (or
        # any other listeners) that the column has been removed
        # Just send a reference to self, and the index of the column that was
        # removed
        self._notify("column_removed", self, indx)
        return self

    def change_attrib(self, col_name, attrib, new_value):
        """
        Change an attribute (in the ``KEYWORD_ATTRIBUTES`` list) of a `Column`.

        Parameters
        ----------
        col_name : str or int
            The column name or index to change

        attrib : str
            The attribute name

        new_value : object
            The new value for the attribute
        """
        setattr(self[col_name], attrib, new_value)

    def change_name(self, col_name, new_name):
        """
        Change a `Column`'s name.

        Parameters
        ----------
        col_name : str
            The current name of the column

        new_name : str
            The new name of the column
        """
        if new_name != col_name and new_name in self.names:
            raise ValueError(f"New name {new_name} already exists.")
        else:
            self.change_attrib(col_name, "name", new_name)

    def change_unit(self, col_name, new_unit):
        """
        Change a `Column`'s unit.

        Parameters
        ----------
        col_name : str or int
            The column name or index

        new_unit : str
            The new unit for the column
        """
        self.change_attrib(col_name, "unit", new_unit)

    def info(self, attrib="all", output=None):
        """
        Get attribute(s) information of the column definition.

        Parameters
        ----------
        attrib : str
            Can be one or more of the attributes listed in
            ``astropy.io.fits.column.KEYWORD_ATTRIBUTES``.  The default is
            ``"all"`` which will print out all attributes.  It forgives plurals
            and blanks.  If there are two or more attribute names, they must be
            separated by comma(s).

        output : file-like, optional
            File-like object to output to.  Outputs to stdout by default.
            If `False`, returns the attributes as a `dict` instead.

        Notes
        -----
        This function doesn't return anything by default; it just prints to
        stdout.
        """
        if output is None:
            output = sys.stdout

        if attrib.strip().lower() in ["all", ""]:
            lst = KEYWORD_ATTRIBUTES
        else:
            lst = attrib.split(",")
            for idx in range(len(lst)):
                lst[idx] = lst[idx].strip().lower()
                if lst[idx][-1] == "s":
                    lst[idx] = list[idx][:-1]

        ret = {}

        for attr in lst:
            if output:
                if attr not in KEYWORD_ATTRIBUTES:
                    output.write(
                        f"'{attr}' is not an attribute of the column definitions.\n"
                    )
                    continue
                output.write(f"{attr}:\n")
                output.write(f"    {getattr(self, attr + 's')}\n")
            else:
                ret[attr] = getattr(self, attr + "s")

        if not output:
            return ret


class _AsciiColDefs(ColDefs):
    """ColDefs implementation for ASCII tables."""

    _padding_byte = " "
    _col_format_cls = _AsciiColumnFormat

    def __init__(self, input, ascii=True):
        super().__init__(input)

        # if the format of an ASCII column has no width, add one
        if not isinstance(input, _AsciiColDefs):
            self._update_field_metrics()
        else:
            for idx, s in enumerate(input.starts):
                self.columns[idx].start = s

            self._spans = input.spans
            self._width = input._width

    @lazyproperty
    def dtype(self):
        dtype = {}

        for j in range(len(self)):
            data_type = "S" + str(self.spans[j])
            dtype[self.names[j]] = (data_type, self.starts[j] - 1)

        return np.dtype(dtype)

    @property
    def spans(self):
        """A list of the widths of each field in the table."""
        return self._spans

    @lazyproperty
    def _recformats(self):
        if len(self) == 1:
            widths = []
        else:
            widths = [y - x for x, y in pairwise(self.starts)]

        # Widths is the width of each field *including* any space between
        # fields; this is so that we can map the fields to string records in a
        # Numpy recarray
        widths.append(self._width - self.starts[-1] + 1)
        return ["S" + str(w) for w in widths]

    def add_col(self, column):
        super().add_col(column)
        self._update_field_metrics()

    def del_col(self, col_name):
        super().del_col(col_name)
        self._update_field_metrics()

    def _update_field_metrics(self):
        """
        Updates the list of the start columns, the list of the widths of each
        field, and the total width of each record in the table.
        """
        spans = [0] * len(self.columns)
        end_col = 0  # Refers to the ASCII text column, not the table col
        for idx, col in enumerate(self.columns):
            width = col.format.width

            # Update the start columns and column span widths taking into
            # account the case that the starting column of a field may not
            # be the column immediately after the previous field
            if not col.start:
                col.start = end_col + 1
            end_col = col.start + width - 1
            spans[idx] = width

        self._spans = spans
        self._width = end_col


# Utilities


class _VLF(np.ndarray):
    """Variable length field object."""

    def __new__(cls, input, dtype="S"):
        """
        Parameters
        ----------
        input
            a sequence of variable-sized elements.
        """
        if dtype == "S":
            try:
                # this handles ['abc'] and [['a','b','c']]
                # equally, beautiful!
                input = [chararray.array(x, itemsize=1) for x in input]
            except Exception:
                raise ValueError(f"Inconsistent input data array: {input}")

        a = np.array(input, dtype=object)
        self = np.ndarray.__new__(cls, shape=(len(input),), buffer=a, dtype=object)
        self.max = 0
        self.element_dtype = dtype
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.max = obj.max
        self.element_dtype = obj.element_dtype

    def __setitem__(self, key, value):
        """
        To make sure the new item has consistent data type to avoid
        misalignment.
        """
        if isinstance(value, np.ndarray) and value.dtype == self.dtype:
            pass
        elif isinstance(value, chararray.chararray) and value.itemsize == 1:
            pass
        elif self.element_dtype == "S":
            value = chararray.array(value, itemsize=1)
        else:
            value = np.array(value, dtype=self.element_dtype)
        np.ndarray.__setitem__(self, key, value)
        nelem = value.shape
        len_value = np.prod(nelem)
        self.max = max(self.max, len_value)

    def tolist(self):
        return [list(item) for item in super().tolist()]


def _get_index(names, key):
    """
    Get the index of the ``key`` in the ``names`` list.

    The ``key`` can be an integer or string.  If integer, it is the index
    in the list.  If string,

        a. Field (column) names are case sensitive: you can have two
           different columns called 'abc' and 'ABC' respectively.

        b. When you *refer* to a field (presumably with the field
           method), it will try to match the exact name first, so in
           the example in (a), field('abc') will get the first field,
           and field('ABC') will get the second field.

        If there is no exact name matched, it will try to match the
        name with case insensitivity.  So, in the last example,
        field('Abc') will cause an exception since there is no unique
        mapping.  If there is a field named "XYZ" and no other field
        name is a case variant of "XYZ", then field('xyz'),
        field('Xyz'), etc. will get this field.
    """
    if _is_int(key):
        indx = int(key)
    elif isinstance(key, str):
        # try to find exact match first
        try:
            indx = names.index(key.rstrip())
        except ValueError:
            # try to match case-insentively,
            _key = key.lower().rstrip()
            names = [n.lower().rstrip() for n in names]
            count = names.count(_key)  # occurrence of _key in names
            if count == 1:
                indx = names.index(_key)
            elif count == 0:
                raise KeyError(f"Key '{key}' does not exist.")
            else:  # multiple match
                raise KeyError(f"Ambiguous key name '{key}'.")
    else:
        raise KeyError(f"Illegal key '{key!r}'.")

    return indx


def _unwrapx(input, output, repeat):
    """
    Unwrap the X format column into a Boolean array.

    Parameters
    ----------
    input
        input ``Uint8`` array of shape (`s`, `nbytes`)

    output
        output Boolean array of shape (`s`, `repeat`)

    repeat
        number of bits
    """
    pow2 = np.array([128, 64, 32, 16, 8, 4, 2, 1], dtype="uint8")
    nbytes = ((repeat - 1) // 8) + 1
    for i in range(nbytes):
        _min = i * 8
        _max = min((i + 1) * 8, repeat)
        for j in range(_min, _max):
            output[..., j] = np.bitwise_and(input[..., i], pow2[j - i * 8])


def _wrapx(input, output, repeat):
    """
    Wrap the X format column Boolean array into an ``UInt8`` array.

    Parameters
    ----------
    input
        input Boolean array of shape (`s`, `repeat`)

    output
        output ``Uint8`` array of shape (`s`, `nbytes`)

    repeat
        number of bits
    """
    output[...] = 0  # reset the output
    nbytes = ((repeat - 1) // 8) + 1
    unused = nbytes * 8 - repeat
    for i in range(nbytes):
        _min = i * 8
        _max = min((i + 1) * 8, repeat)
        for j in range(_min, _max):
            if j != _min:
                np.left_shift(output[..., i], 1, output[..., i])
            np.add(output[..., i], input[..., j], output[..., i])

    # shift the unused bits
    np.left_shift(output[..., i], unused, output[..., i])


def _makep(array, descr_output, format, nrows=None):
    """
    Construct the P (or Q) format column array, both the data descriptors and
    the data.  It returns the output "data" array of data type `dtype`.

    The descriptor location will have a zero offset for all columns
    after this call.  The final offset will be calculated when the file
    is written.

    Parameters
    ----------
    array
        input object array

    descr_output
        output "descriptor" array of data type int32 (for P format arrays) or
        int64 (for Q format arrays)--must be nrows long in its first dimension

    format
        the _FormatP object representing the format of the variable array

    nrows : int, optional
        number of rows to create in the column; defaults to the number of rows
        in the input array
    """
    # TODO: A great deal of this is redundant with FITS_rec._convert_p; see if
    # we can merge the two somehow.

    _offset = 0

    if not nrows:
        nrows = len(array)

    data_output = _VLF([None] * nrows, dtype=format.dtype)

    if format.dtype == "S":
        _nbytes = 1
    else:
        _nbytes = np.array([], dtype=format.dtype).itemsize

    for idx in range(nrows):
        if idx < len(array):
            rowval = array[idx]
        else:
            if format.dtype == "S":
                rowval = " " * data_output.max
            else:
                rowval = [0] * data_output.max
        if format.dtype == "S":
            data_output[idx] = chararray.array(encode_ascii(rowval), itemsize=1)
        else:
            data_output[idx] = np.array(rowval, dtype=format.dtype)

        nelem = data_output[idx].shape
        descr_output[idx, 0] = np.prod(nelem)
        descr_output[idx, 1] = _offset
        _offset += descr_output[idx, 0] * _nbytes

    return data_output


def _parse_tformat(tform):
    """Parse ``TFORMn`` keyword for a binary table into a
    ``(repeat, format, option)`` tuple.
    """
    try:
        (repeat, format, option) = TFORMAT_RE.match(tform.strip()).groups()
    except Exception:
        # TODO: Maybe catch this error use a default type (bytes, maybe?) for
        # unrecognized column types.  As long as we can determine the correct
        # byte width somehow..
        raise VerifyError(f"Format {tform!r} is not recognized.")

    if repeat == "":
        repeat = 1
    else:
        repeat = int(repeat)

    return (repeat, format.upper(), option)


def _parse_ascii_tformat(tform, strict=False):
    """
    Parse the ``TFORMn`` keywords for ASCII tables into a ``(format, width,
    precision)`` tuple (the latter is always zero unless format is one of 'E',
    'F', or 'D').
    """
    match = TFORMAT_ASCII_RE.match(tform.strip())
    if not match:
        raise VerifyError(f"Format {tform!r} is not recognized.")

    # Be flexible on case
    format = match.group("format")
    if format is None:
        # Floating point format
        format = match.group("formatf").upper()
        width = match.group("widthf")
        precision = match.group("precision")
        if width is None or precision is None:
            if strict:
                raise VerifyError(
                    "Format {!r} is not unambiguously an ASCII table format."
                )
            else:
                width = 0 if width is None else width
                precision = 1 if precision is None else precision
    else:
        format = format.upper()
        width = match.group("width")
        if width is None:
            if strict:
                raise VerifyError(
                    "Format {!r} is not unambiguously an ASCII table format."
                )
            else:
                # Just use a default width of 0 if unspecified
                width = 0
        precision = 0

    def convert_int(val):
        msg = (
            "Format {!r} is not valid--field width and decimal precision "
            "must be integers."
        )
        try:
            val = int(val)
        except (ValueError, TypeError):
            raise VerifyError(msg.format(tform))

        return val

    if width and precision:
        # This should only be the case for floating-point formats
        width, precision = convert_int(width), convert_int(precision)
    elif width:
        # Just for integer/string formats; ignore precision
        width = convert_int(width)
    else:
        # For any format, if width was unspecified use the set defaults
        width, precision = ASCII_DEFAULT_WIDTHS[format]

    if width <= 0:
        raise VerifyError(
            f"Format {tform!r} not valid--field width must be a positive integeter."
        )

    if precision >= width:
        raise VerifyError(
            f"Format {tform!r} not valid--the number of decimal digits "
            f"must be less than the format's total width {width}."
        )

    return format, width, precision


def _parse_tdim(tdim):
    """Parse the ``TDIM`` value into a tuple (may return an empty tuple if
    the value ``TDIM`` value is empty or invalid).
    """
    m = tdim and TDIM_RE.match(tdim)
    if m:
        dims = m.group("dims")
        return tuple(int(d.strip()) for d in dims.split(","))[::-1]

    # Ignore any dim values that don't specify a multidimensional column
    return ()


def _scalar_to_format(value):
    """
    Given a scalar value or string, returns the minimum FITS column format
    that can represent that value.  'minimum' is defined by the order given in
    FORMATORDER.
    """
    # First, if value is a string, try to convert to the appropriate scalar
    # value
    for type_ in (int, float, complex):
        try:
            value = type_(value)
            break
        except ValueError:
            continue

    numpy_dtype_str = np.min_scalar_type(value).str
    numpy_dtype_str = numpy_dtype_str[1:]  # Strip endianness

    try:
        fits_format = NUMPY2FITS[numpy_dtype_str]
        return FITSUPCONVERTERS.get(fits_format, fits_format)
    except KeyError:
        return "A" + str(len(value))


def _cmp_recformats(f1, f2):
    """
    Compares two numpy recformats using the ordering given by FORMATORDER.
    """
    if f1[0] == "S" and f2[0] == "S":
        return cmp(int(f1[1:]), int(f2[1:]))
    else:
        f1, f2 = NUMPY2FITS[f1], NUMPY2FITS[f2]
        return cmp(FORMATORDER.index(f1), FORMATORDER.index(f2))


def _convert_fits2record(format):
    """
    Convert FITS format spec to record format spec.
    """
    repeat, dtype, option = _parse_tformat(format)

    if dtype in FITS2NUMPY:
        if dtype == "A":
            output_format = FITS2NUMPY[dtype] + str(repeat)
            # to accommodate both the ASCII table and binary table column
            # format spec, i.e. A7 in ASCII table is the same as 7A in
            # binary table, so both will produce 'S7'.
            # Technically the FITS standard does not allow this but it's a very
            # common mistake
            if format.lstrip()[0] == "A" and option != "":
                # make sure option is integer
                output_format = FITS2NUMPY[dtype] + str(int(option))
        else:
            repeat_str = ""
            if repeat != 1:
                repeat_str = str(repeat)
            output_format = repeat_str + FITS2NUMPY[dtype]

    elif dtype == "X":
        output_format = _FormatX(repeat)
    elif dtype == "P":
        output_format = _FormatP.from_tform(format)
    elif dtype == "Q":
        output_format = _FormatQ.from_tform(format)
    elif dtype == "F":
        output_format = "f8"
    else:
        raise ValueError(f"Illegal format `{format}`.")

    return output_format


def _convert_record2fits(format):
    """
    Convert record format spec to FITS format spec.
    """
    recformat, kind, dtype = _dtype_to_recformat(format)
    shape = dtype.shape
    itemsize = dtype.base.itemsize
    if dtype.char == "U" or (
        dtype.subdtype is not None and dtype.subdtype[0].char == "U"
    ):
        # Unicode dtype--itemsize is 4 times actual ASCII character length,
        # which what matters for FITS column formats
        # Use dtype.base and dtype.subdtype --dtype for multi-dimensional items
        itemsize = itemsize // 4

    option = str(itemsize)

    ndims = len(shape)
    repeat = 1
    if ndims > 0:
        nel = np.array(shape, dtype="i8").prod()
        if nel > 1:
            repeat = nel

    if kind == "S":
        # This is a kludge that will place string arrays into a
        # single field, so at least we won't lose data.  Need to
        # use a TDIM keyword to fix this, declaring as (slength,
        # dim1, dim2, ...)  as mwrfits does

        ntot = int(repeat) * int(option)

        output_format = str(ntot) + "A"
    elif recformat in NUMPY2FITS:  # record format
        if repeat != 1:
            repeat = str(repeat)
        else:
            repeat = ""
        output_format = repeat + NUMPY2FITS[recformat]
    else:
        raise ValueError(f"Illegal format `{format}`.")

    return output_format


def _dtype_to_recformat(dtype):
    """
    Utility function for converting a dtype object or string that instantiates
    a dtype (e.g. 'float32') into one of the two character Numpy format codes
    that have been traditionally used by Astropy.
    """
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    kind = dtype.base.kind

    if kind in ("U", "S"):
        recformat = kind = "S"
    else:
        itemsize = dtype.base.itemsize
        recformat = kind + str(itemsize)

    return recformat, kind, dtype


def _convert_format(format, reverse=False):
    """
    Convert FITS format spec to record format spec.  Do the opposite if
    reverse=True.
    """
    if reverse:
        return _convert_record2fits(format)
    else:
        return _convert_fits2record(format)


def _convert_ascii_format(format, reverse=False):
    """Convert ASCII table format spec to record format spec."""
    if reverse:
        recformat, kind, dtype = _dtype_to_recformat(format)
        itemsize = dtype.itemsize

        if kind == "S":
            return "A" + str(itemsize)
        elif NUMPY2FITS.get(recformat) == "L":
            # Special case for logical/boolean types--for ASCII tables we
            # represent these as single character columns containing 'T' or 'F'
            # (a la the storage format for Logical columns in binary tables)
            return "A1"
        elif kind == "i":
            # Use for the width the maximum required to represent integers
            # of that byte size plus 1 for signs, but use a minimum of the
            # default width (to keep with existing behavior)
            width = 1 + len(str(2 ** (itemsize * 8)))
            width = max(width, ASCII_DEFAULT_WIDTHS["I"][0])
            return "I" + str(width)
        elif kind == "f":
            # This is tricky, but go ahead and use D if float-64, and E
            # if float-32 with their default widths
            if itemsize >= 8:
                format = "D"
            else:
                format = "E"
            width = ".".join(str(w) for w in ASCII_DEFAULT_WIDTHS[format])
            return format + width
        # TODO: There may be reasonable ways to represent other Numpy types so
        # let's see what other possibilities there are besides just 'S', 'i',
        # and 'f'.  If it doesn't have a reasonable ASCII representation then
        # raise an exception
    else:
        format, width, precision = _parse_ascii_tformat(format)

        # This gives a sensible "default" dtype for a given ASCII
        # format code
        recformat = ASCII2NUMPY[format]

        # The following logic is taken from CFITSIO:
        # For integers, if the width <= 4 we can safely use 16-bit ints for all
        # values, if width >= 10 we may need to accommodate 64-bit ints.
        # values [for the non-standard J format code just always force 64-bit]
        if format == "I":
            if width <= 4:
                recformat = "i2"
            elif width > 9:
                recformat = "i8"
        elif format == "A":
            recformat += str(width)

        return recformat


def _parse_tdisp_format(tdisp):
    """
    Parse the ``TDISPn`` keywords for ASCII and binary tables into a
    ``(format, width, precision, exponential)`` tuple (the TDISP values
    for ASCII and binary are identical except for 'Lw',
    which is only present in BINTABLE extensions.

    Parameters
    ----------
    tdisp : str
        TDISPn FITS Header keyword.  Used to specify display formatting.

    Returns
    -------
    formatc: str
        The format characters from TDISPn
    width: str
        The width int value from TDISPn
    precision: str
        The precision int value from TDISPn
    exponential: str
        The exponential int value from TDISPn

    """
    # Use appropriate regex for format type
    tdisp = tdisp.strip()
    fmt_key = (
        tdisp[0]
        if tdisp[0] != "E" or (len(tdisp) > 1 and tdisp[1] not in "NS")
        else tdisp[:2]
    )
    try:
        tdisp_re = TDISP_RE_DICT[fmt_key]
    except KeyError:
        raise VerifyError(f"Format {tdisp} is not recognized.")

    match = tdisp_re.match(tdisp.strip())
    if not match or match.group("formatc") is None:
        raise VerifyError(f"Format {tdisp} is not recognized.")

    formatc = match.group("formatc")
    width = match.group("width")
    precision = None
    exponential = None

    # Some formats have precision and exponential
    if tdisp[0] in ("I", "B", "O", "Z", "F", "E", "G", "D"):
        precision = match.group("precision")
        if precision is None:
            precision = 1
    if tdisp[0] in ("E", "D", "G") and tdisp[1] not in ("N", "S"):
        exponential = match.group("exponential")
        if exponential is None:
            exponential = 1

    # Once parsed, check format dict to do conversion to a formatting string
    return formatc, width, precision, exponential


def _fortran_to_python_format(tdisp):
    """
    Turn the TDISPn fortran format pieces into a final Python format string.
    See the format_type definitions above the TDISP_FMT_DICT. If codes is
    changed to take advantage of the exponential specification, will need to
    add it as another input parameter.

    Parameters
    ----------
    tdisp : str
        TDISPn FITS Header keyword.  Used to specify display formatting.

    Returns
    -------
    format_string: str
        The TDISPn keyword string translated into a Python format string.
    """
    format_type, width, precision, exponential = _parse_tdisp_format(tdisp)

    try:
        fmt = TDISP_FMT_DICT[format_type]
        return fmt.format(width=width, precision=precision)

    except KeyError:
        raise VerifyError(f"Format {format_type} is not recognized.")


def python_to_tdisp(format_string, logical_dtype=False):
    """
    Turn the Python format string to a TDISP FITS compliant format string. Not
    all formats convert. these will cause a Warning and return None.

    Parameters
    ----------
    format_string : str
        TDISPn FITS Header keyword.  Used to specify display formatting.
    logical_dtype : bool
        True is this format type should be a logical type, 'L'. Needs special
        handling.

    Returns
    -------
    tdsip_string: str
        The TDISPn keyword string translated into a Python format string.
    """
    fmt_to_tdisp = {
        "a": "A",
        "s": "A",
        "d": "I",
        "b": "B",
        "o": "O",
        "x": "Z",
        "X": "Z",
        "f": "F",
        "F": "F",
        "g": "G",
        "G": "G",
        "e": "E",
        "E": "E",
    }

    if format_string in [None, "", "{}"]:
        return None

    # Strip out extra format characters that aren't a type or a width/precision
    if format_string[0] == "{" and format_string != "{}":
        fmt_str = format_string.lstrip("{:").rstrip("}")
    elif format_string[0] == "%":
        fmt_str = format_string.lstrip("%")
    else:
        fmt_str = format_string

    precision, sep = "", ""

    # Character format, only translate right aligned, and don't take zero fills
    if fmt_str[-1].isdigit() and fmt_str[0] == ">" and fmt_str[1] != "0":
        ftype = fmt_to_tdisp["a"]
        width = fmt_str[1:]

    elif fmt_str[-1] == "s" and fmt_str != "s":
        ftype = fmt_to_tdisp["a"]
        width = fmt_str[:-1].lstrip("0")

    # Number formats, don't take zero fills
    elif fmt_str[-1].isalpha() and len(fmt_str) > 1 and fmt_str[0] != "0":
        ftype = fmt_to_tdisp[fmt_str[-1]]
        fmt_str = fmt_str[:-1]

        # If format has a "." split out the width and precision
        if "." in fmt_str:
            width, precision = fmt_str.split(".")
            sep = "."
            if width == "":
                key = ftype if ftype != "G" else "F"
                width = str(
                    int(precision)
                    + (ASCII_DEFAULT_WIDTHS[key][0] - ASCII_DEFAULT_WIDTHS[key][1])
                )
        # Otherwise we just have a width
        else:
            width = fmt_str

    else:
        warnings.warn(
            f"Format {format_string} cannot be mapped to the accepted TDISPn "
            "keyword values.  Format will not be moved into TDISPn keyword.",
            AstropyUserWarning,
        )
        return None

    # Catch logical data type, set the format type back to L in this case
    if logical_dtype:
        ftype = "L"

    return ftype + width + sep + precision
