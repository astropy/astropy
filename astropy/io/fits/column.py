# Licensed under a 3-clause BSD style license - see PYFITS.rst

import copy
import operator
import re
import sys
import warnings
import weakref

from functools import reduce

import numpy as np
from numpy import char as chararray

from .card import Card
from .util import pairwise, _is_int, _convert_array, encode_ascii, cmp
from .verify import VerifyError, VerifyWarning

from ...extern.six import string_types, iteritems
from ...utils import lazyproperty, isiterable, indent
from ...utils.compat import ignored


__all__ = ['Column', 'ColDefs', 'Delayed']


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
FITS2NUMPY = {'L': 'i1', 'B': 'u1', 'I': 'i2', 'J': 'i4', 'K': 'i8', 'E': 'f4',
              'D': 'f8', 'C': 'c8', 'M': 'c16', 'A': 'a'}

# the inverse dictionary of the above
NUMPY2FITS = dict([(val, key) for key, val in iteritems(FITS2NUMPY)])
# Normally booleans are represented as ints in pyfits, but if passed in a numpy
# boolean array, that should be supported
NUMPY2FITS['b1'] = 'L'
# Add unsigned types, which will be stored as signed ints with a TZERO card.
NUMPY2FITS['u2'] = 'I'
NUMPY2FITS['u4'] = 'J'
NUMPY2FITS['u8'] = 'K'

# This is the order in which values are converted to FITS types
# Note that only double precision floating point/complex are supported
FORMATORDER = ['L', 'B', 'I', 'J', 'K', 'D', 'M', 'A']

# mapping from ASCII table TFORM data type to numpy data type
# A: Character
# I: Integer (32-bit)
# J: Integer (64-bit; non-standard)
# F: Float (32-bit; fixed decimal notation)
# E: Float (32-bit; exponential notation)
# D: Float (64-bit; exponential notation, always 64-bit by convention)
ASCII2NUMPY = {'A': 'a', 'I': 'i4', 'J': 'i8', 'F': 'f4', 'E': 'f4',
               'D': 'f8'}

# Maps FITS ASCII column format codes to the appropriate Python string
# formatting codes for that type.
ASCII2STR = {'A': 's', 'I': 'd', 'J': 'd', 'F': 'f', 'E': 'E', 'D': 'E'}

# For each ASCII table format code, provides a default width (and decimal
# precision) for when one isn't given explicity in the column format
ASCII_DEFAULT_WIDTHS= {'A': (1, 0), 'I': (10, 0), 'J': (15, 0),
                       'E': (15, 7), 'F': (16, 7), 'D': (25, 17)}




# lists of column/field definition common names and keyword names, make
# sure to preserve the one-to-one correspondence when updating the list(s).
# Use lists, instead of dictionaries so the names can be displayed in a
# preferred order.
KEYWORD_NAMES = ['TTYPE', 'TFORM', 'TUNIT', 'TNULL', 'TSCAL', 'TZERO',
                 'TDISP', 'TBCOL', 'TDIM']
KEYWORD_ATTRIBUTES = ['name', 'format', 'unit', 'null', 'bscale', 'bzero',
                      'disp', 'start', 'dim']
"""This is a list of the attributes that can be set on `Column` objects."""

# TFORMn regular expression
TFORMAT_RE = re.compile(r'(?P<repeat>^[0-9]*)(?P<format>[LXBIJKAEDCMPQ])'
                        r'(?P<option>[!-~]*)', re.I)

# TFORMn for ASCII tables; two different versions depending on whether
# the format is floating-point or not; allows empty values for width
# in which case defaults are used
TFORMAT_ASCII_RE = re.compile(r'(?:(?P<format>[AIJ])(?P<width>[0-9]+)?)|'
                              r'(?:(?P<formatf>[FED])'
                              r'(?:(?P<widthf>[0-9]+)\.'
                              r'(?P<precision>[0-9]+))?)')

# table definition keyword regular expression
TDEF_RE = re.compile(r'(?P<label>^T[A-Z]*)(?P<num>[1-9][0-9 ]*$)')

# table dimension keyword regular expression (fairly flexible with whitespace)
TDIM_RE = re.compile(r'\(\s*(?P<dims>(?:\d+,\s*)+\s*\d+)\s*\)\s*')

# value for ASCII table cell with value = TNULL
# this can be reset by user.
ASCIITNULL = 0

# The default placeholder to use for NULL values in ASCII tables when
# converting from binary to ASCII tables
DEFAULT_ASCII_TNULL = '---'


class Delayed(object):
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
        self = super(_ColumnFormat, cls).__new__(cls, format)
        self.repeat, self.format, self.option = _parse_tformat(format)
        self.format = self.format.upper()
        if self.format in ('P', 'Q'):
            # TODO: There should be a generic factory that returns either
            # _FormatP or _FormatQ as appropriate for a given TFORMn
            if self.format == 'P':
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
            repeat = ''
        else:
            repeat = str(self.repeat)

        return '%s%s%s' % (repeat, self.format, self.option)


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
    possible and may result in a `~.exceptions.ValueError`.
    """

    def __new__(cls, format, strict=False):
        self = super(_AsciiColumnFormat, cls).__new__(cls, format)
        self.format, self.width, self.precision = \
            _parse_ascii_tformat(format, strict)

        # This is to support handling logical (boolean) data from binary tables
        # in an ASCII table
        self._pseudo_logical = False
        return self

    @classmethod
    def from_column_format(cls, format):
        inst = cls.from_recformat(format.recformat)
        # Hack
        if format.format == 'L':
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

        if self.format in ('E', 'F', 'D'):
            return '%s%s.%s' % (self.format, self.width, self.precision)

        return '%s%s' % (self.format, self.width)


class _FormatX(str):
    """For X format in binary tables."""

    def __new__(cls, repeat=1):
        nbytes = ((repeat - 1) // 8) + 1
        # use an array, even if it is only ONE u1 (i.e. use tuple always)
        obj = super(_FormatX, cls).__new__(cls, repr((nbytes,)) + 'u1')
        obj.repeat = repeat
        return obj

    @property
    def tform(self):
        return '%sX' % self.repeat


# TODO: Table column formats need to be verified upon first reading the file;
# as it is, an invalid P format will raise a VerifyError from some deep,
# unexpected place
class _FormatP(str):
    """For P format in variable length table."""

    # As far as I can tell from my reading of the FITS standard, a type code is
    # *required* for P and Q formats; there is no default
    _format_re_template = (r'(?P<repeat>\d+)?%s(?P<dtype>[LXBIJKAEDCM])'
                            '(?:\((?P<max>\d*)\))?')
    _format_code = 'P'
    _format_re = re.compile(_format_re_template % _format_code)
    _descriptor_format = '2i4'

    def __new__(cls, dtype, repeat=None, max=None):
        obj = super(_FormatP, cls).__new__(cls, cls._descriptor_format)
        obj.format = NUMPY2FITS[dtype]
        obj.dtype = dtype
        obj.repeat = repeat
        obj.max = max
        return obj

    @classmethod
    def from_tform(cls, format):
        m = cls._format_re.match(format)
        if not m or m.group('dtype') not in FITS2NUMPY:
            raise VerifyError('Invalid column format: %s' % format)
        repeat = m.group('repeat')
        array_dtype = m.group('dtype')
        max = m.group('max')
        if not max:
            max = None
        return cls(FITS2NUMPY[array_dtype], repeat=repeat, max=max)

    @property
    def tform(self):
        repeat = '' if self.repeat is None else self.repeat
        max = '' if self.max is None else self.max
        return '%s%s%s(%s)' % (repeat, self._format_code, self.format, max)


class _FormatQ(_FormatP):
    """Carries type description of the Q format for variable length arrays.

    The Q format is like the P format but uses 64-bit integers in the array
    descriptors, allowing for heaps stored beyond 2GB into a file.
    """

    _format_code = 'Q'
    _format_re = re.compile(_FormatP._format_re_template % _format_code)
    _descriptor_format = '2i8'


class Column(object):
    """
    Class which contains the definition of one column, e.g.  ``ttype``,
    ``tform``, etc. and the array containing values for the column.
    """

    def __init__(self, name=None, format=None, unit=None, null=None,
                 bscale=None, bzero=None, disp=None, start=None, dim=None,
                 array=None, ascii=None):
        """
        Construct a `Column` by specifying attributes.  All attributes
        except `format` can be optional.

        Parameters
        ----------
        name : str, optional
            column name, corresponding to ``TTYPE`` keyword

        format : str, optional
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
            initialize an ndarray) providing intial data for this column.
            The array will be automatically converted, if possible, to the data
            format of the column.  In the case were non-trivial ``bscale``
            and/or ``bzero`` arguments are given, the values in the array must
            be the *physical* values--that is, the values of column as if the
            scaling has already been applied (the array stored on the column
            object will then be converted back to its storage values).

        ascii : bool, optional
            set `True` if this describes a column for an ASCII table; this
            may be required to disambiguate the column format
        """

        if format is None:
            raise ValueError('Must specify format to construct Column.')

        # any of the input argument (except array) can be a Card or just
        # a number/string
        kwargs = {'ascii': ascii}
        for attr in KEYWORD_ATTRIBUTES:
            value = locals()[attr]  # get the argument's value

            if isinstance(value, Card):
                value = value.value

            kwargs[attr] = value

        valid_kwargs, invalid_kwargs = self._verify_keywords(**kwargs)

        if invalid_kwargs:
            msg = ['The following keyword arguments to Column were invalid:']

            for val in invalid_kwargs.values():
                msg.append(indent(val[1]))

            raise VerifyError('\n'.join(msg))


        for attr in KEYWORD_ATTRIBUTES:
            setattr(self, attr, valid_kwargs.get(attr))

        # TODO: For PyFITS 3.3 try to eliminate the following two special cases
        # for recformat and dim:
        # This is not actually stored as an attribute on columns for some
        # reason
        recformat = valid_kwargs['recformat']

        # The 'dim' keyword's original value is stored in self.dim, while
        # *only* the tuple form is stored in self._dims.
        self._dims = self.dim
        self.dim = dim

        # Zero-length formats are legal in the FITS format, but since they
        # are not supported by numpy we mark columns that use them as
        # "phantom" columns, that are not considered when reading the data
        # as a record array.
        if self.format[0] == '0' or \
           (self.format[-1] == '0' and self.format[-2].isalpha()):
            self._phantom = True
            array = None
        else:
            self._phantom = False

        # Awful hack to use for now to keep track of whether the column holds
        # pseudo-unsigned int data
        self._pseudo_unsigned_ints = False

        # if the column data is not ndarray, make it to be one, i.e.
        # input arrays can be just list or tuple, not required to be ndarray
        # does not include Object array because there is no guarantee
        # the elements in the object array are consistent.
        if not isinstance(array,
                          (np.ndarray, chararray.chararray, Delayed)):
            try:  # try to convert to a ndarray first
                if array is not None:
                    array = np.array(array)
            except:
                try:  # then try to convert it to a strings array
                    itemsize = int(recformat[1:])
                    array = chararray.array(array, itemsize=itemsize)
                except ValueError:
                    # then try variable length array
                    # Note: This includes _FormatQ by inheritance
                    if isinstance(recformat, _FormatP):
                        array = _VLF(array, dtype=recformat.dtype)
                    else:
                        raise ValueError('Data is inconsistent with the '
                                         'format `%s`.' % format)

        array = self._convert_to_valid_data_type(array)

        # We have required (through documentation) that arrays passed in to
        # this constructor are already in their physical values, so we make
        # note of that here
        if isinstance(array, np.ndarray):
            self._physical_values = True
        else:
            self._physical_values = False

        self.array = array

    def __repr__(self):
        text = ''
        for attr in KEYWORD_ATTRIBUTES:
            value = getattr(self, attr)
            if value is not None:
                text += attr + ' = ' + repr(value) + '; '
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

    @lazyproperty
    def dtype(self):
        return np.dtype(_convert_format(self.format))

    def copy(self):
        """
        Return a copy of this `Column`.
        """
        tmp = Column(format='I')  # just use a throw-away format
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
            with ignored(VerifyError):
                # legit recarray format?
                recformat = format
                format = cls.from_recformat(format)

        try:
            # legit FITS format?
            format = cls(format)
            recformat = format.recformat
        except VerifyError:
            raise VerifyError('Illegal format `%s`.' % format)

        return format, recformat

    @classmethod
    def _verify_keywords(cls, name=None, format=None, unit=None, null=None,
                         bscale=None, bzero=None, disp=None, start=None,
                         dim=None, ascii=None):
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

        format, recformat = cls._determine_formats(format, start, dim, ascii)
        valid.update(format=format, recformat=recformat)

        # Currently we don't have any validation for name, unit, bscale, or
        # bzero so include those by default
        # TODO: Add validation for these keywords, obviously
        for k, v in [('name', name), ('unit', unit), ('bscale', bscale),
                     ('bzero', bzero)]:
            if v is not None and v != '':
                valid[k] = v

        # Validate null option
        # Note: Enough code exists that thinks empty strings are sensible
        # inputs for these options that we need to treat '' as None
        if null is not None and null != '':
            msg = None
            if isinstance(format, _AsciiColumnFormat):
                null = str(null)
                if len(null) > format.width:
                    msg = (
                        "ASCII table null option (TNULLn) is longer than "
                        "the column's character width and will be truncated "
                        "(got %r)." % null)
            else:
                if not _is_int(null):
                    # Make this an exception instead of a warning, since any
                    # non-int value is meaningless
                    msg = (
                        'Column null option (TNULLn) must be an integer for '
                        'binary table columns (got %r).  The invalid value '
                        'will be ignored for the purpose of formatting '
                        'the data in this column.' % null)

                tnull_formats = ('B', 'I', 'J', 'K')

                if not (format.format in tnull_formats or
                        (format.format in ('P', 'Q') and
                         format.p_format in tnull_formats)):
                    # TODO: We should also check that TNULLn's integer value
                    # is in the range allowed by the column's format
                    msg = (
                        'Column null option (TNULLn) is invalid for binary '
                        'table columns of type %r (got %r).  The invalid '
                        'value will be ignored for the purpose of formatting '
                        'the data in this column.' % (format, null))

            if msg is None:
                valid['null'] = null
            else:
                invalid['null'] = (null, msg)

        # Validate the disp option
        # TODO: Add full parsing and validation of TDISPn keywords
        if disp is not None and disp != '':
            msg = None
            if not isinstance(disp, string_types):
                msg = (
                    'Column disp option (TDISPn) must be a string (got %r).'
                    'The invalid value will be ignored for the purpose of '
                    'formatting the data in this column.' % disp)

            if (isinstance(format, _AsciiColumnFormat) and
                    disp[0].upper() == 'L'):
                # disp is at least one character long and has the 'L' format
                # which is not recognized for ASCII tables
                msg = (
                    "Column disp option (TDISPn) may not use the 'L' format "
                    "with ASCII table columns.  The invalid value will be "
                    "ignored for the purpose of formatting the data in this "
                    "column.")

            if msg is None:
                valid['disp'] = disp
            else:
                invalid['disp'] = (disp, msg)

        # Validate the start option
        if start is not None and start != '':
            msg = None
            if not isinstance(format, _AsciiColumnFormat):
                # The 'start' option only applies to ASCII columns
                msg = (
                    'Column start option (TBCOLn) is not allowed for binary '
                    'table columns (got %r).  The invalid keyword will be '
                    'ignored for the purpose of formatting the data in this '
                    'column.'% start)
            try:
                start = int(start)
            except (TypeError, ValueError):
                pass

            if not _is_int(start) and start < 1:
                msg = (
                    'Column start option (TBCOLn) must be a positive integer '
                    '(got %r).  The invalid value will be ignored for the '
                    'purpose of formatting the data in this column.' % start)

            if msg is None:
                valid['start'] = start
            else:
                invalid['start'] = (start, msg)

        # Process TDIMn options
        # ASCII table columns can't have a TDIMn keyword associated with it;
        # for now we just issue a warning and ignore it.
        # TODO: This should be checked by the FITS verification code
        if dim is not None and dim != '':
            msg = None
            dims_tuple = tuple()
            # NOTE: If valid, the dim keyword's value in the the valid dict is
            # a tuple, not the original string; if invalid just the original
            # string is returned
            if isinstance(format, _AsciiColumnFormat):
                msg = (
                    'Column dim option (TDIMn) is not allowed for ASCII table '
                    'columns (got %r).  The invalid keyword will be ignored '
                    'for the purpose of formatting this column.' % dim)

            elif isinstance(dim, string_types):
                dims_tuple = _parse_tdim(dim)
            elif isinstance(dim, tuple):
                dims_tuple = dim
            else:
                msg = (
                    "`dim` argument must be a string containing a valid value "
                    "for the TDIMn header keyword associated with this column, "
                    "or a tuple containing the C-order dimensions for the "
                    "column.  The invalid value will be ignored for the purpose "
                    "of formatting this column.")

            if dims_tuple:
                if reduce(operator.mul, dims_tuple) > format.repeat:
                    msg = (
                        "The repeat count of the column format %r for column %r "
                        "is fewer than the number of elements per the TDIM "
                        "argument %r.  The invalid TDIMn value will be ignored "
                        "for the purpose of formatting this column." %
                        (name, format, dim))

            if msg is None:
                valid['dim'] = dims_tuple
            else:
                invalid['dim'] = (dim, msg)

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

        # If the given format string is unabiguously a Numpy dtype or one of
        # the Numpy record format type specifiers supported by PyFITS then that
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
            format, recformat = cls._convert_format(format,
                                                    _AsciiColumnFormat)
        else:
            # The format is already acceptable and unambiguous
            recformat = format.recformat

        return format, recformat

    @classmethod
    def _guess_format(cls, format, start, dim):
        if start and dim:
            # This is impossible; this can't be a valid FITS column
            raise ValueError(
                'Columns cannot have both a start (TCOLn) and dim '
                '(TDIMn) option, since the former is only applies to '
                'ASCII tables, and the latter is only valid for binary '
                'tables.')
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
            with ignored(VerifyError):
                format = _AsciiColumnFormat(format, strict=True)

            # A safe guess which reflects the existing behavior of previous
            # PyFITS versions
            guess_format = _ColumnFormat

        try:
            format, recformat = cls._convert_format(format, guess_format)
        except VerifyError:
            # For whatever reason our guess was wrong (for example if we got
            # just 'F' that's not a valid binary format, but it an ASCII format
            # code albeit with the width/precision ommitted
            guess_format = (_AsciiColumnFormat
                            if guess_format is _ColumnFormat
                            else _ColumnFormat)
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
            if 'P' in format or 'Q' in format:
                return array
            elif 'A' in format:
                if array.dtype.char in 'SU':
                    if dims:
                        # The 'last' dimension (first in the order given
                        # in the TDIMn keyword itself) is the number of
                        # characters in each string
                        fsize = dims[-1]
                    else:
                        fsize = np.dtype(format.recformat).itemsize
                    return chararray.array(array, itemsize=fsize)
                else:
                    return _convert_array(array, np.dtype(format.recformat))
            elif 'L' in format:
                # boolean needs to be scaled back to storage values ('T', 'F')
                if array.dtype == np.dtype('bool'):
                    return np.where(array == False, ord('F'), ord('T'))
                else:
                    return np.where(array == 0, ord('F'), ord('T'))
            elif 'X' in format:
                return _convert_array(array, np.dtype('uint8'))
            else:
                # Preserve byte order of the original array for now; see #77
                # TODO: For some reason we drop the format repeat here; need
                # to investigate why that was and if it's something we can
                # avoid doing...
                new_format = _convert_format(format.format)
                numpy_format = array.dtype.byteorder + new_format

                # Handle arrays passed in as unsigned ints as pseudo-unsigned
                # int arrays; blatantly tacked in here for now--we need columns
                # to have explicit knowledge of whether they treated as
                # pseudo-unsigned
                bzeros = {2: np.uint16(2**15), 4: np.uint32(2**31),
                          8: np.uint64(2**63)}
                if (array.dtype.kind == 'u' and
                        array.dtype.itemsize in bzeros and
                        self.bscale in (1, None, '') and
                        self.bzero == bzeros[array.dtype.itemsize]):
                    # Basically the array is uint, has scale == 1.0, and the
                    # bzero is the appropriate value for a pseudo-unsigned
                    # integer of the input dtype, then go ahead and assume that
                    # uint is assumed
                    numpy_format = numpy_format.replace('i', 'u')
                    self._pseudo_unsigned_ints = True

                return _convert_array(array, np.dtype(numpy_format))


class ColDefs(object):
    """
    Column definitions class.

    It has attributes corresponding to the `Column` attributes
    (e.g. `ColDefs` has the attribute ``names`` while `Column`
    has ``name``). Each attribute in `ColDefs` is a list of
    corresponding attribute values from all `Column` objects.
    """

    _padding_byte = '\x00'
    _col_format_cls = _ColumnFormat

    def __new__(cls, input, tbtype=None, ascii=False):
        if tbtype is not None:
            warnings.warn(
                'The ``tbtype`` argument to `ColDefs` is deprecated as of '
                'Astropy 0.4; instead the appropriate table type should be '
                'inferred from the formats of the supplied columns.  Use the '
                '``ascii=True`` argument to ensure that ASCII table columns '
                'are used.', AstropyDeprecationWarning)
        else:
            tbtype = 'BinTableHDU'  # The old default

        # Backards-compat support
        # TODO: Remove once the tbtype argument is removed entirely
        if tbtype == 'BinTableHDU':
            klass = cls
        elif tbtype == 'TableHDU':
            klass = _AsciiColDefs
        else:
            raise ValueError('Invalid table type: %s.' % tbtype)

        if (hasattr(input, '_columns_type') and
                issubclass(input._columns_type, ColDefs)):
            klass = input._columns_type
        elif (hasattr(input, '_col_format_cls') and
                issubclass(input._col_format_cls, _AsciiColumnFormat)):
            klass = _AsciiColDefs

        if ascii:  # force ASCII if this has been explicitly requested
            klass = _AsciiColDefs

        return object.__new__(klass)

    def __getnewargs__(self):
        return (self._arrays,)

    def __init__(self, input, tbtype=None, ascii=False):
        """
        Parameters
        ----------

        input : sequence of `Column`, `ColDefs`, other
            An existing table HDU, an existing `ColDefs`, or any multi-field
            Numpy array or `numpy.recarray`.

        **(Deprecated)** tbtype : str, optional
            which table HDU, ``"BinTableHDU"`` (default) or
            ``"TableHDU"`` (text table).
            Now ColDefs for a normal (binary) table by default, but converted
            automatically to ASCII table ColDefs in the appropriate contexts
            (namely, when creating an ASCII table).

        ascii : bool
        """

        from .hdu.table import _TableBaseHDU
        from .fitsrec import FITS_rec

        if isinstance(input, ColDefs):
            self._init_from_coldefs(input)
        elif (isinstance(input, FITS_rec) and hasattr(input, '_coldefs') and
                input._coldefs):
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
            raise TypeError('Input to ColDefs must be a table HDU, a list '
                            'of Columns, or a record/field array.')

    def _init_from_coldefs(self, coldefs):
        """Initialize from an existing ColDefs object (just copy the
        columns and convert their formats if necessary).
        """

        self.columns = [self._copy_column(col) for col in coldefs]

    def _init_from_sequence(self, columns):
        for idx, col in enumerate(columns):
            if not isinstance(col, Column):
                raise TypeError(
                    'Element %d in the ColDefs input is not a Column.' % idx)

        self._init_from_coldefs(columns)

    def _init_from_array(self, array):
        self.columns = []
        for idx in range(len(array.dtype)):
            cname = array.dtype.names[idx]
            ftype = array.dtype.fields[cname][0]
            format = self._col_format_cls.from_recformat(ftype)

            # Determine the appropriate dimensions for items in the column
            # (typically just 1D)
            dim = array.dtype[idx].shape[::-1]
            if dim and (len(dim) > 1 or 'A' in format):
                if 'A' in format:
                    # n x m string arrays must include the max string
                    # length in their dimensions (e.g. l x n x m)
                    dim = (array.dtype[idx].base.itemsize,) + dim
                dim = repr(dim).replace(' ', '')
            else:
                dim = None

            # Check for unsigned ints.
            bzero = None
            if 'I' in format and ftype == np.dtype('uint16'):
                bzero = np.uint16(2**15)
            elif 'J' in format and ftype == np.dtype('uint32'):
                bzero = np.uint32(2**31)
            elif 'K' in format and ftype == np.dtype('uint64'):
                bzero = np.uint64(2**63)

            c = Column(name=cname, format=format,
                       array=array.view(np.ndarray)[cname], bzero=bzero,
                       dim=dim)
            self.columns.append(c)

    def _init_from_table(self, table):
        hdr = table._header
        nfields = hdr['TFIELDS']
        self._width = hdr['NAXIS1']
        self._shape = hdr['NAXIS2']

        # go through header keywords to pick out column definition keywords
        # definition dictionaries for each field
        col_keywords = [{} for i in range(nfields)]
        for keyword, value in iteritems(hdr):
            key = TDEF_RE.match(keyword)
            try:
                keyword = key.group('label')
            except:
                continue  # skip if there is no match
            if keyword in KEYWORD_NAMES:
                col = int(key.group('num'))
                if col <= nfields and col > 0:
                    idx = KEYWORD_NAMES.index(keyword)
                    attr = KEYWORD_ATTRIBUTES[idx]
                    if attr == 'format':
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
                    'Invalid keyword for column %d: %s' % (idx + 1, val[1]),
                    VerifyWarning)
            # Special cases for recformat and dim
            # TODO: Try to eliminate the need for these special cases
            del valid_kwargs['recformat']
            if 'dim' in valid_kwargs:
                valid_kwargs['dim'] = kwargs['dim']
            col_keywords[idx] = valid_kwargs

        # data reading will be delayed
        for col in range(nfields):
            col_keywords[col]['array'] = Delayed(table, col)

        # now build the columns
        self.columns = [Column(**attrs) for attrs in col_keywords]
        self._listener = weakref.proxy(table)

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
        new_column.format = self._col_format_cls.from_column_format(
                column.format)

        # Handle a few special cases of column format options that are not
        # compatible between ASCII an binary tables
        # TODO: This is sort of hacked in right now; we really neet
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
            if (new_column.disp is not None and
                    new_column.disp.upper().startswith('L')):
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
        if cname in KEYWORD_ATTRIBUTES and name[-1] == 's':
            attr = []
            for col in self:
                val = getattr(col, cname)
                if val is not None:
                    attr.append(val)
                else:
                    attr.append('')
            return attr
        raise AttributeError(name)

    @lazyproperty
    def dtype(self):
        recformats = [f for idx, f in enumerate(self._recformats)
                      if not self[idx]._phantom]
        formats = ','.join(recformats)
        names = [n for idx, n in enumerate(self.names)
                 if not self[idx]._phantom]
        return np.rec.format_parser(formats, names, None).dtype

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
        x = self.columns[key]
        if _is_int(key):
            return x
        else:
            return ColDefs(x)

    def __len__(self):
        return len(self.columns)

    def __repr__(self):
        rep = 'ColDefs('
        if hasattr(self, 'columns') and self.columns:
            # The hasattr check is mostly just useful in debugging sessions
            # where self.columns may not be defined yet
            rep += '\n    '
            rep += '\n    '.join([repr(c) for c in self.columns])
            rep += '\n'
        rep += ')'
        return rep

    def __add__(self, other, option='left'):
        if isinstance(other, Column):
            b = [other]
        elif isinstance(other, ColDefs):
            b = list(other.columns)
        else:
            raise TypeError('Wrong type of input.')
        if option == 'left':
            tmp = list(self.columns) + b
        else:
            tmp = b + list(self.columns)
        return ColDefs(tmp)

    def __radd__(self, other):
        return self.__add__(other, 'right')

    def __sub__(self, other):
        if not isinstance(other, (list, tuple)):
            other = [other]
        _other = [_get_index(self.names, key) for key in other]
        indx = list(range(len(self)))
        for x in _other:
            indx.remove(x)
        tmp = [self[i] for i in indx]
        return ColDefs(tmp)

    def _update_listener(self):
        if hasattr(self, '_listener'):
            try:
                if self._listener._data_loaded:
                    del self._listener.data
                self._listener.columns = self
            except ReferenceError:
                del self._listener

    def add_col(self, column):
        """
        Append one `Column` to the column definition.
        """

        assert isinstance(column, Column)

        for cname in KEYWORD_ATTRIBUTES:
            attr = getattr(self, cname + 's')
            attr.append(getattr(column, cname))

        self._arrays.append(column.array)
        # Obliterate caches of certain things
        del self.dtype
        del self._recformats
        del self._dims

        self.columns.append(column)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()
        return self

    def del_col(self, col_name):
        """
        Delete (the definition of) one `Column`.

        col_name : str or int
            The column's name or index
        """

        indx = _get_index(self.names, col_name)

        for cname in KEYWORD_ATTRIBUTES:
            attr = getattr(self, cname + 's')
            del attr[indx]

        del self._arrays[indx]
        # Obliterate caches of certain things
        del self.dtype
        del self._recformats
        del self._dims

        del self.columns[indx]

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()
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

        indx = _get_index(self.names, col_name)
        getattr(self, attrib + 's')[indx] = new_value

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

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
            raise ValueError('New name %s already exists.' % new_name)
        else:
            self.change_attrib(col_name, 'name', new_name)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

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

        self.change_attrib(col_name, 'unit', new_unit)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

    def info(self, attrib='all', output=None):
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

        output : file, optional
            File-like object to output to.  Outputs to stdout by default.
            If `False`, returns the attributes as a `dict` instead.

        Notes
        -----
        This function doesn't return anything by default; it just prints to
        stdout.
        """

        if output is None:
            output = sys.stdout

        if attrib.strip().lower() in ['all', '']:
            lst = KEYWORD_ATTRIBUTES
        else:
            lst = attrib.split(',')
            for idx in range(len(lst)):
                lst[idx] = lst[idx].strip().lower()
                if lst[idx][-1] == 's':
                    lst[idx] = list[idx][:-1]

        ret = {}

        for attr in lst:
            if output:
                if attr not in KEYWORD_ATTRIBUTES:
                    output.write("'%s' is not an attribute of the column "
                                 "definitions.\n" % attr)
                    continue
                output.write("%s:\n" % attr)
                output.write('    %s\n' % getattr(self, attr + 's'))
            else:
                ret[attr] = getattr(self, attr + 's')

        if not output:
            return ret


class _AsciiColDefs(ColDefs):
    """ColDefs implementation for ASCII tables."""

    _padding_byte = ' '
    _col_format_cls = _AsciiColumnFormat

    def __init__(self, input, tbtype=None, ascii=True):
        super(_AsciiColDefs, self).__init__(input)

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
        _itemsize = self.spans[-1] + self.starts[-1] - 1
        dtype = {}

        for j in range(len(self)):
            data_type = 'S' + str(self.spans[j])
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
        return ['a' + str(w) for w in widths]

    def add_col(self, column):
        super(_AsciiColDefs, self).add_col(column)
        self._update_field_metrics()

    def del_col(self, col_name):
        super(_AsciiColDefs, self).del_col(col_name)
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


class _VLF(np.ndarray):
    """Variable length field object."""

    def __new__(cls, input, dtype='a'):
        """
        Parameters
        ----------
        input
            a sequence of variable-sized elements.
        """

        if dtype == 'a':
            try:
                # this handles ['abc'] and [['a','b','c']]
                # equally, beautiful!
                input = [chararray.array(x, itemsize=1) for x in input]
            except:
                raise ValueError(
                    'Inconsistent input data array: {0}'.format(input))

        a = np.array(input, dtype=np.object)
        self = np.ndarray.__new__(cls, shape=(len(input),), buffer=a,
                                  dtype=np.object)
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
        elif self.element_dtype == 'a':
            value = chararray.array(value, itemsize=1)
        else:
            value = np.array(value, dtype=self.element_dtype)
        np.ndarray.__setitem__(self, key, value)
        self.max = max(self.max, len(value))


def _get_index(names, key):
    """
    Get the index of the `key` in the `names` list.

    The `key` can be an integer or string.  If integer, it is the index
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
    elif isinstance(key, string_types):
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
                raise KeyError("Key '%s' does not exist." % key)
            else:              # multiple match
                raise KeyError("Ambiguous key name '%s'." % key)
    else:
        raise KeyError("Illegal key '%s'." % repr(key))

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

    pow2 = np.array([128, 64, 32, 16, 8, 4, 2, 1], dtype='uint8')
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
    n = min(len(array), nrows)

    data_output = _VLF([None] * nrows, dtype=format.dtype)

    if format.dtype == 'a':
        _nbytes = 1
    else:
        _nbytes = np.array([], dtype=format.dtype).itemsize

    for idx in range(nrows):
        if idx < len(array):
            rowval = array[idx]
        else:
            if format.dtype == 'a':
                rowval = ' ' * data_output.max
            else:
                rowval = [0] * data_output.max
        if format.dtype == 'a':
            data_output[idx] = chararray.array(encode_ascii(rowval),
                                               itemsize=1)
        else:
            data_output[idx] = np.array(rowval, dtype=format.dtype)

        descr_output[idx, 0] = len(data_output[idx])
        descr_output[idx, 1] = _offset
        _offset += len(data_output[idx]) * _nbytes

    return data_output


def _parse_tformat(tform):
    """Parse ``TFORMn`` keyword for a binary table into a
    ``(repeat, format, option)`` tuple.
    """

    try:
        (repeat, format, option) = TFORMAT_RE.match(tform.strip()).groups()
    except:
        # TODO: Maybe catch this error use a default type (bytes, maybe?) for
        # unrecognized column types.  As long as we can determine the correct
        # byte width somehow..
        raise VerifyError('Format %r is not recognized.' % tform)

    if repeat == '':
        repeat = 1
    else:
        repeat = int(repeat)

    return (repeat, format.upper(), option)


def _parse_ascii_tformat(tform, strict=False):
    """Parse the ``TFORMn`` keywords for ASCII tables into a
    ``(format, width, precision)`` tuple (the latter is zero unless
    width is one of 'E', 'F', or 'D').
    """

    match = TFORMAT_ASCII_RE.match(tform.strip())
    if not match:
        raise VerifyError('Format %r is not recognized.' % tform)

    # Be flexible on case
    format = match.group('format')
    if format is None:
        # Floating point format
        format = match.group('formatf').upper()
        width = match.group('widthf')
        precision = match.group('precision')
        if width is None or precision is None:
            if strict:
                raise VerifyError('Format %r is not unambiguously an ASCII '
                                  'table format.')
            else:
                width = 0 if width is None else width
                precision = 1 if precision is None else precision
    else:
        format = format.upper()
        width = match.group('width')
        if width is None:
            if strict:
                raise VerifyError('Format %r is not unambiguously an ASCII '
                                  'table format.')
            else:
                # Just use a default width of 0 if unspecified
                width = 0
        precision = 0

    def convert_int(val):
        msg = ('Format %r is not valid--field width and decimal precision '
               'must be positive integers.')
        try:
            val = int(val)
        except (ValueError, TypeError):
            raise VerifyError(msg % tform)

        if val <= 0:
            raise VerifyError(msg % tform)

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

    if precision >= width:
        raise VerifyError("Format %r not valid--the number of decimal digits "
                          "must be less than the format's total width %s." &
                          (tform, width))

    return format, width, precision


def _parse_tdim(tdim):
    """Parse the ``TDIM`` value into a tuple (may return an empty tuple if
    the value ``TDIM`` value is empty or invalid).
    """

    m = tdim and TDIM_RE.match(tdim)
    if m:
        dims = m.group('dims')
        return tuple(int(d.strip()) for d in dims.split(','))[::-1]

    # Ignore any dim values that don't specify a multidimensional column
    return tuple()


def _scalar_to_format(value):
    """
    Given a scalar value or string, returns the minimum FITS column format
    that can represent that value.  'minimum' is defined by the order given in
    FORMATORDER.
    """

    # TODO: Numpy 1.6 and up has a min_scalar_type() function that can handle
    # this; in the meantime we have to use our own implementation (which for
    # now is pretty naive)

    # First, if value is a string, try to convert to the appropriate scalar
    # value
    for type_ in (int, float, complex):
        try:
            value = type_(value)
            break
        except ValueError:
            continue

    if isinstance(value, int) and value in (0, 1):
        # Could be a boolean
        return 'L'
    elif isinstance(value, int):
        for char in ('B', 'I', 'J', 'K'):
            type_ = np.dtype(FITS2NUMPY[char]).type
            if type_(value) == value:
                return char
    elif isinstance(value, float):
        # For now just assume double precision
        return 'D'
    elif isinstance(value, complex):
        return 'M'
    else:
        return 'A' + str(len(value))


def _cmp_recformats(f1, f2):
    """
    Compares two numpy recformats using the ordering given by FORMATORDER.
    """

    if f1[0] == 'a' and f2[0] == 'a':
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
        if dtype == 'A':
            output_format = FITS2NUMPY[dtype] + str(repeat)
            # to accomodate both the ASCII table and binary table column
            # format spec, i.e. A7 in ASCII table is the same as 7A in
            # binary table, so both will produce 'a7'.
            # Technically the FITS standard does not allow this but it's a very
            # common mistake
            if format.lstrip()[0] == 'A' and option != '':
                # make sure option is integer
                output_format = FITS2NUMPY[dtype] + str(int(option))
        else:
            repeat_str = ''
            if repeat != 1:
                repeat_str = str(repeat)
            output_format = repeat_str + FITS2NUMPY[dtype]

    elif dtype == 'X':
        output_format = _FormatX(repeat)
    elif dtype == 'P':
        output_format = _FormatP.from_tform(format)
    elif dtype == 'Q':
        output_format = _FormatQ.from_tform(format)
    elif dtype == 'F':
        output_format = 'f8'
    else:
        raise ValueError('Illegal format %s.' % format)

    return output_format


def _convert_record2fits(format):
    """
    Convert record format spec to FITS format spec.
    """

    recformat, kind, dtype = _dtype_to_recformat(format)
    shape = dtype.shape
    option = str(dtype.base.itemsize)

    ndims = len(shape)
    repeat = 1
    if ndims > 0:
        nel = np.array(shape, dtype='i8').prod()
        if nel > 1:
            repeat = nel

    if kind == 'a':
        # This is a kludge that will place string arrays into a
        # single field, so at least we won't lose data.  Need to
        # use a TDIM keyword to fix this, declaring as (slength,
        # dim1, dim2, ...)  as mwrfits does

        ntot = int(repeat) * int(option)

        output_format = str(ntot) + 'A'
    elif recformat in NUMPY2FITS:  # record format
        if repeat != 1:
            repeat = str(repeat)
        else:
            repeat = ''
        output_format = repeat + NUMPY2FITS[recformat]
    else:
        raise ValueError('Illegal format %s.' % format)

    return output_format


def _dtype_to_recformat(dtype):
    """
    Utility function for converting a dtype object or string that instantiates
    a dtype (e.g. 'float32') into one of the two character Numpy format codes
    that have been traditionally used by PyFITS.

    In particular, use of 'a' to refer to character data is long since
    deprecated in Numpy, but PyFITS remains heavily invested in its use
    (something to try to get away from sooner rather than later).
    """

    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    kind = dtype.base.kind
    itemsize = dtype.base.itemsize
    recformat = kind + str(itemsize)

    if kind in ('U', 'S'):
        recformat = kind = 'a'

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

        if kind == 'a':
            return 'A' + str(itemsize)
        elif NUMPY2FITS.get(recformat) == 'L':
            # Special case for logical/boolean types--for ASCII tables we
            # represent these as single character columns containing 'T' or 'F'
            # (a la the storage format for Logical columns in binary tables)
            return 'A1'
        elif kind == 'i':
            # Use for the width the maximum required to represent integers
            # of that byte size plus 1 for signs, but use a minumum of the
            # default width (to keep with existing behavior)
            width = 1 + len(str(2 ** (itemsize * 8)))
            width = max(width, ASCII_DEFAULT_WIDTHS['I'][0])
            return 'I' + str(width)
        elif kind == 'f':
            # This is tricky, but go ahead and use D if float-64, and E
            # if float-32 with their default widths
            if itemsize >= 8:
                format = 'D'
            else:
                format = 'E'
            width = '.'.join(str(w) for w in ASCII_DEFAULT_WIDTHS[format])
            return format + width
        # TODO: There may be reasonable ways to represent other Numpy types so
        # let's see what other possibilities there are besides just 'a', 'i',
        # and 'f'.  If it doesn't have a reasonable ASCII representation then
        # raise an exception
    else:
        format, width, precision = _parse_ascii_tformat(format)

        # This gives a sensible "default" dtype for a given ASCII
        # format code
        recformat = ASCII2NUMPY[format]

        # The following logic is taken from CFITSIO:
        # For integers, if the width <= 4 we can safely use 16-bit ints for all
        # values [for the non-standard J format code just always force 64-bit]
        if format == 'I' and width <= 4:
            recformat = 'i2'
        elif format == 'F' and width > 7:
            # 32-bit floats (the default) may not be accurate enough to support
            # all values that can fit in this field, so upgrade to 64-bit
            recformat = 'f8'
        elif format == 'E' and precision > 6:
            # Again upgrade to a 64-bit int if we require greater decimal
            # precision
            recformat = 'f8'
        elif format == 'A':
            recformat += str(width)

        return recformat
