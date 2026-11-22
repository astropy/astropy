# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module handles the conversion of various VOTABLE datatypes
to/from TABLEDATA_ and BINARY_ formats.
"""

# STDLIB
import re
import struct
import sys
from itertools import repeat
from math import prod

# THIRD-PARTY
import numpy as np
from numpy import ma

# ASTROPY
from astropy.utils.xml.writer import xml_escape_cdata

# LOCAL
from .exceptions import (
    E01,
    E02,
    E03,
    E04,
    E05,
    E06,
    E24,
    W01,
    W30,
    W31,
    W39,
    W46,
    W47,
    W49,
    W51,
    W55,
    W60,
    W61,
    W62,
    vo_raise,
    vo_warn,
    warn_or_raise,
)
from .fast_converters import (
    fast_binparse_bit,
    fast_binparse_bool,
    fast_binparse_double,
    fast_binparse_float,
    fast_binparse_int,
    fast_binparse_long,
    fast_binparse_short,
    fast_binparse_ubyte,
)
from .util import version_compare

__all__ = ["Converter", "get_converter", "table_column_to_votable_datatype"]


pedantic_array_splitter = re.compile(r" +")
array_splitter = re.compile(r"\s+|(?:\s*,\s*)")
"""
A regex to handle splitting values on either whitespace or commas.

SPEC: Usage of commas is not actually allowed by the spec, but many
files in the wild use them.
"""

_zero_int = b"\0\0\0\0"
_empty_bytes = b""
_zero_byte = b"\0"


if sys.byteorder == "little":

    def _ensure_bigendian(x):
        if x.dtype.byteorder != ">":
            return x.byteswap()
        return x

else:

    def _ensure_bigendian(x):
        if x.dtype.byteorder == "<":
            return x.byteswap()
        return x


def _make_masked_array(data, mask):
    """
    Masked arrays of zero length that also have a mask of zero length
    cause problems in Numpy (at least in 1.6.2).  This function
    creates a masked array from data and a mask, unless it is zero
    length.
    """
    # np.ma doesn't like setting mask to []
    if len(data):
        return ma.array(np.array(data), mask=np.array(mask, dtype="bool"))
    else:
        return ma.array(np.array(data))


def bitarray_to_bool(data, length):
    """
    Converts a bit array (a string of bits in a bytes object) to a
    boolean Numpy array.

    Parameters
    ----------
    data : bytes
        The bit array.  The most significant byte is read first.

    length : int
        The number of bits to read.  The least significant bits in the
        data bytes beyond length will be ignored.

    Returns
    -------
    array : numpy bool array
    """
    return np.unpackbits(
        np.frombuffer(data, dtype=np.uint8), count=length if length >= 0 else None
    ).astype(bool, copy=False)


def bool_to_bitarray(value):
    """
    Converts a numpy boolean array to a bit array (a string of bits in
    a bytes object).

    Parameters
    ----------
    value : numpy bool array

    Returns
    -------
    bit_array : bytes
        The first value in the input array will be the most
        significant bit in the result.  The length will be `floor((N +
        7) / 8)` where `N` is the length of `value`.
    """
    bytes = np.packbits(np.ma.filled(value, fill_value=False).flat).tobytes()
    return struct.pack(f"{len(bytes)}B", *bytes)


class Converter:
    """
    The base class for all converters.  Each subclass handles
    converting a specific VOTABLE data type to/from the TABLEDATA_ and
    BINARY_ on-disk representations.

    Parameters
    ----------
    field : `~astropy.io.votable.tree.Field`
        object describing the datatype

    config : dict
        The parser configuration dictionary

    pos : tuple
        The position in the XML file where the FIELD object was
        found.  Used for error messages.

    """

    def __init__(self, field, config=None, pos=None):
        pass

    @staticmethod
    def _parse_length(read):
        return struct.unpack(">I", read(4))[0]

    @staticmethod
    def _write_length(length):
        return struct.pack(">I", int(length))

    def supports_empty_values(self, config):
        """
        Returns True when the field can be completely empty.
        """
        return config.get("version_1_3_or_later")

    def parse(self, value, config=None, pos=None):
        """
        Convert the string *value* from the TABLEDATA_ format into an
        object with the correct native in-memory datatype and mask flag.

        Parameters
        ----------
        value : str
            value in TABLEDATA format
        config : dict, optional
            The parser configuration dictionary
        pos : tuple, optional
            The position in the XML file where the value was found.  Used
            for error messages.

        Returns
        -------
        native : tuple
            A two-element tuple of: value, mask.
            The value as a Numpy array or scalar, and *mask* is True
            if the value is missing.
        """
        raise NotImplementedError("This datatype must implement a 'parse' method.")

    def parse_scalar(self, value, config=None, pos=None):
        """
        Parse a single scalar of the underlying type of the converter.
        For non-array converters, this is equivalent to parse.  For
        array converters, this is used to parse a single
        element of the array.

        Parameters
        ----------
        value : str
            value in TABLEDATA format
        config : dict, optional
            The parser configuration dictionary
        pos : tuple, optional
            The position in the XML file where the value was found.  Used
            for error messages.

        Returns
        -------
        native : (2,) tuple
            (value, mask)
            The value as a Numpy array or scalar, and *mask* is True
            if the value is missing.
        """
        return self.parse(value, config, pos)

    def output(self, value, mask):
        """
        Convert the object *value* (in the native in-memory datatype)
        to a unicode string suitable for serializing in the TABLEDATA_
        format.

        Parameters
        ----------
        value
            The value, the native type corresponding to this converter

        mask : bool
            If `True`, will return the string representation of a
            masked value.

        Returns
        -------
        tabledata_repr : unicode
        """
        raise NotImplementedError("This datatype must implement a 'output' method.")

    def binparse(self, read):
        """
        Reads some number of bytes from the BINARY_ format
        representation by calling the function *read*, and returns the
        native in-memory object representation for the datatype
        handled by *self*.

        Parameters
        ----------
        read : function
            A function that given a number of bytes, returns a byte
            string.

        Returns
        -------
        native : (2,) tuple
            (value, mask). The value as a Numpy array or scalar, and *mask* is
            True if the value is missing.
        """
        raise NotImplementedError("This datatype must implement a 'binparse' method.")

    def binoutput(self, value, mask):
        """
        Convert the object *value* in the native in-memory datatype to
        a string of bytes suitable for serialization in the BINARY_
        format.

        Parameters
        ----------
        value
            The value, the native type corresponding to this converter

        mask : bool
            If `True`, will return the string representation of a
            masked value.

        Returns
        -------
        bytes : bytes
            The binary representation of the value, suitable for
            serialization in the BINARY_ format.
        """
        raise NotImplementedError("This datatype must implement a 'binoutput' method.")


def _string_length(arraysize):
    """
    The string length part of a ``char`` or ``unicodeChar`` arraysize: "10"
    for "10x4", and "10*" for "10*".

    Only the first dimension is the length of a string.  The others give
    the shape of an array of strings, which `StringArray` and
    `StringVarArray` handle.
    """
    return arraysize.split("x")[0]


def _string_bound(arraysize):
    """
    The string length of an arraysize as ``(n, bounded)``: ``(10, False)``
    for "10x4", ``(10, True)`` for "10*", and ``(0, False)`` for "*" or None.
    """
    if arraysize is None:
        return 0, False
    length = _string_length(arraysize)
    n = length.removesuffix("*")
    if not n.isdigit():
        return 0, False
    return int(n), length.endswith("*")


def _shape_suffix(shape):
    """VOTable lists the fastest varying dimension first, numpy last."""
    return "".join(f"x{n}" for n in reversed(shape))


def _truncate_utf8(data, size):
    """Cut to at most ``size`` bytes without splitting a UTF-8 sequence."""
    if len(data) <= size:
        return data
    end = size
    while end > 0 and data[end] & 0xC0 == 0x80:
        end -= 1
    return data[:end]


def _write_config(config):
    if config.get("verify") == "ignore":
        return {**config, "verify": "warn"}
    return config


class Char(Converter):
    """
    Handles the char datatype.

    From VOTable 1.6 this is UTF-8, with arraysize counting bytes and the
    optional ``width`` giving the maximum number of characters. Earlier
    versions only allow ASCII, but binary data is read and written as UTF-8
    for every version since many older files contain it anyway.

    Missing values are not handled for string or unicode types.
    """

    default = _empty_bytes

    def __init__(self, field, config=None, pos=None):
        if config is None:
            config = {}
        self._config = config
        self._write_config = _write_config(config)

        Converter.__init__(self, field, config, pos)

        self.is_utf8_version = self._config.get("version_1_6_or_later", False)
        self.field_name = field.name

        if field.arraysize is None:
            vo_warn(W47, (), config, pos)
            field.arraysize = "1"

        self._width_limit = None

        if field.arraysize == "*":
            self.format = "O"
            self.binparse = self._binparse_var
            self.binoutput = self._binoutput_var
            self.arraysize = "*"
        else:
            # Check if this is a bounded variable-length field
            length = _string_length(field.arraysize)
            is_variable = length.endswith("*")
            numeric_part = length.removesuffix("*")
            try:
                self.arraysize = int(numeric_part)
            except ValueError:
                vo_raise(E01, (numeric_part, "char", field.ID), config)

            n_chars = self.arraysize
            if self.is_utf8_version and field.width is not None:
                # A string never has more characters than UTF-8 bytes.
                n_chars = min(field.width, self.arraysize)
                if n_chars < self.arraysize:
                    self._width_limit = n_chars
            self.format = f"U{n_chars:d}"

            # For bounded variable-length fields use the variable methods
            if is_variable:
                self.binparse = self._binparse_var
                self.binoutput = self._binoutput_var
            else:
                self.binparse = self._binparse_fixed
                self.binoutput = self._binoutput_fixed

    def supports_empty_values(self, config):
        return True

    def parse(self, value, config=None, pos=None):
        if self.arraysize != "*":
            length = len(value.encode("utf-8")) if self.is_utf8_version else len(value)
            if length > self.arraysize:
                vo_warn(W46, ("char", self.arraysize), config, pos)
        if not self.is_utf8_version and not value.isascii():
            vo_warn(W55, (self.field_name, value), config, pos)
        self._check_width(value, config, pos)
        return value, False

    def _check_width(self, value, config=None, pos=None):
        if self._width_limit is not None and len(value) > self._width_limit:
            vo_warn(W61, (self.field_name, self._width_limit, value), config, pos)

    def output(self, value, mask):
        if mask:
            return ""

        if not isinstance(value, str):
            value = value.decode("utf-8")

        if not self.is_utf8_version and not value.isascii():
            warn_or_raise(E24, E24, (value, self.field_name), self._write_config)

        if self.arraysize != "*":
            if self.is_utf8_version:
                encoded = value.encode("utf-8", errors="replace")
                if len(encoded) > self.arraysize:
                    vo_warn(W46, ("char", self.arraysize), self._write_config)
                    value = _truncate_utf8(encoded, self.arraysize).decode("utf-8")
            elif len(value) > self.arraysize:
                vo_warn(W46, ("char", self.arraysize), self._write_config)
                value = value[: self.arraysize]

        return xml_escape_cdata(value)

    def _decode(self, data):
        try:
            value = data.decode("utf-8")
        except UnicodeDecodeError:
            vo_warn(W60, (self.field_name, data), self._config)
            return data.decode("utf-8", errors="replace")
        if not self.is_utf8_version and not value.isascii():
            vo_warn(W55, (self.field_name, value), self._config)
        return value

    def _encode(self, value):
        if not self.is_utf8_version and not value.isascii():
            warn_or_raise(E24, E24, (value, self.field_name), self._write_config)
        return value.encode("utf-8", errors="replace")

    def _binparse_var(self, read):
        length = self._parse_length(read)
        if self.arraysize != "*" and length > self.arraysize:
            vo_warn(W46, ("char", self.arraysize), None, None)

        value = self._decode(read(length))
        self._check_width(value, self._config)
        return value, False

    def _binparse_fixed(self, read):
        s = read(self.arraysize)
        end = s.find(_zero_byte)
        if end != -1:
            s = s[:end]
        value = self._decode(s)
        self._check_width(value, self._config)
        return value, False

    def _binoutput_var(self, value, mask):
        if mask or value is None or value == "":
            return _zero_int
        if isinstance(value, str):
            value = self._encode(value)

        if self.arraysize != "*" and len(value) > self.arraysize:
            vo_warn(W46, ("char", self.arraysize), self._write_config)
            value = _truncate_utf8(value, self.arraysize)

        return self._write_length(len(value)) + value

    def _binoutput_fixed(self, value, mask):
        if mask or value is None:
            value_bytes = _empty_bytes
        elif isinstance(value, str):
            value_bytes = self._encode(value)
        else:
            value_bytes = value

        if len(value_bytes) > self.arraysize:
            vo_warn(W46, ("char", self.arraysize), self._write_config)
            value_bytes = _truncate_utf8(value_bytes, self.arraysize)
        return value_bytes.ljust(self.arraysize, _zero_byte)


class UnicodeChar(Converter):
    """
    Handles the unicodeChar data type. UTF-16-BE.

    The unicodeChar datatype is deprecated from VOTable 1.6, which uses UTF-8
    ``char`` instead. Files using it can still be read and written.

    Missing values are not handled for string or unicode types.
    """

    default = ""

    def __init__(self, field, config=None, pos=None):
        if config is None:
            config = {}
        self._config = config
        self._write_config = _write_config(config)
        self.field_name = field.name
        Converter.__init__(self, field, config, pos)

        if field.arraysize is None:
            vo_warn(W47, (), config, pos)
            field.arraysize = "1"

        if field.arraysize == "*":
            self.format = "O"
            self.binparse = self._binparse_var
            self.binoutput = self._binoutput_var
            self.arraysize = "*"
        else:
            # Check if this is a bounded variable-length field
            length = _string_length(field.arraysize)
            is_variable = length.endswith("*")
            numeric_part = length.removesuffix("*")
            try:
                self.arraysize = int(numeric_part)
            except ValueError:
                vo_raise(E01, (numeric_part, "unicode", field.ID), config)

            self.format = f"U{self.arraysize:d}"

            if is_variable:
                self.binparse = self._binparse_var
                self.binoutput = self._binoutput_var
            else:
                self.binparse = self._binparse_fixed
                self.binoutput = self._binoutput_fixed

            self._struct_format = f">{self.arraysize * 2:d}s"

    def parse(self, value, config=None, pos=None):
        if self.arraysize != "*" and len(value) > self.arraysize:
            vo_warn(W46, ("unicodeChar", self.arraysize), config, pos)
        return value, False

    def _check_bmp(self, value):
        if not value.isascii() and max(value) > "\uffff":
            vo_warn(W62, (self.field_name, value), self._write_config)

    def _encode(self, value, n_units):
        self._check_bmp(value)
        encoded = value.encode("utf_16_be")
        if n_units is not None and len(encoded) > 2 * n_units:
            vo_warn(W46, ("unicodeChar", n_units), self._write_config)
            encoded = encoded[: 2 * n_units]
            # Don't leave half of a surrogate pair.
            if encoded and 0xD8 <= encoded[-2] <= 0xDB:
                encoded = encoded[:-2]
        return encoded

    def output(self, value, mask):
        if mask:
            return ""
        value = str(value)
        if self.arraysize == "*":
            self._check_bmp(value)
        else:
            value = self._encode(value, self.arraysize).decode("utf_16_be")
        return xml_escape_cdata(value)

    def _binparse_var(self, read):
        length = self._parse_length(read)

        if self.arraysize != "*" and length > self.arraysize:
            vo_warn(W46, ("unicodeChar", self.arraysize), None, None)

        return read(length * 2).decode("utf_16_be"), False

    def _binparse_fixed(self, read):
        s = struct.unpack(self._struct_format, read(self.arraysize * 2))[0]
        s = s.decode("utf_16_be")
        end = s.find("\0")
        if end != -1:
            return s[:end], False
        return s, False

    def _binoutput_var(self, value, mask):
        if mask or value is None or value == "":
            return _zero_int

        n_units = None if self.arraysize == "*" else self.arraysize
        encoded = self._encode(value, n_units)

        return self._write_length(len(encoded) // 2) + encoded

    def _binoutput_fixed(self, value, mask):
        if mask:
            value = ""
        encoded = self._encode(value, self.arraysize)
        return encoded.ljust(2 * self.arraysize, _zero_byte)


class Array(Converter):
    """
    Handles both fixed and variable-lengths arrays.
    """

    def __init__(self, field, config=None, pos=None):
        if config is None:
            config = {}
        Converter.__init__(self, field, config, pos)
        if config.get("verify", "ignore") == "exception":
            self._splitter = self._splitter_pedantic
        else:
            self._splitter = self._splitter_lax

    def parse_scalar(self, value, config=None, pos=0):
        return self._base.parse_scalar(value, config, pos)

    @staticmethod
    def _splitter_pedantic(value, config=None, pos=None):
        return pedantic_array_splitter.split(value)

    @staticmethod
    def _splitter_lax(value, config=None, pos=None):
        if "," in value:
            vo_warn(W01, (), config, pos)
        return array_splitter.split(value)


class VarArray(Array):
    """
    Handles variable lengths arrays (i.e. where *arraysize* is '*').
    """

    format = "O"

    def __init__(self, field, base, arraysize, config=None, pos=None):
        Array.__init__(self, field, config)

        self._base = base
        self.default = np.array([], dtype=self._base.format)

    def output(self, value, mask):
        output = self._base.output
        # *mask* is for the cell as a whole; the mask of each element is
        # in the mask of *value*, as binoutput reads it.
        mask = np.logical_or(mask, ma.getmaskarray(value))
        result = [output(x, m) for x, m in np.broadcast(value, mask)]
        return " ".join(result)

    def binparse(self, read):
        length = self._parse_length(read)

        result = []
        result_mask = []
        binparse = self._base.binparse
        for i in range(length):
            val, mask = binparse(read)
            result.append(val)
            result_mask.append(mask)

        return _make_masked_array(result, result_mask), False

    def binoutput(self, value, mask):
        if value is None or len(value) == 0:
            return _zero_int

        length = len(value)
        result = [self._write_length(length)]
        binoutput = self._base.binoutput
        for x, m in zip(value, value.mask):
            result.append(binoutput(x, m))
        return _empty_bytes.join(result)


class ArrayVarArray(VarArray):
    """
    Handles an array of variable-length arrays, i.e. where *arraysize*
    ends in '*'.
    """

    def parse(self, value, config=None, pos=None):
        if value.strip() == "":
            return ma.array([]), False

        parts = self._splitter(value, config, pos)
        items = self._base._items
        parse_parts = self._base.parse_parts
        if len(parts) % items != 0:
            vo_raise(E02, (items, len(parts)), config, pos)
        result = []
        result_mask = []
        for i in range(0, len(parts), items):
            value, mask = parse_parts(parts[i : i + items], config, pos)
            result.append(value)
            result_mask.append(mask)

        return _make_masked_array(result, result_mask), False


class _StringArrayMixin:
    """
    Reading and writing arrays of fixed-length strings as TABLEDATA text.

    The strings are stored end to end with no separator, each padded with
    blanks to the string length.  ``self._scalar`` is the `Char` or
    `UnicodeChar` converter for one string.
    """

    @property
    def _utf8(self):
        # From VOTable 1.6 the string length of a char array is in UTF-8 bytes
        return getattr(self._scalar, "is_utf8_version", False)

    def _split(self, value):
        # The padding of the last string never arrives: the XML parser
        # strips trailing whitespace from the cell.
        width = self._scalar.arraysize
        if self._utf8:
            data = value.encode("utf-8")
            parts = (
                self._scalar._decode(data[i : i + width])
                for i in range(0, len(data), width)
            )
            return [x.rstrip(" ") for x in parts]
        return [value[i : i + width].rstrip(" ") for i in range(0, len(value), width)]

    def _parse_strings(self, parts, shape, config, pos):
        parse = self._scalar.parse
        result = [parse(x, config, pos)[0] for x in parts]
        return np.array(result, dtype=self._scalar.format).reshape(shape)

    def _flat(self, value, mask):
        """Pairs of (string, mask) for every string of *value*."""
        # Lists of Python objects are much faster to loop over than arrays
        strings = np.asarray(value).ravel().tolist()
        if np.ndim(mask) == 0:
            return zip(strings, repeat(bool(mask)))
        return zip(strings, np.asarray(mask).ravel().tolist(), strict=True)

    def _binoutput_strings(self, value, mask):
        binoutput = self._scalar.binoutput
        return _empty_bytes.join(binoutput(x, m) for x, m in self._flat(value, mask))

    def output(self, value, mask):
        width = self._scalar.arraysize
        output = self._scalar.output
        result = []
        for x, m in self._flat(value, mask):
            if m:
                x = ""
            elif isinstance(x, bytes):
                x = x.decode("utf-8")
            # Pad before escaping, so that "&amp;" counts as one character
            if self._utf8:
                encoded = x.encode("utf-8", errors="replace")
                if len(encoded) > width:
                    vo_warn(W46, (self._datatype, width), None, None)
                    encoded = _truncate_utf8(encoded, width)
                padded = encoded.ljust(width, b" ").decode("utf-8")
            else:
                if len(x) > width:
                    vo_warn(W46, (self._datatype, width), None, None)
                padded = f"{x[:width]:<{width}}"
            result.append(output(padded, False))
        return "".join(result)


class StringArray(_StringArrayMixin, Array):
    """
    Handles a fixed-size array of fixed-length strings, such as
    ``arraysize="10x4"``: 4 strings of 10 characters.

    Each string is converted by *base*, the `Char` or `UnicodeChar`
    converter of the field.
    """

    def __init__(self, field, base, arraysize, config=None, pos=None):
        Array.__init__(self, field, config, pos)

        self._base = self._scalar = base
        self._datatype = field.datatype
        self._arraysize = arraysize
        self.format = f"{tuple(arraysize)}{base.format}"
        self._items = prod(arraysize)
        self.default = np.zeros(arraysize, dtype=base.format)

    def parse(self, value, config=None, pos=None):
        parts = self._split(value)
        if len(parts) > self._items:
            warn_or_raise(E02, E02, (self._items, len(parts)), config, pos)
        # Too few is not an error: trailing empty strings are all padding,
        # which the XML parser strips.
        parts = (parts + [""] * self._items)[: self._items]
        return self.parse_parts(parts, config, pos)

    def parse_parts(self, parts, config=None, pos=None):
        result = self._parse_strings(parts, self._arraysize, config, pos)
        return result, np.zeros(self._arraysize, dtype="bool")

    def binparse(self, read):
        binparse = self._base.binparse
        result = [binparse(read)[0] for _ in range(self._items)]
        result = np.array(result, dtype=self._base.format).reshape(self._arraysize)
        return result, np.zeros(self._arraysize, dtype="bool")

    def binoutput(self, value, mask):
        return self._binoutput_strings(value, mask)


class StringVarArray(_StringArrayMixin, VarArray):
    """
    Handles a variable number of fixed-length strings, such as
    ``arraysize="10x*"``, or of fixed-size arrays of them, such as
    ``arraysize="10x4x*"``.

    *base* is the `Char` or `UnicodeChar` converter of the field in the
    first case, and a `StringArray` in the second.
    """

    def __init__(self, field, base, arraysize, config=None, pos=None):
        VarArray.__init__(self, field, base, arraysize, config, pos)

        self._scalar = getattr(base, "_base", base)
        self._datatype = field.datatype
        self._shape = (-1, *arraysize)
        self._items = prod(arraysize)
        # The BINARY formats store the length of the array in primitives:
        # characters, or UTF-8 bytes for char from VOTable 1.6
        self._chars = self._items * self._scalar.arraysize
        self.default = np.zeros((0, *arraysize), dtype=self._scalar.format)

    def parse(self, value, config=None, pos=None):
        parts = self._split(value)
        # As in StringArray.parse, trailing empty strings arrive stripped.
        parts += [""] * (-len(parts) % self._items)
        return ma.array(self._parse_strings(parts, self._shape, config, pos)), False

    def binparse(self, read):
        length = self._parse_length(read) // self._chars
        binparse = self._base.binparse
        result = [binparse(read)[0] for _ in range(length)]
        result = np.array(result, dtype=self._scalar.format).reshape(self._shape)
        return ma.array(result), False

    def output(self, value, mask):
        # *mask* is for the cell as a whole.  Which strings are masked is
        # in the mask of *value*.
        strings_mask = ma.getmask(value)
        if strings_mask is not ma.nomask:
            mask = mask | strings_mask
        return super().output(value, mask)

    def binoutput(self, value, mask):
        if value is None or len(value) == 0:
            return _zero_int

        length = self._write_length(len(value) * self._chars)
        return length + self._binoutput_strings(value, ma.getmaskarray(value))


class ScalarVarArray(VarArray):
    """
    Handles a variable-length array of numeric scalars.
    """

    def parse(self, value, config=None, pos=None):
        if value.strip() == "":
            return ma.array([]), False

        parts = self._splitter(value, config, pos)

        parse = self._base.parse
        result = []
        result_mask = []
        for x in parts:
            value, mask = parse(x, config, pos)
            result.append(value)
            result_mask.append(mask)

        return _make_masked_array(result, result_mask), False


class NumericArray(Array):
    """
    Handles a fixed-length array of numeric scalars.
    """

    vararray_type = ArrayVarArray

    def __init__(self, field, base, arraysize, config=None, pos=None):
        Array.__init__(self, field, config, pos)

        self._base = base
        self._arraysize = arraysize
        self.format = f"{tuple(arraysize)}{base.format}"
        self._items = prod(arraysize)
        self._memsize = np.dtype(self.format).itemsize
        self._bigendian_format = ">" + self.format

        self.default = np.empty(arraysize, dtype=self._base.format)
        self.default[...] = self._base.default

    def parse(self, value, config=None, pos=None):
        if config is None:
            config = {}
        elif config["version_1_3_or_later"] and value == "":
            return np.zeros(self._arraysize, dtype=self._base.format), True
        parts = self._splitter(value, config, pos)
        if len(parts) != self._items:
            warn_or_raise(E02, E02, (self._items, len(parts)), config, pos)
        if config.get("verify", "ignore") == "exception":
            return self.parse_parts(parts, config, pos)
        else:
            if len(parts) == self._items:
                pass
            elif len(parts) > self._items:
                parts = parts[: self._items]
            else:
                parts = parts + ([self._base.default] * (self._items - len(parts)))
            return self.parse_parts(parts, config, pos)

    def parse_parts(self, parts, config=None, pos=None):
        base_parse = self._base.parse
        result = []
        result_mask = []
        for x in parts:
            value, mask = base_parse(x, config, pos)
            result.append(value)
            result_mask.append(mask)
        result = np.array(result, dtype=self._base.format).reshape(self._arraysize)
        result_mask = np.array(result_mask, dtype="bool").reshape(self._arraysize)
        return result, result_mask

    def output(self, value, mask):
        base_output = self._base.output
        value = np.asarray(value)
        mask = np.asarray(mask)
        if mask.size <= 1:
            func = np.broadcast
        else:  # When mask is already array but value is scalar, this prevents broadcast
            func = zip
        return " ".join(base_output(x, m) for x, m in func(value.flat, mask.flat))

    def binparse(self, read):
        result = np.frombuffer(read(self._memsize), dtype=self._bigendian_format)[0]
        result_mask = self._base.is_null(result)
        return result, result_mask

    def binoutput(self, value, mask):
        filtered = self._base.filter_array(value, mask)
        filtered = _ensure_bigendian(filtered)
        return filtered.tobytes()


class Numeric(Converter):
    """
    The base class for all numeric data types.
    """

    array_type = NumericArray
    vararray_type = ScalarVarArray
    null = None

    def __init__(self, field, config=None, pos=None):
        Converter.__init__(self, field, config, pos)

        self._memsize = np.dtype(self.format).itemsize
        self._bigendian_format = ">" + self.format
        if field.values.null is not None:
            self.null = np.asarray(field.values.null, dtype=self.format)
            self.default = self.null
            self.is_null = self._is_null
        else:
            self.is_null = np.isnan

    def binparse(self, read):
        result = np.frombuffer(read(self._memsize), dtype=self._bigendian_format)
        return result[0], self.is_null(result[0])

    def _is_null(self, value):
        return value == self.null


class FloatingPoint(Numeric):
    """
    The base class for floating-point datatypes.
    """

    default = np.nan

    def __init__(self, field, config=None, pos=None):
        if config is None:
            config = {}

        Numeric.__init__(self, field, config, pos)

        if self.null is not None:
            self._fast_null_param = float(self.null)
            self._fast_has_null = True
        else:
            self._fast_null_param = 0.0
            self._fast_has_null = False

        precision = field.precision
        width = field.width

        if precision is None:
            format_parts = ["{!s:>"]
        else:
            format_parts = ["{:"]

        if width is not None:
            format_parts.append(str(width))

        if precision is not None:
            if precision.startswith("E"):
                format_parts.append(f".{int(precision[1:]):d}g")
            elif precision.startswith("F"):
                format_parts.append(f".{int(precision[1:]):d}f")
            else:
                format_parts.append(f".{int(precision):d}f")

        format_parts.append("}")

        self._output_format = "".join(format_parts)

        self.nan = np.array(np.nan, self.format)

        if self.null is None:
            self._null_output = "NaN"
            self._null_binoutput = self.binoutput(self.nan, False)
            self.filter_array = self._filter_nan
        else:
            self._null_output = self.output(np.asarray(self.null), False)
            self._null_binoutput = self.binoutput(np.asarray(self.null), False)
            self.filter_array = self._filter_null

        if config.get("verify", "ignore") == "exception":
            self.parse = self._parse_pedantic
        else:
            self.parse = self._parse_permissive

    def supports_empty_values(self, config):
        return True

    def _parse_pedantic(self, value, config=None, pos=None):
        if value.strip() == "":
            return self.null, True
        f = float(value)
        return f, self.is_null(f)

    def _parse_permissive(self, value, config=None, pos=None):
        try:
            f = float(value)
            return f, self.is_null(f)
        except ValueError:
            # IRSA VOTables use the word 'null' to specify empty values,
            # but this is not defined in the VOTable spec.
            if value.strip() != "":
                vo_warn(W30, value, config, pos)
            return self.null, True

    @property
    def output_format(self):
        return self._output_format

    def output(self, value, mask):
        if mask:
            return self._null_output
        if np.isfinite(value):
            if not np.isscalar(value):
                value = value.dtype.type(value)
            result = self._output_format.format(value)
            if result.startswith("array"):
                raise RuntimeError()
            if self._output_format[2] == "s" and result.endswith(".0"):
                result = result[:-2]
            return result
        elif np.isnan(value):
            return "NaN"
        elif np.isposinf(value):
            return "+InF"
        elif np.isneginf(value):
            return "-InF"
        # Should never raise
        vo_raise(f"Invalid floating point value '{value}'")

    def binoutput(self, value, mask):
        if mask:
            return self._null_binoutput

        value = _ensure_bigendian(value)
        return value.tobytes()

    def _filter_nan(self, value, mask):
        return np.where(mask, np.nan, value)

    def _filter_null(self, value, mask):
        return np.where(mask, self.null, value)


class Double(FloatingPoint):
    """
    Handles the double datatype.  Double-precision IEEE floating-point.
    """

    format = "f8"

    def binparse(self, read):
        data = read(8)
        return fast_binparse_double(data, 0, self._fast_null_param, self._fast_has_null)


class Float(FloatingPoint):
    """
    Handles the float datatype.  Single-precision IEEE floating-point.
    """

    format = "f4"

    def binparse(self, read):
        data = read(4)
        return fast_binparse_float(data, 0, self._fast_null_param, self._fast_has_null)


class Integer(Numeric):
    """
    The base class for all the integral datatypes.
    """

    default = 0

    def __init__(self, field, config=None, pos=None):
        Numeric.__init__(self, field, config, pos)

        if self.null is not None:
            self._fast_null_param = int(self.null)
            self._fast_has_null = True
        else:
            self._fast_null_param = 0
            self._fast_has_null = False

    def parse(self, value, config=None, pos=None):
        if config is None:
            config = {}
        mask = False
        if isinstance(value, str):
            value = value.lower()
            if value == "":
                if config["version_1_3_or_later"]:
                    mask = True
                else:
                    warn_or_raise(W49, W49, (), config, pos)
                if self.null is not None:
                    value = self.null
                else:
                    value = self.default
            elif value == "nan":
                mask = True
                if self.null is None:
                    warn_or_raise(W31, W31, (), config, pos)
                    value = self.default
                else:
                    value = self.null
            elif value.startswith("0x"):
                value = int(value[2:], 16)
            else:
                value = int(value, 10)
        else:
            value = int(value)
        if self.null is not None and value == self.null:
            mask = True

        if value < self.val_range[0]:
            warn_or_raise(W51, W51, (value, self.bit_size), config, pos)
            value = self.val_range[0]
        elif value > self.val_range[1]:
            warn_or_raise(W51, W51, (value, self.bit_size), config, pos)
            value = self.val_range[1]

        return value, mask

    def output(self, value, mask):
        if mask:
            if self.null is None:
                warn_or_raise(W31, W31)
                return "NaN"
            return str(self.null)
        return str(value)

    def binoutput(self, value, mask):
        if mask:
            if self.null is None:
                vo_raise(W31)
            else:
                value = self.null

        value = _ensure_bigendian(value)
        return value.tobytes()

    def filter_array(self, value, mask):
        if np.any(mask):
            if self.null is not None:
                return np.where(mask, self.null, value)
            else:
                vo_raise(W31)
        return value


class UnsignedByte(Integer):
    """
    Handles the unsignedByte datatype.  Unsigned 8-bit integer.
    """

    format = "u1"
    val_range = (0, 255)
    bit_size = "8-bit unsigned"

    def binparse(self, read):
        data = read(1)
        return fast_binparse_ubyte(data, 0, self._fast_null_param, self._fast_has_null)


class Short(Integer):
    """
    Handles the short datatype.  Signed 16-bit integer.
    """

    format = "i2"
    val_range = (-32768, 32767)
    bit_size = "16-bit"

    def binparse(self, read):
        data = read(2)
        return fast_binparse_short(data, 0, self._fast_null_param, self._fast_has_null)


class Int(Integer):
    """
    Handles the int datatype.  Signed 32-bit integer.
    """

    format = "i4"
    val_range = (-2147483648, 2147483647)
    bit_size = "32-bit"

    def binparse(self, read):
        data = read(4)
        return fast_binparse_int(data, 0, self._fast_null_param, self._fast_has_null)


class Long(Integer):
    """
    Handles the long datatype.  Signed 64-bit integer.
    """

    format = "i8"
    val_range = (-9223372036854775808, 9223372036854775807)
    bit_size = "64-bit"

    def binparse(self, read):
        data = read(8)
        return fast_binparse_long(data, 0, self._fast_null_param, self._fast_has_null)


class ComplexArrayVarArray(VarArray):
    """
    Handles an array of variable-length arrays of complex numbers.
    """

    def parse(self, value, config=None, pos=None):
        if value.strip() == "":
            return ma.array([]), True

        parts = self._splitter(value, config, pos)
        items = self._base._items
        parse_parts = self._base.parse_parts
        if len(parts) % items != 0:
            vo_raise(E02, (items, len(parts)), config, pos)
        result = []
        result_mask = []
        for i in range(0, len(parts), items):
            value, mask = parse_parts(parts[i : i + items], config, pos)
            result.append(value)
            result_mask.append(mask)

        return _make_masked_array(result, result_mask), False


class ComplexVarArray(VarArray):
    """
    Handles a variable-length array of complex numbers.
    """

    def parse(self, value, config=None, pos=None):
        if value.strip() == "":
            return ma.array([]), True

        parts = self._splitter(value, config, pos)
        parse_parts = self._base.parse_parts
        result = []
        result_mask = []
        for i in range(0, len(parts), 2):
            value = [float(x) for x in parts[i : i + 2]]
            value, mask = parse_parts(value, config, pos)
            result.append(value)
            result_mask.append(mask)

        return (
            _make_masked_array(np.array(result, dtype=self._base.format), result_mask),
            False,
        )


class ComplexArray(NumericArray):
    """
    Handles a fixed-size array of complex numbers.
    """

    vararray_type = ComplexArrayVarArray

    def __init__(self, field, base, arraysize, config=None, pos=None):
        NumericArray.__init__(self, field, base, arraysize, config, pos)
        self._items *= 2

    def parse(self, value, config=None, pos=None):
        parts = self._splitter(value, config, pos)
        if parts == [""]:
            parts = []
        return self.parse_parts(parts, config, pos)

    def parse_parts(self, parts, config=None, pos=None):
        if len(parts) != self._items:
            vo_raise(E02, (self._items, len(parts)), config, pos)
        base_parse = self._base.parse_parts
        result = []
        result_mask = []
        for i in range(0, self._items, 2):
            value = [float(x) for x in parts[i : i + 2]]
            value, mask = base_parse(value, config, pos)
            result.append(value)
            result_mask.append(mask)
        result = np.array(result, dtype=self._base.format).reshape(self._arraysize)
        result_mask = np.array(result_mask, dtype="bool").reshape(self._arraysize)
        return result, result_mask


class Complex(FloatingPoint, Array):
    """
    The base class for complex numbers.
    """

    array_type = ComplexArray
    vararray_type = ComplexVarArray
    default = np.nan

    def __init__(self, field, config=None, pos=None):
        FloatingPoint.__init__(self, field, config, pos)
        Array.__init__(self, field, config, pos)

    def parse(self, value, config=None, pos=None):
        stripped = value.strip()
        if stripped == "" or stripped.lower() == "nan":
            return np.nan, True
        splitter = self._splitter
        parts = [float(x) for x in splitter(value, config, pos)]
        if len(parts) != 2:
            vo_raise(E03, (value,), config, pos)
        return self.parse_parts(parts, config, pos)

    _parse_permissive = parse
    _parse_pedantic = parse

    def parse_parts(self, parts, config=None, pos=None):
        value = complex(*parts)
        return value, self.is_null(value)

    def output(self, value, mask):
        if mask:
            if self.null is None:
                return "NaN"
            else:
                value = self.null
        real = self._output_format.format(float(value.real))
        imag = self._output_format.format(float(value.imag))
        if self._output_format[2] == "s":
            real = real.removesuffix(".0")
            imag = imag.removesuffix(".0")
        return real + " " + imag


class FloatComplex(Complex):
    """
    Handle floatComplex datatype.  Pair of single-precision IEEE
    floating-point numbers.
    """

    format = "c8"


class DoubleComplex(Complex):
    """
    Handle doubleComplex datatype.  Pair of double-precision IEEE
    floating-point numbers.
    """

    format = "c16"


class BitArray(NumericArray):
    """
    Handles an array of bits.
    """

    vararray_type = ArrayVarArray

    def __init__(self, field, base, arraysize, config=None, pos=None):
        NumericArray.__init__(self, field, base, arraysize, config, pos)

        self._bytes = ((self._items - 1) // 8) + 1

    @staticmethod
    def _splitter_pedantic(value, config=None, pos=None):
        return list(re.sub(r"\s", "", value))

    @staticmethod
    def _splitter_lax(value, config=None, pos=None):
        if "," in value:
            vo_warn(W01, (), config, pos)
        return list(re.sub(r"\s|,", "", value))

    def output(self, value, mask):
        if np.any(mask):
            vo_warn(W39)
        value = np.asarray(value)
        mapping = {False: "0", True: "1"}
        return "".join(mapping[x] for x in value.flat)

    def binparse(self, read):
        data = read(self._bytes)
        result = bitarray_to_bool(data, self._items)
        result = result.reshape(self._arraysize)
        result_mask = np.zeros(self._arraysize, dtype="b1")
        return result, result_mask

    def binoutput(self, value, mask):
        if np.any(mask):
            vo_warn(W39)
        return bool_to_bitarray(value)


class Bit(Converter):
    """
    Handles the bit datatype.
    """

    format = "b1"
    array_type = BitArray
    vararray_type = ScalarVarArray
    default = False
    binary_one = b"\x08"
    binary_zero = b"\0"

    def parse(self, value, config=None, pos=None):
        if config is None:
            config = {}
        mapping = {"1": True, "0": False}
        if value is False or value.strip() == "":
            if not config["version_1_3_or_later"]:
                warn_or_raise(W49, W49, (), config, pos)
            return False, True
        else:
            try:
                return mapping[value], False
            except KeyError:
                vo_raise(E04, (value,), config, pos)

    def output(self, value, mask):
        if mask:
            vo_warn(W39)

        if value:
            return "1"
        else:
            return "0"

    def binparse(self, read):
        data = read(1)
        return fast_binparse_bit(data)

    def binoutput(self, value, mask):
        if mask:
            vo_warn(W39)

        if value:
            return self.binary_one
        return self.binary_zero


class BooleanArray(NumericArray):
    """
    Handles an array of boolean values.
    """

    vararray_type = ArrayVarArray

    def binparse(self, read):
        data = read(self._items)
        binparse = self._base.binparse_value
        result = []
        result_mask = []
        for char in data:
            value, mask = binparse(char)
            result.append(value)
            result_mask.append(mask)
        result = np.array(result, dtype="b1").reshape(self._arraysize)
        result_mask = np.array(result_mask, dtype="b1").reshape(self._arraysize)
        return result, result_mask

    def binoutput(self, value, mask):
        binoutput = self._base.binoutput
        value = np.asarray(value)
        mask = np.asarray(mask)
        result = [binoutput(x, m) for x, m in np.broadcast(value.flat, mask.flat)]
        return _empty_bytes.join(result)


class Boolean(Converter):
    """
    Handles the boolean datatype.
    """

    format = "b1"
    array_type = BooleanArray
    vararray_type = ScalarVarArray
    default = False
    binary_question_mark = b"?"
    binary_true = b"T"
    binary_false = b"F"

    def parse(self, value, config=None, pos=None):
        if value == "":
            return False, True
        if value is False:
            return False, True
        mapping = {
            "TRUE": (True, False),
            "FALSE": (False, False),
            "1": (True, False),
            "0": (False, False),
            "T": (True, False),
            "F": (False, False),
            "\0": (False, True),
            " ": (False, True),
            "?": (False, True),
            "": (False, True),
        }
        try:
            return mapping[value.upper()]
        except KeyError:
            vo_raise(E05, (value,), config, pos)

    def output(self, value, mask):
        if mask:
            return "?"
        if value:
            return "T"
        return "F"

    def binparse(self, read):
        data = read(1)
        return fast_binparse_bool(data)

    _binparse_mapping = {
        ord("T"): (True, False),
        ord("t"): (True, False),
        ord("1"): (True, False),
        ord("F"): (False, False),
        ord("f"): (False, False),
        ord("0"): (False, False),
        ord("\0"): (False, True),
        ord(" "): (False, True),
        ord("?"): (False, True),
    }

    def binparse_value(self, value):
        try:
            return self._binparse_mapping[value]
        except KeyError:
            vo_raise(E05, (value,))

    def binoutput(self, value, mask):
        if mask:
            return self.binary_question_mark
        if value:
            return self.binary_true
        return self.binary_false


converter_mapping = {
    "double": Double,
    "float": Float,
    "bit": Bit,
    "boolean": Boolean,
    "unsignedByte": UnsignedByte,
    "short": Short,
    "int": Int,
    "long": Long,
    "floatComplex": FloatComplex,
    "doubleComplex": DoubleComplex,
    "char": Char,
    "unicodeChar": UnicodeChar,
}


def get_converter(field, config=None, pos=None):
    """
    Get an appropriate converter instance for a given field.

    Parameters
    ----------
    field : astropy.io.votable.tree.Field

    config : dict, optional
        Parser configuration dictionary

    pos : tuple
        Position in the input XML file.  Used for error messages.

    Returns
    -------
    converter : astropy.io.votable.converters.Converter
    """
    if config is None:
        config = {}

    if field.datatype not in converter_mapping:
        vo_raise(E06, (field.datatype, field.ID), config)

    cls = converter_mapping[field.datatype]
    converter = cls(field, config, pos)

    arraysize = field.arraysize

    # With numeric datatypes, special things need to happen for
    # arrays.
    if field.datatype not in ("char", "unicodeChar") and arraysize is not None:
        if arraysize[-1] == "*":
            arraysize = arraysize[:-1]
            last_x = arraysize.rfind("x")
            if last_x == -1:
                arraysize = ""
            else:
                arraysize = arraysize[:last_x]
            fixed = False
        else:
            fixed = True

        if arraysize != "":
            arraysize = [int(x) for x in arraysize.split("x")]
            arraysize.reverse()
        else:
            arraysize = []

        if arraysize != []:
            converter = converter.array_type(field, converter, arraysize, config)

        if not fixed:
            converter = converter.vararray_type(field, converter, arraysize, config)

    if field.datatype in ("char", "unicodeChar") and "x" in (arraysize or ""):
        # The first dimension is the string length, already handled by the
        # scalar converter
        dims = arraysize.split("x")[1:]
        fixed = dims[-1] != "*"
        if not fixed:
            dims.pop()
        shape = [int(x) if x.isdecimal() else 0 for x in reversed(dims)]
        if 0 in (converter.arraysize, *shape):
            vo_raise(E01, (arraysize, field.datatype, field.ID), config, pos)

        if shape:
            converter = StringArray(field, converter, shape, config)
        if not fixed:
            converter = StringVarArray(field, converter, shape, config)

    return converter


numpy_dtype_to_field_mapping = {
    np.float64().dtype.num: "double",
    np.float32().dtype.num: "float",
    np.bool_().dtype.num: "bit",
    np.uint8().dtype.num: "unsignedByte",
    np.int16().dtype.num: "short",
    np.int32().dtype.num: "int",
    np.int64().dtype.num: "long",
    np.complex64().dtype.num: "floatComplex",
    np.complex128().dtype.num: "doubleComplex",
    np.str_().dtype.num: "unicodeChar",
    np.bytes_().dtype.num: "char",
}


def _all_matching_dtype(column):
    first_dtype = False
    first_shape = ()
    for x in column:
        if not isinstance(x, np.ndarray) or len(x) == 0:
            continue

        if first_dtype is False:
            first_dtype = x.dtype
            first_shape = x.shape[1:]
        elif first_dtype != x.dtype:
            return False, ()
        elif first_shape != x.shape[1:]:
            first_shape = ()
    return first_dtype, first_shape


def numpy_to_votable_dtype(dtype, shape):
    """
    Converts a numpy dtype and shape to a dictionary of attributes for
    a VOTable FIELD element and correspond to that type.

    Parameters
    ----------
    dtype : Numpy dtype instance

    shape : tuple

    Returns
    -------
    attributes : dict
        A dict containing 'datatype' and 'arraysize' keys that can be
        set on a VOTable FIELD element.
    """
    if dtype.num not in numpy_dtype_to_field_mapping:
        raise TypeError(f"{dtype!r} can not be represented in VOTable")

    if dtype.char in "SU":
        # The string length comes first, then the shape with the fastest
        # varying dimension first, as get_converter reads it.
        if dtype.char == "S":
            datatype, length = "char", dtype.itemsize
        else:
            datatype, length = "unicodeChar", dtype.itemsize // 4
        arraysize = "x".join(str(x) for x in (length, *reversed(shape)))
        return {"datatype": datatype, "arraysize": arraysize}
    else:
        result = {"datatype": numpy_dtype_to_field_mapping[dtype.num]}
        if len(shape):
            result["arraysize"] = "x".join(str(x) for x in shape)

        return result


def _utf8_char_datatype(column, min_arraysize=0):
    n_chars = column.dtype.itemsize // 4
    data = np.asarray(column)
    if data.size:
        n_bytes = np.char.str_len(np.char.encode(data, "utf-8"))
        max_bytes = int(n_bytes.max())
        is_ascii = np.array_equal(n_bytes, np.char.str_len(data))
    else:
        max_bytes = 0
        is_ascii = True

    arraysize = max(n_chars, max_bytes, min_arraysize, 1)
    result = {"datatype": "char", "arraysize": str(arraysize)}
    # Readers may assume ASCII when width equals arraysize (VOTable 1.6 4.2).
    # None also replaces any stale width from the column metadata.
    result["width"] = n_chars if n_chars and (is_ascii or n_chars < arraysize) else None
    return result


def _max_utf8_bytes(column):
    """The longest UTF-8 encoded string in *column*, also inside arrays."""
    n = 0
    for x in column:
        for s in [x] if isinstance(x, str) else np.ravel(x):
            if isinstance(s, str):
                n = max(n, len(s.encode("utf-8")))
    return n


def _all_strings(column):
    return all(
        isinstance(x, (str, bytes)) or x is np.ma.masked or x is None for x in column
    )


def table_column_to_votable_datatype(column, *, version=None):
    """
    Given a `astropy.table.Column` instance, returns the attributes
    necessary to create a VOTable FIELD element that corresponds to
    the type of the column.

    This necessarily must perform some heuristics to determine the
    type of variable length arrays fields, since they are not directly
    supported by Numpy.

    If the column has dtype of "object", it performs the following
    tests:

       - If all elements are byte or unicode strings, it creates a
         variable-length byte or unicode field, respectively.

       - If all elements are numpy arrays of the same dtype and with a
         consistent shape in all but the first dimension, it creates a
         variable length array of fixed sized arrays.  If the dtypes
         match, but the shapes do not, a variable length array is
         created.

    If the dtype of the input is not understood, it sets the data type
    to the most inclusive: a variable length unicodeChar array.

    For VOTable 1.6 and later, string columns are written as UTF-8 ``char``,
    also when they were read from a ``unicodeChar`` FIELD. Fixed-width
    columns get an arraysize in bytes and a ``width`` in characters.

    Parameters
    ----------
    column : `astropy.table.Column` instance

    version : str, optional
        The VOTable version being written. If not given, the pre-1.6
        mapping is used.

    Returns
    -------
    attributes : dict
        A dict containing 'datatype' and 'arraysize' keys, and possibly
        'width', that can be set on a VOTable FIELD element.
    """
    utf8_char = version is not None and version_compare(version, "1.6") >= 0
    votable_string_dtype = None
    max_length = None
    original_arraysize = None

    if column.info.meta is not None:
        votable_string_dtype = column.info.meta.get("_votable_string_dtype")
        # Check if we have stored the original arraysize with bounds
        # If so, extract the max length
        original_arraysize = column.info.meta.get("_votable_arraysize")
        if (
            original_arraysize is not None
            and original_arraysize.endswith("*")
            and not original_arraysize == "*"
        ):
            max_length = original_arraysize[:-1]

    if column.dtype.char == "O":
        arraysize = "*"

        # If max length is stored use it to create a bounded var-length array
        if max_length is not None:
            arraysize = f"{max_length}*"

        if utf8_char and votable_string_dtype == "unicodeChar":
            # The unicodeChar bound counted characters, char counts bytes.
            if max_length is not None:
                n, bounded = _string_bound(original_arraysize)
                rest = original_arraysize[len(_string_length(original_arraysize)) :]
                n = max(n, _max_utf8_bytes(column))
                arraysize = f"{n}{'*' if bounded else ''}{rest}"
            return {"datatype": "char", "arraysize": arraysize}
        elif votable_string_dtype is not None:
            return {"datatype": votable_string_dtype, "arraysize": arraysize}
        elif isinstance(column[0], np.ndarray):
            dtype, shape = _all_matching_dtype(column)
            if dtype is not False:
                result = numpy_to_votable_dtype(dtype, shape)
                if "arraysize" not in result:
                    result["arraysize"] = arraysize
                else:
                    result["arraysize"] += "*"
                return result

        if utf8_char and _all_strings(column):
            return {"datatype": "char", "arraysize": arraysize}

        # All bets are off, do the most generic thing
        return {"datatype": "unicodeChar", "arraysize": arraysize}

    if utf8_char and column.dtype.char == "U":
        n, bounded = _string_bound(original_arraysize)
        result = _utf8_char_datatype(column, n)
        shape = column.shape[1:]
        if shape:
            result["arraysize"] += _shape_suffix(shape)
        elif bounded:
            result["arraysize"] += "*"
        return result

    # For fixed size string columns, datatype here will be unicodeChar,
    # but honor the original FIELD datatype if present.
    result = numpy_to_votable_dtype(column.dtype, column.shape[1:])
    if result["datatype"] == "unicodeChar" and votable_string_dtype == "char":
        result["datatype"] = "char"

    # Keep the original string length unless the column has grown wider
    if original_arraysize is not None and "arraysize" in result:
        n, bounded = _string_bound(original_arraysize)
        derived = _string_length(result["arraysize"])
        if n and derived.isdigit():
            rest = result["arraysize"][len(derived) :]
            n = max(n, int(derived))
            result["arraysize"] = f"{n}{'*' if bounded and not rest else ''}{rest}"
        else:
            result["arraysize"] = original_arraysize

    return result
