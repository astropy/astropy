# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

"""
Fast binary converters for VOTable Binary2 parsing.

This module replaces the slow python binary parsing in astropy's VOTable
implementation with more performant Cython code.
"""

import numpy as np
import sys
cimport numpy as cnp
cimport cython
from libc.string cimport memcpy
from libc.stdlib cimport malloc, free

cnp.import_array()

ctypedef cnp.float64_t float64_t
ctypedef cnp.float32_t float32_t
ctypedef cnp.int64_t int64_t
ctypedef cnp.int32_t int32_t
ctypedef cnp.int16_t int16_t
ctypedef cnp.uint8_t uint8_t
ctypedef cnp.uint8_t bool_t

cdef bint NEED_BYTESWAP = sys.byteorder == 'little'

def fast_bitarray_to_bool(data, int length):
    """Convert packed bits to bool array"""
    cdef const unsigned char[::1] byte_data

    if hasattr(data, 'tobytes'):
        byte_data = memoryview(data.tobytes())
    elif isinstance(data, (bytes, bytearray)):
        byte_data = memoryview(data)
    else:
        byte_data = data

    cdef cnp.ndarray[bool_t, ndim=1] results = np.empty(length, dtype=np.bool_)
    cdef int i = 0
    cdef int byte_idx = 0
    cdef int bit_no
    cdef unsigned char byte_val
    cdef bool_t bit_val

    while i < length and byte_idx < byte_data.shape[0]:
        byte_val = byte_data[byte_idx]
        for bit_no in range(7, -1, -1):
            if i >= length:
                break
            bit_val = (byte_val & (1 << bit_no)) != 0
            results[i] = bit_val
            i += 1
        byte_idx += 1

    return results

def fast_bool_to_bitarray(value):
    """
    Cython implementation of bool_to_bitarray.
    """
    if np.isscalar(value):
        flat_array = np.array([bool(value)], dtype=np.bool_)
    else:
        value_array = np.asarray(value, dtype=np.bool_)
        flat_array = value_array.ravel()

    cdef int length = len(flat_array)
    cdef int num_bytes = (length + 7) // 8
    cdef cnp.ndarray[unsigned char, ndim=1] result = np.zeros(num_bytes, dtype=np.uint8)

    cdef int i
    cdef int byte_idx
    cdef int bit_no

    for i in range(length):
        byte_idx = i // 8
        bit_no = 7 - (i % 8)
        if flat_array[i]:
            result[byte_idx] |= (1 << bit_no)

    return result.tobytes()

cdef inline void swap_bytes_64(unsigned char* data) noexcept nogil:
    """64-bit byte swap."""
    cdef unsigned char temp
    temp = data[0]; data[0] = data[7]; data[7] = temp
    temp = data[1]; data[1] = data[6]; data[6] = temp
    temp = data[2]; data[2] = data[5]; data[5] = temp
    temp = data[3]; data[3] = data[4]; data[4] = temp

cdef inline void swap_bytes_32(unsigned char* data) noexcept nogil:
    """32-bit byte swap."""
    cdef unsigned char t
    t = data[0]; data[0] = data[3]; data[3] = t
    t = data[1]; data[1] = data[2]; data[2] = t

cdef inline void swap_bytes_16(unsigned char* data) noexcept nogil:
    """16-bit byte swap."""
    cdef unsigned char temp
    temp = data[0]; data[0] = data[1]; data[1] = temp


def fast_binparse_double(const unsigned char[::1] data, int offset=0):
    """Parse 8-byte double from big-endian VOTable data."""
    cdef float64_t value
    cdef unsigned char temp_data[8]
    cdef int i

    for i in range(8):
        temp_data[i] = data[offset + i]

    if NEED_BYTESWAP:
        swap_bytes_64(temp_data)

    value = (<float64_t*>temp_data)[0]

    cdef bint is_null = value != value

    return value, is_null


def fast_binparse_float(const unsigned char[::1] data, int offset=0):
    """Parse 4-byte float."""
    cdef float32_t value
    cdef unsigned char temp_data[4]
    cdef int i

    for i in range(4):
        temp_data[i] = data[offset + i]

    if NEED_BYTESWAP:
        swap_bytes_32(temp_data)

    value = (<float32_t*>temp_data)[0]

    cdef bint is_null = value != value

    return value, is_null

def fast_binparse_long(const unsigned char[::1] data, int offset=0):
    """Parse 8-byte signed integer."""
    cdef int64_t value
    cdef unsigned char temp_data[8]
    cdef int i

    for i in range(8):
        temp_data[i] = data[offset + i]

    if NEED_BYTESWAP:
        swap_bytes_64(temp_data)

    value = (<int64_t*>temp_data)[0]

    return value, False

def fast_binparse_int(const unsigned char[::1] data, int offset=0):
    """Parse 4-byte signed int."""
    cdef int32_t value
    cdef unsigned char temp_data[4]
    cdef int i

    for i in range(4):
        temp_data[i] = data[offset + i]

    if NEED_BYTESWAP:
        swap_bytes_32(temp_data)

    value = (<int32_t*>temp_data)[0]

    return value, False

def fast_binparse_short(const unsigned char[::1] data, int offset=0):
    """Parse 2-byte signed short."""
    cdef int16_t value
    cdef unsigned char temp_data[2]
    cdef int i

    for i in range(2):
        temp_data[i] = data[offset + i]

    if NEED_BYTESWAP:
        swap_bytes_16(temp_data)

    value = (<int16_t*>temp_data)[0]

    return value, False

def fast_binparse_ubyte(const unsigned char[::1] data, int offset=0):
    """Parse unsigned byte. """
    cdef uint8_t value = data[offset]
    return value, False

def fast_binparse_bool(const unsigned char[::1] data, int offset=0):
    """Parse boolean from character representation."""
    cdef unsigned char byte_val = data[offset]
    cdef bint value = False
    cdef bint is_null = False

    if byte_val == ord(b'T') or byte_val == ord(b't') or byte_val == ord(b'1'):
        value = True
    elif byte_val == ord(b'F') or byte_val == ord(b'f') or byte_val == ord(b'0'):
        value = False
    elif byte_val == 0 or byte_val == ord(b' ') or byte_val == ord(b'?'):
        value = False
        is_null = True
    else:
        value = False
        is_null = True

    return value, is_null

def fast_binparse_bit(const unsigned char[::1] data, int offset=0):
    """Parse single bit from VOTable bit field"""
    cdef unsigned char byte_val = data[offset]
    cdef bint value = (byte_val & 0x8) != 0
    return value, False

# Mapping of datatype strings to according parsing function
cdef dict FAST_BINPARSE_MAP = {
    'double': fast_binparse_double,
    'float': fast_binparse_float,
    'long': fast_binparse_long,
    'int': fast_binparse_int,
    'short': fast_binparse_short,
    'unsignedByte': fast_binparse_ubyte,
    'boolean': fast_binparse_bool,
    'bit': fast_binparse_bit,
}

def get_fast_binparse_func(datatype):
    """Get the binparse function for a given datatype."""
    return FAST_BINPARSE_MAP.get(datatype)
