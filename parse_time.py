# -*- coding: utf-8 -*-
""" Example of wrapping a C library function that accepts a C double array as
    input using the numpy.ctypeslib. """
import time

import numpy as np
import numpy.ctypeslib as npct
from ctypes import c_int

import erfa

# input type for the parse_iso_times function
# must be a double array, with single dimension that is contiguous
array_1d_char = npct.ndpointer(dtype=np.uint8, ndim=1, flags='C_CONTIGUOUS')
array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.intc, ndim=1, flags='C_CONTIGUOUS')

# load the library, using numpy mechanisms
libpt = npct.load_library("libparse_time.so", ".")

# setup the return types and argument types
libpt.parse_iso_times.restype = c_int
libpt.parse_iso_times.argtypes = [array_1d_char, c_int, c_int,
                                  array_1d_int, array_1d_int, array_1d_int,
                                  array_1d_int, array_1d_int, array_1d_double]
libpt.check_unicode.restype = c_int
libpt.check_unicode.argtypes = [array_1d_char, c_int]

n_times = 1000000
# This fails as expected:
#   val1 = np.array(['2020-01-01 1ᛦ:13:14.4324'] * n_times)
val1 = np.array(['2020-01-01 12:13:14.4324'] * n_times)
t0 = time.time()
char_size = 4 if val1.dtype.kind == 'U' else 1
val1_str_len = int(val1.dtype.itemsize // char_size)
chars = val1.ravel().view(np.uint8)
if char_size == 4:
    status = libpt.check_unicode(chars, len(chars) // 4)
    if status < 0:
        raise ValueError('input is not pure ASCII')
    chars = chars[::4]
chars = np.ascontiguousarray(chars)

year = np.zeros(n_times, dtype=np.intc)
month = np.zeros(n_times, dtype=np.intc)
day = np.zeros(n_times, dtype=np.intc)
hour = np.zeros(n_times, dtype=np.intc)
minute = np.zeros(n_times, dtype=np.intc)
second = np.zeros(n_times, dtype=np.double)

status = libpt.parse_iso_times(chars, n_times, val1_str_len,
                               year, month, day, hour, minute, second)
time_per_parse = (time.time() - t0) / n_times / 1e-9
print(f'parse time = {time_per_parse:.0f} nsec')

jd1, jd2 = erfa.dtf2d(b'UTC', year, month, day, hour, minute, second)
time_per_parse = (time.time() - t0) / n_times / 1e-9
print(f'total time = {time_per_parse:.0f} nsec')
