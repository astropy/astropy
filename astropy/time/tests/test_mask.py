# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
import numpy as np

from .. import Time

allclose_sec = functools.partial(np.allclose, rtol=2. ** -52,
                                 atol=2. ** -52 * 24 * 3600)  # 20 ps atol
is_masked = np.ma.is_masked


def test_simple():
    t = Time([1, 2, 3], format='cxcsec')
    assert t.masked is False

    # Before masking, format output is not a masked array (it is an ndarray
    # like always)
    assert not isinstance(t.value, np.ma.MaskedArray)
    assert not isinstance(t.unix, np.ma.MaskedArray)

    t[2] = np.ma.masked
    assert t.masked is True
    assert np.all(t.mask == [False, False, True])
    assert allclose_sec(t.value[:2], [1, 2])
    assert is_masked(t.value[2])
    assert is_masked(t[2].value)

    # After masking format output is a masked array
    assert isinstance(t.value, np.ma.MaskedArray)
    assert isinstance(t.unix, np.ma.MaskedArray)
    # Todo : test all formats


def test_str():
    t = Time(['2000:001', '2000:002'])
    t[1] = np.ma.masked
    assert str(t) == "['2000:001:00:00:00.000' --]"
    assert repr(t) == "<Time object: scale='utc' format='yday' value=['2000:001:00:00:00.000' --]>"
    assert repr(t.iso) == "masked_array(['2000-01-01 00:00:00.000' --])"

    # Assign value to unmask
    t[1] = '2000:111'
    assert str(t) == "['2000:001:00:00:00.000' '2000:111:00:00:00.000']"
    assert t.masked is False


def test_transform():
    t = Time(['2000:001', '2000:002'])
    t[1] = np.ma.masked

    # Change scale (this tests the ERFA machinery with masking as well)
    t_ut1 = t.ut1
    assert is_masked(t_ut1.value[1])
    assert not is_masked(t_ut1.value[0])
    assert np.all(t_ut1.mask == [False, True])
    
    # Change format
    t_unix = t.unix
    assert is_masked(t_unix[1])
    assert not is_masked(t_unix[0])
    assert np.all(t_unix.mask == [False, True])
    
