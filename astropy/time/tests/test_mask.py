# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
import numpy as np

from ...utils.compat import NUMPY_LT_1_14
from ...tests.helper import pytest
from .. import Time

allclose_sec = functools.partial(np.allclose, rtol=2. ** -52,
                                 atol=2. ** -52 * 24 * 3600)  # 20 ps atol
is_masked = np.ma.is_masked


def test_simple():
    t = Time([1, 2, 3], format='cxcsec')
    assert t.masked is False
    assert np.all(t.mask == [False, False, False])

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


def test_scalar_init():
    t = Time('2000:001')
    assert t.masked is False
    assert t.mask == np.array(False)


def test_mask_not_writeable():
    t = Time('2000:001')
    with pytest.raises(AttributeError) as err:
        t.mask = True
    assert "can't set attribute" in str(err)

    t = Time(['2000:001'])
    with pytest.raises(ValueError) as err:
        t.mask[0] = True
    assert "assignment destination is read-only" in str(err)


def test_str():
    t = Time(['2000:001', '2000:002'])
    t[1] = np.ma.masked
    assert str(t) == "['2000:001:00:00:00.000' --]"
    assert repr(t) == "<Time object: scale='utc' format='yday' value=['2000:001:00:00:00.000' --]>"

    if NUMPY_LT_1_14:
        expected = ["masked_array(data = ['2000-01-01 00:00:00.000' --],",
                    "             mask = [False  True],",
                    "       fill_value = N/A)"]
    else:
        expected = ["masked_array(data=['2000-01-01 00:00:00.000', --],",
                    '             mask=[False,  True],',
                    "       fill_value='N/A',",
                    "            dtype='<U23')"]
    assert repr(t.iso).splitlines() == expected

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


def test_masked_input():
    v0 = np.ma.MaskedArray([[1, 2], [3, 4]])  # No masked elements
    v1 = np.ma.MaskedArray([[1, 2], [3, 4]], mask=[[True, False], [False, False]])
    v2 = np.ma.MaskedArray([[10, 20], [30, 40]], mask=[[False, False], [False, True]])

    # Init from various combinations of masked arrays
    t = Time(v0, format='cxcsec')
    assert np.ma.allclose(t.value, v0)
    assert np.all(t.mask == [[False, False], [False, False]])
    assert t.masked is False

    t = Time(v1, format='cxcsec')
    assert np.ma.allclose(t.value, v1)
    assert np.all(t.mask == v1.mask)
    assert np.all(t.value.mask == v1.mask)
    assert t.masked is True

    t = Time(v1, v2, format='cxcsec')
    assert np.ma.allclose(t.value, v1 + v2)
    assert np.all(t.mask == (v1 + v2).mask)
    assert t.masked is True

    t = Time(v0, v1, format='cxcsec')
    assert np.ma.allclose(t.value, v0 + v1)
    assert np.all(t.mask == (v0 + v1).mask)
    assert t.masked is True

    t = Time(0, v2, format='cxcsec')
    assert np.ma.allclose(t.value, v2)
    assert np.all(t.mask == v2.mask)
    assert t.masked is True

    # Init from a string masked array
    t_iso = t.iso
    t2 = Time(t_iso)
    assert np.all(t2.value == t_iso)
    assert np.all(t2.mask == v2.mask)
    assert t2.masked is True
