from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect

import numpy as np

from ...tests.helper import catch_warnings, pytest
from ... import units as u
from ...io import fits
from ...extern import six

from ..nddata import NDData
from ..io import NDIOMixin


# Define minimal class that uses the I/O mixin
class NDDataIO(NDIOMixin, NDData):
    pass


def test_simple_write_read(tmpdir):
    path = tmpdir.join('tmp.fits')
    data_in = np.array([1, 2, 3])
    ndd = NDDataIO(data_in)
    ndd.write(path.strpath)
    # Does the data read directly by io.fits match?
    with fits.open(path.strpath) as f:
        np.testing.assert_array_equal(data_in, f[0].data)
    # Does the data read by our reader match?
    ndd2 = NDDataIO.read(path.strpath)
    np.testing.assert_array_equal(data_in, ndd2.data)


def test_write_read_with_meta(tmpdir):
    path = tmpdir.join('tmp.fits')
    data_in = np.array([1, 2, 3])
    meta_in = {"one": 1, "two": 2}
    ndd = NDDataIO(data_in, meta=meta_in)
    ndd.write(path.strpath)
    # Does the data read directly by io.fits contain the meta?
    # [Note that it will also contain other metadata.]
    with fits.open(path.strpath) as f:
        for k, v in six.iteritems(meta_in):
            print(k, v)
            print(f[0].header)
            assert f[0].header[k] == v

    ndd2 = NDDataIO.read(path.strpath)
    # For now reading from a FITS file gives you a FITS header as meta.
    # Maybe it shouldn't be that way, but for now it is...
    assert isinstance(ndd2.meta, fits.Header)
