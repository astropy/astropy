from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import shutil
import tempfile
import os
import time

import numpy as np

from ... import NDData, StdDevUncertainty, UnknownUncertainty
from ...mixins.ndio import NDIOMixin, read_from_fits, write_to_fits
from .... import units as u


# Define minimal class that uses the I/O mixin
class NDDataIO(NDIOMixin, NDData):
    pass


def test_simple_write_read(tmpdir):
    ndd = NDDataIO([1, 2, 3])
    ndd.read
    ndd.write


class TestIOFunctions(object):

    def setup_class(self):
        self.out_dir = tempfile.mkdtemp()
        self.filename = 'file{0}.fits'

    def teardown_class(self):
        tries = 3
        while tries:
            try:
                shutil.rmtree(self.out_dir)
                break
            except OSError:
                # Probably couldn't delete the file because for whatever
                # reason a handle to it is still open/hasn't been
                # garbage-collected
                time.sleep(0.5)
                tries -= 1

    def temp(self):
        """ Returns the full path to a file in the test temp dir."""
        handle, filename = tempfile.mkstemp(dir=self.out_dir)
        return filename

    def test_nddata_write_read_data_only(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd2.mask is None
        assert ndd2.uncertainty is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_mask_boolean(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]], dtype=bool))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        np.testing.assert_array_equal(ndd1.mask, ndd2.mask)
        assert ndd1.mask.dtype == ndd1.mask.dtype

        assert ndd2.uncertainty is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_mask_not_boolean(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]]))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        np.testing.assert_array_equal(ndd1.mask, ndd2.mask)
        assert ndd1.mask.dtype == ndd1.mask.dtype

        assert ndd2.uncertainty is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_uncertainty_unknown(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=np.random.random((3, 3)))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd1.uncertainty.uncertainty_type == ndd2.uncertainty.uncertainty_type
        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        # only compare kind of dtype because "float64 != >f8 "
        assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd2.mask is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_uncertainty_unknown_explicit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=UnknownUncertainty(np.random.random((3, 3))))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd1.uncertainty.uncertainty_type == ndd2.uncertainty.uncertainty_type
        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        # only compare kind of dtype because "float64 != >f8 "
        assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd2.mask is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_uncertainty_stddev_explicit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=StdDevUncertainty(np.random.random((3, 3))))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd1.uncertainty.uncertainty_type == ndd2.uncertainty.uncertainty_type
        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        # only compare kind of dtype because "float64 != >f8 "
        assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd2.mask is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_uncertainty_with_unit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=UnknownUncertainty(np.random.random((3, 3)), 'm'))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd1.uncertainty.uncertainty_type == ndd2.uncertainty.uncertainty_type
        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        # only compare kind of dtype because "float64 != >f8 "
        assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd2.mask is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_uncertainty_with_same_unit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)), unit='m',
                      uncertainty=UnknownUncertainty(np.random.random((3, 3)), 'm'))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd1.uncertainty.uncertainty_type == ndd2.uncertainty.uncertainty_type
        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        # only compare kind of dtype because "float64 != >f8 "
        assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd2.mask is None
        assert ndd2.unit == u.m

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_unit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)), unit='m')
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd2.uncertainty is None

        assert ndd2.mask is None
        assert ndd2.unit == u.m

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_unit_dimensionless(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)), unit='')
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd2.uncertainty is None

        assert ndd2.mask is None
        assert ndd2.unit == u.dimensionless_unscaled

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_complex_unit(self):
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)), unit='m^2 / s / kg')
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert ndd2.uncertainty is None

        assert ndd2.mask is None
        assert ndd2.unit == u.m * u.m / u.s / u.kg

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None
        assert ndd2.meta

    def test_nddata_write_read_meta(self):
        filename = str(self.temp())
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        ndd1 = NDData(np.ones((3, 3)), meta=meta)
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)

        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype == ndd1.data.dtype

        assert all(ndd1.meta[key] == ndd2.meta[key] for key in ndd1.meta)

        assert ndd2.uncertainty is None

        assert ndd2.mask is None
        assert ndd2.unit is None

        # A basic header is written when using fits to read or write files,
        # this will also create some WCS
        assert ndd2.wcs is not None

    def test_nddata_write_read_wcs(self):
        anotherfile = str(self.temp())
        filename = str(self.temp())
        ndd1 = NDData(np.ones((3, 3)))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)
        # Need another round trip to compare "new" wcs
        write_to_fits(ndd2, anotherfile)
        ndd3 = read_from_fits(anotherfile)

        np.testing.assert_array_equal(ndd2.data, ndd3.data)
        assert ndd2.data.dtype == ndd3.data.dtype

        assert ndd3.wcs.wcs.compare(ndd2.wcs.wcs)

        assert ndd3.meta

        assert ndd3.uncertainty is None

        assert ndd3.mask is None
        assert ndd3.unit is None

    def test_nddata_write_read_wcs_slicing(self):
        anotherfile = str(self.temp())
        filename = str(self.temp())
        ndd1 = NDData(np.ones((10, 10)))
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)
        # Need another round trip to compare "new" wcs.
        # Slice data and wcs to have something real to compare but keep header
        # and see if it is updated correctly
        ndd2tmp = NDData(ndd2.data[2:5, 4:8],
                         wcs=ndd2.wcs[2:5, 4:8], meta=ndd2.meta)
        write_to_fits(ndd2tmp, anotherfile)
        ndd3 = read_from_fits(anotherfile)

        np.testing.assert_array_equal(ndd2.data[2:5, 4:8], ndd3.data)
        assert ndd2.data.dtype == ndd3.data.dtype

        assert ndd3.wcs.wcs.compare(ndd2tmp.wcs.wcs)
        assert not ndd3.wcs.wcs.compare(ndd2.wcs.wcs)

        assert ndd3.meta != ndd2.meta
        assert ndd3.meta != ndd2tmp.meta

        assert ndd3.uncertainty is None

        assert ndd3.mask is None
        assert ndd3.unit is None

    def test_nddata_write_read_meta_wcs(self):
        filename = str(self.temp())
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        ndd1 = NDData(np.ones((10, 10)), meta=meta)
        write_to_fits(ndd1, filename)
        ndd2 = read_from_fits(filename)
        # Need another round trip to compare "new" wcs.
        # Slice data and wcs to have something real to compare but keep header
        # and see if it is updated correctly
        ndd2tmp = NDData(ndd2.data[2:5, 4:8],
                         wcs=ndd2.wcs[2:5, 4:8], meta=ndd2.meta)

        anotherfile = str(self.temp())
        write_to_fits(ndd2tmp, anotherfile)
        ndd3 = read_from_fits(anotherfile)

        np.testing.assert_array_equal(ndd2.data[2:5, 4:8], ndd3.data)
        assert ndd2.data.dtype == ndd3.data.dtype

        assert ndd3.wcs.wcs.compare(ndd2tmp.wcs.wcs)
        assert not ndd3.wcs.wcs.compare(ndd2.wcs.wcs)

        assert ndd3.meta != ndd2.meta
        assert ndd3.meta != ndd2tmp.meta
        assert all(ndd1.meta[key] == ndd2.meta[key] for key in ndd1.meta)
        assert all(ndd1.meta[key] == ndd3.meta[key] for key in ndd1.meta)

        assert ndd3.uncertainty is None

        assert ndd3.mask is None
        assert ndd3.unit is None

    def test_ndiomixin_read(self):
        filename = str(self.temp())

        data = np.ones((10, 10), )
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        unit = 'adu'
        mask = np.random.random((10, 10)) > 0.5
        uncertainty = UnknownUncertainty(np.ones((5, 5)))

        ndd1 = NDData(data, uncertainty=uncertainty, unit=unit, meta=meta, mask=mask)
        write_to_fits(ndd1, filename)
        ndd2 = NDDataIO.read(filename, format='simple_fits')

        assert isinstance(ndd2, NDDataIO)

        np.testing.assert_array_equal(ndd2.data, ndd1.data)

        np.testing.assert_array_equal(ndd2.mask, ndd1.mask)

        assert ndd1.uncertainty.__class__ == ndd2.uncertainty.__class__

        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd2.uncertainty.array)
        assert ndd1.uncertainty.unit == ndd2.uncertainty.unit

        assert ndd1.unit == ndd2.unit
        assert all(ndd1.meta[key] == ndd2.meta[key] for key in ndd1.meta)

        anotherfile = str(self.temp())
        ndd2.write(anotherfile, format='simple_fits')
        ndd3 = NDDataIO.read(anotherfile, format='simple_fits')

        assert isinstance(ndd3, NDDataIO)

        np.testing.assert_array_equal(ndd3.data, ndd1.data)

        np.testing.assert_array_equal(ndd3.mask, ndd1.mask)

        assert ndd1.uncertainty.__class__ == ndd3.uncertainty.__class__

        np.testing.assert_array_equal(ndd1.uncertainty.array, ndd3.uncertainty.array)
        assert ndd1.uncertainty.unit == ndd3.uncertainty.unit

        assert ndd1.unit == ndd3.unit
        assert all(ndd1.meta[key] == ndd3.meta[key] for key in ndd1.meta)
