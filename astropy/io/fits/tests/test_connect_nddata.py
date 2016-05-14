from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ....nddata import NDData, StdDevUncertainty, UnknownUncertainty, NDIOMixin
from ..connect import write_nddata_fits, read_nddata_fits

import numpy as np


# Define minimal class that uses the I/O mixin
class NDDataIO(NDIOMixin, NDData):
    pass


class TestIOFunctions(object):
    counter = 0
    filename = 'file{0}.fits'

    def temp(self, tmpdir):
        self.counter += 1
        return str(tmpdir.join(self.filename.format(self.counter)))

    def compare_nddata(self, ndd1, ndd2, compare_meta=True):
        # Compare if the data is equal:
        np.testing.assert_array_equal(ndd1.data, ndd2.data)
        assert ndd1.data.dtype.kind == ndd1.data.dtype.kind

        # Compare if mask is equal
        if ndd1.mask is not None:
            np.testing.assert_array_equal(ndd1.mask, ndd2.mask)
            assert ndd1.mask.dtype.kind == ndd1.mask.dtype.kind
        else:
            assert ndd1.mask == ndd2.mask

        # Compare if uncertainty is equal
        if ndd1.uncertainty is not None:
            assert ndd1.uncertainty.__class__ == ndd2.uncertainty.__class__
            np.testing.assert_array_equal(ndd1.uncertainty.array,
                                          ndd2.uncertainty.array)
            assert ndd1.uncertainty.array.dtype.kind == ndd2.uncertainty.array.dtype.kind
            assert ndd1.uncertainty.unit == ndd2.uncertainty.unit
        else:
            assert ndd1.uncertainty == ndd2.uncertainty

        # Units equal?
        assert ndd1.unit == ndd2.unit

        # WCS equal?
        if ndd1.wcs is not None:
            assert ndd1.wcs.wcs.compare(ndd2.wcs.wcs)
        else:
            pass

        # Compare meta only if comparison makes sense and only for those
        # keywords that were present in the original (because WCS operations
        # may alter/insert meta attributes).
        if compare_meta:
            assert all(ndd1.meta[key] == ndd2.meta[key] for key in ndd1.meta)

    def test_nddata_write_read_data_only(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_data_only_with_dtype(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename, dtype=np.int32)

        assert ndd2.data.dtype == np.int32

    def test_nddata_write_read_mask_boolean(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]], dtype=bool))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_mask_not_boolean(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]]))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_uncertainty_unknown(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=np.random.random((3, 3)))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_uncertainty_unknown_explicit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=UnknownUncertainty(np.random.random((3, 3))))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_uncertainty_stddev_explicit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=StdDevUncertainty(np.random.random((3, 3))))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_uncertainty_with_unit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)),
                      uncertainty=UnknownUncertainty(np.random.random((3, 3)), 'm'))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_uncertainty_with_same_unit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)), unit='m',
                      uncertainty=UnknownUncertainty(np.random.random((3, 3)), 'm'))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_unit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)), unit='m')

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_unit_dimensionless(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)), unit='')

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_complex_unit(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)), unit='m^2 / s / kg')

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_unit_uppercase(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)))
        ndd1.meta['BUNIT'] = 'ADU'
        # ADU cannot be parsed but it can be parsed if it tries lowercase. We
        # need to change 'kw_unit' during writing though otherwise it will be
        # deleted.

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename, kw_unit='blub')
        ndd2 = read_nddata_fits(filename)

        # It couldn't be parsed so it will have no unit
        # TODO: Catch info-message here
        ndd3 = NDData(ndd1, unit=None)

        self.compare_nddata(ndd3, ndd2)

    def test_nddata_write_read_unit_deletes_keyword(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)))
        ndd1.meta['BUNIT'] = 'adu'

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)  # this should delete BUNIT keyword!
        ndd2 = read_nddata_fits(filename)

        # Since we had no unit
        ndd3 = NDData(ndd1, unit=None)
        del ndd3.meta['BUNIT']

        self.compare_nddata(ndd3, ndd2)

    def test_nddata_write_read_meta(self, tmpdir):
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        ndd1 = NDData(np.ones((3, 3)), meta=meta)

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        self.compare_nddata(ndd1, ndd2)

    def test_nddata_write_read_wcs(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)))

        # Write and read to generate basic wcs
        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        # Need another round trip to compare "new" wcs
        anotherfile = str(self.temp(tmpdir))
        write_nddata_fits(ndd2, anotherfile)
        ndd3 = read_nddata_fits(anotherfile)

        self.compare_nddata(ndd1, ndd2)
        self.compare_nddata(ndd2, ndd3)

    def test_nddata_write_read_wcs_no_toheader(self, tmpdir):
        ndd1 = NDData(np.ones((3, 3)), wcs=5)  # int has no "to_header"-method

        # Write and read to generate basic wcs
        filename = str(self.temp(tmpdir))
        # Writing would fail if the missing method wasn't catched.
        write_nddata_fits(ndd1, filename)

    def test_nddata_write_read_wcs_slicing(self, tmpdir):
        ndd1 = NDData(np.ones((10, 10)))

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        # Need another round trip to compare "new" wcs.
        # Slice data and wcs to have something real to compare but keep header
        # and see if it is updated correctly
        ndd2tmp = NDData(ndd2.data[2:5, 4:8],
                         wcs=ndd2.wcs[2:5, 4:8], meta=ndd2.meta)

        anotherfile = str(self.temp(tmpdir))
        write_nddata_fits(ndd2tmp, anotherfile)
        ndd3 = read_nddata_fits(anotherfile)

        self.compare_nddata(ndd1, ndd2)
        self.compare_nddata(ndd2tmp, ndd3, False)

        # Extra tests:
        assert not ndd3.wcs.wcs.compare(ndd2.wcs.wcs)
        assert ndd3.meta != ndd2.meta
        assert ndd3.meta != ndd2tmp.meta
        np.testing.assert_array_equal(ndd2.data[2:5, 4:8], ndd3.data)

    def test_nddata_write_read_meta_wcs(self, tmpdir):
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        ndd1 = NDData(np.ones((10, 10)), meta=meta)

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = read_nddata_fits(filename)

        ndd2tmp = NDData(ndd2.data[2:5, 4:8],
                         wcs=ndd2.wcs[2:5, 4:8], meta=ndd2.meta)

        anotherfile = str(self.temp(tmpdir))
        write_nddata_fits(ndd2tmp, anotherfile)
        ndd3 = read_nddata_fits(anotherfile)

        self.compare_nddata(ndd1, ndd2)
        self.compare_nddata(ndd2tmp, ndd3, False)

        # Extra tests:
        np.testing.assert_array_equal(ndd2.data[2:5, 4:8], ndd3.data)
        assert ndd3.wcs.wcs.compare(ndd2tmp.wcs.wcs)
        assert not ndd3.wcs.wcs.compare(ndd2.wcs.wcs)
        assert all(ndd1.meta[key] == ndd2.meta[key] for key in ndd1.meta)
        assert all(ndd1.meta[key] == ndd3.meta[key] for key in ndd1.meta)

    def test_ndiomixin_read(self, tmpdir):
        data = np.ones((10, 10))
        meta = dict([(j, i) for i, j in enumerate('abcdefghijklmnopqrstuvwxyz')])
        unit = 'adu'
        mask = np.random.random((10, 10)) > 0.5
        uncertainty = UnknownUncertainty(np.ones((5, 5)))
        ndd1 = NDData(data, uncertainty=uncertainty, unit=unit, meta=meta, mask=mask)

        filename = str(self.temp(tmpdir))
        write_nddata_fits(ndd1, filename)
        ndd2 = NDDataIO.read(filename, format='simple_fits')

        anotherfile = str(self.temp(tmpdir))
        ndd2.write(anotherfile, format='simple_fits')
        ndd3 = NDDataIO.read(anotherfile, format='simple_fits')

        self.compare_nddata(ndd1, ndd2)
        self.compare_nddata(ndd1, ndd3)

        # Extra tests:
        assert isinstance(ndd2, NDDataIO)
        assert isinstance(ndd3, NDDataIO)

    # TODO: Add one test for NDDataRef and NDDataArray!
