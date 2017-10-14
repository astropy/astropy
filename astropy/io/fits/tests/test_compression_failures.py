# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys
import warnings

import pytest
import numpy as np

from ....io import fits
from ..compression import compress_hdu

from . import FitsTestCase


MAX_INT = np.iinfo(np.intc).max
MAX_LONG = np.iinfo(np.long).max
MAX_LONGLONG = np.iinfo(np.longlong).max


class TestCompressionFunction(FitsTestCase):
    def test_wrong_argument_number(self):
        with pytest.raises(TypeError):
            compress_hdu(1, 2)

    def test_unknown_compression_type(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header['ZCMPTYPE'] = 'fun'
        with pytest.raises(ValueError) as exc:
            compress_hdu(hdu)
        assert 'Unrecognized compression type: fun' in str(exc)

    def test_zbitpix_unknown(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header['ZBITPIX'] = 13
        with pytest.raises(ValueError) as exc:
            compress_hdu(hdu)
        assert 'Invalid value for BITPIX: 13' in str(exc)

    def test_data_none(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu.data = None
        with pytest.raises(TypeError) as exc:
            compress_hdu(hdu)
        assert 'CompImageHDU.data must be a numpy.ndarray' in str(exc)

    def test_missing_internal_header(self):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        del hdu._header
        with pytest.raises(AttributeError) as exc:
            compress_hdu(hdu)
        assert '_header' in str(exc)

    @pytest.mark.parametrize('kw', ['ZNAXIS', 'ZBITPIX'])
    def test_header_missing_keyword(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        del hdu._header[kw]
        with pytest.raises(KeyError) as exc:
            compress_hdu(hdu)
        assert kw in str(exc)

    @pytest.mark.parametrize('kw', ['ZNAXIS', 'ZVAL1'])
    def test_header_value_int_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = MAX_INT + 1
        with pytest.raises(OverflowError):
            compress_hdu(hdu)

    @pytest.mark.parametrize('kw', ['ZTILE1'])
    def test_header_value_long_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = MAX_LONG + 1
        with pytest.raises(OverflowError):
            compress_hdu(hdu)

    @pytest.mark.parametrize('kw', ['NAXIS1', 'NAXIS2', 'TNULL1', 'PCOUNT', 'THEAP'])
    def test_header_value_longlong_overflow(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = MAX_LONGLONG + 1
        with pytest.raises(OverflowError):
            compress_hdu(hdu)

    @pytest.mark.parametrize('kw', ['NAXIS1', 'NAXIS2', 'TFIELDS', 'PCOUNT'])
    def test_header_value_negative(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = -1
        with pytest.raises(ValueError) as exc:
            compress_hdu(hdu)
        assert '{} should not be negative.'.format(kw) in str(exc)

    @pytest.mark.parametrize(
        ('kw', 'limit'),
        [('ZNAXIS', 999),
         ('TFIELDS', 999)])
    def test_header_value_exceeds_custom_limit(self, kw, limit):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = limit + 1
        with pytest.raises(ValueError) as exc:
            compress_hdu(hdu)
        assert kw in str(exc)

    @pytest.mark.parametrize('kw', ['TTYPE1', 'TFORM1', 'ZCMPTYPE'])
    def test_header_value_no_string(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = 1
        with pytest.raises(TypeError):
            compress_hdu(hdu)

    @pytest.mark.parametrize('kw', ['TZERO1', 'TSCAL1'])
    def test_header_value_no_double(self, kw):
        hdu = fits.CompImageHDU(np.ones((10, 10)))
        hdu._header[kw] = '1'
        with pytest.raises(TypeError):
            compress_hdu(hdu)
