# Licensed under a 3-clause BSD style license - see PYFITS.rst

import platform

import numpy as np

from ....io import fits
from . import FitsTestCase
from ....tests.helper import pytest


class TestUintFunctions(FitsTestCase):
    def test_uint16(self):
        hdu = fits.PrimaryHDU(np.array([-3, -2, -1, 0, 1, 2, 3]))
        hdu.scale('int16', '', bzero=2 ** 15)
        hdu.writeto(self.temp('tempfile.fits'))
        hdul = fits.open(self.temp('tempfile.fits'), uint=True)
        assert hdul[0].data.dtype == np.uint16
        assert np.all(hdul[0].data ==
                      np.array([(2 ** 16) - 3, (2 ** 16) - 2, (2 ** 16) - 1,
                                0, 1, 2, 3],
                               dtype=np.uint16))
        hdul.writeto(self.temp('tempfile1.fits'))
        hdul1 = fits.open(self.temp('tempfile1.fits'), uint16=True)
        assert (hdul[0].data == hdul1[0].data).all()
        assert hdul[0].section[:1].dtype.name == 'uint16'
        assert (hdul[0].section[:1] == hdul[0].data[:1]).all()
        hdul.close()
        hdul1.close()

    def test_uint32(self):
        hdu = fits.PrimaryHDU(np.array([-3, -2, -1, 0, 1, 2, 3]))
        hdu.scale('int32', '', bzero=2 ** 31)
        hdu.writeto(self.temp('tempfile.fits'))
        hdul = fits.open(self.temp('tempfile.fits'), uint=True)
        assert hdul[0].data.dtype == np.uint32
        assert np.all(hdul[0].data ==
                      np.array([(2 ** 32) - 3, (2 ** 32) - 2, (2 ** 32) - 1,
                                0, 1, 2, 3], dtype=np.uint32))

        hdul.writeto(self.temp('tempfile1.fits'))
        hdul1 = fits.open(self.temp('tempfile1.fits'), uint=True)

        assert (hdul[0].data == hdul1[0].data).all()
        assert hdul[0].section[:1].dtype.name == 'uint32'
        assert (hdul[0].section[:1] == hdul[0].data[:1]).all()
        hdul.close()
        hdul1.close()

    def test_uint64(self):
        if platform.architecture()[0] == '64bit':
            hdu = fits.PrimaryHDU(np.array([-3, -2, -1, 0, 1, 2, 3]))
            hdu.scale('int64', '', bzero=2 ** 63)
            hdu.writeto(self.temp('tempfile.fits'))
            hdul = fits.open(self.temp('tempfile.fits'), uint=True)
            assert hdul[0].data.dtype == np.uint64
            assert np.all(hdul[0].data ==
                          np.array([(2 ** 64) - 3, (2 ** 64) - 2,
                                    (2 ** 64) - 1, 0, 1, 2, 3],
                                   dtype=np.uint64))
            hdul.writeto(self.temp('tempfile1.fits'))
            hdul1 = fits.open(self.temp('tempfile1.fits'), uint=True)
            assert (hdul[0].data == hdul1[0].data).all()
            assert hdul[0].section[:1].dtype.name == 'uint64'
            assert (hdul[0].section[:1] == hdul[0].data[:1]).all()
            hdul.close()
            hdul1.close()

    def test_uint_compressed(self):
        hdu = fits.CompImageHDU(np.array([-3, -2, -1, 0, 1, 2, 3]))
        hdu.scale('int32', '', bzero=2 ** 31)
        hdu.writeto(self.temp('temp.fits'))
        with fits.open(self.temp('temp.fits'), uint=True) as hdul:
            assert hdul[1].data.dtype == np.uint32
            assert (hdul[1].data ==
                    np.array([(2 ** 32) - 3, (2 ** 32) - 2, (2 ** 32) - 1, 0,
                              1, 2, 3], dtype=np.uint32)).all()
            hdul.writeto(self.temp('temp2.fits'))
            with fits.open(self.temp('temp2.fits'), uint=True) as hdul2:
                assert (hdul[1].data == hdul2[1].data).all()
                # TODO: Enable these lines if CompImageHDUs ever grow .section
                # support
                #assert_equal(hdul[1].section[:1].dtype.name, 'uint32')
                #assert_true((hdul[1].section[:1] == hdul[1].data[:1]).all())

    @pytest.mark.parametrize(('dtype',),[('u2',),('u4',),('u8',)])
    def test_uint_columns(self,dtype):
        #
        # Construct array
        #
        dtype_map = {'u2':np.uint16,'u4':np.uint32,'u8':np.uint64}
        itype_map = {'u2':np.int16,'u4':np.int32,'u8':np.int64}
        format_map = {'u2':'I','u4':'J','u8':'K'}
        bits = 8*int(dtype[1])
        bzero = dtype_map[dtype](2**(bits-1))
        one = dtype_map[dtype](1)
        u0 = np.arange(bits+1,dtype=dtype_map[dtype])
        u = 2**u0 - one
        if bits == 64:
            u[63] = bzero - one
            u[64] = u[63] + u[63] + one
        uu = (u - bzero).view(itype_map[dtype])
        #
        # Construct a table from explicit column
        #
        col = fits.Column(name=dtype,array=u,format=format_map[dtype],bzero=bzero)
        table = fits.new_table([col])
        assert (table.data[dtype] == u).all()
        assert (table.data.base[dtype] == uu).all()
        hdu0 = fits.PrimaryHDU()
        hdulist = fits.HDUList([hdu0,table])
        hdulist.writeto(self.temp('tempfile.fits'))
        #
        # Test write of unsigned int
        #
        del hdulist
        with fits.open(self.temp('tempfile.fits')) as hdulist2:
            hdudata = hdulist2[1].data
            assert (hdudata[dtype] == u).all()
            assert (hdudata[dtype].dtype == dtype_map[dtype])
            assert (hdudata.base[dtype] == uu).all()
        #
        # Construct recarray then write out that.
        #
        v = u.view(dtype=[(dtype,dtype_map[dtype])])
        fits.writeto(self.temp('tempfile2.fits'),v)
        with fits.open(self.temp('tempfile2.fits')) as hdulist3:
            hdudata3 = hdulist3[1].data
            assert (hdudata3.base[dtype] == table.data.base[dtype]).all()
            assert (hdudata3[dtype] == table.data[dtype]).all()
            assert (hdudata3[dtype] == u).all()
