# Licensed under a 3-clause BSD style license - see PYFITS.rst

import platform

import numpy as np

from ....io import fits
from . import FitsTestCase


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

    def test_uint_columns(self):
        #
        # Construct arrays
        #
        bzero16 = np.uint16(2**15)
        bzero32 = np.uint32(2**31)
        bzero64 = np.uint64(2**63)
        one64 = np.uint64(1)
        i0 = 2**np.arange(17,dtype=np.uint16) - 1
        i = np.concatenate((i0,i0[::-1],i0,i0[::-1]))[0:65]
        ii = np.array(np.array(i,dtype=np.int32) - np.int32(2**15),dtype=np.int16)
        j0 = 2**np.arange(33,dtype=np.uint32) - 1
        j = np.concatenate((j0,j0[::-1]))[0:65]
        jj = np.array(np.array(j,dtype=np.int64) - np.int64(2**31),dtype=np.int32)
        k0 = np.arange(65,dtype=np.uint64)
        k = np.zeros(k0.shape,dtype=k0.dtype)
        k[0:63] = 2**k0[0:63] - 1
        k[63] = bzero64 - one64
        k[64] = k[63] + k[63] + one64
        kk = (k - bzero64).view(np.int64)
        #
        # Construct a table from explicit columns
        #
        col1 = fits.Column(name='uint16',array=i,format='I',bzero=bzero16)
        col2 = fits.Column(name='uint32',array=j,format='J',bzero=bzero32)
        col3 = fits.Column(name='uint64',array=k,format='K',bzero=bzero64)
        table = fits.new_table([col1,col2,col3])
        assert (table.data['uint16'] == i).all()
        assert (table.data['uint32'] == j).all()
        assert (table.data['uint64'] == k).all()
        assert (table.data.base['uint16'] == ii).all()
        assert (table.data.base['uint32'] == jj).all()
        assert (table.data.base['uint64'] == kk).all()
        hdu0 = fits.PrimaryHDU()
        hdulist = fits.HDUList([hdu0,table])
        hdulist.writeto(self.temp('tempfile.fits'))
        #
        # Test write of unsigned int
        #
        del hdulist
        with fits.open(self.temp('tempfile.fits')) as hdulist2:
            hdudata = hdulist2[1].data
            assert (hdudata['uint16'] == i).all()
            assert (hdudata['uint16'].dtype == np.uint16)
            assert (hdudata['uint32'] == j).all()
            assert (hdudata['uint32'].dtype == np.uint32)
            assert (hdudata['uint64'] == k).all()
            assert (hdudata['uint64'].dtype == np.uint64)
            assert (hdudata.base['uint16'] == ii).all()
            assert (hdudata.base['uint32'] == jj).all()
            assert (hdudata.base['uint64'] == kk).all()
        #
        # Construc recarray then write out that.
        #
        uint_array = np.array(zip(i,j,k),dtype=[('uint16',np.uint16),('uint32',np.uint32),('uint64',np.uint64)])
        fits.writeto(self.temp('tempfile2.fits'),uint_array)
        with fits.open(self.temp('tempfile2.fits')) as hdulist3:
            hdudata3 = hdulist3[1].data
            assert (hdudata3.base['uint16'] == table.data.base['uint16']).all()
            assert (hdudata3['uint16'] == table.data['uint16']).all()
            assert (hdudata3['uint16'] == i).all()
            assert (hdudata3.base['uint32'] == table.data.base['uint32']).all()
            assert (hdudata3['uint32'] == table.data['uint32']).all()
            assert (hdudata3['uint32'] == j).all()
            assert (hdudata3.base['uint64'] == table.data.base['uint64']).all()
            assert (hdudata3['uint64'] == table.data['uint64']).all()
            assert (hdudata3['uint64'] == k).all()
