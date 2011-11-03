from __future__ import division

import platform

import numpy as np

import pyfits
from pyfits.tests import PyfitsTestCase

from nose.tools import assert_equal

class TestUintFunctions(PyfitsTestCase):
    def test_uint16(self):
        hdu = pyfits.PrimaryHDU(np.array([-3,-2,-1,0,1,2,3]))
        hdu.scale('int16', '', bzero=2**15)
        hdu.writeto(self.temp('tempfile.fits'))
        hdul = pyfits.open(self.temp('tempfile.fits'), uint=True)
        assert_equal(hdul[0].data.dtype, np.uint16)
        assert_equal(np.all(hdul[0].data ==
                         np.array([(2**16)-3, (2**16)-2, (2**16)-1, 0, 1, 2, 3],
                                  dtype=np.uint16)), True)
        hdul.writeto(self.temp('tempfile1.fits'))
        hdul1 = pyfits.open(self.temp('tempfile1.fits'), uint16=True)
        assert_equal(np.all(hdul[0].data == hdul1[0].data), True)
        assert_equal(hdul[0].section[:1].dtype, np.dtype('uint16'))
        hdul.close()
        hdul1.close()

    def test_uint32(self):
        hdu = pyfits.PrimaryHDU(np.array([-3, -2, -1, 0, 1, 2, 3]))
        hdu.scale('int32', '', bzero=2**31)
        hdu.writeto(self.temp('tempfile.fits'))
        hdul = pyfits.open(self.temp('tempfile.fits'), uint=True)
        assert_equal(hdul[0].data.dtype, np.uint32)
        assert_equal(np.all(hdul[0].data ==
                         np.array([(2**32)-3, (2**32)-2, (2**32)-1, 0, 1, 2, 3],
                         dtype=np.uint32)), True)
        hdul.writeto(self.temp('tempfile1.fits'))
        hdul1 = pyfits.open(self.temp('tempfile1.fits'), uint=True)
        assert_equal(np.all(hdul[0].data == hdul1[0].data), True)
        assert_equal(hdul[0].section[:1].dtype, np.dtype('uint32'))
        hdul.close()
        hdul1.close()

    def test_uint64(self):
        if platform.architecture()[0] == '64bit':
            hdu = pyfits.PrimaryHDU(np.array([-3,-2,-1,0,1,2,3]))
            hdu.scale('int64', '', bzero=2**63)
            hdu.writeto(self.temp('tempfile.fits'))
            hdul = pyfits.open(self.temp('tempfile.fits'), uint=True)
            assert_equal(hdul[0].data.dtype, np.uint64)
            assert_equal(np.all(hdul[0].data ==
                             np.array([(2**64)-3,(2**64)-2,(2**64)-1,0,1,2,3],
                             dtype=np.uint64)), True)
            hdul.writeto(self.temp('tempfile1.fits'))
            hdul1 = pyfits.open(self.temp('tempfile1.fits'),uint=True)
            assert_equal(np.all(hdul[0].data == hdul1[0].data), True)
            assert_equal(hdul[0].section[:1].dtype, np.dtype('uint64'))
            hdul.close()
            hdul1.close()
