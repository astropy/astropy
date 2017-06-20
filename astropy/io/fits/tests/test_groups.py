# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import with_statement

import os
import time

import pytest
import numpy as np

from . import FitsTestCase
from .test_table import comparerecords
from ....io import fits


class TestGroupsFunctions(FitsTestCase):
    def test_open(self):
        with fits.open(self.data('random_groups.fits')) as hdul:
            assert isinstance(hdul[0], fits.GroupsHDU)
            naxes = (3, 1, 128, 1, 1)
            parameters = ['UU', 'VV', 'WW', 'BASELINE', 'DATE']
            info = [(0, 'PRIMARY', 1, 'GroupsHDU', 147, naxes, 'float32',
                     '3 Groups  5 Parameters')]
            assert hdul.info(output=False) == info

            ghdu = hdul[0]
            assert ghdu.parnames == parameters
            assert list(ghdu.data.dtype.names) == parameters + ['DATA']

            assert isinstance(ghdu.data, fits.GroupData)
            # The data should be equal to the number of groups
            assert ghdu.header['GCOUNT'] == len(ghdu.data)
            assert ghdu.data.data.shape == (len(ghdu.data),) + naxes[::-1]
            assert ghdu.data.parnames == parameters

            assert isinstance(ghdu.data[0], fits.Group)
            assert len(ghdu.data[0]) == len(parameters) + 1
            assert ghdu.data[0].data.shape == naxes[::-1]
            assert ghdu.data[0].parnames == parameters

    def test_open_groups_in_update_mode(self):
        """
        Test that opening a file containing a groups HDU in update mode and
        then immediately closing it does not result in any unnecessary file
        modifications.

        Similar to
        test_image.TestImageFunctions.test_open_scaled_in_update_mode().
        """

        # Copy the original file before making any possible changes to it
        self.copy_file('random_groups.fits')
        mtime = os.stat(self.temp('random_groups.fits')).st_mtime

        time.sleep(1)

        fits.open(self.temp('random_groups.fits'), mode='update',
                  memmap=False).close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp('random_groups.fits')).st_mtime

    def test_random_groups_data_update(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/3730 and
        for https://github.com/spacetelescope/PyFITS/issues/102
        """

        self.copy_file('random_groups.fits')
        with fits.open(self.temp('random_groups.fits'), mode='update') as h:
            h[0].data['UU'] = 0.42

        with fits.open(self.temp('random_groups.fits'), mode='update') as h:
            assert np.all(h[0].data['UU'] == 0.42)

    def test_parnames_round_trip(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/130

        Ensures that opening a random groups file in update mode or writing it
        to a new file does not cause any change to the parameter names.
        """

        # Because this test tries to update the random_groups.fits file, let's
        # make a copy of it first (so that the file doesn't actually get
        # modified in the off chance that the test fails
        self.copy_file('random_groups.fits')

        parameters = ['UU', 'VV', 'WW', 'BASELINE', 'DATE']
        with fits.open(self.temp('random_groups.fits'), mode='update') as h:
            assert h[0].parnames == parameters
            h.flush()
        # Open again just in read-only mode to ensure the parnames didn't
        # change
        with fits.open(self.temp('random_groups.fits')) as h:
            assert h[0].parnames == parameters
            h.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert h[0].parnames == parameters

    def test_groupdata_slice(self):
        """
        A simple test to ensure that slicing GroupData returns a new, smaller
        GroupData object, as is the case with a normal FITS_rec.  This is a
        regression test for an as-of-yet unreported issue where slicing
        GroupData returned a single Group record.
        """

        with fits.open(self.data('random_groups.fits')) as hdul:
            s = hdul[0].data[1:]
            assert isinstance(s, fits.GroupData)
            assert len(s) == 2
            assert hdul[0].data.parnames == s.parnames

    def test_group_slice(self):
        """
        Tests basic slicing a single group record.
        """

        # A very basic slice test
        with fits.open(self.data('random_groups.fits')) as hdul:
            g = hdul[0].data[0]
            s = g[2:4]
            assert len(s) == 2
            assert s[0] == g[2]
            assert s[-1] == g[-3]
            s = g[::-1]
            assert len(s) == 6
            assert (s[0] == g[-1]).all()
            assert s[-1] == g[0]
            s = g[::2]
            assert len(s) == 3
            assert s[0] == g[0]
            assert s[1] == g[2]
            assert s[2] == g[4]

    def test_create_groupdata(self):
        """
        Basic test for creating GroupData from scratch.
        """

        imdata = np.arange(100.0)
        imdata.shape = (10, 1, 1, 2, 5)
        pdata1 = np.arange(10, dtype=np.float32) + 0.1
        pdata2 = 42.0
        x = fits.hdu.groups.GroupData(imdata, parnames=['abc', 'xyz'],
                                      pardata=[pdata1, pdata2], bitpix=-32)
        assert x.parnames == ['abc', 'xyz']
        assert (x.par('abc') == pdata1).all()
        assert (x.par('xyz') == ([pdata2] * len(x))).all()
        assert (x.data == imdata).all()

        # Test putting the data into a GroupsHDU and round-tripping it
        ghdu = fits.GroupsHDU(data=x)
        ghdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            hdr = h[0].header
            assert hdr['GCOUNT'] == 10
            assert hdr['PCOUNT'] == 2
            assert hdr['NAXIS'] == 5
            assert hdr['NAXIS1'] == 0
            assert hdr['NAXIS2'] == 5
            assert hdr['NAXIS3'] == 2
            assert hdr['NAXIS4'] == 1
            assert hdr['NAXIS5'] == 1
            assert h[0].data.parnames == ['abc', 'xyz']
            assert comparerecords(h[0].data, x)

    def test_duplicate_parameter(self):
        """
        Tests support for multiple parameters of the same name, and ensures
        that the data in duplicate parameters are returned as a single summed
        value.
        """

        imdata = np.arange(100.0)
        imdata.shape = (10, 1, 1, 2, 5)
        pdata1 = np.arange(10, dtype=np.float32) + 1
        pdata2 = 42.0
        x = fits.hdu.groups.GroupData(imdata, parnames=['abc', 'xyz', 'abc'],
                                      pardata=[pdata1, pdata2, pdata1],
                                      bitpix=-32)

        assert x.parnames == ['abc', 'xyz', 'abc']
        assert (x.par('abc') == pdata1 * 2).all()
        assert x[0].par('abc') == 2

        # Test setting a parameter
        x[0].setpar(0, 2)
        assert x[0].par('abc') == 3
        pytest.raises(ValueError, x[0].setpar, 'abc', 2)
        x[0].setpar('abc', (2, 3))
        assert x[0].par('abc') == 5
        assert x.par('abc')[0] == 5
        assert (x.par('abc')[1:] == pdata1[1:] * 2).all()

        # Test round-trip
        ghdu = fits.GroupsHDU(data=x)
        ghdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            hdr = h[0].header
            assert hdr['PCOUNT'] == 3
            assert hdr['PTYPE1'] == 'abc'
            assert hdr['PTYPE2'] == 'xyz'
            assert hdr['PTYPE3'] == 'abc'
            assert x.parnames == ['abc', 'xyz', 'abc']
            assert x.dtype.names == ('abc', 'xyz', '_abc', 'DATA')
            assert x.par('abc')[0] == 5
            assert (x.par('abc')[1:] == pdata1[1:] * 2).all()
