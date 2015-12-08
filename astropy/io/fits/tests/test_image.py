# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division, with_statement

import math
import os
import time
import warnings

import numpy as np

from ....io import fits
from ....utils.exceptions import (AstropyDeprecationWarning,
                                  AstropyPendingDeprecationWarning)
from ....tests.helper import pytest, raises, catch_warnings, ignore_warnings
from ..hdu.compressed import SUBTRACTIVE_DITHER_1, DITHER_SEED_CHECKSUM
from .test_table import comparerecords

from . import FitsTestCase


class TestImageFunctions(FitsTestCase):
    def test_constructor_name_arg(self):
        """Like the test of the same name in test_table.py"""

        hdu = fits.ImageHDU()
        assert hdu.name == ''
        assert 'EXTNAME' not in hdu.header
        hdu.name = 'FOO'
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

        # Passing name to constructor
        hdu = fits.ImageHDU(name='FOO')
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

        # And overriding a header with a different extname
        hdr = fits.Header()
        hdr['EXTNAME'] = 'EVENTS'
        hdu = fits.ImageHDU(header=hdr, name='FOO')
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

    def test_constructor_copies_header(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/153

        Ensure that a header from one HDU is copied when used to initialize new
        HDU.
        """

        ifd = fits.HDUList(fits.PrimaryHDU())
        phdr = ifd[0].header
        phdr['FILENAME'] = 'labq01i3q_rawtag.fits'

        primary_hdu = fits.PrimaryHDU(header=phdr)
        ofd = fits.HDUList(primary_hdu)
        ofd[0].header['FILENAME'] = 'labq01i3q_flt.fits'

        # Original header should be unchanged
        assert phdr['FILENAME'] == 'labq01i3q_rawtag.fits'

    @raises(ValueError)
    def test_open(self):
        # The function "open" reads a FITS file into an HDUList object.  There
        # are three modes to open: "readonly" (the default), "append", and
        # "update".

        # Open a file read-only (the default mode), the content of the FITS
        # file are read into memory.
        r = fits.open(self.data('test0.fits'))  # readonly

        # data parts are latent instantiation, so if we close the HDUList
        # without touching data, data can not be accessed.
        r.close()
        r[1].data[:2, :2]

    def test_open_2(self):
        r = fits.open(self.data('test0.fits'))

        info = ([(0, 'PRIMARY', 'PrimaryHDU', 138, (), '', '')] +
                [(x, 'SCI', 'ImageHDU', 61, (40, 40), 'int16', '')
                 for x in range(1, 5)])

        try:
            assert r.info(output=False) == info
        finally:
            r.close()

    def test_primary_with_extname(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/151

        Tests that the EXTNAME keyword works with Primary HDUs as well, and
        interacts properly with the .name attribute.  For convenience
        hdulist['PRIMARY'] will still refer to the first HDU even if it has an
        EXTNAME not equal to 'PRIMARY'.
        """

        prihdr = fits.Header([('EXTNAME', 'XPRIMARY'), ('EXTVER', 1)])
        hdul = fits.HDUList([fits.PrimaryHDU(header=prihdr)])
        assert 'EXTNAME' in hdul[0].header
        assert hdul[0].name == 'XPRIMARY'
        assert hdul[0].name == hdul[0].header['EXTNAME']

        info = [(0, 'XPRIMARY', 'PrimaryHDU', 5, (), '', '')]
        assert hdul.info(output=False) == info

        assert hdul['PRIMARY'] is hdul['XPRIMARY']
        assert hdul['PRIMARY'] is hdul[('XPRIMARY', 1)]

        hdul[0].name = 'XPRIMARY2'
        assert hdul[0].header['EXTNAME'] == 'XPRIMARY2'

        hdul.writeto(self.temp('test.fits'))
        with fits.open(self.temp('test.fits')) as hdul:
            assert hdul[0].name == 'XPRIMARY2'

    @ignore_warnings(AstropyDeprecationWarning)
    def test_io_manipulation(self):
        # This legacy test also tests numerous deprecated interfaces for
        # backwards compatibility

        # Get a keyword value.  An extension can be referred by name or by
        # number.  Both extension and keyword names are case insensitive.
        with fits.open(self.data('test0.fits')) as r:
            assert r['primary'].header['naxis'] == 0
            assert r[0].header['naxis'] == 0

            # If there are more than one extension with the same EXTNAME value,
            # the EXTVER can be used (as the second argument) to distinguish
            # the extension.
            assert r['sci', 1].header['detector'] == 1

            # append (using "update()") a new card
            r[0].header['xxx'] = 1.234e56

            assert (str(r[0].header.ascard[-3:]) ==
                    "EXPFLAG = 'NORMAL            ' / Exposure interruption indicator                \n"
                    "FILENAME= 'vtest3.fits'        / File name                                      \n"
                    "XXX     =            1.234E+56                                                  ")

            # rename a keyword
            r[0].header.rename_key('filename', 'fname')
            pytest.raises(ValueError, r[0].header.rename_key, 'fname',
                          'history')

            pytest.raises(ValueError, r[0].header.rename_key, 'fname',
                          'simple')
            r[0].header.rename_key('fname', 'filename')

            # get a subsection of data
            assert np.array_equal(r[2].data[:3, :3],
                                  np.array([[349, 349, 348],
                                            [349, 349, 347],
                                            [347, 350, 349]], dtype=np.int16))

            # We can create a new FITS file by opening a new file with "append"
            # mode.
            with fits.open(self.temp('test_new.fits'), mode='append') as n:
                # Append the primary header and the 2nd extension to the new
                # file.
                n.append(r[0])
                n.append(r[2])

                # The flush method will write the current HDUList object back
                # to the newly created file on disk.  The HDUList is still open
                # and can be further operated.
                n.flush()
                assert n[1].data[1, 1] == 349

                # modify a data point
                n[1].data[1, 1] = 99

                # When the file is closed, the most recent additions of
                # extension(s) since last flush() will be appended, but any HDU
                # already existed at the last flush will not be modified
            del n

            # If an existing file is opened with "append" mode, like the
            # readonly mode, the HDU's will be read into the HDUList which can
            # be modified in memory but can not be written back to the original
            # file.  A file opened with append mode can only add new HDU's.
            os.rename(self.temp('test_new.fits'),
                      self.temp('test_append.fits'))

            with fits.open(self.temp('test_append.fits'), mode='append') as a:

                # The above change did not take effect since this was made
                # after the flush().
                assert a[1].data[1, 1] == 349
                a.append(r[1])
            del a

            # When changes are made to an HDUList which was opened with
            # "update" mode, they will be written back to the original file
            # when a flush/close is called.
            os.rename(self.temp('test_append.fits'),
                      self.temp('test_update.fits'))

            with fits.open(self.temp('test_update.fits'), mode='update') as u:

                # When the changes do not alter the size structures of the
                # original (or since last flush) HDUList, the changes are
                # written back "in place".
                assert u[0].header['rootname'] == 'U2EQ0201T'
                u[0].header['rootname'] = 'abc'
                assert u[1].data[1, 1] == 349
                u[1].data[1, 1] = 99
                u.flush()

                # If the changes affect the size structure, e.g. adding or
                # deleting HDU(s), header was expanded or reduced beyond
                # existing number of blocks (2880 bytes in each block), or
                # change the data size, the HDUList is written to a temporary
                # file, the original file is deleted, and the temporary file is
                # renamed to the original file name and reopened in the update
                # mode.  To a user, these two kinds of updating writeback seem
                # to be the same, unless the optional argument in flush or
                # close is set to 1.
                del u[2]
                u.flush()

                # the write method in HDUList class writes the current HDUList,
                # with all changes made up to now, to a new file.  This method
                # works the same disregard the mode the HDUList was opened
                # with.
                u.append(r[3])
                u.writeto(self.temp('test_new.fits'))
            del u

        # Another useful new HDUList method is readall.  It will "touch" the
        # data parts in all HDUs, so even if the HDUList is closed, we can
        # still operate on the data.
        with fits.open(self.data('test0.fits')) as r:
            r.readall()
            assert r[1].data[1, 1] == 315

        # create an HDU with data only
        data = np.ones((3, 5), dtype=np.float32)
        hdu = fits.ImageHDU(data=data, name='SCI')
        assert np.array_equal(hdu.data,
                              np.array([[1.,  1.,  1.,  1.,  1.],
                                        [1.,  1.,  1.,  1.,  1.],
                                        [1.,  1.,  1.,  1.,  1.]],
                                       dtype=np.float32))

        # create an HDU with header and data
        # notice that the header has the right NAXIS's since it is constructed
        # with ImageHDU
        hdu2 = fits.ImageHDU(header=r[1].header, data=np.array([1, 2],
                             dtype='int32'))

        assert (str(hdu2.header.ascard[1:5]) ==
            "BITPIX  =                   32 / array data type                                \n"
            "NAXIS   =                    1 / number of array dimensions                     \n"
            "NAXIS1  =                    2                                                  \n"
            "PCOUNT  =                    0 / number of parameters                           ")

    def test_memory_mapping(self):
        # memory mapping
        f1 = fits.open(self.data('test0.fits'), memmap=1)
        f1.close()

    def test_verification_on_output(self):
        # verification on output
        # make a defect HDUList first
        x = fits.ImageHDU()
        hdu = fits.HDUList(x)  # HDUList can take a list or one single HDU
        with catch_warnings() as w:
            hdu.verify()
        text = "HDUList's 0th element is not a primary HDU."
        assert len(w) == 3
        assert text in str(w[1].message)

        with catch_warnings() as w:
            hdu.writeto(self.temp('test_new2.fits'), 'fix')
        text = ("HDUList's 0th element is not a primary HDU.  "
                "Fixed by inserting one as 0th HDU.")
        assert len(w) == 3
        assert text in str(w[1].message)

    def test_section(self):
        # section testing
        fs = fits.open(self.data('arange.fits'))
        assert np.array_equal(fs[0].section[3, 2, 5], 357)
        assert np.array_equal(
            fs[0].section[3, 2, :],
            np.array([352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362]))
        assert np.array_equal(fs[0].section[3, 2, 4:],
                              np.array([356, 357, 358, 359, 360, 361, 362]))
        assert np.array_equal(fs[0].section[3, 2, :8],
                              np.array([352, 353, 354, 355, 356, 357, 358, 359]))
        assert np.array_equal(fs[0].section[3, 2, -8:8],
                              np.array([355, 356, 357, 358, 359]))
        assert np.array_equal(
            fs[0].section[3, 2:5, :],
            np.array([[352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362],
                      [363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373],
                      [374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384]]))

        assert np.array_equal(fs[0].section[3, :, :][:3, :3],
                              np.array([[330, 331, 332],
                                        [341, 342, 343],
                                        [352, 353, 354]]))

        dat = fs[0].data
        assert np.array_equal(fs[0].section[3, 2:5, :8], dat[3, 2:5, :8])
        assert np.array_equal(fs[0].section[3, 2:5, 3], dat[3, 2:5, 3])

        assert np.array_equal(fs[0].section[3:6, :, :][:3, :3, :3],
                              np.array([[[330, 331, 332],
                                         [341, 342, 343],
                                         [352, 353, 354]],
                                        [[440, 441, 442],
                                         [451, 452, 453],
                                         [462, 463, 464]],
                                        [[550, 551, 552],
                                         [561, 562, 563],
                                         [572, 573, 574]]]))

        assert np.array_equal(fs[0].section[:, :, :][:3, :2, :2],
                              np.array([[[0,   1],
                                         [11,  12]],
                                        [[110, 111],
                                         [121, 122]],
                                        [[220, 221],
                                         [231, 232]]]))

        assert np.array_equal(fs[0].section[:, 2, :], dat[:, 2, :])
        assert np.array_equal(fs[0].section[:, 2:5, :], dat[:, 2:5, :])
        assert np.array_equal(fs[0].section[3:6, 3, :], dat[3:6, 3, :])
        assert np.array_equal(fs[0].section[3:6, 3:7, :], dat[3:6, 3:7, :])

        assert np.array_equal(fs[0].section[:, ::2], dat[:, ::2])
        assert np.array_equal(fs[0].section[:, [1, 2, 4], 3],
                              dat[:, [1, 2, 4], 3])
        assert np.array_equal(
            fs[0].section[:, np.array([True, False, True]), :],
            dat[:, np.array([True, False, True]), :])

        assert np.array_equal(
            fs[0].section[3:6, 3, :, ...], dat[3:6, 3, :, ...])
        assert np.array_equal(fs[0].section[..., ::2], dat[..., ::2])
        assert np.array_equal(fs[0].section[..., [1, 2, 4], 3],
                              dat[..., [1, 2, 4], 3])

    def test_section_data_single(self):
        a = np.array([1])
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        sec = hdul[0].section
        dat = hdul[0].data
        assert np.array_equal(sec[0], dat[0])
        assert np.array_equal(sec[...], dat[...])
        assert np.array_equal(sec[..., 0], dat[..., 0])
        assert np.array_equal(sec[0, ...], dat[0, ...])

    def test_section_data_square(self):
        a = np.arange(4).reshape((2, 2))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert (d.section[:, :] == dat[:, :]).all()
        assert (d.section[0, :] == dat[0, :]).all()
        assert (d.section[1, :] == dat[1, :]).all()
        assert (d.section[:, 0] == dat[:, 0]).all()
        assert (d.section[:, 1] == dat[:, 1]).all()
        assert (d.section[0, 0] == dat[0, 0]).all()
        assert (d.section[0, 1] == dat[0, 1]).all()
        assert (d.section[1, 0] == dat[1, 0]).all()
        assert (d.section[1, 1] == dat[1, 1]).all()
        assert (d.section[0:1, 0:1] == dat[0:1, 0:1]).all()
        assert (d.section[0:2, 0:1] == dat[0:2, 0:1]).all()
        assert (d.section[0:1, 0:2] == dat[0:1, 0:2]).all()
        assert (d.section[0:2, 0:2] == dat[0:2, 0:2]).all()

    def test_section_data_cube(self):
        a = np.arange(18).reshape((2, 3, 3))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data

        # TODO: Generate these perumtions instead of having them all written
        # out, yeesh!
        assert (d.section[:, :, :] == dat[:, :, :]).all()
        assert (d.section[:, :] == dat[:, :]).all()
        assert (d.section[:] == dat[:]).all()
        assert (d.section[0, :, :] == dat[0, :, :]).all()
        assert (d.section[1, :, :] == dat[1, :, :]).all()
        assert (d.section[0, 0, :] == dat[0, 0, :]).all()
        assert (d.section[0, 1, :] == dat[0, 1, :]).all()
        assert (d.section[0, 2, :] == dat[0, 2, :]).all()
        assert (d.section[1, 0, :] == dat[1, 0, :]).all()
        assert (d.section[1, 1, :] == dat[1, 1, :]).all()
        assert (d.section[1, 2, :] == dat[1, 2, :]).all()
        assert (d.section[0, 0, 0] == dat[0, 0, 0]).all()
        assert (d.section[0, 0, 1] == dat[0, 0, 1]).all()
        assert (d.section[0, 0, 2] == dat[0, 0, 2]).all()
        assert (d.section[0, 1, 0] == dat[0, 1, 0]).all()
        assert (d.section[0, 1, 1] == dat[0, 1, 1]).all()
        assert (d.section[0, 1, 2] == dat[0, 1, 2]).all()
        assert (d.section[0, 2, 0] == dat[0, 2, 0]).all()
        assert (d.section[0, 2, 1] == dat[0, 2, 1]).all()
        assert (d.section[0, 2, 2] == dat[0, 2, 2]).all()
        assert (d.section[1, 0, 0] == dat[1, 0, 0]).all()
        assert (d.section[1, 0, 1] == dat[1, 0, 1]).all()
        assert (d.section[1, 0, 2] == dat[1, 0, 2]).all()
        assert (d.section[1, 1, 0] == dat[1, 1, 0]).all()
        assert (d.section[1, 1, 1] == dat[1, 1, 1]).all()
        assert (d.section[1, 1, 2] == dat[1, 1, 2]).all()
        assert (d.section[1, 2, 0] == dat[1, 2, 0]).all()
        assert (d.section[1, 2, 1] == dat[1, 2, 1]).all()
        assert (d.section[1, 2, 2] == dat[1, 2, 2]).all()
        assert (d.section[:, 0, 0] == dat[:, 0, 0]).all()
        assert (d.section[:, 0, 1] == dat[:, 0, 1]).all()
        assert (d.section[:, 0, 2] == dat[:, 0, 2]).all()
        assert (d.section[:, 1, 0] == dat[:, 1, 0]).all()
        assert (d.section[:, 1, 1] == dat[:, 1, 1]).all()
        assert (d.section[:, 1, 2] == dat[:, 1, 2]).all()
        assert (d.section[:, 2, 0] == dat[:, 2, 0]).all()
        assert (d.section[:, 2, 1] == dat[:, 2, 1]).all()
        assert (d.section[:, 2, 2] == dat[:, 2, 2]).all()
        assert (d.section[0, :, 0] == dat[0, :, 0]).all()
        assert (d.section[0, :, 1] == dat[0, :, 1]).all()
        assert (d.section[0, :, 2] == dat[0, :, 2]).all()
        assert (d.section[1, :, 0] == dat[1, :, 0]).all()
        assert (d.section[1, :, 1] == dat[1, :, 1]).all()
        assert (d.section[1, :, 2] == dat[1, :, 2]).all()
        assert (d.section[:, :, 0] == dat[:, :, 0]).all()
        assert (d.section[:, :, 1] == dat[:, :, 1]).all()
        assert (d.section[:, :, 2] == dat[:, :, 2]).all()
        assert (d.section[:, 0, :] == dat[:, 0, :]).all()
        assert (d.section[:, 1, :] == dat[:, 1, :]).all()
        assert (d.section[:, 2, :] == dat[:, 2, :]).all()

        assert (d.section[:, :, 0:1] == dat[:, :, 0:1]).all()
        assert (d.section[:, :, 0:2] == dat[:, :, 0:2]).all()
        assert (d.section[:, :, 0:3] == dat[:, :, 0:3]).all()
        assert (d.section[:, :, 1:2] == dat[:, :, 1:2]).all()
        assert (d.section[:, :, 1:3] == dat[:, :, 1:3]).all()
        assert (d.section[:, :, 2:3] == dat[:, :, 2:3]).all()
        assert (d.section[0:1, 0:1, 0:1] == dat[0:1, 0:1, 0:1]).all()
        assert (d.section[0:1, 0:1, 0:2] == dat[0:1, 0:1, 0:2]).all()
        assert (d.section[0:1, 0:1, 0:3] == dat[0:1, 0:1, 0:3]).all()
        assert (d.section[0:1, 0:1, 1:2] == dat[0:1, 0:1, 1:2]).all()
        assert (d.section[0:1, 0:1, 1:3] == dat[0:1, 0:1, 1:3]).all()
        assert (d.section[0:1, 0:1, 2:3] == dat[0:1, 0:1, 2:3]).all()
        assert (d.section[0:1, 0:2, 0:1] == dat[0:1, 0:2, 0:1]).all()
        assert (d.section[0:1, 0:2, 0:2] == dat[0:1, 0:2, 0:2]).all()
        assert (d.section[0:1, 0:2, 0:3] == dat[0:1, 0:2, 0:3]).all()
        assert (d.section[0:1, 0:2, 1:2] == dat[0:1, 0:2, 1:2]).all()
        assert (d.section[0:1, 0:2, 1:3] == dat[0:1, 0:2, 1:3]).all()
        assert (d.section[0:1, 0:2, 2:3] == dat[0:1, 0:2, 2:3]).all()
        assert (d.section[0:1, 0:3, 0:1] == dat[0:1, 0:3, 0:1]).all()
        assert (d.section[0:1, 0:3, 0:2] == dat[0:1, 0:3, 0:2]).all()
        assert (d.section[0:1, 0:3, 0:3] == dat[0:1, 0:3, 0:3]).all()
        assert (d.section[0:1, 0:3, 1:2] == dat[0:1, 0:3, 1:2]).all()
        assert (d.section[0:1, 0:3, 1:3] == dat[0:1, 0:3, 1:3]).all()
        assert (d.section[0:1, 0:3, 2:3] == dat[0:1, 0:3, 2:3]).all()
        assert (d.section[0:1, 1:2, 0:1] == dat[0:1, 1:2, 0:1]).all()
        assert (d.section[0:1, 1:2, 0:2] == dat[0:1, 1:2, 0:2]).all()
        assert (d.section[0:1, 1:2, 0:3] == dat[0:1, 1:2, 0:3]).all()
        assert (d.section[0:1, 1:2, 1:2] == dat[0:1, 1:2, 1:2]).all()
        assert (d.section[0:1, 1:2, 1:3] == dat[0:1, 1:2, 1:3]).all()
        assert (d.section[0:1, 1:2, 2:3] == dat[0:1, 1:2, 2:3]).all()
        assert (d.section[0:1, 1:3, 0:1] == dat[0:1, 1:3, 0:1]).all()
        assert (d.section[0:1, 1:3, 0:2] == dat[0:1, 1:3, 0:2]).all()
        assert (d.section[0:1, 1:3, 0:3] == dat[0:1, 1:3, 0:3]).all()
        assert (d.section[0:1, 1:3, 1:2] == dat[0:1, 1:3, 1:2]).all()
        assert (d.section[0:1, 1:3, 1:3] == dat[0:1, 1:3, 1:3]).all()
        assert (d.section[0:1, 1:3, 2:3] == dat[0:1, 1:3, 2:3]).all()
        assert (d.section[1:2, 0:1, 0:1] == dat[1:2, 0:1, 0:1]).all()
        assert (d.section[1:2, 0:1, 0:2] == dat[1:2, 0:1, 0:2]).all()
        assert (d.section[1:2, 0:1, 0:3] == dat[1:2, 0:1, 0:3]).all()
        assert (d.section[1:2, 0:1, 1:2] == dat[1:2, 0:1, 1:2]).all()
        assert (d.section[1:2, 0:1, 1:3] == dat[1:2, 0:1, 1:3]).all()
        assert (d.section[1:2, 0:1, 2:3] == dat[1:2, 0:1, 2:3]).all()
        assert (d.section[1:2, 0:2, 0:1] == dat[1:2, 0:2, 0:1]).all()
        assert (d.section[1:2, 0:2, 0:2] == dat[1:2, 0:2, 0:2]).all()
        assert (d.section[1:2, 0:2, 0:3] == dat[1:2, 0:2, 0:3]).all()
        assert (d.section[1:2, 0:2, 1:2] == dat[1:2, 0:2, 1:2]).all()
        assert (d.section[1:2, 0:2, 1:3] == dat[1:2, 0:2, 1:3]).all()
        assert (d.section[1:2, 0:2, 2:3] == dat[1:2, 0:2, 2:3]).all()
        assert (d.section[1:2, 0:3, 0:1] == dat[1:2, 0:3, 0:1]).all()
        assert (d.section[1:2, 0:3, 0:2] == dat[1:2, 0:3, 0:2]).all()
        assert (d.section[1:2, 0:3, 0:3] == dat[1:2, 0:3, 0:3]).all()
        assert (d.section[1:2, 0:3, 1:2] == dat[1:2, 0:3, 1:2]).all()
        assert (d.section[1:2, 0:3, 1:3] == dat[1:2, 0:3, 1:3]).all()
        assert (d.section[1:2, 0:3, 2:3] == dat[1:2, 0:3, 2:3]).all()
        assert (d.section[1:2, 1:2, 0:1] == dat[1:2, 1:2, 0:1]).all()
        assert (d.section[1:2, 1:2, 0:2] == dat[1:2, 1:2, 0:2]).all()
        assert (d.section[1:2, 1:2, 0:3] == dat[1:2, 1:2, 0:3]).all()
        assert (d.section[1:2, 1:2, 1:2] == dat[1:2, 1:2, 1:2]).all()
        assert (d.section[1:2, 1:2, 1:3] == dat[1:2, 1:2, 1:3]).all()
        assert (d.section[1:2, 1:2, 2:3] == dat[1:2, 1:2, 2:3]).all()
        assert (d.section[1:2, 1:3, 0:1] == dat[1:2, 1:3, 0:1]).all()
        assert (d.section[1:2, 1:3, 0:2] == dat[1:2, 1:3, 0:2]).all()
        assert (d.section[1:2, 1:3, 0:3] == dat[1:2, 1:3, 0:3]).all()
        assert (d.section[1:2, 1:3, 1:2] == dat[1:2, 1:3, 1:2]).all()
        assert (d.section[1:2, 1:3, 1:3] == dat[1:2, 1:3, 1:3]).all()
        assert (d.section[1:2, 1:3, 2:3] == dat[1:2, 1:3, 2:3]).all()

    def test_section_data_four(self):
        a = np.arange(256).reshape((4, 4, 4, 4))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert (d.section[:, :, :, :] == dat[:, :, :, :]).all()
        assert (d.section[:, :, :] == dat[:, :, :]).all()
        assert (d.section[:, :] == dat[:, :]).all()
        assert (d.section[:] == dat[:]).all()
        assert (d.section[0, :, :, :] == dat[0, :, :, :]).all()
        assert (d.section[0, :, 0, :] == dat[0, :, 0, :]).all()
        assert (d.section[:, :, 0, :] == dat[:, :, 0, :]).all()
        assert (d.section[:, 1, 0, :] == dat[:, 1, 0, :]).all()
        assert (d.section[:, :, :, 1] == dat[:, :, :, 1]).all()

    def test_section_data_scaled(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/143

        This is like test_section_data_square but uses a file containing scaled
        image data, to test that sections can work correctly with scaled data.
        """

        hdul = fits.open(self.data('scale.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert (d.section[:, :] == dat[:, :]).all()
        assert (d.section[0, :] == dat[0, :]).all()
        assert (d.section[1, :] == dat[1, :]).all()
        assert (d.section[:, 0] == dat[:, 0]).all()
        assert (d.section[:, 1] == dat[:, 1]).all()
        assert (d.section[0, 0] == dat[0, 0]).all()
        assert (d.section[0, 1] == dat[0, 1]).all()
        assert (d.section[1, 0] == dat[1, 0]).all()
        assert (d.section[1, 1] == dat[1, 1]).all()
        assert (d.section[0:1, 0:1] == dat[0:1, 0:1]).all()
        assert (d.section[0:2, 0:1] == dat[0:2, 0:1]).all()
        assert (d.section[0:1, 0:2] == dat[0:1, 0:2]).all()
        assert (d.section[0:2, 0:2] == dat[0:2, 0:2]).all()

        # Test without having accessed the full data first
        hdul = fits.open(self.data('scale.fits'))
        d = hdul[0]
        assert (d.section[:, :] == dat[:, :]).all()
        assert (d.section[0, :] == dat[0, :]).all()
        assert (d.section[1, :] == dat[1, :]).all()
        assert (d.section[:, 0] == dat[:, 0]).all()
        assert (d.section[:, 1] == dat[:, 1]).all()
        assert (d.section[0, 0] == dat[0, 0]).all()
        assert (d.section[0, 1] == dat[0, 1]).all()
        assert (d.section[1, 0] == dat[1, 0]).all()
        assert (d.section[1, 1] == dat[1, 1]).all()
        assert (d.section[0:1, 0:1] == dat[0:1, 0:1]).all()
        assert (d.section[0:2, 0:1] == dat[0:2, 0:1]).all()
        assert (d.section[0:1, 0:2] == dat[0:1, 0:2]).all()
        assert (d.section[0:2, 0:2] == dat[0:2, 0:2]).all()
        assert not d._data_loaded

    def test_do_not_scale_image_data(self):
        hdul = fits.open(self.data('scale.fits'), do_not_scale_image_data=True)
        assert hdul[0].data.dtype == np.dtype('>i2')
        hdul = fits.open(self.data('scale.fits'))
        assert hdul[0].data.dtype == np.dtype('float32')

    def test_append_uint_data(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/56
        (BZERO and BSCALE added in the wrong location when appending scaled
        data)
        """

        fits.writeto(self.temp('test_new.fits'), data=np.array([],
                     dtype='uint8'))
        d = np.zeros([100, 100]).astype('uint16')
        fits.append(self.temp('test_new.fits'), data=d)
        f = fits.open(self.temp('test_new.fits'), uint=True)
        assert f[1].data.dtype == 'uint16'

    def test_uint_header_consistency(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2305

        This ensures that an HDU containing unsigned integer data always has
        the apppriate BZERO value in its header.
        """

        for int_size in (16, 32, 64):
            # Just make an array of some unsigned ints that wouldn't fit in a
            # signed int array of the same bit width
            max_uint = (2 ** int_size) - 1
            if int_size == 64:
                # Otherwise may get an overflow error, at least on Python 2
                max_uint = np.uint64(int_size)

            dtype = 'uint%d' % int_size
            arr = np.empty(100, dtype=dtype)
            arr.fill(max_uint)
            arr -= np.arange(100, dtype=dtype)

            uint_hdu = fits.PrimaryHDU(data=arr)
            assert np.all(uint_hdu.data == arr)
            assert uint_hdu.data.dtype.name == 'uint%d' % int_size
            assert 'BZERO' in uint_hdu.header
            assert uint_hdu.header['BZERO'] == (2 ** (int_size - 1))

            filename = 'uint%d.fits' % int_size
            uint_hdu.writeto(self.temp(filename))

            with fits.open(self.temp(filename), uint=True) as hdul:
                new_uint_hdu = hdul[0]
                assert np.all(new_uint_hdu.data == arr)
                assert new_uint_hdu.data.dtype.name == 'uint%d' % int_size
                assert 'BZERO' in new_uint_hdu.header
                assert new_uint_hdu.header['BZERO'] == (2 ** (int_size - 1))

    def test_blanks(self):
        """Test image data with blank spots in it (which should show up as
        NaNs in the data array.
        """

        arr = np.zeros((10, 10), dtype=np.int32)
        # One row will be blanks
        arr[1] = 999
        hdu = fits.ImageHDU(data=arr)
        hdu.header['BLANK'] = 999
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        assert np.isnan(hdul[1].data[1]).all()

    def test_invalid_blanks(self):
        """
        Test that invalid use of the BLANK keyword leads to an appropriate
        warning, and that the BLANK keyword is ignored when returning the
        HDU data.

        Regression test for https://github.com/astropy/astropy/issues/3865
        """

        arr = np.arange(5, dtype=np.float64)
        hdu = fits.PrimaryHDU(data=arr)
        hdu.header['BLANK'] = 2

        with catch_warnings() as w:
            hdu.writeto(self.temp('test_new.fits'))
            # Allow the HDU to be written, but there should be a warning
            # when writing a header with BLANK when then data is not
            # int
            assert len(w) == 1
            assert "Invalid 'BLANK' keyword in header" in str(w[0].message)

        # Should also get a warning when opening the file, and the BLANK
        # value should not be applied
        with catch_warnings() as w:
            with fits.open(self.temp('test_new.fits')) as h:
                assert len(w) == 1
                assert "Invalid 'BLANK' keyword in header" in str(w[0].message)
                assert np.all(arr == h[0].data)

    def test_scale_back_with_blanks(self):
        """
        Test that when auto-rescaling integer data with "blank" values (where
        the blanks are replaced by NaN in the float data), that the "BLANK"
        keyword is removed from the header.

        Further, test that when using the ``scale_back=True`` option the blank
        values are restored properly.

        Regression test for https://github.com/astropy/astropy/issues/3865
        """

        # Make the sample file
        arr = np.arange(5, dtype=np.int32)
        hdu = fits.PrimaryHDU(data=arr)
        hdu.scale('int16', bscale=1.23)

        # Creating data that uses BLANK is currently kludgy--a separate issue
        # TODO: Rewrite this test when scaling with blank support is better
        # supported

        # Let's just add a value to the data that should be converted to NaN
        # when it is read back in:
        hdu.data[0] = 9999
        hdu.header['BLANK'] = 9999
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            data = hdul[0].data
            assert np.isnan(data[0])
            hdul.writeto(self.temp('test2.fits'))

        # Now reopen the newly written file.  It should not have a 'BLANK'
        # keyword
        with catch_warnings() as w:
            with fits.open(self.temp('test2.fits')) as hdul2:
                assert len(w) == 0
                assert 'BLANK' not in hdul2[0].header
                data = hdul2[0].data
                assert np.isnan(data[0])

        # Finally, test that scale_back keeps the BLANKs correctly
        with fits.open(self.temp('test.fits'), scale_back=True,
                       mode='update') as hdul3:
            data = hdul3[0].data
            assert np.isnan(data[0])

        with fits.open(self.temp('test.fits'),
                       do_not_scale_image_data=True) as hdul4:
            assert hdul4[0].header['BLANK'] == 9999
            assert hdul4[0].header['BSCALE'] == 1.23
            assert hdul4[0].data[0] == 9999

    def test_bzero_with_floats(self):
        """Test use of the BZERO keyword in an image HDU containing float
        data.
        """

        arr = np.zeros((10, 10)) - 1
        hdu = fits.ImageHDU(data=arr)
        hdu.header['BZERO'] = 1.0
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        arr += 1
        assert (hdul[1].data == arr).all()

    def test_rewriting_large_scaled_image(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/84 and
        https://aeon.stsci.edu/ssb/trac/pyfits/ticket/101
        """

        hdul = fits.open(self.data('fixed-1890.fits'))
        orig_data = hdul[0].data
        with ignore_warnings():
            hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[0].data == orig_data).all()
        hdul.close()

        # Just as before, but this time don't touch hdul[0].data before writing
        # back out--this is the case that failed in
        # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/84
        hdul = fits.open(self.data('fixed-1890.fits'))
        with ignore_warnings():
            hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[0].data == orig_data).all()
        hdul.close()

        # Test opening/closing/reopening a scaled file in update mode
        hdul = fits.open(self.data('fixed-1890.fits'),
                         do_not_scale_image_data=True)
        hdul.writeto(self.temp('test_new.fits'), clobber=True,
                     output_verify='silentfix')
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        orig_data = hdul[0].data
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'), mode='update')
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[0].data == orig_data).all()
        hdul = fits.open(self.temp('test_new.fits'))
        hdul.close()

    def test_image_update_header(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/105

        Replacing the original header to an image HDU and saving should update
        the NAXISn keywords appropriately and save the image data correctly.
        """

        # Copy the original file before saving to it
        self.copy_file('test0.fits')
        with fits.open(self.temp('test0.fits'), mode='update') as hdul:
            orig_data = hdul[1].data.copy()
            hdr_copy = hdul[1].header.copy()
            del hdr_copy['NAXIS*']
            hdul[1].header = hdr_copy

        with fits.open(self.temp('test0.fits')) as hdul:
            assert (orig_data == hdul[1].data).all()

    def test_open_scaled_in_update_mode(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/119
        (Don't update scaled image data if the data is not read)

        This ensures that merely opening and closing a file containing scaled
        image data does not cause any change to the data (or the header).
        Changes should only occur if the data is accessed.
        """

        # Copy the original file before making any possible changes to it
        self.copy_file('scale.fits')
        mtime = os.stat(self.temp('scale.fits')).st_mtime

        time.sleep(1)

        fits.open(self.temp('scale.fits'), mode='update').close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp('scale.fits')).st_mtime

        # Insert a slight delay to ensure the mtime does change when the file
        # is changed
        time.sleep(1)

        hdul = fits.open(self.temp('scale.fits'), 'update')
        orig_data = hdul[0].data
        hdul.close()

        # Now the file should be updated with the rescaled data
        assert mtime != os.stat(self.temp('scale.fits')).st_mtime
        hdul = fits.open(self.temp('scale.fits'), mode='update')
        assert hdul[0].data.dtype == np.dtype('>f4')
        assert hdul[0].header['BITPIX'] == -32
        assert 'BZERO' not in hdul[0].header
        assert 'BSCALE' not in hdul[0].header
        assert (orig_data == hdul[0].data).all()

        # Try reshaping the data, then closing and reopening the file; let's
        # see if all the changes are preseved properly
        hdul[0].data.shape = (42, 10)
        hdul.close()

        hdul = fits.open(self.temp('scale.fits'))
        assert hdul[0].shape == (42, 10)
        assert hdul[0].data.dtype == np.dtype('>f4')
        assert hdul[0].header['BITPIX'] == -32
        assert 'BZERO' not in hdul[0].header
        assert 'BSCALE' not in hdul[0].header

    def test_scale_back(self):
        """A simple test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/120

        The scale_back feature for image HDUs.
        """

        self.copy_file('scale.fits')
        with fits.open(self.temp('scale.fits'), mode='update',
                       scale_back=True) as hdul:
            orig_bitpix = hdul[0].header['BITPIX']
            orig_bzero = hdul[0].header['BZERO']
            orig_bscale = hdul[0].header['BSCALE']
            orig_data = hdul[0].data.copy()
            hdul[0].data[0] = 0

        with fits.open(self.temp('scale.fits'),
                       do_not_scale_image_data=True) as hdul:
            assert hdul[0].header['BITPIX'] == orig_bitpix
            assert hdul[0].header['BZERO'] == orig_bzero
            assert hdul[0].header['BSCALE'] == orig_bscale

            zero_point = int(math.floor(-orig_bzero / orig_bscale))
            assert (hdul[0].data[0] == zero_point).all()

        with fits.open(self.temp('scale.fits')) as hdul:
            assert (hdul[0].data[1:] == orig_data[1:]).all()

    def test_image_none(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/27
        """

        with fits.open(self.data('test0.fits')) as h:
            h[1].data
            h[1].data = None
            h[1].writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert h[1].data is None
            assert h[1].header['NAXIS'] == 0
            assert 'NAXIS1' not in h[1].header
            assert 'NAXIS2' not in h[1].header

    def test_invalid_blank(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2711

        If the BLANK keyword contains an invalid value it should be ignored for
        any calculations (though a warning should be issued).
        """

        data = np.arange(100, dtype=np.float64)
        hdu = fits.PrimaryHDU(data)
        hdu.header['BLANK'] = 'nan'
        hdu.writeto(self.temp('test.fits'))

        with catch_warnings() as w:
            with fits.open(self.temp('test.fits')) as hdul:
                assert np.all(hdul[0].data == data)

        assert len(w) == 2
        msg = "Invalid value for 'BLANK' keyword in header"
        assert msg in str(w[0].message)
        msg = "Invalid 'BLANK' keyword"
        assert msg in str(w[1].message)

    def test_scaled_image_fromfile(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2710
        """

        # Make some sample data
        a = np.arange(100, dtype=np.float32)

        hdu = fits.PrimaryHDU(data=a.copy())
        hdu.scale(bscale=1.1)
        hdu.writeto(self.temp('test.fits'))

        with open(self.temp('test.fits'), 'rb') as f:
            file_data = f.read()

        hdul = fits.HDUList.fromstring(file_data)
        assert np.allclose(hdul[0].data, a)


class TestCompressedImage(FitsTestCase):
    def test_empty(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2595
        """

        hdu = fits.CompImageHDU()
        assert hdu.data is None
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits'), mode='update') as hdul:
            assert len(hdul) == 2
            assert isinstance(hdul[1], fits.CompImageHDU)
            assert hdul[1].data is None

            # Now test replacing the empty data with an array and see what
            # happens
            hdul[1].data = np.arange(100, dtype=np.int32)

        with fits.open(self.temp('test.fits')) as hdul:
            assert len(hdul) == 2
            assert isinstance(hdul[1], fits.CompImageHDU)
            assert np.all(hdul[1].data == np.arange(100, dtype=np.int32))

    @pytest.mark.parametrize(
        ('data', 'compression_type', 'quantize_level', 'byte_order'),
        sum([[(np.zeros((2, 10, 10), dtype=np.float32), 'RICE_1', 16, bo),
              (np.zeros((2, 10, 10), dtype=np.float32), 'GZIP_1', -0.01, bo),
              (np.zeros((100, 100)) + 1, 'HCOMPRESS_1', 16, bo)]
             for bo in ('<', '>')], []))
    def test_comp_image(self, data, compression_type, quantize_level,
                        byte_order):
        data = data.newbyteorder(byte_order)
        primary_hdu = fits.PrimaryHDU()
        ofd = fits.HDUList(primary_hdu)
        chdu = fits.CompImageHDU(data, name='SCI',
                                 compressionType=compression_type,
                                 quantizeLevel=quantize_level)
        ofd.append(chdu)
        ofd.writeto(self.temp('test_new.fits'), clobber=True)
        ofd.close()
        with fits.open(self.temp('test_new.fits')) as fd:
            assert (fd[1].data == data).all()
            assert fd[1].header['NAXIS'] == chdu.header['NAXIS']
            assert fd[1].header['NAXIS1'] == chdu.header['NAXIS1']
            assert fd[1].header['NAXIS2'] == chdu.header['NAXIS2']
            assert fd[1].header['BITPIX'] == chdu.header['BITPIX']

    @ignore_warnings(AstropyPendingDeprecationWarning)
    def test_comp_image_hcompression_1_invalid_data(self):
        """
        Tests compression with the HCOMPRESS_1 algorithm with data that is
        not 2D and has a non-2D tile size.
        """

        pytest.raises(ValueError, fits.CompImageHDU,
                      np.zeros((2, 10, 10), dtype=np.float32), name='SCI',
                      compressionType='HCOMPRESS_1', quantizeLevel=16,
                      tileSize=[2, 10, 10])

    @ignore_warnings(AstropyPendingDeprecationWarning)
    def test_comp_image_hcompress_image_stack(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/171

        Tests that data containing more than two dimensions can be
        compressed with HCOMPRESS_1 so long as the user-supplied tile size can
        be flattened to two dimensions.
        """

        cube = np.arange(300, dtype=np.float32).reshape((3, 10, 10))
        hdu = fits.CompImageHDU(data=cube, name='SCI',
                                compressionType='HCOMPRESS_1',
                                quantizeLevel=16, tileSize=[5, 5, 1])
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert (hdul['SCI'].data == cube).all()

    def test_subtractive_dither_seed(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/32

        Ensure that when floating point data is compressed with the
        SUBTRACTIVE_DITHER_1 quantization method that the correct ZDITHER0 seed
        is added to the header, and that the data can be correctly
        decompressed.
        """

        array = np.arange(100.0).reshape(10, 10)
        csum = (array[0].view('uint8').sum() % 10000) + 1
        hdu = fits.CompImageHDU(data=array,
                                quantize_method=SUBTRACTIVE_DITHER_1,
                                dither_seed=DITHER_SEED_CHECKSUM)
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert isinstance(hdul[1], fits.CompImageHDU)
            assert 'ZQUANTIZ' in hdul[1]._header
            assert hdul[1]._header['ZQUANTIZ'] == 'SUBTRACTIVE_DITHER_1'
            assert 'ZDITHER0' in hdul[1]._header
            assert hdul[1]._header['ZDITHER0'] == csum
            assert np.all(hdul[1].data == array)

    def test_disable_image_compression(self):
        with catch_warnings():
            # No warnings should be displayed in this case
            warnings.simplefilter('error')
            with fits.open(self.data('comp.fits'),
                           disable_image_compression=True) as hdul:
                # The compressed image HDU should show up as a BinTableHDU, but
                # *not* a CompImageHDU
                assert isinstance(hdul[1], fits.BinTableHDU)
                assert not isinstance(hdul[1], fits.CompImageHDU)

        with fits.open(self.data('comp.fits')) as hdul:
            assert isinstance(hdul[1], fits.CompImageHDU)

    def test_open_comp_image_in_update_mode(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/167

        Similar to test_open_scaled_in_update_mode(), but specifically for
        compressed images.
        """

        # Copy the original file before making any possible changes to it
        self.copy_file('comp.fits')
        mtime = os.stat(self.temp('comp.fits')).st_mtime

        time.sleep(1)

        fits.open(self.temp('comp.fits'), mode='update').close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp('comp.fits')).st_mtime

    def test_open_scaled_in_update_mode_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 2

        Identical to test_open_scaled_in_update_mode() but with a compressed
        version of the scaled image.
        """

        # Copy+compress the original file before making any possible changes to
        # it
        with fits.open(self.data('scale.fits'),
                       do_not_scale_image_data=True) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data,
                                     header=hdul[0].header)
            chdu.writeto(self.temp('scale.fits'))
        mtime = os.stat(self.temp('scale.fits')).st_mtime

        time.sleep(1)

        fits.open(self.temp('scale.fits'), mode='update').close()

        # Ensure that no changes were made to the file merely by immediately
        # opening and closing it.
        assert mtime == os.stat(self.temp('scale.fits')).st_mtime

        # Insert a slight delay to ensure the mtime does change when the file
        # is changed
        time.sleep(1)

        hdul = fits.open(self.temp('scale.fits'), 'update')
        hdul[1].data
        hdul.close()

        # Now the file should be updated with the rescaled data
        assert mtime != os.stat(self.temp('scale.fits')).st_mtime
        hdul = fits.open(self.temp('scale.fits'), mode='update')
        assert hdul[1].data.dtype == np.dtype('float32')
        assert hdul[1].header['BITPIX'] == -32
        assert 'BZERO' not in hdul[1].header
        assert 'BSCALE' not in hdul[1].header

        # Try reshaping the data, then closing and reopening the file; let's
        # see if all the changes are preseved properly
        hdul[1].data.shape = (42, 10)
        hdul.close()

        hdul = fits.open(self.temp('scale.fits'))
        assert hdul[1].shape == (42, 10)
        assert hdul[1].data.dtype == np.dtype('float32')
        assert hdul[1].header['BITPIX'] == -32
        assert 'BZERO' not in hdul[1].header
        assert 'BSCALE' not in hdul[1].header

    def test_write_comp_hdu_direct_from_existing(self):
        with fits.open(self.data('comp.fits')) as hdul:
            hdul[1].writeto(self.temp('test.fits'))

        with fits.open(self.data('comp.fits')) as hdul1:
            with fits.open(self.temp('test.fits')) as hdul2:
                assert np.all(hdul1[1].data == hdul2[1].data)
                assert comparerecords(hdul1[1].compressed_data,
                                      hdul2[1].compressed_data)

    def test_rewriting_large_scaled_image_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 1

        Identical to test_rewriting_large_scaled_image() but with a compressed
        image.
        """

        with fits.open(self.data('fixed-1890.fits'),
                       do_not_scale_image_data=True) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data,
                                     header=hdul[0].header)
            chdu.writeto(self.temp('fixed-1890-z.fits'))

        hdul = fits.open(self.temp('fixed-1890-z.fits'))
        orig_data = hdul[1].data
        with ignore_warnings():
            hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[1].data == orig_data).all()
        hdul.close()

        # Just as before, but this time don't touch hdul[0].data before writing
        # back out--this is the case that failed in
        # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/84
        hdul = fits.open(self.temp('fixed-1890-z.fits'))
        with ignore_warnings():
            hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[1].data == orig_data).all()
        hdul.close()

        # Test opening/closing/reopening a scaled file in update mode
        hdul = fits.open(self.temp('fixed-1890-z.fits'),
                         do_not_scale_image_data=True)
        hdul.writeto(self.temp('test_new.fits'), clobber=True,
                     output_verify='silentfix')
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        orig_data = hdul[1].data
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'), mode='update')
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[1].data == orig_data).all()
        hdul = fits.open(self.temp('test_new.fits'))
        hdul.close()

    def test_scale_back_compressed(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/88 3

        Identical to test_scale_back() but uses a compressed image.
        """

        # Create a compressed version of the scaled image
        with fits.open(self.data('scale.fits'),
                       do_not_scale_image_data=True) as hdul:
            chdu = fits.CompImageHDU(data=hdul[0].data,
                                     header=hdul[0].header)
            chdu.writeto(self.temp('scale.fits'))

        with fits.open(self.temp('scale.fits'), mode='update',
                       scale_back=True) as hdul:
            orig_bitpix = hdul[1].header['BITPIX']
            orig_bzero = hdul[1].header['BZERO']
            orig_bscale = hdul[1].header['BSCALE']
            orig_data = hdul[1].data.copy()
            hdul[1].data[0] = 0

        with fits.open(self.temp('scale.fits'),
                       do_not_scale_image_data=True) as hdul:
            assert hdul[1].header['BITPIX'] == orig_bitpix
            assert hdul[1].header['BZERO'] == orig_bzero
            assert hdul[1].header['BSCALE'] == orig_bscale

            zero_point = int(math.floor(-orig_bzero / orig_bscale))
            assert (hdul[1].data[0] == zero_point).all()

        with fits.open(self.temp('scale.fits')) as hdul:
            assert (hdul[1].data[1:] == orig_data[1:]).all()
            # Extra test to ensure that after everything the data is still the
            # same as in the original uncompressed version of the image
            with fits.open(self.data('scale.fits')) as hdul2:
                # Recall we made the same modification to the data in hdul
                # above
                hdul2[0].data[0] = 0
                assert (hdul[1].data == hdul2[0].data).all()

    def test_lossless_gzip_compression(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/198"""

        noise = np.random.normal(size=(1000, 1000))

        chdu1 = fits.CompImageHDU(data=noise, compressionType='GZIP_1')
        # First make a test image with lossy compression and make sure it
        # wasn't compressed perfectly.  This shouldn't happen ever, but just to
        # make sure the test non-trivial.
        chdu1.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as h:
            assert np.abs(noise - h[1].data).max() > 0.0

        del h

        chdu2 = fits.CompImageHDU(data=noise, compressionType='GZIP_1',
                                  quantizeLevel=0.0)  # No quantization
        with ignore_warnings():
            chdu2.writeto(self.temp('test.fits'), clobber=True)

        with fits.open(self.temp('test.fits')) as h:
            assert (noise == h[1].data).all()

    def test_compression_column_tforms(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/199"""

        # Some interestingly tiled data so that some of it is quantized and
        # some of it ends up just getting gzip-compressed
        data2 = ((np.arange(1, 8, dtype=np.float32) * 10)[:, np.newaxis] +
                 np.arange(1, 7))
        np.random.seed(1337)
        data1 = np.random.uniform(size=(6 * 4, 7 * 4))
        data1[:data2.shape[0], :data2.shape[1]] = data2
        chdu = fits.CompImageHDU(data1, compressionType='RICE_1',
                                 tileSize=(6, 7))
        chdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits'),
                       disable_image_compression=True) as h:
            assert h[1].header['TFORM1'] == '1PB(30)'
            assert h[1].header['TFORM2'] == '1PB(359)'

    def test_compression_update_header(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/23
        """

        self.copy_file('comp.fits')
        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            assert isinstance(hdul[1], fits.CompImageHDU)
            hdul[1].header['test1'] = 'test'
            hdul[1]._header['test2'] = 'test2'

        with fits.open(self.temp('comp.fits')) as hdul:
            assert 'test1' in hdul[1].header
            assert hdul[1].header['test1'] == 'test'
            assert 'test2' in hdul[1].header
            assert  hdul[1].header['test2'] == 'test2'

        # Test update via index now:
        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            hdr = hdul[1].header
            hdr[hdr.index('TEST1')] = 'foo'

        with fits.open(self.temp('comp.fits')) as hdul:
            assert hdul[1].header['TEST1'] == 'foo'

        # Test slice updates
        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            hdul[1].header['TEST*'] = 'qux'

        with fits.open(self.temp('comp.fits')) as hdul:
            assert list(hdul[1].header['TEST*'].values()) == ['qux', 'qux']

        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            hdr = hdul[1].header
            idx = hdr.index('TEST1')
            hdr[idx:idx + 2] = 'bar'

        with fits.open(self.temp('comp.fits')) as hdul:
            assert list(hdul[1].header['TEST*'].values()) == ['bar', 'bar']

        # Test updating a specific COMMENT card duplicate
        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            hdul[1].header[('COMMENT', 1)] = 'I am fire. I am death!'

        with fits.open(self.temp('comp.fits')) as hdul:
            assert hdul[1].header['COMMENT'][1] == 'I am fire. I am death!'
            assert hdul[1]._header['COMMENT'][1] == 'I am fire. I am death!'

        # Test deleting by keyword and by slice
        with fits.open(self.temp('comp.fits'), mode='update') as hdul:
            hdr = hdul[1].header
            del hdr['COMMENT']
            idx = hdr.index('TEST1')
            del hdr[idx:idx + 2]

        with fits.open(self.temp('comp.fits')) as hdul:
            assert 'COMMENT' not in hdul[1].header
            assert 'COMMENT' not in hdul[1]._header
            assert 'TEST1' not in hdul[1].header
            assert 'TEST1' not in hdul[1]._header
            assert 'TEST2' not in hdul[1].header
            assert 'TEST2' not in hdul[1]._header

    def test_compression_update_header_with_reserved(self):
        """
        Ensure that setting reserved keywords related to the table data
        structure on CompImageHDU image headers fails.
        """

        def test_set_keyword(hdr, keyword, value):
            with catch_warnings() as w:
                hdr[keyword] = value
                assert len(w) == 1
                assert str(w[0].message).startswith(
                        "Keyword %r is reserved" % keyword)
                assert keyword not in hdr

        with fits.open(self.data('comp.fits')) as hdul:
            hdr = hdul[1].header
            test_set_keyword(hdr, 'TFIELDS', 8)
            test_set_keyword(hdr, 'TTYPE1', 'Foo')
            test_set_keyword(hdr, 'ZCMPTYPE', 'ASDF')
            test_set_keyword(hdr, 'ZVAL1', 'Foo')

    def test_compression_header_append(self):
        with fits.open(self.data('comp.fits')) as hdul:
            imghdr = hdul[1].header
            tblhdr = hdul[1]._header
            with catch_warnings() as w:
                imghdr.append('TFIELDS')
                assert len(w) == 1
                assert 'TFIELDS' not in imghdr

            imghdr.append(('FOO', 'bar', 'qux'), end=True)
            assert 'FOO' in imghdr
            assert imghdr[-1] == 'bar'
            assert 'FOO' in tblhdr
            assert tblhdr[-1] == 'bar'

            imghdr.append(('CHECKSUM', 'abcd1234'))
            assert 'CHECKSUM' in imghdr
            assert imghdr['CHECKSUM'] == 'abcd1234'
            assert 'CHECKSUM' not in tblhdr
            assert 'ZHECKSUM' in tblhdr
            assert tblhdr['ZHECKSUM'] == 'abcd1234'

    def test_compression_header_insert(self):
        with fits.open(self.data('comp.fits')) as hdul:
            imghdr = hdul[1].header
            tblhdr = hdul[1]._header
            # First try inserting a restricted keyword
            with catch_warnings() as w:
                imghdr.insert(1000, 'TFIELDS')
                assert len(w) == 1
                assert 'TFIELDS' not in imghdr
                assert tblhdr.count('TFIELDS') == 1

            # First try keyword-relative insert
            imghdr.insert('TELESCOP', ('OBSERVER', 'Phil Plait'))
            assert 'OBSERVER' in imghdr
            assert imghdr.index('OBSERVER') == imghdr.index('TELESCOP') - 1
            assert 'OBSERVER' in tblhdr
            assert tblhdr.index('OBSERVER') == tblhdr.index('TELESCOP') - 1

            # Next let's see if an index-relative insert winds up being
            # sensible
            idx = imghdr.index('OBSERVER')
            imghdr.insert('OBSERVER', ('FOO',))
            assert 'FOO' in imghdr
            assert imghdr.index('FOO') == idx
            assert 'FOO' in tblhdr
            assert tblhdr.index('FOO') == tblhdr.index('OBSERVER') - 1

    def test_compression_header_set_before_after(self):
        with fits.open(self.data('comp.fits')) as hdul:
            imghdr = hdul[1].header
            tblhdr = hdul[1]._header

            with catch_warnings() as w:
                imghdr.set('ZBITPIX', 77, 'asdf', after='XTENSION')
                assert len(w) == 1
                assert 'ZBITPIX' not in imghdr
                assert tblhdr.count('ZBITPIX') == 1
                assert tblhdr['ZBITPIX'] != 77

            # Move GCOUNT before PCOUNT (not that there's any reason you'd
            # *want* to do that, but it's just a test...)
            imghdr.set('GCOUNT', 99, before='PCOUNT')
            assert imghdr.index('GCOUNT') == imghdr.index('PCOUNT') - 1
            assert imghdr['GCOUNT'] == 99
            assert tblhdr.index('ZGCOUNT') == tblhdr.index('ZPCOUNT') - 1
            assert tblhdr['ZGCOUNT'] == 99
            assert tblhdr.index('PCOUNT') == 5
            assert tblhdr.index('GCOUNT') == 6
            assert tblhdr['GCOUNT'] == 1

            imghdr.set('GCOUNT', 2, after='PCOUNT')
            assert imghdr.index('GCOUNT') == imghdr.index('PCOUNT') + 1
            assert imghdr['GCOUNT'] == 2
            assert tblhdr.index('ZGCOUNT') == tblhdr.index('ZPCOUNT') + 1
            assert tblhdr['ZGCOUNT'] == 2
            assert tblhdr.index('PCOUNT') == 5
            assert tblhdr.index('GCOUNT') == 6
            assert tblhdr['GCOUNT'] == 1

    def test_compression_header_append_commentary(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2363
        """

        hdu = fits.CompImageHDU(np.array([0], dtype=np.int32))
        hdu.header['COMMENT'] = 'hello world'
        assert hdu.header['COMMENT'] == ['hello world']
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert hdul[1].header['COMMENT'] == ['hello world']

    def test_compression_with_gzip_column(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/71
        """

        arr = np.zeros((2, 7000), dtype='float32')

        # The first row (which will be the first compressed tile) has a very
        # wide range of values that will be difficult to quantize, and should
        # result in use of a GZIP_COMPRESSED_DATA column
        arr[0] = np.linspace(0, 1, 7000)
        arr[1] = np.random.normal(size=7000)

        hdu = fits.CompImageHDU(data=arr)
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            comp_hdu = hdul[1]

            # GZIP-compressed tile should compare exactly
            assert np.all(comp_hdu.data[0] == arr[0])
            # The second tile uses lossy compression and may be somewhat off,
            # so we don't bother comparing it exactly

    def test_duplicate_compression_header_keywords(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/2750

        Tests that the fake header (for the compressed image) can still be read
        even if the real header contained a duplicate ZTENSION keyword (the
        issue applies to any keyword specific to the compression convention,
        however).
        """

        arr = np.arange(100, dtype=np.int32)
        hdu = fits.CompImageHDU(data=arr)

        header = hdu._header
        # append the duplicate keyword
        hdu._header.append(('ZTENSION', 'IMAGE'))
        hdu.writeto(self.temp('test.fits'))

        with fits.open(self.temp('test.fits')) as hdul:
            assert header == hdul[1]._header
            # There's no good reason to have a duplicate keyword, but
            # technically it isn't invalid either :/
            assert hdul[1]._header.count('ZTENSION') == 2
