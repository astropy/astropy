# Licensed under a 3-clause BSD style license - see PYFITS.rst

import glob
import io
import os
import platform
import sys
import copy
import subprocess

import pytest
import numpy as np

from astropy.io.fits.verify import VerifyError
from astropy.io import fits
from astropy.tests.helper import raises, catch_warnings, ignore_warnings
from astropy.utils.exceptions import AstropyUserWarning, AstropyDeprecationWarning

from . import FitsTestCase


class TestHDUListFunctions(FitsTestCase):
    def test_update_name(self):
        hdul = fits.open(self.data('o4sp040b0_raw.fits'))
        hdul[4].name = 'Jim'
        hdul[4].ver = 9
        assert hdul[('JIM', 9)].header['extname'] == 'JIM'

    def test_hdu_file_bytes(self):
        hdul = fits.open(self.data('checksum.fits'))
        res = hdul[0].filebytes()
        assert res == 11520
        res = hdul[1].filebytes()
        assert res == 8640

    def test_hdulist_file_info(self):
        hdul = fits.open(self.data('checksum.fits'))
        res = hdul.fileinfo(0)

        def test_fileinfo(**kwargs):
            assert res['datSpan'] == kwargs.get('datSpan', 2880)
            assert res['resized'] == kwargs.get('resized', False)
            assert res['filename'] == self.data('checksum.fits')
            assert res['datLoc'] == kwargs.get('datLoc', 8640)
            assert res['hdrLoc'] == kwargs.get('hdrLoc', 0)
            assert res['filemode'] == 'readonly'

        res = hdul.fileinfo(1)
        test_fileinfo(datLoc=17280, hdrLoc=11520)

        hdu = fits.ImageHDU(data=hdul[0].data)
        hdul.insert(1, hdu)

        res = hdul.fileinfo(0)
        test_fileinfo(resized=True)

        res = hdul.fileinfo(1)
        test_fileinfo(datSpan=None, resized=True, datLoc=None, hdrLoc=None)

        res = hdul.fileinfo(2)
        test_fileinfo(resized=1, datLoc=17280, hdrLoc=11520)

    def test_create_from_multiple_primary(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/145

        Ensure that a validation error occurs when saving an HDUList containing
        multiple PrimaryHDUs.
        """

        hdul = fits.HDUList([fits.PrimaryHDU(), fits.PrimaryHDU()])
        pytest.raises(VerifyError, hdul.writeto, self.temp('temp.fits'),
                      output_verify='exception')

    def test_append_primary_to_empty_list(self):
        # Tests appending a Simple PrimaryHDU to an empty HDUList.
        hdul = fits.HDUList()
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_extension_to_empty_list(self):
        """Tests appending a Simple ImageHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (100,), 'int32', '')]
        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_table_extension_to_empty_list(self):
        """Tests appending a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = fits.HDUList()
        hdul1 = fits.open(self.data('tb.fits'))
        hdul.append(hdul1[1])
        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (), '', ''),
                (1, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_groupshdu_to_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 1, 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_primary_to_non_empty_list(self):
        """Tests appending a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('arange.fits'))
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 1, 'ImageHDU', 6, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_extension_to_non_empty_list(self):
        """Tests appending a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('tb.fits'))
        hdul.append(hdul[1])

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 11, (), '', ''),
                (1, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    @raises(ValueError)
    def test_append_groupshdu_to_non_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        hdu = fits.GroupsHDU()
        hdul.append(hdu)

    def test_insert_primary_to_empty_list(self):
        """Tests inserting a Simple PrimaryHDU to an empty HDUList."""
        hdul = fits.HDUList()
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_extension_to_empty_list(self):
        """Tests inserting a Simple ImageHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_table_extension_to_empty_list(self):
        """Tests inserting a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = fits.HDUList()
        hdul1 = fits.open(self.data('tb.fits'))
        hdul.insert(0, hdul1[1])

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (), '', ''),
                (1, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_groupshdu_to_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 1, 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_primary_to_non_empty_list(self):
        """Tests inserting a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('arange.fits'))
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(1, hdu)

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 1, 'ImageHDU', 6, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_extension_to_non_empty_list(self):
        """Tests inserting a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('tb.fits'))
        hdul.insert(1, hdul[1])

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 11, (), '', ''),
                (1, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_groupshdu_to_non_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)
        hdu = fits.GroupsHDU()

        with pytest.raises(ValueError):
            hdul.insert(1, hdu)

        info = [(0, 'PRIMARY', 1, 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters'),
                (1, '', 1, 'ImageHDU', 6, (100,), 'int32', '')]

        hdul.insert(0, hdu)

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    @raises(ValueError)
    def test_insert_groupshdu_to_begin_of_hdulist_with_groupshdu(self):
        """
        Tests inserting a Simple GroupsHDU to the beginning of an HDUList
        that that already contains a GroupsHDU.
        """

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.insert(0, hdu)
        hdul.insert(0, hdu)

    def test_insert_extension_to_primary_in_non_empty_list(self):
        # Tests inserting a Simple ExtensionHDU to a non-empty HDUList.
        hdul = fits.open(self.data('tb.fits'))
        hdul.insert(0, hdul[1])

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 4, (), '', ''),
                (1, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 1, 'ImageHDU', 12, (), '', ''),
                (3, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_image_extension_to_primary_in_non_empty_list(self):
        """
        Tests inserting a Simple Image ExtensionHDU to a non-empty HDUList
        as the primary HDU.
        """

        hdul = fits.open(self.data('tb.fits'))
        hdu = fits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', ''),
                (1, '', 1, 'ImageHDU', 12, (), '', ''),
                (2, '', 1, 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_filename(self):
        """Tests the HDUList filename method."""

        hdul = fits.open(self.data('tb.fits'))
        name = hdul.filename()
        assert name == self.data('tb.fits')

    def test_file_like(self):
        """
        Tests the use of a file like object with no tell or seek methods
        in HDUList.writeto(), HDULIST.flush() or astropy.io.fits.writeto()
        """

        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul = fits.HDUList()
        hdul.append(hdu)
        tmpfile = open(self.temp('tmpfile.fits'), 'wb')
        hdul.writeto(tmpfile)
        tmpfile.close()

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_file_like_2(self):
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        tmpfile = open(self.temp('tmpfile.fits'), 'wb')
        hdul = fits.open(tmpfile, mode='ostream')
        hdul.append(hdu)
        hdul.flush()
        tmpfile.close()
        hdul.close()

        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_file_like_3(self):
        tmpfile = open(self.temp('tmpfile.fits'), 'wb')
        fits.writeto(tmpfile, np.arange(100, dtype=np.int32))
        tmpfile.close()
        info = [(0, 'PRIMARY', 1, 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_shallow_copy(self):
        """
        Tests that `HDUList.__copy__()` and `HDUList.copy()` return a
        shallow copy (regression test for #7211).
        """

        n = np.arange(10.0)
        primary_hdu = fits.PrimaryHDU(n)
        hdu = fits.ImageHDU(n)
        hdul = fits.HDUList([primary_hdu, hdu])

        for hdulcopy in (hdul.copy(), copy.copy(hdul)):
            assert isinstance(hdulcopy, fits.HDUList)
            assert hdulcopy is not hdul
            assert hdulcopy[0] is hdul[0]
            assert hdulcopy[1] is hdul[1]

    def test_deep_copy(self):
        """
        Tests that `HDUList.__deepcopy__()` returns a deep copy.
        """

        n = np.arange(10.0)
        primary_hdu = fits.PrimaryHDU(n)
        hdu = fits.ImageHDU(n)
        hdul = fits.HDUList([primary_hdu, hdu])

        hdulcopy = copy.deepcopy(hdul)

        assert isinstance(hdulcopy, fits.HDUList)
        assert hdulcopy is not hdul

        for index in range(len(hdul)):
            assert hdulcopy[index] is not hdul[index]
            assert hdulcopy[index].header == hdul[index].header
            np.testing.assert_array_equal(hdulcopy[index].data, hdul[index].data)

    def test_new_hdu_extname(self):
        """
        Tests that new extension HDUs that are added to an HDUList can be
        properly indexed by their EXTNAME/EXTVER (regression test for
        ticket:48).
        """

        f = fits.open(self.data('test0.fits'))
        hdul = fits.HDUList()
        hdul.append(f[0].copy())
        hdu = fits.ImageHDU(header=f[1].header)
        hdul.append(hdu)

        assert hdul[1].header['EXTNAME'] == 'SCI'
        assert hdul[1].header['EXTVER'] == 1
        assert hdul.index_of(('SCI', 1)) == 1
        assert hdul.index_of(hdu) == len(hdul) - 1

    def test_update_filelike(self):
        """Test opening a file-like object in update mode and resizing the
        HDU.
        """

        sf = io.BytesIO()
        arr = np.zeros((100, 100))
        hdu = fits.PrimaryHDU(data=arr)
        hdu.writeto(sf)

        sf.seek(0)
        arr = np.zeros((200, 200))
        hdul = fits.open(sf, mode='update')
        hdul[0].data = arr
        hdul.flush()

        sf.seek(0)
        hdul = fits.open(sf)
        assert len(hdul) == 1
        assert (hdul[0].data == arr).all()

    def test_flush_readonly(self):
        """Test flushing changes to a file opened in a read only mode."""

        oldmtime = os.stat(self.data('test0.fits')).st_mtime
        hdul = fits.open(self.data('test0.fits'))
        hdul[0].header['FOO'] = 'BAR'
        with catch_warnings(AstropyUserWarning) as w:
            hdul.flush()
        assert len(w) == 1
        assert 'mode is not supported' in str(w[0].message)
        assert oldmtime == os.stat(self.data('test0.fits')).st_mtime

    def test_fix_extend_keyword(self):
        hdul = fits.HDUList()
        hdul.append(fits.PrimaryHDU())
        hdul.append(fits.ImageHDU())
        del hdul[0].header['EXTEND']
        hdul.verify('silentfix')

        assert 'EXTEND' in hdul[0].header
        assert hdul[0].header['EXTEND'] is True

    def test_fix_malformed_naxisj(self):
        """
        Tests that malformed NAXISj values are fixed sensibly.
        """

        hdu = fits.open(self.data('arange.fits'))

        # Malform NAXISj header data
        hdu[0].header['NAXIS1'] = 11.0
        hdu[0].header['NAXIS2'] = '10.0'
        hdu[0].header['NAXIS3'] = '7'

        # Axes cache needs to be malformed as well
        hdu[0]._axes = [11.0, '10.0', '7']

        # Perform verification including the fix
        hdu.verify('silentfix')

        # Check that malformed data was converted
        assert hdu[0].header['NAXIS1'] == 11
        assert hdu[0].header['NAXIS2'] == 10
        assert hdu[0].header['NAXIS3'] == 7

    def test_fix_wellformed_naxisj(self):
        """
        Tests that wellformed NAXISj values are not modified.
        """

        hdu = fits.open(self.data('arange.fits'))

        # Fake new NAXISj header data
        hdu[0].header['NAXIS1'] = 768
        hdu[0].header['NAXIS2'] = 64
        hdu[0].header['NAXIS3'] = 8

        # Axes cache needs to be faked as well
        hdu[0]._axes = [768, 64, 8]

        # Perform verification including the fix
        hdu.verify('silentfix')

        # Check that malformed data was converted
        assert hdu[0].header['NAXIS1'] == 768
        assert hdu[0].header['NAXIS2'] == 64
        assert hdu[0].header['NAXIS3'] == 8

    def test_new_hdulist_extend_keyword(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/114

        Tests that adding a PrimaryHDU to a new HDUList object updates the
        EXTEND keyword on that HDU.
        """

        h0 = fits.Header()
        hdu = fits.PrimaryHDU(header=h0)
        sci = fits.ImageHDU(data=np.array(10))
        image = fits.HDUList([hdu, sci])
        image.writeto(self.temp('temp.fits'))
        assert 'EXTEND' in hdu.header
        assert hdu.header['EXTEND'] is True

    def test_replace_memmaped_array(self):
        # Copy the original before we modify it
        hdul = fits.open(self.data('test0.fits'))
        hdul.writeto(self.temp('temp.fits'))

        hdul = fits.open(self.temp('temp.fits'), mode='update', memmap=True)
        old_data = hdul[1].data.copy()
        hdul[1].data = hdul[1].data + 1
        hdul.close()
        hdul = fits.open(self.temp('temp.fits'), memmap=True)
        assert ((old_data + 1) == hdul[1].data).all()

    def test_open_file_with_end_padding(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/106

        Open files with end padding bytes.
        """

        hdul = fits.open(self.data('test0.fits'),
                         do_not_scale_image_data=True)
        info = hdul.info(output=False)
        hdul.writeto(self.temp('temp.fits'))
        with open(self.temp('temp.fits'), 'ab') as f:
            f.seek(0, os.SEEK_END)
            f.write(b'\0' * 2880)
        with ignore_warnings():
            assert info == fits.info(self.temp('temp.fits'), output=False,
                                     do_not_scale_image_data=True)

    def test_open_file_with_bad_header_padding(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/136

        Open files with nulls for header block padding instead of spaces.
        """

        a = np.arange(100).reshape(10, 10)
        hdu = fits.PrimaryHDU(data=a)
        hdu.writeto(self.temp('temp.fits'))

        # Figure out where the header padding begins and fill it with nulls
        end_card_pos = str(hdu.header).index('END' + ' ' * 77)
        padding_start = end_card_pos + 80
        padding_len = 2880 - padding_start
        with open(self.temp('temp.fits'), 'r+b') as f:
            f.seek(padding_start)
            f.write('\0'.encode('ascii') * padding_len)

        with catch_warnings(AstropyUserWarning) as w:
            with fits.open(self.temp('temp.fits')) as hdul:
                assert (hdul[0].data == a).all()
        assert ('contains null bytes instead of spaces' in
                str(w[0].message))
        assert len(w) == 1
        assert len(hdul) == 1
        assert str(hdul[0].header) == str(hdu.header)

    def test_update_with_truncated_header(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/148

        Test that saving an update where the header is shorter than the
        original header doesn't leave a stump from the old header in the file.
        """

        data = np.arange(100)
        hdu = fits.PrimaryHDU(data=data)
        idx = 1
        while len(hdu.header) < 34:
            hdu.header['TEST{}'.format(idx)] = idx
            idx += 1
        hdu.writeto(self.temp('temp.fits'), checksum=True)

        with fits.open(self.temp('temp.fits'), mode='update') as hdul:
            # Modify the header, forcing it to be rewritten
            hdul[0].header['TEST1'] = 2

        with fits.open(self.temp('temp.fits')) as hdul:
            assert (hdul[0].data == data).all()

    @pytest.mark.xfail(platform.system() == 'Windows',
                       reason='https://github.com/astropy/astropy/issues/5797')
    def test_update_resized_header(self):
        """
        Test saving updates to a file where the header is one block smaller
        than before, and in the case where the heade ris one block larger than
        before.
        """

        data = np.arange(100)
        hdu = fits.PrimaryHDU(data=data)
        idx = 1
        while len(str(hdu.header)) <= 2880:
            hdu.header['TEST{}'.format(idx)] = idx
            idx += 1
        orig_header = hdu.header.copy()
        hdu.writeto(self.temp('temp.fits'))

        with fits.open(self.temp('temp.fits'), mode='update') as hdul:
            while len(str(hdul[0].header)) > 2880:
                del hdul[0].header[-1]

        with fits.open(self.temp('temp.fits')) as hdul:
            assert hdul[0].header == orig_header[:-1]
            assert (hdul[0].data == data).all()

        with fits.open(self.temp('temp.fits'), mode='update') as hdul:
            idx = 101
            while len(str(hdul[0].header)) <= 2880 * 2:
                hdul[0].header['TEST{}'.format(idx)] = idx
                idx += 1
            # Touch something in the data too so that it has to be rewritten
            hdul[0].data[0] = 27

        with fits.open(self.temp('temp.fits')) as hdul:
            assert hdul[0].header[:-37] == orig_header[:-1]
            assert hdul[0].data[0] == 27
            assert (hdul[0].data[1:] == data[1:]).all()

    def test_update_resized_header2(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/150

        This is similar to test_update_resized_header, but specifically tests a
        case of multiple consecutive flush() calls on the same HDUList object,
        where each flush() requires a resize.
        """

        data1 = np.arange(100)
        data2 = np.arange(100) + 100
        phdu = fits.PrimaryHDU(data=data1)
        hdu = fits.ImageHDU(data=data2)

        phdu.writeto(self.temp('temp.fits'))

        with fits.open(self.temp('temp.fits'), mode='append') as hdul:
            hdul.append(hdu)

        with fits.open(self.temp('temp.fits'), mode='update') as hdul:
            idx = 1
            while len(str(hdul[0].header)) <= 2880 * 2:
                hdul[0].header['TEST{}'.format(idx)] = idx
                idx += 1
            hdul.flush()
            hdul.append(hdu)

        with fits.open(self.temp('temp.fits')) as hdul:
            assert (hdul[0].data == data1).all()
            assert hdul[1].header == hdu.header
            assert (hdul[1].data == data2).all()
            assert (hdul[2].data == data2).all()

    @ignore_warnings()
    def test_hdul_fromstring(self):
        """
        Test creating the HDUList structure in memory from a string containing
        an entire FITS file.  This is similar to test_hdu_fromstring but for an
        entire multi-extension FITS file at once.
        """

        # Tests HDUList.fromstring for all of Astropy's built in test files
        def test_fromstring(filename):
            with fits.open(filename) as hdul:
                orig_info = hdul.info(output=False)
                with open(filename, 'rb') as f:
                    dat = f.read()

                hdul2 = fits.HDUList.fromstring(dat)

                assert orig_info == hdul2.info(output=False)
                for idx in range(len(hdul)):
                    assert hdul[idx].header == hdul2[idx].header
                    if hdul[idx].data is None or hdul2[idx].data is None:
                        assert hdul[idx].data == hdul2[idx].data
                    elif (hdul[idx].data.dtype.fields and
                          hdul2[idx].data.dtype.fields):
                        # Compare tables
                        for n in hdul[idx].data.names:
                            c1 = hdul[idx].data[n]
                            c2 = hdul2[idx].data[n]
                            assert (c1 == c2).all()
                    elif (any(dim == 0 for dim in hdul[idx].data.shape) or
                          any(dim == 0 for dim in hdul2[idx].data.shape)):
                        # For some reason some combinations of Python and Numpy
                        # on Windows result in MemoryErrors when trying to work
                        # on memmap arrays with more than one dimension but
                        # some dimensions of size zero, so include a special
                        # case for that
                        return hdul[idx].data.shape == hdul2[idx].data.shape
                    else:
                        np.testing.assert_array_equal(hdul[idx].data,
                                                      hdul2[idx].data)

        for filename in glob.glob(os.path.join(self.data_dir, '*.fits')):
            if sys.platform == 'win32' and filename == 'zerowidth.fits':
                # Running this test on this file causes a crash in some
                # versions of Numpy on Windows.  See ticket:
                # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/174
                continue
            elif filename.endswith('variable_length_table.fits'):
                # Comparing variable length arrays is non-trivial and thus
                # skipped at this point.
                # TODO: That's probably possible, so one could make it work.
                continue
            test_fromstring(filename)

        # Test that creating an HDUList from something silly raises a TypeError
        pytest.raises(TypeError, fits.HDUList.fromstring, ['a', 'b', 'c'])

    def test_save_backup(self):
        """Test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/121

        Save backup of file before flushing changes.
        """

        self.copy_file('scale.fits')

        with ignore_warnings():
            with fits.open(self.temp('scale.fits'), mode='update',
                           save_backup=True) as hdul:
                # Make some changes to the original file to force its header
                # and data to be rewritten
                hdul[0].header['TEST'] = 'TEST'
                hdul[0].data[0] = 0

        assert os.path.exists(self.temp('scale.fits.bak'))
        with fits.open(self.data('scale.fits'),
                       do_not_scale_image_data=True) as hdul1:
            with fits.open(self.temp('scale.fits.bak'),
                           do_not_scale_image_data=True) as hdul2:
                assert hdul1[0].header == hdul2[0].header
                assert (hdul1[0].data == hdul2[0].data).all()

        with ignore_warnings():
            with fits.open(self.temp('scale.fits'), mode='update',
                           save_backup=True) as hdul:
                # One more time to see if multiple backups are made
                hdul[0].header['TEST2'] = 'TEST'
                hdul[0].data[0] = 1

        assert os.path.exists(self.temp('scale.fits.bak'))
        assert os.path.exists(self.temp('scale.fits.bak.1'))

    def test_replace_mmap_data(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/25

        Replacing the mmap'd data of one file with mmap'd data from a
        different file should work.
        """

        arr_a = np.arange(10)
        arr_b = arr_a * 2

        def test(mmap_a, mmap_b):
            hdu_a = fits.PrimaryHDU(data=arr_a)
            hdu_a.writeto(self.temp('test_a.fits'), overwrite=True)
            hdu_b = fits.PrimaryHDU(data=arr_b)
            hdu_b.writeto(self.temp('test_b.fits'), overwrite=True)

            with fits.open(self.temp('test_a.fits'), mode='update',
                           memmap=mmap_a) as hdul_a:
                with fits.open(self.temp('test_b.fits'),
                               memmap=mmap_b) as hdul_b:
                    hdul_a[0].data = hdul_b[0].data

            with fits.open(self.temp('test_a.fits')) as hdul_a:
                assert np.all(hdul_a[0].data == arr_b)

        with ignore_warnings():
            test(True, True)

            # Repeat the same test but this time don't mmap A
            test(False, True)

            # Finally, without mmaping B
            test(True, False)

    def test_replace_mmap_data_2(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/25

        Replacing the mmap'd data of one file with mmap'd data from a
        different file should work.  Like test_replace_mmap_data but with
        table data instead of image data.
        """

        arr_a = np.arange(10)
        arr_b = arr_a * 2

        def test(mmap_a, mmap_b):
            col_a = fits.Column(name='a', format='J', array=arr_a)
            col_b = fits.Column(name='b', format='J', array=arr_b)
            hdu_a = fits.BinTableHDU.from_columns([col_a])
            hdu_a.writeto(self.temp('test_a.fits'), overwrite=True)
            hdu_b = fits.BinTableHDU.from_columns([col_b])
            hdu_b.writeto(self.temp('test_b.fits'), overwrite=True)

            with fits.open(self.temp('test_a.fits'), mode='update',
                           memmap=mmap_a) as hdul_a:
                with fits.open(self.temp('test_b.fits'),
                               memmap=mmap_b) as hdul_b:
                    hdul_a[1].data = hdul_b[1].data

            with fits.open(self.temp('test_a.fits')) as hdul_a:
                assert 'b' in hdul_a[1].columns.names
                assert 'a' not in hdul_a[1].columns.names
                assert np.all(hdul_a[1].data['b'] == arr_b)

        with ignore_warnings():
            test(True, True)

            # Repeat the same test but this time don't mmap A
            test(False, True)

            # Finally, without mmaping B
            test(True, False)

    def test_extname_in_hdulist(self):
        """
        Tests to make sure that the 'in' operator works.

        Regression test for https://github.com/astropy/astropy/issues/3060
        """
        hdulist = fits.open(self.data('o4sp040b0_raw.fits'))
        hdulist.append(fits.ImageHDU(name='a'))

        assert 'a' in hdulist
        assert 'A' in hdulist
        assert ('a', 1) in hdulist
        assert ('A', 1) in hdulist
        assert 'b' not in hdulist
        assert ('a', 2) not in hdulist
        assert ('b', 1) not in hdulist
        assert ('b', 2) not in hdulist
        assert hdulist[0] in hdulist
        assert fits.ImageHDU() not in hdulist

    def test_overwrite_vs_clobber(self):
        hdulist = fits.HDUList([fits.PrimaryHDU()])
        hdulist.writeto(self.temp('test_overwrite.fits'))
        hdulist.writeto(self.temp('test_overwrite.fits'), overwrite=True)
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            hdulist.writeto(self.temp('test_overwrite.fits'), clobber=True)
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) == '"clobber" was '
                    'deprecated in version 2.0 and will be removed in a '
                    'future version. Use argument "overwrite" instead.')

    def test_invalid_hdu_key_in_contains(self):
        """
        Make sure invalid keys in the 'in' operator return False.
        Regression test for https://github.com/astropy/astropy/issues/5583
        """
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist.append(fits.ImageHDU())
        hdulist.append(fits.ImageHDU())

        # A more or less random assortment of things which are not valid keys.
        bad_keys = [None, 3.5, {}]

        for key in bad_keys:
            assert not (key in hdulist)

    def test_iteration_of_lazy_loaded_hdulist(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5585
        """
        hdulist = fits.HDUList(fits.PrimaryHDU())
        hdulist.append(fits.ImageHDU(name='SCI'))
        hdulist.append(fits.ImageHDU(name='SCI'))
        hdulist.append(fits.ImageHDU(name='nada'))
        hdulist.append(fits.ImageHDU(name='SCI'))

        filename = self.temp('many_extension.fits')
        hdulist.writeto(filename)
        f = fits.open(filename)

        # Check that all extensions are read if f is not sliced
        all_exts = [ext for ext in f]
        assert len(all_exts) == 5

        # Reload the file to ensure we are still lazy loading
        f.close()
        f = fits.open(filename)

        # Try a simple slice with no conditional on the ext. This is essentially
        # the reported failure.
        all_exts_but_zero = [ext for ext in f[1:]]
        assert len(all_exts_but_zero) == 4

        # Reload the file to ensure we are still lazy loading
        f.close()
        f = fits.open(filename)

        # Check whether behavior is proper if the upper end of the slice is not
        # omitted.
        read_exts = [ext for ext in f[1:4] if ext.header['EXTNAME'] == 'SCI']
        assert len(read_exts) == 2

    def test_proper_error_raised_on_non_fits_file_with_unicode(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5594

        The failure shows up when (in python 3+) you try to open a file
        with unicode content that is not actually a FITS file. See:
        https://github.com/astropy/astropy/issues/5594#issuecomment-266583218
        """
        import codecs
        filename = self.temp('not-fits-with-unicode.fits')
        with codecs.open(filename, mode='w', encoding='utf=8') as f:
            f.write(u'Ce\xe7i ne marche pas')

        # This should raise an OSError because there is no end card.
        with pytest.raises(OSError):
            with pytest.warns(AstropyUserWarning, match='non-ASCII characters '
                              'are present in the FITS file header'):
                fits.open(filename)

    def test_no_resource_warning_raised_on_non_fits_file(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/6168

        The ResourceWarning shows up when (in python 3+) you try to
        open a non-FITS file when using a filename.
        """

        # To avoid creating the file multiple times the tests are
        # all included in one test file. See the discussion to the
        # PR at https://github.com/astropy/astropy/issues/6168
        #
        filename = self.temp('not-fits.fits')
        with open(filename, mode='w') as f:
            f.write('# header line\n')
            f.write('0.1 0.2\n')

        # Opening the file should raise an OSError however the file
        # is opened (there are two distinct code paths, depending on
        # whether ignore_missing_end is True or False).
        #
        # Explicit tests are added to make sure the file handle is not
        # closed when passed in to fits.open. In this case the ResourceWarning
        # was not raised, but a check is still included.
        #
        with catch_warnings(ResourceWarning) as ws:

            # Make sure that files opened by the user are not closed
            with open(filename, mode='rb') as f:
                with pytest.raises(OSError):
                    fits.open(f, ignore_missing_end=False)

                assert not f.closed

            with open(filename, mode='rb') as f:
                with pytest.raises(OSError):
                    fits.open(f, ignore_missing_end=True)

                assert not f.closed

            with pytest.raises(OSError):
                fits.open(filename, ignore_missing_end=False)

            with pytest.raises(OSError):
                fits.open(filename, ignore_missing_end=True)

        assert len(ws) == 0

    def test_pop_with_lazy_load(self):
        filename = self.data('checksum.fits')
        hdul = fits.open(filename)
        # Try popping the hdulist before doing anything else. This makes sure
        # that https://github.com/astropy/astropy/issues/7185 is fixed.
        hdu = hdul.pop()
        assert len(hdul) == 1

        # Read the file again and try popping from the beginning
        hdul2 = fits.open(filename)
        hdu2 = hdul2.pop(0)
        assert len(hdul2) == 1

        # Just a sanity check
        hdul3 = fits.open(filename)
        assert len(hdul3) == 2
        assert hdul3[0].header == hdu2.header
        assert hdul3[1].header == hdu.header

    def test_pop_extname(self):
        hdul = fits.open(self.data('o4sp040b0_raw.fits'))
        assert len(hdul) == 7
        hdu1 = hdul[1]
        hdu4 = hdul[4]
        hdu_popped = hdul.pop(('SCI', 2))
        assert len(hdul) == 6
        assert hdu_popped is hdu4
        hdu_popped = hdul.pop('SCI')
        assert len(hdul) == 5
        assert hdu_popped is hdu1

    def test_write_hdulist_to_stream(self):
        """
        Unit test for https://github.com/astropy/astropy/issues/7435
        Ensure that an HDUList can be written to a stream in Python 2
        """
        data = np.array([[1,2,3],[4,5,6]])
        hdu = fits.PrimaryHDU(data)
        hdulist = fits.HDUList([hdu])

        with open(self.temp('test.fits'), 'wb') as fout:
            p = subprocess.Popen(["cat"], stdin=subprocess.PIPE, stdout=fout)
            hdulist.writeto(p.stdin)
