# Licensed under a 3-clause BSD style license - see PYFITS.rst

import glob
import io
import os
import sys

import numpy as np

from ..verify import VerifyError
from ....io import fits
from ....tests.helper import pytest, raises, catch_warnings, ignore_warnings
from ....utils.exceptions import AstropyUserWarning

from . import FitsTestCase


class TestHDUListFunctions(FitsTestCase):
    @ignore_warnings()
    def test_update_name(self):
        hdul = fits.open(self.data('o4sp040b0_raw.fits'))
        hdul[4].update_ext_name('Jim', "added by Jim")
        hdul[4].update_ext_version(9, "added by Jim")
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
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_extension_to_empty_list(self):
        """Tests appending a Simple ImageHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (100,), 'int32', '')]
        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_table_extension_to_empty_list(self):
        """Tests appending a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = fits.HDUList()
        hdul1 = fits.open(self.data('tb.fits'))
        hdul.append(hdul1[1])
        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), '', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_groupshdu_to_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_primary_to_non_empty_list(self):
        """Tests appending a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('arange.fits'))
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_extension_to_non_empty_list(self):
        """Tests appending a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('tb.fits'))
        hdul.append(hdul[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), '', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_extension_to_empty_list(self):
        """Tests inserting a Simple ImageHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_table_extension_to_empty_list(self):
        """Tests inserting a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = fits.HDUList()
        hdul1 = fits.open(self.data('tb.fits'))
        hdul.insert(0, hdul1[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), '', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_groupshdu_to_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_primary_to_non_empty_list(self):
        """Tests inserting a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('arange.fits'))
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(1, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_extension_to_non_empty_list(self):
        """Tests inserting a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = fits.open(self.data('tb.fits'))
        hdul.insert(1, hdul[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), '', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

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

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), '',
                 '1 Groups  0 Parameters'),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), '', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'ImageHDU', 12, (), '', ''),
                (3, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', ''),
                (1, '', 'ImageHDU', 12, (), '', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_file_like_2(self):
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        tmpfile = open(self.temp('tmpfile.fits'), 'wb')
        hdul = fits.open(tmpfile, mode='ostream')
        hdul.append(hdu)
        hdul.flush()
        tmpfile.close()
        hdul.close()

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_file_like_3(self):
        tmpfile = open(self.temp('tmpfile.fits'), 'wb')
        fits.writeto(tmpfile, np.arange(100, dtype=np.int32))
        tmpfile.close()
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_new_hdu_extname(self):
        """
        Tests that new extension HDUs that are added to an HDUList can be
        properly indexed by their EXTNAME/EXTVER (regression test for
        ticket:48).
        """

        f = fits.open(self.data('test0.fits'))
        hdul = fits.HDUList()
        hdul.append(f[0].copy())
        hdul.append(fits.ImageHDU(header=f[1].header))

        assert hdul[1].header['EXTNAME'] == 'SCI'
        assert hdul[1].header['EXTVER'] == 1
        assert hdul.index_of(('SCI', 1)) == 1

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
        assert hdul[0].header['EXTEND'] == True

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
        assert hdu.header['EXTEND'] == True

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

        a = np.arange(100).reshape((10, 10))
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
            hdu.header['TEST%d' % idx] = idx
            idx += 1
        hdu.writeto(self.temp('temp.fits'), checksum=True)

        with fits.open(self.temp('temp.fits'), mode='update') as hdul:
            # Modify the header, forcing it to be rewritten
            hdul[0].header['TEST1'] = 2

        with fits.open(self.temp('temp.fits')) as hdul:
            assert (hdul[0].data == data).all()

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
            hdu.header['TEST%d' % idx] = idx
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
                hdul[0].header['TEST%d' % idx] = idx
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
                hdul[0].header['TEST%d' % idx] = idx
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

        # Tests HDUList.fromstring for all of PyFITS' built in test files
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
                # versions of Numpy on Windows.  See PyFITS ticket
                # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/174
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
            hdu_a.writeto(self.temp('test_a.fits'), clobber=True)
            hdu_b = fits.PrimaryHDU(data=arr_b)
            hdu_b.writeto(self.temp('test_b.fits'), clobber=True)

            hdul_a = fits.open(self.temp('test_a.fits'), mode='update',
                               memmap=mmap_a)
            hdul_b = fits.open(self.temp('test_b.fits'), memmap=mmap_b)
            hdul_a[0].data = hdul_b[0].data
            hdul_a.close()
            hdul_b.close()

            hdul_a = fits.open(self.temp('test_a.fits'))

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
            hdu_a.writeto(self.temp('test_a.fits'), clobber=True)
            hdu_b = fits.BinTableHDU.from_columns([col_b])
            hdu_b.writeto(self.temp('test_b.fits'), clobber=True)

            hdul_a = fits.open(self.temp('test_a.fits'), mode='update',
                               memmap=mmap_a)
            hdul_b = fits.open(self.temp('test_b.fits'), memmap=mmap_b)
            hdul_a[1].data = hdul_b[1].data
            hdul_a.close()
            hdul_b.close()

            hdul_a = fits.open(self.temp('test_a.fits'))

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

        hdulist = fits.HDUList()
        hdulist.append(fits.ImageHDU(name='a'))

        assert 'a' in hdulist
        assert 'A' in hdulist
        assert ('a', 1) in hdulist
        assert ('A', 1) in hdulist
        assert 'b' not in hdulist
        assert ('a', 2) not in hdulist
        assert ('b', 1) not in hdulist
        assert ('b', 2) not in hdulist
