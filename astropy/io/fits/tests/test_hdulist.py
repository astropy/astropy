import io
import os
import warnings

import numpy as np

from ....io import fits
from ....tests.helper import pytest, raises

from . import FitsTestCase


class TestHDUListFunctions(FitsTestCase):
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
            assert res['resized'] == kwargs.get('resized', 0)
            assert res['filename'] == self.data('checksum.fits')
            assert res['datLoc'] == kwargs.get('datLoc', 8640)
            assert res['hdrLoc'] == kwargs.get('hdrLoc', 0)
            assert res['filemode'] == 'copyonwrite'

        res = hdul.fileinfo(1)
        test_fileinfo(datLoc=17280, hdrLoc=11520)

        hdu = fits.ImageHDU(data=hdul[0].data)
        hdul.insert(1, hdu)

        res = hdul.fileinfo(0)
        test_fileinfo(resized=1)

        res = hdul.fileinfo(1)
        test_fileinfo(datSpan=None, resized=1, datLoc=None, hdrLoc=None)

        res = hdul.fileinfo(2)
        test_fileinfo(resized=1, datLoc=17280, hdrLoc=11520)

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

        assert fits.info(self.temp('test-append.fits'), output=False) ==info

    def test_append_table_extension_to_empty_list(self):
        """Tests appending a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = fits.HDUList()
        hdul1 = fits.open(self.data('tb.fits'))
        hdul.append(hdul1[1])
        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-append.fits'))

        assert fits.info(self.temp('test-append.fits'), output=False) == info

    def test_append_groupshdu_to_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), 'int16', ''),
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
        hdu = fits.ImageHDU(np.arange(100,dtype=np.int32))
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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert hdul.info(output=False) == info

        hdul.writeto(self.temp('test-insert.fits'))

        assert fits.info(self.temp('test-insert.fits'), output=False) == info

    def test_insert_groupshdu_to_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = fits.HDUList()
        hdu = fits.GroupsHDU()
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters')]

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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), 'int16', ''),
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

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters'),
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

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'ImageHDU', 12, (), 'uint8', ''),
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
                (1, '', 'ImageHDU', 12, (), 'uint8', ''),
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
        tmpfile = open(self.temp('tmpfile.fits'), 'w')
        hdul.writeto(tmpfile)
        tmpfile.close()

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert fits.info(self.temp('tmpfile.fits'), output=False) == info

    def test_file_like_2(self):
        hdu = fits.PrimaryHDU(np.arange(100, dtype=np.int32))
        tmpfile = open(self.temp('tmpfile.fits'), 'w')
        hdul = fits.open(tmpfile, mode='ostream')
        hdul.append(hdu)
        hdul.flush()
        tmpfile.close()

        hdul2 = fits.open(self.temp('tmpfile.fits'))
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert hdul2.info(output=False) == info

    def test_file_like_3(self):

        tmpfile = open(self.temp('tmpfile.fits'), 'w')
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
        hdul[0].header.update('FOO', 'BAR')
        with warnings.catch_warnings(record=True) as w:
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
