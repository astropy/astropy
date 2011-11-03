from __future__ import division # confidence high
from __future__ import with_statement

import os

import numpy as np

import pyfits
from pyfits.util import BytesIO
from pyfits.tests import PyfitsTestCase
from pyfits.tests.util import catch_warnings

from nose.tools import assert_equal, assert_raises, assert_true


class TestHDUListFunctions(PyfitsTestCase):
    def test_update_name(self):
        hdul = pyfits.open(self.data('o4sp040b0_raw.fits'))
        hdul[4].update_ext_name('Jim', "added by Jim")
        hdul[4].update_ext_version(9, "added by Jim")
        assert_equal(hdul[('JIM', 9)].header['extname'], 'JIM')

    def test_hdu_file_bytes(self):
        hdul = pyfits.open(self.data('checksum.fits'))
        res = hdul[0].filebytes()
        assert_equal(res, 11520)
        res = hdul[1].filebytes()
        assert_equal(res, 8640)

    def test_hdulist_file_info(self):
        hdul = pyfits.open(self.data('checksum.fits'))
        res = hdul.fileinfo(0)

        def test_fileinfo(**kwargs):
            assert_equal(res['datSpan'], kwargs.get('datSpan', 2880))
            assert_equal(res['resized'], kwargs.get('resized', 0))
            assert_equal(res['filename'], self.data('checksum.fits'))
            assert_equal(res['datLoc'], kwargs.get('datLoc', 8640))
            assert_equal(res['hdrLoc'], kwargs.get('hdrLoc', 0))
            assert_equal(res['filemode'], 'copyonwrite')

        res = hdul.fileinfo(1)
        test_fileinfo(datLoc=17280, hdrLoc=11520)

        hdu = pyfits.ImageHDU(data=hdul[0].data)
        hdul.insert(1, hdu)

        res = hdul.fileinfo(0)
        test_fileinfo(resized=1)

        res = hdul.fileinfo(1)
        test_fileinfo(datSpan=None, resized=1, datLoc=None, hdrLoc=None)

        res = hdul.fileinfo(2)
        test_fileinfo(resized=1, datLoc=17280, hdrLoc=11520)

    def test_append_primary_to_empty_list(self):
        # Tests appending a Simple PrimaryHDU to an empty HDUList.
        hdul = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_extension_to_empty_list(self):
        """Tests appending a Simple ImageHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (100,), 'int32', '')]
        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_table_extension_to_empty_list(self):
        """Tests appending a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = pyfits.HDUList()
        hdul1 = pyfits.open(self.data('tb.fits'))
        hdul.append(hdul1[1])
        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_groupshdu_to_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.GroupsHDU()
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_primary_to_non_empty_list(self):
        """Tests appending a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = pyfits.open(self.data('arange.fits'))
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_extension_to_non_empty_list(self):
        """Tests appending a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = pyfits.open(self.data('tb.fits'))
        hdul.append(hdul[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), 'int16', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-append.fits'))

        assert_equal(pyfits.info(self.temp('test-append.fits'), output=False),
                     info)

    def test_append_groupshdu_to_non_empty_list(self):
        """Tests appending a Simple GroupsHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.append(hdu)
        hdu = pyfits.GroupsHDU()

        assert_raises(ValueError, hdul.append, hdu)

    def test_insert_primary_to_empty_list(self):
        """Tests inserting a Simple PrimaryHDU to an empty HDUList."""
        hdul = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_extension_to_empty_list(self):
        """Tests inserting a Simple ImageHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.ImageHDU(np.arange(100,dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (100,), 'int32', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_table_extension_to_empty_list(self):
        """Tests inserting a Simple Table ExtensionHDU to a empty HDUList."""

        hdul = pyfits.HDUList()
        hdul1 = pyfits.open(self.data('tb.fits'))
        hdul.insert(0, hdul1[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_groupshdu_to_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.GroupsHDU()
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_primary_to_non_empty_list(self):
        """Tests inserting a Simple PrimaryHDU to a non-empty HDUList."""

        hdul = pyfits.open(self.data('arange.fits'))
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(1, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 7, (11, 10, 7), 'int32', ''),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_extension_to_non_empty_list(self):
        """Tests inserting a Simple ExtensionHDU to a non-empty HDUList."""

        hdul = pyfits.open(self.data('tb.fits'))
        hdul.insert(1, hdul[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 11, (), 'int16', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_groupshdu_to_non_empty_list(self):
        """Tests inserting a Simple GroupsHDU to an empty HDUList."""

        hdul = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)
        hdu = pyfits.GroupsHDU()

        assert_raises(ValueError, hdul.insert, hdul, 1, hdu)

        info = [(0, 'PRIMARY', 'GroupsHDU', 8, (), 'uint8',
                 '   1 Groups  0 Parameters'),
                (1, '', 'ImageHDU', 6, (100,), 'int32', '')]

        hdul.insert(0, hdu)

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_groupshdu_to_begin_of_hdulist_with_groupshdu(self):
        """
        Tests inserting a Simple GroupsHDU to the beginning of an HDUList
        that that already contains a GroupsHDU.
        """

        hdul = pyfits.HDUList()
        hdu = pyfits.GroupsHDU()
        hdul.insert(0, hdu)

        assert_raises(ValueError, hdul.insert, hdul, 0, hdu)

    def test_insert_extension_to_primary_in_non_empty_list(self):
        # Tests inserting a Simple ExtensionHDU to a non-empty HDUList.
        hdul = pyfits.open(self.data('tb.fits'))
        hdul.insert(0, hdul[1])

        info = [(0, 'PRIMARY', 'PrimaryHDU', 4, (), 'uint8', ''),
                (1, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', ''),
                (2, '', 'ImageHDU', 12, (), 'uint8', ''),
                (3, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_insert_image_extension_to_primary_in_non_empty_list(self):
        """
        Tests inserting a Simple Image ExtensionHDU to a non-empty HDUList
        as the primary HDU.
        """

        hdul = pyfits.open(self.data('tb.fits'))
        hdu = pyfits.ImageHDU(np.arange(100, dtype=np.int32))
        hdul.insert(0, hdu)

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', ''),
                (1, '', 'ImageHDU', 12, (), 'uint8', ''),
                (2, '', 'BinTableHDU', 24, '2R x 4C', '[1J, 3A, 1E, 1L]', '')]

        assert_equal(hdul.info(output=False), info)

        hdul.writeto(self.temp('test-insert.fits'))

        assert_equal(pyfits.info(self.temp('test-insert.fits'), output=False),
                     info)

    def test_filename(self):
        """Tests the HDUList filename method."""

        hdul = pyfits.open(self.data('tb.fits'))
        name = hdul.filename()
        assert_equal(name, self.data('tb.fits'))

    def test_file_like(self):
        """
        Tests the use of a file like object with no tell or seek methods
        in HDUList.writeto(), HDULIST.flush() or pyfits.writeto()
        """

        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        hdul = pyfits.HDUList()
        hdul.append(hdu)
        tmpfile = open(self.temp('tmpfile.fits'), 'w')
        hdul.writeto(tmpfile)
        tmpfile.close()

        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]

        assert_equal(pyfits.info(self.temp('tmpfile.fits'), output=False),
                     info)

    def test_file_like_2(self):
        hdu = pyfits.PrimaryHDU(np.arange(100, dtype=np.int32))
        tmpfile = open(self.temp('tmpfile.fits'), 'w')
        hdul = pyfits.open(tmpfile, mode='ostream')
        hdul.append(hdu)
        hdul.flush()
        tmpfile.close()

        hdul2 = pyfits.open(self.temp('tmpfile.fits'))
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert_equal(hdul2.info(output=False), info)

    def test_file_like_3(self):

        tmpfile = open(self.temp('tmpfile.fits'), 'w')
        pyfits.writeto(tmpfile, np.arange(100, dtype=np.int32))
        tmpfile.close()
        info = [(0, 'PRIMARY', 'PrimaryHDU', 5, (100,), 'int32', '')]
        assert_equal(pyfits.info(self.temp('tmpfile.fits'), output=False),
                     info)

    def test_new_hdu_extname(self):
        """
        Tests that new extension HDUs that are added to an HDUList can be
        properly indexed by their EXTNAME/EXTVER (regression test for
        ticket:48).
        """

        f = pyfits.open(self.data('test0.fits'))
        hdul = pyfits.HDUList()
        hdul.append(f[0].copy())
        hdul.append(pyfits.ImageHDU(header=f[1].header))

        assert_equal(hdul[1].header['EXTNAME'], 'SCI')
        assert_equal(hdul[1].header['EXTVER'], 1)
        assert_equal(hdul.index_of(('SCI', 1)), 1)

    def test_update_filelike(self):
        """Test opening a file-like object in update mode and resizing the
        HDU.
        """

        sf = BytesIO()
        arr = np.zeros((100, 100))
        hdu = pyfits.PrimaryHDU(data=arr)
        hdu.writeto(sf)

        sf.seek(0)
        arr = np.zeros((200, 200))
        hdul = pyfits.open(sf, mode='update')
        hdul[0].data = arr
        hdul.flush()

        sf.seek(0)
        hdul = pyfits.open(sf)
        assert_equal(len(hdul), 1)
        assert_true((hdul[0].data == arr).all())

    def test_flush_readonly(self):
        """Test flushing changes to a file opened in a read only mode."""

        oldmtime = os.stat(self.data('test0.fits')).st_mtime
        hdul = pyfits.open(self.data('test0.fits'))
        hdul[0].header.update('FOO', 'BAR')
        with catch_warnings(record=True) as w:
            hdul.flush()
            assert_equal(len(w), 1)
            assert_true('mode is not supported' in str(w[0].message))
        assert_equal(oldmtime, os.stat(self.data('test0.fits')).st_mtime)

    def test_fix_extend_keyword(self):
        hdul = pyfits.HDUList()
        hdul.append(pyfits.PrimaryHDU())
        hdul.append(pyfits.ImageHDU())
        del hdul[0].header['EXTEND']
        hdul.verify('silentfix')

        assert_true('EXTEND' in hdul[0].header)
        assert_equal(hdul[0].header['EXTEND'], True)
