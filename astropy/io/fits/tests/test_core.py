from __future__ import division # confidence high
from __future__ import with_statement

import gzip
import os
import warnings
import zipfile

import numpy as np

import pyfits
from pyfits.util import BytesIO
from pyfits.tests import PyfitsTestCase
from pyfits.tests.util import catch_warnings

from nose.tools import assert_equal, assert_raises, assert_true, assert_false


class TestCore(PyfitsTestCase):
    def test_with_statement(self):
        with pyfits.open(self.data('ascii.fits')) as f:
            pass

    def test_naxisj_check(self):
        hdulist = pyfits.open(self.data('o4sp040b0_raw.fits'))

        hdulist[1].header.update("NAXIS3", 500)

        assert 'NAXIS3' in hdulist[1].header
        hdulist.verify('fix')
        assert 'NAXIS3' not in hdulist[1].header

    def test_byteswap(self):
        p = pyfits.PrimaryHDU()
        l = pyfits.HDUList()

        n = np.zeros(3, dtype='i2')
        n[0] = 1
        n[1] = 60000
        n[2] = 2

        c = pyfits.Column(name='foo', format='i2', bscale=1, bzero=32768,
                          array=n)
        t = pyfits.new_table([c])

        l.append(p)
        l.append(t)

        l.writeto(self.temp('test.fits'), clobber=True)

        p = pyfits.open(self.temp('test.fits'))
        assert_equal(p[1].data[1]['foo'], 60000.0)

    def test_add_del_columns(self):
        p = pyfits.ColDefs([])
        p.add_col(pyfits.Column(name='FOO', format='3J'))
        p.add_col(pyfits.Column(name='BAR', format='1I'))
        assert_equal(p.names, ['FOO', 'BAR'])
        p.del_col('FOO')
        assert_equal(p.names, ['BAR'])

    def test_add_del_columns2(self):
        hdulist = pyfits.open(self.data('tb.fits'))
        table = hdulist[1]
        assert_equal(table.data.dtype.names, ('c1', 'c2', 'c3', 'c4'))
        assert_equal(table.columns.names, ['c1', 'c2', 'c3', 'c4'])
        #old_data = table.data.base.copy().view(pyfits.FITS_rec)
        table.columns.del_col('c1')
        assert_equal(table.data.dtype.names, ('c2', 'c3', 'c4'))
        assert_equal(table.columns.names, ['c2', 'c3', 'c4'])

        #for idx in range(len(old_data)):
        #    assert_equal(np.all(old_data[idx][1:], table.data[idx]))

        table.columns.del_col('c3')
        assert_equal(table.data.dtype.names, ('c2', 'c4'))
        assert_equal(table.columns.names, ['c2', 'c4'])

        #for idx in range(len(old_data)):
        #    assert_equal(np.all(old_data[idx][2:], table.data[idx]))

        table.columns.add_col(pyfits.Column('foo', '3J'))
        assert_equal(table.data.dtype.names, ('c2', 'c4', 'foo'))
        assert_equal(table.columns.names, ['c2', 'c4', 'foo'])

        hdulist.writeto(self.temp('test.fits'), clobber=True)
        hdulist = pyfits.open(self.temp('test.fits'))
        table = hdulist[1]
        assert_equal(table.data.dtype.names, ('c1', 'c2', 'c3'))
        assert_equal(table.columns.names, ['c1', 'c2', 'c3'])

    def test_update_header_card(self):
        """A very basic test for the Header.update method--I'd like to add a
        few more cases to this at some point.
        """

        header = pyfits.Header()
        comment = 'number of bits per data pixel'
        header.update('BITPIX', 16, comment)
        assert 'BITPIX' in header
        assert_equal(header['BITPIX'], 16)
        assert_equal(header.ascard['BITPIX'].comment, comment)

        header.update('BITPIX', 32, savecomment=True)
        # Make sure the value has been updated, but the comment was preserved
        assert_equal(header['BITPIX'], 32)
        assert_equal(header.ascard['BITPIX'].comment, comment)

        # The comment should still be preserved--savecomment only takes effect if
        # a new comment is also specified
        header.update('BITPIX', 16)
        assert_equal(header.ascard['BITPIX'].comment, comment)
        header.update('BITPIX', 16, 'foobarbaz', savecomment=True)
        assert_equal(header.ascard['BITPIX'].comment, comment)

    def test_set_card_value(self):
        """Similar to test_update_header_card(), but tests the the
        `header['FOO'] = 'bar'` method of updating card values.
        """

        header = pyfits.Header()
        comment = 'number of bits per data pixel'
        card = pyfits.Card.fromstring('BITPIX  = 32 / %s' % comment)
        header.ascard.append(card)

        header['BITPIX'] = 32

        assert 'BITPIX' in header
        assert_equal(header['BITPIX'], 32)
        assert_equal(header.ascard['BITPIX'].key, 'BITPIX')
        assert_equal(header.ascard['BITPIX'].value, 32)
        assert_equal(header.ascard['BITPIX'].comment, comment)

    def test_uint(self):
        hdulist_f = pyfits.open(self.data('o4sp040b0_raw.fits'))
        hdulist_i = pyfits.open(self.data('o4sp040b0_raw.fits'), uint=True)

        assert_equal(hdulist_f[1].data.dtype, np.float32)
        assert_equal(hdulist_i[1].data.dtype, np.uint16)
        assert np.all(hdulist_f[1].data == hdulist_i[1].data)

    def test_fix_missing_card_append(self):
        hdu = pyfits.ImageHDU()
        errs = hdu.req_cards('TESTKW', None, None, 'foo', 'silentfix', [])
        assert_equal(len(errs), 1)
        assert 'TESTKW' in hdu.header
        assert_equal(hdu.header['TESTKW'], 'foo')
        assert_equal(hdu.header.ascard[-1].key, 'TESTKW')

    def test_fix_invalid_keyword_value(self):
        hdu = pyfits.ImageHDU()
        hdu.header.update('TESTKW', 'foo')
        errs = hdu.req_cards('TESTKW', None,
                             lambda v: v == 'foo', 'foo', 'ignore', [])
        assert_equal(len(errs), 0)

        # Now try a test that will fail, and ensure that an error will be
        # raised in 'exception' mode
        errs = hdu.req_cards('TESTKW', None, lambda v: v == 'bar', 'bar',
                             'exception', [])
        assert_equal(len(errs), 1)
        assert_equal(errs[0], "'TESTKW' card has invalid value 'foo'.")

        # See if fixing will work
        hdu.req_cards('TESTKW', None, lambda v: v == 'bar', 'bar', 'silentfix',
                      [])
        assert_equal(hdu.header['TESTKW'], 'bar')

    def test_unfixable_missing_card(self):
        class TestHDU(pyfits.hdu.base.NonstandardExtHDU):
            def _verify(self, option='warn'):
                errs = super(TestHDU, self)._verify(option)
                hdu.req_cards('TESTKW', None, None, None, 'fix', errs)
                return errs

        hdu = TestHDU(header=pyfits.Header())
        assert_raises(pyfits.VerifyError, hdu.verify, 'fix')

    def test_exception_on_verification_error(self):
        hdu = pyfits.ImageHDU()
        del hdu.header['NAXIS']
        assert_raises(pyfits.VerifyError, hdu.verify, 'exception')

    def test_ignore_verification_error(self):
        hdu = pyfits.ImageHDU()
        # The default here would be to issue a warning; ensure that no warnings
        # or exceptions are raised
        with catch_warnings():
            warnings.simplefilter('error')
            del hdu.header['NAXIS']
            try:
                hdu.verify('ignore')
            except Exception, e:
                self.fail('An exception occurred when the verification error '
                          'should have been ignored: %s' % e)
        # Make sure the error wasn't fixed either, silently or otherwise
        assert_false('NAXIS' in hdu.header)

    def test_unrecognized_verify_option(self):
        hdu = pyfits.ImageHDU()
        assert_raises(ValueError, hdu.verify, 'foobarbaz')


class TestFileFunctions(PyfitsTestCase):
    """Tests various basic I/O operations, specifically in the
    pyfits.file._File class.
    """

    def test_open_gzipped(self):
        assert_equal(len(pyfits.open(self._make_gzip_file())), 5)

    def test_open_gzipped_writeable(self):
        """Opening gzipped files in a writeable mode should fail."""

        gf = self._make_gzip_file()
        assert_raises(IOError, pyfits.open, gf, 'update')
        assert_raises(IOError, pyfits.open, gf, 'append')

    def test_open_zipped(self):
        assert_equal(len(pyfits.open(self._make_zip_file())), 5)

    def test_open_zipped_writeable(self):
        """Opening zipped files in a writeable mode should fail."""

        zf = self._make_zip_file()
        assert_raises(IOError, pyfits.open, zf, 'update')
        assert_raises(IOError, pyfits.open, zf, 'append')

    def test_open_multipe_member_zipfile(self):
        """Opening zip files containing more than one member files should fail
        as there's no obvious way to specify which file is the FITS file to
        read.
        """

        zfile = zipfile.ZipFile(self.temp('test0.zip'), 'w')
        zfile.write(self.data('test0.fits'))
        zfile.writestr('foo', 'bar')
        zfile.close()

        assert_raises(IOError, pyfits.open, zfile.filename)

    def test_read_open_file(self):
        """Read from an existing file object."""

        with open(self.data('test0.fits'), 'rb') as f:
            assert_equal(len(pyfits.open(f)), 5)

    def test_read_closed_file(self):
        """Read from an existing file object that's been closed."""

        f = open(self.data('test0.fits'), 'rb')
        f.close()
        assert_equal(len(pyfits.open(f)), 5)

    def test_read_open_gzip_file(self):
        """Read from an open gzip file object."""

        gf = gzip.GzipFile(self._make_gzip_file())
        try:
            assert_equal(len(pyfits.open(gf)), 5)
        finally:
            gf.close()

    def test_read_file_like_object(self):
        """Test reading a FITS file from a file-like object."""

        filelike = BytesIO()
        with open(self.data('test0.fits'), 'rb') as f:
            filelike.write(f.read())
        filelike.seek(0)
        assert_equal(len(pyfits.open(filelike)), 5)

    def test_updated_file_permissions(self):
        """
        Regression test for #79.  Tests that when a FITS file is modified in
        update mode, the file permissions are preserved.
        """

        filename = self.temp('test.fits')
        hdul = [pyfits.PrimaryHDU(), pyfits.ImageHDU()]
        hdul = pyfits.HDUList(hdul)
        hdul.writeto(filename)

        old_mode = os.stat(filename).st_mode

        hdul = pyfits.open(filename, mode='update')
        hdul.insert(1, pyfits.ImageHDU())
        hdul.flush()
        hdul.close()

        assert_equal(old_mode, os.stat(filename).st_mode)

    def _make_gzip_file(self):
        gzfile = self.temp('test0.fits.gz')
        with open(self.data('test0.fits'), 'rb') as f:
            gz = gzip.open(gzfile, 'wb')
            gz.write(f.read())
            gz.close()

        return gzfile

    def _make_zip_file(self, mode='copyonwrite'):
        zfile = zipfile.ZipFile(self.temp('test0.fits.zip'), 'w')
        zfile.write(self.data('test0.fits'))
        zfile.close()

        return zfile.filename


class TestStreamingFunctions(PyfitsTestCase):
    """Test functionality of the StreamingHDU class."""

    def test_streaming_hdu(self):
        shdu = self._make_streaming_hdu(self.temp('new.fits'))
        assert_true(isinstance(shdu.size, int))
        assert_equal(shdu.size, 100)

    def test_streaming_hdu_file_wrong_mode(self):
        """Test that streaming an HDU to a file opened in the wrong mode
        fails as expected.
        """

        with open(self.temp('new.fits'), 'wb') as f:
            header = pyfits.Header()
            assert_raises(ValueError, pyfits.StreamingHDU, f, header)

    def test_streaming_hdu_write_file(self):
        """Test streaming an HDU to an open file object."""

        arr = np.zeros((5, 5), dtype=np.int32)
        with open(self.temp('new.fits'), 'ab+') as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)
            assert_true(shdu.writecomplete)
            assert_equal(shdu.size, 100)
        hdul = pyfits.open(self.temp('new.fits'))
        assert_equal(len(hdul), 1)
        assert_true((hdul[0].data == arr).all())

    def test_streaming_hdu_write_file_like(self):
        """Test streaming an HDU to an open file-like object."""

        arr = np.zeros((5, 5), dtype=np.int32)
        # The file-like object underlying a StreamingHDU must be in binary mode
        sf = BytesIO()
        shdu = self._make_streaming_hdu(sf)
        shdu.write(arr)
        assert_true(shdu.writecomplete)
        assert_equal(shdu.size, 100)

        sf.seek(0)
        hdul = pyfits.open(sf)
        assert_equal(len(hdul), 1)
        assert_true((hdul[0].data == arr).all())

    def test_streaming_hdu_append_extension(self):
        arr = np.zeros((5, 5), dtype=np.int32)
        with open(self.temp('new.fits'), 'ab+') as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)
        # Doing this again should update the file with an extension
        with open(self.temp('new.fits'), 'ab+') as f:
            shdu = self._make_streaming_hdu(f)
            shdu.write(arr)

    # TODO: This test is temporarily borrowed from the header-refactoring
    # branch as a regression test for ticket #69; it can be removed when the
    # header-refactoring branch is merged into trunk
    def test_update_comment(self):
        hdul = pyfits.open(self.data('arange.fits'))
        hdul[0].header.update('FOO', 'BAR', 'BAZ')
        hdul.writeto(self.temp('test.fits'))

        hdul = pyfits.open(self.temp('test.fits'), mode='update')
        hdul[0].header.ascard['FOO'].comment = 'QUX'
        hdul.close()

        hdul = pyfits.open(self.temp('test.fits'))
        assert_equal(hdul[0].header.ascard['FOO'].comment, 'QUX')

    def _make_streaming_hdu(self, fileobj):
        hd = pyfits.Header()
        hd.update('SIMPLE', True, 'conforms to FITS standard')
        hd.update('BITPIX', 32, 'array data type')
        hd.update('NAXIS', 2, 'number of array dimensions')
        hd.update('NAXIS1', 5)
        hd.update('NAXIS2', 5)
        hd.update('EXTEND', True)
        return pyfits.StreamingHDU(fileobj, hd)
