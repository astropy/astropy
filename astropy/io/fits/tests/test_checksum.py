# Licensed under a 3-clause BSD style license - see PYFITS.rst

import sys
import warnings
import sys

import numpy as np

from .test_table import comparerecords
from ..hdu.base import _ValidHDU
from ....io import fits
from ....tests.helper import pytest

from . import FitsTestCase


class TestChecksumFunctions(FitsTestCase):

    # All checksums have been verified against CFITSIO
    def setup(self):
        super(TestChecksumFunctions, self).setup()
        self._oldfilters = warnings.filters[:]
        warnings.filterwarnings(
            'error',
            message='Checksum verification failed')
        warnings.filterwarnings(
            'error',
            message='Datasum verification failed')

        # Monkey-patch the _get_timestamp method so that the checksum
        # timestamps (and hence the checksum themselves) are always the same
        self._old_get_timestamp = _ValidHDU._get_timestamp
        _ValidHDU._get_timestamp = lambda self: '2013-12-20T13:36:10'

    def teardown(self):
        super(TestChecksumFunctions, self).teardown()
        warnings.filters = self._oldfilters
        _ValidHDU._get_timestamp = self._old_get_timestamp

    def test_sample_file(self):
        hdul = fits.open(self.data('checksum.fits'), checksum=True)
        hdul.close()

    def test_image_create(self):
        n = np.arange(100, dtype=np.int64)
        hdu = fits.PrimaryHDU(n)
        hdu.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert (hdu.data == hdul[0].data).all()
            assert 'CHECKSUM' in hdul[0].header
            assert 'DATASUM' in hdul[0].header

            if not sys.platform.startswith('win32'):
                # The checksum ends up being different on Windows, possibly due
                # to slight floating point differences
                assert hdul[0].header['CHECKSUM'] == 'ZHMkeGKjZGKjbGKj'
                assert hdul[0].header['DATASUM'] == '4950'

    def test_nonstandard_checksum(self):
        hdu = fits.PrimaryHDU(np.arange(10.0 ** 6, dtype=np.float64))
        hdu.writeto(self.temp('tmp.fits'), clobber=True,
                    checksum='nonstandard')
        del hdu
        with fits.open(self.temp('tmp.fits'), checksum='nonstandard') as hdul:
            assert 'CHECKSUM' in hdul[0].header
            assert 'DATASUM' in hdul[0].header

            if not sys.platform.startswith('win32'):
                # The checksum ends up being different on Windows, possibly due
                # to slight floating point differences
                assert hdul[0].header['CHECKSUM'] == 'jD4Am942jC48j948'
                assert hdul[0].header['DATASUM'] == '4164005614'

    def test_scaled_data(self):
        with fits.open(self.data('scale.fits')) as hdul:
            orig_data = hdul[0].data.copy()
            hdul[0].scale('int16', 'old')
            hdul.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
            with fits.open(self.temp('tmp.fits'), checksum=True) as hdul1:
                assert (hdul1[0].data == orig_data).all()
                assert 'CHECKSUM' in hdul1[0].header
                assert hdul1[0].header['CHECKSUM'] == 'cUmaeUjZcUjacUjW'
                assert 'DATASUM' in hdul1[0].header
                assert hdul1[0].header['DATASUM'] == '1891563534'

    def test_uint16_data(self):
        checksums = [
            ('aDcXaCcXaCcXaCcX', '0'), ('84DJ83CH83CH83CH', '1746888714'),
            ('VhqQWZoQVfoQVZoQ', '0'), ('4cPp5aOn4aOn4aOn', '0'),
            ('97ZYI7WV97WVG7WV', '1756785133'), ('UhqdUZnbUfnbUZnb', '0'),
            ('4cQJ5aN94aNG4aN9', '0')]
        with fits.open(self.data('o4sp040b0_raw.fits'), uint=True) as hdul:
            hdul.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
            with fits.open(self.temp('tmp.fits'), uint=True,
                           checksum=True) as hdul1:
                for idx, (hdu_a, hdu_b) in enumerate(zip(hdul, hdul1)):
                    if hdu_a.data is None or hdu_b.data is None:
                        assert hdu_a.data is hdu_b.data
                    else:
                        assert (hdu_a.data == hdu_b.data).all()

                    assert 'CHECKSUM' in hdul[idx].header
                    assert hdul[idx].header['CHECKSUM'] == checksums[idx][0]
                    assert 'DATASUM' in hdul[idx].header
                    assert hdul[idx].header['DATASUM'] == checksums[idx][1]

    def test_groups_hdu_data(self):
        imdata = np.arange(100.0)
        imdata.shape = (10, 1, 1, 2, 5)
        pdata1 = np.arange(10) + 0.1
        pdata2 = 42
        x = fits.hdu.groups.GroupData(imdata, parnames=[str('abc'), str('xyz')],
                                      pardata=[pdata1, pdata2], bitpix=-32)
        hdu = fits.GroupsHDU(x)
        hdu.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert comparerecords(hdul[0].data, hdu.data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == '3eDQAZDO4dDOAZDO'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '2797758084'

    def test_binary_table_data(self):
        a1 = np.array(['NGC1001', 'NGC1002', 'NGC1003'])
        a2 = np.array([11.1, 12.3, 15.2])
        col1 = fits.Column(name='target', format='20A', array=a1)
        col2 = fits.Column(name='V_mag', format='E', array=a2)
        cols = fits.ColDefs([col1, col2])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert comparerecords(tbhdu.data, hdul[1].data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'
            assert 'CHECKSUM' in hdul[1].header
            assert hdul[1].header['CHECKSUM'] == 'aD1Oa90MaC0Ma90M'
            assert 'DATASUM' in hdul[1].header
            assert hdul[1].header['DATASUM'] == '1062205743'

    def test_variable_length_table_data(self):
        c1 = fits.Column(name='var', format='PJ()',
                         array=np.array([[45.0, 56], np.array([11, 12, 13])],
                         'O'))
        c2 = fits.Column(name='xyz', format='2I', array=[[11, 3], [12, 4]])
        tbhdu = fits.BinTableHDU.from_columns([c1, c2])
        tbhdu.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert comparerecords(tbhdu.data, hdul[1].data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'
            assert 'CHECKSUM' in hdul[1].header
            assert hdul[1].header['CHECKSUM'] == 'YIGoaIEmZIEmaIEm'
            assert 'DATASUM' in hdul[1].header
            assert hdul[1].header['DATASUM'] == '1507485'

    def test_ascii_table_data(self):
        a1 = np.array(['abc', 'def'])
        r1 = np.array([11.0, 12.0])
        c1 = fits.Column(name='abc', format='A3', array=a1)
        # This column used to be E format, but the single-precision float lost
        # too much precision when scaling so it was changed to a D
        c2 = fits.Column(name='def', format='D', array=r1, bscale=2.3,
                         bzero=0.6)
        c3 = fits.Column(name='t1', format='I', array=[91, 92, 93])
        x = fits.ColDefs([c1, c2, c3])
        hdu = fits.TableHDU.from_columns(x)
        hdu.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert comparerecords(hdu.data, hdul[1].data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'

            if not sys.platform.startswith('win32'):
                # The checksum ends up being different on Windows, possibly due
                # to slight floating point differences
                assert 'CHECKSUM' in hdul[1].header
                assert hdul[1].header['CHECKSUM'] == '51IDA1G981GCA1G9'
                assert 'DATASUM' in hdul[1].header
                assert hdul[1].header['DATASUM'] == '1948208413'

    def test_compressed_image_data(self):
        with fits.open(self.data('comp.fits')) as h1:
            h1.writeto(self.temp('tmp.fits'), clobber=True, checksum=True)
            with fits.open(self.temp('tmp.fits'), checksum=True) as h2:
                assert np.all(h1[1].data == h2[1].data)
                assert 'CHECKSUM' in h2[0].header
                assert h2[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
                assert 'DATASUM' in h2[0].header
                assert h2[0].header['DATASUM'] == '0'
                assert 'CHECKSUM' in h2[1].header
                assert h2[1].header['CHECKSUM'] == 'ZeAbdb8aZbAabb7a'
                assert 'DATASUM' in h2[1].header
                assert h2[1].header['DATASUM'] == '113055149'

    def test_compressed_image_data_int16(self):
        n = np.arange(100, dtype='int16')
        hdu = fits.ImageHDU(n)
        comp_hdu = fits.CompImageHDU(hdu.data, hdu.header)
        comp_hdu.writeto(self.temp('tmp.fits'), checksum=True)
        hdu.writeto(self.temp('uncomp.fits'), checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert np.all(hdul[1].data == comp_hdu.data)
            assert np.all(hdul[1].data == hdu.data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'

            assert 'CHECKSUM' in hdul[1].header
            assert hdul[1]._header['CHECKSUM'] == 'J5cCJ5c9J5cAJ5c9'
            assert 'DATASUM' in hdul[1].header
            assert hdul[1]._header['DATASUM'] == '2453673070'
            assert 'CHECKSUM' in hdul[1].header

            with fits.open(self.temp('uncomp.fits'), checksum=True) as hdul2:
                header_comp = hdul[1]._header
                header_uncomp = hdul2[1].header
                assert 'ZHECKSUM' in header_comp
                assert 'CHECKSUM' in header_uncomp
                assert header_uncomp['CHECKSUM'] == 'ZE94eE91ZE91bE91'
                assert header_comp['ZHECKSUM'] == header_uncomp['CHECKSUM']
                assert 'ZDATASUM' in header_comp
                assert 'DATASUM' in header_uncomp
                assert header_uncomp['DATASUM'] == '160565700'
                assert header_comp['ZDATASUM'] == header_uncomp['DATASUM']

    def test_compressed_image_data_float32(self):
        n = np.arange(100, dtype='float32')
        hdu = fits.ImageHDU(n)
        comp_hdu = fits.CompImageHDU(hdu.data, hdu.header)
        comp_hdu.writeto(self.temp('tmp.fits'), checksum=True)
        hdu.writeto(self.temp('uncomp.fits'), checksum=True)
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            assert np.all(hdul[1].data == comp_hdu.data)
            assert np.all(hdul[1].data == hdu.data)
            assert 'CHECKSUM' in hdul[0].header
            assert hdul[0].header['CHECKSUM'] == 'D8iBD6ZAD6fAD6ZA'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'

            assert 'CHECKSUM' in hdul[1].header
            assert 'DATASUM' in hdul[1].header

            if not sys.platform.startswith('win32'):
                # The checksum ends up being different on Windows, possibly due
                # to slight floating point differences
                assert hdul[1]._header['CHECKSUM'] == 'eATIf3SHe9SHe9SH'
                assert hdul[1]._header['DATASUM'] == '1277667818'

            with fits.open(self.temp('uncomp.fits'), checksum=True) as hdul2:
                header_comp = hdul[1]._header
                header_uncomp = hdul2[1].header
                assert 'ZHECKSUM' in header_comp
                assert 'CHECKSUM' in header_uncomp
                assert header_uncomp['CHECKSUM'] == 'Cgr5FZo2Cdo2CZo2'
                assert header_comp['ZHECKSUM'] == header_uncomp['CHECKSUM']
                assert 'ZDATASUM' in header_comp
                assert 'DATASUM' in header_uncomp
                assert header_uncomp['DATASUM'] == '2393636889'
                assert header_comp['ZDATASUM'] == header_uncomp['DATASUM']

    def test_open_with_no_keywords(self):
        hdul = fits.open(self.data('arange.fits'), checksum=True)
        hdul.close()

    def test_append(self):
        hdul = fits.open(self.data('tb.fits'))
        hdul.writeto(self.temp('tmp.fits'), clobber=True)
        n = np.arange(100)
        fits.append(self.temp('tmp.fits'), n, checksum=True)
        hdul.close()
        hdul = fits.open(self.temp('tmp.fits'), checksum=True)
        assert hdul[0]._checksum is None
        hdul.close()

    def test_writeto_convenience(self):
        n = np.arange(100)
        fits.writeto(self.temp('tmp.fits'), n, clobber=True, checksum=True)
        hdul = fits.open(self.temp('tmp.fits'), checksum=True)
        self._check_checksums(hdul[0])
        hdul.close()

    def test_hdu_writeto(self):
        n = np.arange(100, dtype='int16')
        hdu = fits.ImageHDU(n)
        hdu.writeto(self.temp('tmp.fits'), checksum=True)
        hdul = fits.open(self.temp('tmp.fits'), checksum=True)
        self._check_checksums(hdul[0])
        hdul.close()

    def test_hdu_writeto_existing(self):
        """
        Tests that when using writeto with checksum=True, a checksum and
        datasum are added to HDUs that did not previously have one.

        Regression test for https://github.com/spacetelescope/PyFITS/issues/8
        """

        with fits.open(self.data('tb.fits')) as hdul:
            hdul.writeto(self.temp('test.fits'), checksum=True)

        with fits.open(self.temp('test.fits')) as hdul:
            assert 'CHECKSUM' in hdul[0].header
            # These checksums were verified against CFITSIO
            assert hdul[0].header['CHECKSUM'] == '7UgqATfo7TfoATfo'
            assert 'DATASUM' in hdul[0].header
            assert hdul[0].header['DATASUM'] == '0'
            assert 'CHECKSUM' in hdul[1].header
            assert hdul[1].header['CHECKSUM'] == '99daD8bX98baA8bU'
            assert 'DATASUM' in hdul[1].header
            assert hdul[1].header['DATASUM'] == '1829680925'

    def test_datasum_only(self):
        n = np.arange(100, dtype='int16')
        hdu = fits.ImageHDU(n)
        hdu.writeto(self.temp('tmp.fits'), clobber=True, checksum='datasum')
        with fits.open(self.temp('tmp.fits'), checksum=True) as hdul:
            if not (hasattr(hdul[0], '_datasum') and hdul[0]._datasum):
                pytest.fail(msg='Missing DATASUM keyword')

            if not (hasattr(hdul[0], '_checksum') and not hdul[0]._checksum):
                pytest.fail(msg='Non-empty CHECKSUM keyword')

            if not (hasattr(hdul[0], '_datasum_comment') and
                    hdul[0]._datasum_comment):
                pytest.fail(msg='Missing DATASUM Card comment')

            if not (hasattr(hdul[0], '_checksum_comment') and
                    not hdul[0]._checksum_comment):
                pytest.fail(msg='Non-empty CHECKSUM Card comment')

    def test_open_update_mode_preserve_checksum(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/148 where
        checksums are being removed from headers when a file is opened in
        update mode, even though no changes were made to the file.
        """

        self.copy_file('checksum.fits')

        with fits.open(self.temp('checksum.fits')) as hdul:
            data = hdul[1].data.copy()

        hdul = fits.open(self.temp('checksum.fits'), mode='update')
        hdul.close()

        with fits.open(self.temp('checksum.fits')) as hdul:
            assert 'CHECKSUM' in hdul[1].header
            assert 'DATASUM' in hdul[1].header
            assert comparerecords(data, hdul[1].data)

    def test_open_update_mode_update_checksum(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/148, part
        2.  This ensures that if a file contains a checksum, the checksum is
        updated when changes are saved to the file, even if the file was opened
        with the default of checksum=False.

        An existing checksum and/or datasum are only stripped if the file is
        opened with checksum='remove'.
        """

        self.copy_file('checksum.fits')
        with fits.open(self.temp('checksum.fits')) as hdul:
            header = hdul[1].header.copy()
            data = hdul[1].data.copy()

        with fits.open(self.temp('checksum.fits'), mode='update') as hdul:
            hdul[1].header['FOO'] = 'BAR'
            hdul[1].data[0]['TIME'] = 42

        with fits.open(self.temp('checksum.fits')) as hdul:
            header2 = hdul[1].header
            data2 = hdul[1].data
            assert header2[:-3] == header[:-2]
            assert 'CHECKSUM' in header2
            assert 'DATASUM' in header2
            assert header2['FOO'] == 'BAR'
            assert (data2['TIME'][1:] == data['TIME'][1:]).all()
            assert data2['TIME'][0] == 42

        with fits.open(self.temp('checksum.fits'), mode='update',
                       checksum='remove') as hdul:
            pass

        with fits.open(self.temp('checksum.fits')) as hdul:
            header2 = hdul[1].header
            data2 = hdul[1].data
            assert header2[:-1] == header[:-2]
            assert 'CHECKSUM' not in header2
            assert 'DATASUM' not in header2
            assert header2['FOO'] == 'BAR'
            assert (data2['TIME'][1:] == data['TIME'][1:]).all()
            assert data2['TIME'][0] == 42

    def _check_checksums(self, hdu):
        if not (hasattr(hdu, '_datasum') and hdu._datasum):
            pytest.fail(msg='Missing DATASUM keyword')

        if not (hasattr(hdu, '_checksum') and hdu._checksum):
            pytest.fail(msg='Missing CHECKSUM keyword')

        if not (hasattr(hdu, '_datasum_comment') and hdu._datasum_comment):
            pytest.fail(msg='Missing DATASUM Card comment')

        if not (hasattr(hdu, '_checksum_comment') and hdu._checksum_comment):
            pytest.fail(msg='Missing CHECKSUM Card comment')
