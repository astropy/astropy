# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division, with_statement

import os
import shutil
import warnings

import pytest
import numpy as np

from ....extern import six  # pylint: disable=W0611
from ....io import fits
from ....table import Table
from .. import printdiff
from ....tests.helper import catch_warnings

from . import FitsTestCase


class TestConvenience(FitsTestCase):

    @pytest.mark.skipif('six.PY2')
    def test_resource_warning(self):
        warnings.simplefilter('always', ResourceWarning)
        with catch_warnings() as w:
            data = fits.getdata(self.data('test0.fits'))
            assert len(w) == 0

        with catch_warnings() as w:
            header = fits.getheader(self.data('test0.fits'))
            assert len(w) == 0

    def test_fileobj_not_closed(self):
        """
        Tests that file-like objects are not closed after being passed
        to convenience functions.

        Regression test for https://github.com/astropy/astropy/issues/5063
        """

        f = open(self.data('test0.fits'), 'rb')
        data = fits.getdata(f)
        assert not f.closed

        f.seek(0)
        header = fits.getheader(f)
        assert not f.closed

        f.close()  # Close it now

    def test_table_to_hdu(self):
        table = Table([[1, 2, 3], ['a', 'b', 'c'], [2.3, 4.5, 6.7]],
                      names=['a', 'b', 'c'], dtype=['i', 'U1', 'f'])
        table['a'].unit = 'm/s'
        table['b'].unit = 'not-a-unit'

        with catch_warnings() as w:
            hdu = fits.table_to_hdu(table)
            assert len(w) == 1
            assert str(w[0].message).startswith("'not-a-unit' did not parse as"
                                                " fits unit")

        # Check that TUNITn cards appear in the correct order
        # (https://github.com/astropy/astropy/pull/5720)
        assert hdu.header.index('TUNIT1') < hdu.header.index('TTYPE2')

        assert isinstance(hdu, fits.BinTableHDU)
        filename = self.temp('test_table_to_hdu.fits')
        hdu.writeto(filename, overwrite=True)

    def test_table_to_hdu_convert_comment_convention(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/6079
        """
        table = Table([[1, 2, 3], ['a', 'b', 'c'], [2.3, 4.5, 6.7]],
                      names=['a', 'b', 'c'], dtype=['i', 'U1', 'f'])
        table.meta['comments'] = ['This', 'is', 'a', 'comment']
        hdu = fits.table_to_hdu(table)

        assert hdu.header.get('comment') == ['This', 'is', 'a', 'comment']
        with pytest.raises(ValueError):
            hdu.header.index('comments')

    def test_table_writeto_header(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5988
        """
        data = np.zeros((5, ), dtype=[('x', np.float), ('y', np.int)])
        h_in = fits.Header()
        h_in['ANSWER'] = (42.0, 'LTU&E')
        filename = self.temp('tabhdr42.fits')
        fits.writeto(filename, data=data, header=h_in, overwrite=True)
        h_out = fits.getheader(filename, ext=1)
        assert h_out['ANSWER'] == 42

    def test_image_extension_update_header(self):
        """
        Test that _makehdu correctly includes the header. For example in the
        fits.update convenience function.
        """
        filename = self.temp('twoextension.fits')

        hdus = [fits.PrimaryHDU(np.zeros((10, 10))),
                fits.ImageHDU(np.zeros((10, 10)))]

        fits.HDUList(hdus).writeto(filename)

        fits.update(filename,
                    np.zeros((10, 10)),
                    header=fits.Header([('WHAT', 100)]),
                    ext=1)
        h_out = fits.getheader(filename, ext=1)
        assert h_out['WHAT'] == 100

    def test_printdiff(self):
        """
        Test that FITSDiff can run the different inputs without crashing.
        """

        # Testing different string input options
        assert printdiff(self.data('arange.fits'),
                         self.data('blank.fits')) is None
        assert printdiff(self.data('arange.fits'),
                         self.data('blank.fits'), ext=0) is None
        assert printdiff(self.data('o4sp040b0_raw.fits'),
                         self.data('o4sp040b0_raw.fits'),
                         extname='sci') is None

        # This may seem weird, but check printdiff to see, need to test
        # incorrect second file
        with pytest.raises(IOError):
            printdiff('o4sp040b0_raw.fits', 'fakefile.fits', extname='sci')

        # Test HDU object inputs
        with fits.open(self.data('stddata.fits'), mode='readonly') as in1:
            with fits.open(self.data('checksum.fits'), mode='readonly') as in2:

                assert printdiff(in1[0], in2[0]) is None

                with pytest.raises(ValueError):
                    printdiff(in1[0], in2[0], ext=0)

                assert printdiff(in1, in2) is None

                with pytest.raises(NotImplementedError):
                    printdiff(in1, in2, 0)

    def test_tabledump(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/6937
        """
        # copy fits file to the temp directory
        filename = self.data('tb.fits')
        temp_filename = self.temp('tb.fits')
        shutil.copyfile(filename, temp_filename)

        # test without datafile
        fits.tabledump(temp_filename)
        assert os.path.isfile(self.temp('tb_1.txt'))

        # test with datafile
        fits.tabledump(temp_filename, datafile=self.temp('test_tb.txt'))
        assert os.path.isfile(self.temp('test_tb.txt'))

    @pytest.mark.parametrize('mode', ['wb', 'wb+', 'ab', 'ab+'])
    def test_append_filehandle(self, tmpdir, mode):

        append_file = tmpdir.join('append.fits')
        with append_file.open(mode) as handle:
            fits.append(filename=handle, data=np.ones((4, 4)))
