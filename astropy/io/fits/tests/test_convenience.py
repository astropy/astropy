# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division, with_statement

import warnings

import numpy as np

from ....extern import six  # pylint: disable=W0611
from ....io import fits
from ....table import Table
from ....tests.helper import pytest, catch_warnings

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

    def test_table_to_hdu(self):
        table = Table([[1, 2, 3], ['a', 'b', 'c'], [2.3, 4.5, 6.7]],
                      names=['a', 'b', 'c'], dtype=['i', 'U1', 'f'])
        table['a'].unit = 'm/s'
        table['b'].unit = 'not-a-unit'

        hdu = fits.table_to_hdu(table)

        # Check that TUNITn cards appear in the correct order
        # (https://github.com/astropy/astropy/pull/5720)
        assert hdu.header.index('TUNIT1') < hdu.header.index('TTYPE2')

        assert isinstance(hdu, fits.BinTableHDU)
        filename = self.temp('test_table_to_hdu.fits')
        hdu.writeto(filename, overwrite=True)

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
