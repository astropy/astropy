# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division, with_statement

import warnings

from ....extern import six  # pylint: disable=W0611
from ....io import fits
from ....table import Table
from .. import printdiff
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

    def test_table_to_hdu(self, tmpdir):
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
        filename = str(tmpdir.join('test_table_to_hdu.fits'))
        hdu.writeto(filename, overwrite=True)

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
