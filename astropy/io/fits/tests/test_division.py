import numpy as np

from ...io import fits
from . import FitsTestCase
from .util import CaptureStdout


class TestDivisionFunctions(FitsTestCase):
    """Test code units that rely on correct integer division."""

    def test_rec_from_string(self):
        t1 = fits.open(self.data('tb.fits'))
        s = t1[1].data.tostring()
        a1 = np.rec.array(
                s,
                dtype=np.dtype([('c1', '>i4'), ('c2', '|S3'),
                                ('c3', '>f4'), ('c4', '|i1')]))

    def test_card_ncards(self):
        c1 = fits.Card('temp', 80.0, 'temperature')
        assert type(c1._ncards()) == type(1)

    def test_card_with_continue(self):
        h = fits.PrimaryHDU()
        with CaptureStdout() as f:
            h.header.update('abc', 'abcdefg'*20)
            assert f.getvalue() == ''

    def test_valid_hdu_size(self):
        t1 = fits.open(self.data('tb.fits'))
        assert type(t1[1].size) == type(1)

    def test_hdu_get_size(self):
        with CaptureStdout() as f:
            t1 = fits.open(self.data('tb.fits'))
            assert f.getvalue() == ''

    def test_section(self):
        # section testing
        fs = fits.open(self.data('arange.fits'))
        with CaptureStdout() as f:
            assert_equal(fs[0].section[3,2,5], np.array([357]))
            assert f.getvalue() == ''
