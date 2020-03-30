# Licensed under a 3-clause BSD style license - see PYFITS.rst

import numpy as np

from astropy.io import fits
from . import FitsTestCase
from astropy.tests.helper import catch_warnings


class TestDivisionFunctions(FitsTestCase):
    """Test code units that rely on correct integer division."""

    def test_rec_from_string(self):
        with fits.open(self.data('tb.fits')) as t1:
            s = t1[1].data.tobytes()
        np.rec.array(
            s,
            dtype=np.dtype([('c1', '>i4'), ('c2', '|S3'),
                            ('c3', '>f4'), ('c4', '|i1')]),
            shape=len(s) // 12)

    def test_card_with_continue(self):
        h = fits.PrimaryHDU()
        with catch_warnings() as w:
            h.header['abc'] = 'abcdefg' * 20
        assert len(w) == 0

    def test_valid_hdu_size(self):
        with fits.open(self.data('tb.fits')) as t1:
            assert type(t1[1].size) is type(1)  # noqa

    def test_hdu_get_size(self):
        with catch_warnings() as w:
            with fits.open(self.data('tb.fits')) as _:
                pass
        assert len(w) == 0

    def test_section(self, capsys):
        # section testing
        fs = fits.open(self.data('arange.fits'))
        with catch_warnings() as w:
            assert np.all(fs[0].section[3, 2, 5] == np.array([357]))
            assert len(w) == 0
