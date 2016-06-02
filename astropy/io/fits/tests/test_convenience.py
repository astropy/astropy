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

    def test_table_to_hdu(self, tmpdir):
        table = Table([[1, 2, 3], ['a', 'b', 'c'], [2.3, 4.5, 6.7]],
                      names=['a', 'b', 'c'], dtype=['i', 'U1', 'f'])
        hdu = fits.table_to_hdu(table)
        assert isinstance(hdu, fits.BinTableHDU)
        filename = str(tmpdir.join('test_table_to_hdu.fits'))
        hdu.writeto(filename, clobber=True)
