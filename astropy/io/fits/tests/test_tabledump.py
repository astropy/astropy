# Licensed under a 3-clause BSD style license - see PYFITS.rst


import os
import warnings
import unittest
import pytest
import numpy as np

from ....io import fits
from ....table import Table
from .. import printdiff
from ....tests.helper import catch_warnings

from . import FitsTestCase

class TestConvenience(FitsTestCase):

    def test_tabledump(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/6937
        """
        # test without datafile
        filename = self.data('tb.fits')
        fits.tabledump(filename)
        assert os.path.isfile(self.data('tb_1.txt'))
        os.remove(self.data('tb_1.txt'))

        # test with datafile
        fits.tabledump(filename, datafile=self.temp('test_tb.txt'))
        assert os.path.isfile(self.temp('test_tb.txt'))

        # test with cdfile, hfile
        """
        test for https://github.com//astropy/astropy/issues/7015
        """
        fits.tabledump(filename, datafile=self.temp('test_tb.txt'), \
        hfile=self.temp('hfile.txt'), cdfile=self.temp('cdfile.txt'),\
        overwrite=True )
        assert os.path.isfile(self.temp('hfile.txt'))
        assert os.path.isfile(self.temp('cdfile.txt'))
        hdump = fits.tableload(self.temp('test_tb.txt'),\
        hfile=self.temp('hfile.txt'), cdfile=self.temp('cdfile.txt'))
        assert os.path.isfile(self.temp('test_tb.txt'))
        hdump.writeto(self.temp('tableb_tb.fits'))
        assert printdiff(self.temp('tableb_tb.fits'), filename) is None
