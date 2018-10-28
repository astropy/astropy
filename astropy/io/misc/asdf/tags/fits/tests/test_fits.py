# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

import pytest
import numpy as np

from astropy.io import fits

asdf = pytest.importorskip('asdf', minversion='2.0.0.dev0')
from asdf.tests import helpers


def test_complex_structure(tmpdir):
    with fits.open(os.path.join(
            os.path.dirname(__file__), 'data', 'complex.fits'), memmap=False) as hdulist:
        tree = {
            'fits': hdulist
            }

        helpers.assert_roundtrip_tree(tree, tmpdir)


def test_fits_table(tmpdir):
    a = np.array(
        [(0, 1), (2, 3)],
        dtype=[(str('A'), int), (str('B'), int)])

    h = fits.HDUList()
    h.append(fits.BinTableHDU.from_columns(a))
    tree = {'fits': h}

    def check_yaml(content):
        assert b'!core/table' in content

    helpers.assert_roundtrip_tree(tree, tmpdir, raw_yaml_check_func=check_yaml)
