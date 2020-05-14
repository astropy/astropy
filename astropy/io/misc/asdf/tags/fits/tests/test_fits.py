# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

import pytest
import numpy as np

from astropy.io import fits
from astropy import __minimum_asdf_version__

asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)
from asdf.tests import helpers

from astropy.io.misc.asdf.tags.tests.helpers import run_schema_example_test


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
        dtype=[('A', int), ('B', int)])

    h = fits.HDUList()
    h.append(fits.BinTableHDU.from_columns(a))
    tree = {'fits': h}

    def check_yaml(content):
        assert b'!<tag:astropy.org:astropy/table/table-1.0.0>' in content

    helpers.assert_roundtrip_tree(tree, tmpdir, raw_yaml_check_func=check_yaml)


def test_backwards_compat():
    """
    Make sure that we can continue to read FITS HDUs that use the schema from
    the ASDF Standard.

    This test uses the examples in the fits schema from the ASDF Standard,
    since these make no reference to Astropy's own fits definition.
    """

    def check(asdffile):
        assert isinstance(asdffile['example'], fits.HDUList)

    run_schema_example_test('stsci.edu', 'asdf', 'fits/fits', '1.0.0', check)
