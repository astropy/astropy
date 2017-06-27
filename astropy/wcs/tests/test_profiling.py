# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import pytest
import numpy as np

from ...utils.data import get_pkg_data_filenames, get_pkg_data_contents
from ...utils.misc import NumpyRNGContext

from ... import wcs

# hdr_map_file_list = list(get_pkg_data_filenames("maps", pattern="*.hdr"))

# use the base name of the file, because everything we yield
# will show up in the test name in the pandokia report
hdr_map_file_list = [os.path.basename(fname) for fname in get_pkg_data_filenames("maps", pattern="*.hdr")]

# Checking the number of files before reading them in.
# OLD COMMENTS:
# AFTER we tested with every file that we found, check to see that we
# actually have the list we expect.  If N=0, we will not have performed
# any tests at all.  If N < n_data_files, we are missing some files,
# so we will have skipped some tests.  Without this check, both cases
# happen silently!


def test_read_map_files():
    # how many map files we expect to see
    n_map_files = 28

    assert len(hdr_map_file_list) == n_map_files, (
           "test_read_map_files has wrong number data files: found {}, expected "
           " {}".format(len(hdr_map_file_list), n_map_files))


@pytest.mark.parametrize("filename", hdr_map_file_list)
def test_map(filename):
        header = get_pkg_data_contents(os.path.join("maps", filename))
        wcsobj = wcs.WCS(header)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 12, wcsobj.wcs.naxis)
            world = wcsobj.wcs_pix2world(x, 1)
            pix = wcsobj.wcs_world2pix(x, 1)


hdr_spec_file_list = [os.path.basename(fname) for fname in get_pkg_data_filenames("spectra", pattern="*.hdr")]


def test_read_spec_files():
    # how many spec files expected
    n_spec_files = 6

    assert len(hdr_spec_file_list) == n_spec_files, (
            "test_spectra has wrong number data files: found {}, expected "
            " {}".format(len(hdr_spec_file_list), n_spec_files))
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


@pytest.mark.parametrize("filename", hdr_spec_file_list)
def test_spectrum(filename):
    header = get_pkg_data_contents(os.path.join("spectra", filename))
    wcsobj = wcs.WCS(header)
    with NumpyRNGContext(123456789):
        x = np.random.rand(2 ** 16, wcsobj.wcs.naxis)
        world = wcsobj.wcs_pix2world(x, 1)
        pix = wcsobj.wcs_world2pix(x, 1)
