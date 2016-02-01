# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

import numpy as np

from ...utils.data import get_pkg_data_filenames, get_pkg_data_contents
from ...utils.misc import NumpyRNGContext

from ... import wcs


def test_maps():
    def test_map(filename):
        header = get_pkg_data_contents(os.path.join("maps", filename))
        wcsobj = wcs.WCS(header)

        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 12, wcsobj.wcs.naxis)
            world = wcsobj.wcs_pix2world(x, 1)
            pix = wcsobj.wcs_world2pix(x, 1)

    hdr_file_list = list(get_pkg_data_filenames("maps", pattern="*.hdr"))

    # actually perform a test for each one
    for filename in hdr_file_list:

        # use the base name of the file, because everything we yield
        # will show up in the test name in the pandokia report
        filename = os.path.basename(filename)

        # yield a function name and parameters to make a generated test
        yield test_map, filename

    # AFTER we tested with every file that we found, check to see that we
    # actually have the list we expect.  If N=0, we will not have performed
    # any tests at all.  If N < n_data_files, we are missing some files,
    # so we will have skipped some tests.  Without this check, both cases
    # happen silently!

    # how many do we expect to see?
    n_data_files = 28

    if len(hdr_file_list) != n_data_files:
        assert False, (
            "test_maps has wrong number data files: found %d, expected "
            " %d" % (
                len(hdr_file_list), n_data_files))
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


def test_spectra():
    def test_spectrum(filename):
        header = get_pkg_data_contents(os.path.join("spectra", filename))
        wcsobj = wcs.WCS(header)
        with NumpyRNGContext(123456789):
            x = np.random.rand(2 ** 16, wcsobj.wcs.naxis)
            world = wcsobj.wcs_pix2world(x, 1)
            pix = wcsobj.wcs_world2pix(x, 1)

    hdr_file_list = list(get_pkg_data_filenames("spectra", pattern="*.hdr"))

    # actually perform a test for each one
    for filename in hdr_file_list:

        # use the base name of the file, because everything we yield
        # will show up in the test name in the pandokia report
        filename = os.path.basename(filename)

        # yield a function name and parameters to make a generated test
        yield test_spectrum, filename

    # AFTER we tested with every file that we found, check to see that we
    # actually have the list we expect.  If N=0, we will not have performed
    # any tests at all.  If N < n_data_files, we are missing some files,
    # so we will have skipped some tests.  Without this check, both cases
    # happen silently!

    # how many do we expect to see?
    n_data_files = 6

    if len(hdr_file_list) != n_data_files:
        assert False, (
            "test_spectra has wrong number data files: found %d, expected "
            " %d" % (
                len(hdr_file_list), n_data_files))
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.
