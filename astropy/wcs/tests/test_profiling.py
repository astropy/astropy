import glob
import os
import sys

import numpy as np
from numpy.testing import assert_array_almost_equal

from astropy import wcs

ROOT_DIR = None


def setup_module():
    global ROOT_DIR

    # do not use __file__ here - we want to find the data files that
    # belong to the astropy.wcs that we are testing, even if we are
    # not running this test from the installed copy of this file.  Use
    # astropy.wcs.__file__

    ROOT_DIR = os.path.join(os.path.dirname(wcs.__file__), 'tests')


def test_maps():
    def test_map(filename):
        filename = os.path.join(ROOT_DIR, "maps", filename)

        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()
        wcsobj = wcs.WCS(header)

        x = np.random.rand(2 ** 16, wcsobj.wcs.naxis)
        world = wcsobj.wcs_pix2sky(x, 1)
        pix = wcsobj.wcs_sky2pix(x, 1)

    # get the list of the hdr files that we want to test
    hdr_file_list = [x for x in glob.glob(
        os.path.join(ROOT_DIR, "maps", "*.hdr"))]

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
            " %d, looking in %s" % (
                len(hdr_file_list), n_data_files, ROOT_DIR))
        # b.t.w.  If this assert happens, nose reports one more test
        # than it would have otherwise.


def test_spectra():
    def test_spectrum(filename):
        filename = os.path.join(ROOT_DIR, "spectra", filename)

        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()
        wcsobj = wcs.WCS(header)

        x = np.random.rand(2 ** 16, wcsobj.wcs.naxis)
        world = wcsobj.wcs_pix2sky(x, 1)
        pix = wcsobj.wcs_sky2pix(x, 1)

    # get the list of the hdr files that we want to test
    hdr_file_list = [x for x in glob.glob(
        os.path.join(ROOT_DIR, "spectra", "*.hdr"))]

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
            " %d, looking in %s" % (
                len(hdr_file_list), n_data_files, ROOT_DIR))
        # b.t.w.  If this assert happens, nose reports one more test
        # than it would have otherwise.
