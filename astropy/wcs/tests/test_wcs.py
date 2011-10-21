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
    # belong to the wcs that we are testing, even if we are not running
    # this test from the installed copy of this file.  Use wcs.__file__
    ROOT_DIR = os.path.dirname(wcs.__file__) + "/tests"


# test_maps() is a generator
def test_maps():

    # test_map() is the function that is called to perform the generated test
    def test_map(filename):

        # the test parameter is the base name of the file to use; find
        # the file in the installed wcs test directory
        filename = os.path.join(ROOT_DIR, "maps", filename)

        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()

        wcsobj = wcs.WCS(header)

        world = wcsobj.wcs_pix2sky([[97, 97]], 1)

        assert_array_almost_equal(world, [[285.0, -66.25]], decimal=1)

        pix = wcsobj.wcs_sky2pix([[285.0, -66.25]], 1)

        assert_array_almost_equal(pix, [[97, 97]], decimal=0)

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
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


# test_spectra() is a generator
def test_spectra():

    # test_spectrum() is the function that is called to perform the
    # generated test
    def test_spectrum(filename):

        # the test parameter is the base name of the file to use; find
        # the file in the installed wcs test directory
        filename = os.path.join(ROOT_DIR, "spectra", filename)

        fd = open(filename, 'rb')
        header = fd.read()
        fd.close()

        wcsobj = wcs.WCS(header)

        all = wcs.find_all_wcs(header)
        assert len(all) == 9

    # get the list of the hdr files that we want to test
    hdr_file_list = [x for x in glob.glob(os.path.join(
        ROOT_DIR, "spectra", "*.hdr"))]

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
        # b.t.w.  If this assert happens, py.test reports one more test
        # than it would have otherwise.


def test_units():
    u = wcs.UnitConverter("log(MHz)", "ln(Hz)")
    print(u.convert([1, 2, 3, 4]))

basic_units = "m s g rad sr K A mol cd".split()
derived_units = "Hz J W V N Pa C Ohm ohm S F Wb T H lm lx".split()
add_all_units = "eV Jy R G barn".split()
add_sup_units = "a yr pc bit byte Byte".split()
add_sub_units = "mag".split()
general_units = \
              "deg arcmin arcsec mas d h min erg Ry u D DEGREE DEGREES".split()
astro_units = "Angstrom angstrom AU lyr beam solRad solMass solLum Sun".split()
device_units = "adu bin chan count ct photon ph pixel pix voxel".split()

sub_prefixes = "y z a f p n u m c d".split()
sup_prefixes = "da h k M G T P E Z Y".split()


def test_all_units():
    def test_self(x):
        # x appears in the test name. If we would have had an ambiguous
        # test name, we had -xxx added to the unit name.  Remove it if
        # necessary.
        if '-' in x:
            x = x.split('-')[0]

        # here is the test:
        try:
            u = wcs.UnitConverter(x, x)
        except ValueError:
            e = sys.exc_info()[1]
            if str(e).startswith("ERROR 12 in wcsutrne") and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    # list of all the units to test
    all = sorted(basic_units + derived_units + add_all_units + add_sup_units
            + add_sub_units + general_units + astro_units + device_units)

    # Pandokia has non-case-sensitve test names; since the unit name is
    # showing up in the test name, we want to disambiguate any name collisions.
    # Here is a list of all the lower-cased unit name names.
    all_lower = [x.lower() for x in all]

    # here are serial numbers to use to disambiguate
    unique_tags = {}

    for unit in all:
        # disambiguate the test name, if necessary
        l_unit = unit.lower()
        if unit != l_unit and l_unit in all_lower:
            n = unique_tags.get(l_unit, 1)
            unique_tags[n] = n + 1

            # the test will tear off the part after the '-'
            unit = '%s-%d' % (unit, n)

        # perform the test
        yield test_self, unit


def test_unit_prefixes():
    def test_self(x, p):
        unit = p + x
        try:
            u = wcs.UnitConverter(unit, unit)
        except ValueError:
            e = sys.exc_info()[1]
            if str(e) == "Potentially unsafe translation" and \
                    x in ("S", "H", "D"):
                return
            else:
                raise
        assert u.scale == 1.0
        assert u.offset == 0.0
        assert u.power == 1.0

    for unit in (basic_units + derived_units + add_all_units):
        for prefix in (sub_prefixes + sup_prefixes):
            yield test_self, unit, prefix

    for unit in add_sup_units:
        for prefix in sup_prefixes:
            yield test_self, unit, prefix

    for unit in add_sub_units:
        for prefix in sub_prefixes:
            yield test_self, unit, prefix
