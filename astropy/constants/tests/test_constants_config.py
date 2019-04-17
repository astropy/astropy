# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import os
import pytest
import tempfile

from astropy.tests.helper import catch_warnings
from astropy.config import paths


def write_test_config(dir, physical_constants=None,
                      astronomical_constants=None):
    """Write configuration items to a test directory

    Parameters:
    -----------
    dir: `str`
        path to directory to write configuration (in astropy subdirectory)
    physical_constants: None or `str`
        value of config item for physical constants. If None, comment out.
    astronomical_constants: `str`
        value of config item for astronomical constants. If None, comment out.

    """
    cfgdir = os.path.join(dir, 'astropy')
    os.makedirs(os.path.join(cfgdir), exist_ok=True)
    cfgpath = os.path.join(cfgdir, 'astropy.cfg')
    with open(cfgpath, 'w') as cfgfile:
        cfgfile.write("# -*- coding: utf-8 -*-\n\n")
        cfgfile.write("[constants]\n\n")
        if isinstance(physical_constants, str):
            cfgfile.write("physical_constants = '{}'\n"
                          .format(physical_constants))
        else:
            cfgfile.write("# physical_constants = 'codata2014'\n")
        if isinstance(astronomical_constants, str):
            cfgfile.write("astronomical_constants = '{}'\n"
                          .format(astronomical_constants))
        else:
            cfgfile.write("# astronomical_constants = 'iau2015'\n")


def test_prior_config():
    with tempfile.TemporaryDirectory() as tmpdirname:
        write_test_config(tmpdirname, physical_constants='codata2010',
                          astronomical_constants='iau2012')
        with paths.set_temp_config(tmpdirname):
            from astropy.constants import conf

            with catch_warnings():
                conf.reload()

                import astropy.constants as const
                importlib.reload(const)

            assert conf.physical_constants == 'codata2010'
            h = const.h
            # check that the value is the CODATA2010 value
            assert abs(h.value - 6.62606957e-34) < 1e-43
            assert abs(h.si.value - 6.62606957e-34) < 1e-43
            assert abs(h.cgs.value - 6.62606957e-27) < 1e-36

            from astropy.constants.codata2014 import h as h_2014
            # Check it is different from the current value
            assert abs(h.value - h_2014.value) > 4e-42

            assert conf.astronomical_constants == 'iau2012'
            R_earth = const.R_earth
            assert R_earth.value == 6.378136e6

            from astropy.constants.iau2015 import R_earth as R_earth_2015
            # Check it is different from the current value
            assert abs(R_earth.value - R_earth_2015.value) > 10.0

        # Test setting by version
        write_test_config(tmpdirname, physical_constants='astropyconst13',
                          astronomical_constants='astropyconst13')
        with paths.set_temp_config(tmpdirname):
            from astropy.constants import conf

            with catch_warnings():
                conf.reload()

                import astropy.constants as const
                importlib.reload(const)

            assert conf.physical_constants == 'astropyconst13'
            h = const.h
            # check that the value is the CODATA2010 value
            assert abs(h.value - 6.62606957e-34) < 1e-43
            assert abs(h.si.value - 6.62606957e-34) < 1e-43
            assert abs(h.cgs.value - 6.62606957e-27) < 1e-36

            from astropy.constants.codata2014 import h as h_2014
            # Check it is different from the current value
            assert abs(h.value - h_2014.value) > 4e-42

            assert conf.astronomical_constants == 'astropyconst13'
            R_earth = const.R_earth
            assert R_earth.value == 6.378136e6

            from astropy.constants.iau2015 import R_earth as R_earth_2015
            # Check it is different from the current value
            assert abs(R_earth.value - R_earth_2015.value) > 10.0

    # reset state of constants (in part to prevent failures of later tests)
    with catch_warnings():
        conf.reload()
        importlib.reload(const)
    assert conf.physical_constants == 'codata2014'
    h = const.h
    assert abs(h.value - h_2014.value) < 4e-42

    assert conf.astronomical_constants == 'iau2015'
    R_earth = const.R_earth
    assert abs(R_earth.value - R_earth_2015.value) < 0.01

    # check current constants are consistent with units
    import astropy.units as u
    assert (abs((const.M_earth / u.M_earth).to(u.dimensionless_unscaled)
            - 1) < 1.e-6)


def assert_config_outputs(physical_in, physical_out,
                          astronomical_in, astronomical_out):
    """Write inputs to temporary config and assert the outputs

    Parameters:
    -----------
    physical_in: `str` or None
        input physical constants version
    physical_out: `str`
        output physical constants version
    astronomical_in: `str` or None
        input astronomical constants version
    astronomical_out: `str`
        output astronomical constants version
    """
    import astropy.constants as const
    from astropy.constants import conf
    with tempfile.TemporaryDirectory() as tmpdirname:
        write_test_config(tmpdirname, physical_constants=physical_in,
                          astronomical_constants=astronomical_in)
        with paths.set_temp_config(tmpdirname):
            with catch_warnings():
                conf.reload()
            importlib.reload(const)
            assert conf.physical_constants == physical_out
            assert conf.astronomical_constants == astronomical_out


def test_invalid_config():
    """Test invalid config items"""
    with pytest.raises(ValueError):
        assert_config_outputs('cooldata2014', 'codata2014',
                              'iau2012', 'iau2012')

    with pytest.raises(ValueError):
        assert_config_outputs('codata2010', 'codata2010',
                              'allen1976', 'iau2015')


def test_valid_config():
    """Test valid config items"""

    assert_config_outputs('codata2014', 'codata2014',
                          'iau2015', 'iau2015')

    assert_config_outputs('codata2010', 'codata2010',
                          'iau2015', 'iau2015')

    assert_config_outputs('codata2014', 'codata2014',
                          'iau2012', 'iau2012')

    assert_config_outputs('codata2010', 'codata2010',
                          'iau2012', 'iau2012')

    assert_config_outputs('astropyconst13', 'astropyconst13',
                          'astropyconst20', 'astropyconst20')

    assert_config_outputs('astropyconst20', 'astropyconst20',
                          'astropyconst20', 'astropyconst20')

    assert_config_outputs('astropyconst20', 'astropyconst20',
                          'astropyconst13', 'astropyconst13')

    # Check that commenting out values gives the defaults
    assert_config_outputs(None, 'codata2014',
                          None, 'iau2015')
