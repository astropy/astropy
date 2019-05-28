# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import pytest
import sys


from astropy import physical_constants, astronomical_constants


def test_invalid_config():
    """Test invalid config items"""
    if 'astropy.units' in sys.modules:
        del sys.modules['astropy.units']
    if 'astropy.constants' in sys.modules:
        del sys.modules['astropy.constants']
    with pytest.raises(ValueError):
        physical_constants.set('cooldata2014')

    with pytest.raises(ValueError):
        astronomical_constants.set('iau2006')


"""
def test_prior_config():
    with physical_constants.set('codata2010'):
        assert physical_constants.get() == 'codata2010'
        with astronomical_constants.set('iau2012'):

            import astropy.constants as const
            importlib.reload(const)

            h = const.h
            # check that the value is the CODATA2010 value
            assert abs(h.value - 6.62606957e-34) < 1e-43
            assert abs(h.si.value - 6.62606957e-34) < 1e-43
            assert abs(h.cgs.value - 6.62606957e-27) < 1e-36

            from astropy.constants.codata2014 import h as h_2014
            # Check it is different from the current value
            assert abs(h.value - h_2014.value) > 4e-42

            assert astronomical_constants.get() == 'iau2012'
            R_earth = const.R_earth
            assert R_earth.value == 6.378136e6

            from astropy.constants.iau2015 import R_earth as R_earth_2015
            # Check it is different from the current value
            assert abs(R_earth.value - R_earth_2015.value) > 10.0

            # check current constants are consistent with units
            import astropy.units as u
            assert (abs((const.M_earth / u.M_earth)
                        .to(u.dimensionless_unscaled) - 1) < 1.e-6)
"""


def test_previously_imported():
    import astropy.units

    with pytest.raises(RuntimeError):
        physical_constants.set('codata2018')
