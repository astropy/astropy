import pytest

import astropy.io.fits as fits
import astropy.units as u
from astropy.utils.data import get_pkg_data_filename


@pytest.mark.remote_data
def test_unit_aliases_helper():
    with fits.open(get_pkg_data_filename('allsky/allsky_rosat.fits')) as hdul:
        with u.unit_aliases({"counts": u.count}):
            assert u.Unit(hdul[0].header["BUNIT"]) == u.Unit("10**(-6) count/s")


@pytest.mark.remote_data
def test_unit_aliases_registry():
    registry = u.get_current_unit_registry()
    with fits.open(get_pkg_data_filename('allsky/allsky_rosat.fits')) as hdul:
        with registry.unit_aliases({"counts": u.count}):
            assert u.Unit(hdul[0].header["BUNIT"]) == u.Unit("10**(-6) count/s")
