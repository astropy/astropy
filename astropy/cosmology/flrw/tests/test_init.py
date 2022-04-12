# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.__init__.py`."""

##############################################################################
# IMPORTS

import pytest

from astropy.utils import resolve_name
from astropy.utils.exceptions import AstropyDeprecationWarning

##############################################################################
# TESTS
##############################################################################


@pytest.mark.parametrize(
    "attr",
    [
        "H0units_to_invs",
        "a_B_c2",
        "critdens_const",
        "kB_evK",
        "quad",
        "radian_in_arcmin",
        "radian_in_arcsec",
        "sec_to_Gyr",
        "ellipkinc",
        "hyp2f1",
    ],
)
def test_deprecated_private_variables(attr):
    """Test deprecation warnings are raised for private variables."""
    with pytest.warns(AstropyDeprecationWarning):
        resolve_name("astropy", "cosmology", "flrw", attr)


def test_getattr_error_attr_not_found():
    """Test getattr raises error for DNE."""
    with pytest.raises(ImportError):
        from astropy.cosmology.flrw import this_is_not_a_variable  # noqa: F401
