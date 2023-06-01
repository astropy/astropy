# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Testing :mod:`astropy.cosmology.flrw.__init__.py`."""

##############################################################################
# IMPORTS

import pytest

##############################################################################
# TESTS
##############################################################################


def test_getattr_error_attr_not_found():
    """Test getattr raises error for DNE."""
    with pytest.raises(ImportError):
        from astropy.cosmology.flrw import this_is_not_a_variable  # noqa: F401
