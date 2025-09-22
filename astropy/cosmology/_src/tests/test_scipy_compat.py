# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.cosmology._src.flrw.base import quad
from astropy.utils.compat.optional_deps import HAS_SCIPY


@pytest.mark.skipif(HAS_SCIPY, reason="scipy is installed")
def test_optional_deps_functions():
    """Test stand-in functions when optional dependencies not installed."""
    with pytest.raises(ModuleNotFoundError, match="No module named 'scipy.integrate'"):
        quad()
