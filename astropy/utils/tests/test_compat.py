# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Just deprecation tests, the rest is tested throughout astropy.
"""

import pytest

from astropy.utils.exceptions import AstropyPendingDeprecationWarning


def test_copy_if_needed_deprecation():
    with pytest.warns(AstropyPendingDeprecationWarning, match="COPY_IF_NEEDED"):
        from astropy.utils.compat.numpycompat import COPY_IF_NEEDED
    assert COPY_IF_NEEDED is None

    with pytest.warns(AstropyPendingDeprecationWarning, match="COPY_IF_NEEDED"):
        from astropy.utils.compat import COPY_IF_NEEDED
    assert COPY_IF_NEEDED is None
