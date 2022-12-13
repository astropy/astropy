# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.coordinates.tests.helper import skycoord_equal as _skycoord_equal
from astropy.utils.decorators import deprecated

__all__ = ["skycoord_equal"]


@deprecated("5.1", alternative="astropy.coordinates.tests.helper.skycoord_equal")
def skycoord_equal(sc1, sc2):
    return _skycoord_equal(sc1, sc2)
