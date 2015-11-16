# Licensed under a 3-clause BSD style license - see LICENSE.rst
import functools
import itertools

import numpy as np

from ...tests.helper import (pytest, quantity_allclose as allclose,
                             assert_quantity_allclose as assert_allclose)
from .. import Time

def test_corrections():
    from .. import coordinates as coord, units as u
    wht  = coord.EarthLocation.of_site('lapalma')
    star = coord.SkyCoord("08:08:08 +32:00:00",distance=120*u.pc,unit=(u.hour,u.degree),frame='icrs')
    t = Time("2013-02-02T23:00")
    t.location = wht
    hval = t.heliocentric_correction(star)
    bval = t.barycentric_correction(star)
    assert allclose(hval.sec,461.43037870502235)
    assert allclose(bval.sec,460.58538779827836)