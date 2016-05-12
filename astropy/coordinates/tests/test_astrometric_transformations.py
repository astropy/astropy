# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ... import units as u
from ..distances import Distance
from ..builtin_frames import (ICRS, AstrometricICRS)
from .. import SkyCoord
from ...tests.helper import (pytest, quantity_allclose as allclose,
                             assert_quantity_allclose as assert_allclose)


def test_astrometric_unit():
    """Make sure it works with skycoord too."""
    
    origin = ICRS(ra=45*u.deg, dec=90*u.deg)
    astrometric_frame = AstrometricICRS(origin=origin)
    skycoord = SkyCoord([0, 45, 90], [0, 45, 90], frame=ICRS, unit=u.deg) 
    actual = skycoord.transform_to(astrometric_frame)
    
    actual_xyz = actual.cartesian.xyz
    expected = SkyCoord([-45, 0, 45], [-45, 0, 45], frame=astrometric_frame, unit=u.deg) 
    expected_xyz = expected.cartesian.xyz

    assert_allclose(actual_xyz, expected_xyz)