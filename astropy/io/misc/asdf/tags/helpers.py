# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import numpy as np


__all__ = []


def skycoord_equal(sc1, sc2):
    """SkyCoord equality useful for testing and ASDF serialization
    """
    if not sc1.is_equivalent_frame(sc2):
        return False
    if sc1.representation_type is not sc2.representation_type:
        return False
    if sc1.shape != sc2.shape:
        return False  # Maybe raise ValueError corresponding to future numpy behavior
    eq = np.ones(shape=sc1.shape, dtype=bool)
    for comp in sc1.data.components:
        eq &= getattr(sc1.data, comp) == getattr(sc2.data, comp)
    return np.all(eq)
