# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy import testing as npt
from ...tests.helper import pytest

from ... import units as u


"""
These are the tests for coordinate matching.  

Note that this requires scipy.
"""

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching_function():
    from .. import ICRS
    from ..matching import match_coordinates_3d
    #this only uses match_coordinates_3d because that's the actual implementation

    cmatch = ICRS([4, 2.1], [0, 0], unit=(u.degree, u.degree))
    ccatalog = ICRS([1, 2, 3, 4], [0, 0, 0, 0], unit=(u.degree, u.degree))

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog)
    npt.assert_array_equal(idx, [3, 1])
    npt.assert_array_almost_equal(d2d.degree, [0, 0.1])
    assert d3d.value[0] == 0

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog, nthneighbor=2)
    assert np.all(idx == 2)
    npt.assert_array_almost_equal(d2d.degree, [1, 0.9])
    npt.assert_array_less(d3d.value, 0.02)


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching_function_3d_and_sky():
    from .. import ICRS
    from ..matching import match_coordinates_3d, match_coordinates_sky

    cmatch = ICRS([4, 2.1], [0, 0], unit=(u.degree, u.degree), distance=[1, 5] * u.kpc)
    ccatalog = ICRS([1, 2, 3, 4], [0, 0, 0, 0], unit=(u.degree, u.degree), distance=[1, 1, 1, 5] * u.kpc)

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog)
    npt.assert_array_equal(idx, [2, 3])

    npt.assert_array_almost_equal(d2d.degree, [1, 1.9])
    npt.assert_array_less(d3d.value, [0.02, .2])
    assert d3d.value[1] > 0.02

    idx, d2d, d3d = match_coordinates_sky(cmatch, ccatalog)
    npt.assert_array_equal(idx, [2, 1])

    npt.assert_array_almost_equal(d2d.degree, [1, 0.1])
    npt.assert_array_less(d3d.value, [0.02, .2])
    assert d3d.value[1] < 0.02


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_kdtree_storage():
    from .. import ICRS
    from ..matching import match_coordinates_3d

    cmatch = ICRS([4, 2.1], [0, 0], unit=(u.degree, u.degree))
    ccatalog = ICRS([1, 2, 3, 4], [0, 0, 0, 0], unit=(u.degree, u.degree))

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog, storekdtree=False)
    assert not hasattr(ccatalog, '_kdtree')

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog, storekdtree=True)
    assert hasattr(ccatalog, '_kdtree')

    assert not hasattr(ccatalog, 'tilsitcheese')
    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog, storekdtree='tilsitcheese')
    assert hasattr(ccatalog, 'tilsitcheese')
    assert not hasattr(cmatch, 'tilsitcheese')


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching_method():
    from .. import ICRS
    from ...utils import NumpyRNGContext
    from ..matching import match_coordinates_3d, match_coordinates_sky

    with NumpyRNGContext(987654321):
        cmatch = ICRS(np.random.rand(20) * 360., np.random.rand(20) * 180. - 90., unit=(u.degree, u.degree))
        ccatalog = ICRS(np.random.rand(100) * 360., np.random.rand(100) * 180. - 90., unit=(u.degree, u.degree))

    idx1, d2d1, d3d1 = cmatch.match_to_catalog_3d(ccatalog)
    idx2, d2d2, d3d2 = match_coordinates_3d(cmatch, ccatalog)

    npt.assert_array_equal(idx1, idx2)
    npt.assert_array_equal(d2d1, d2d2)
    npt.assert_array_equal(d3d1, d3d2)


    #should be the same as above because there's no distance, but just make sure this method works
    idx1, d2d1, d3d1 = cmatch.match_to_catalog_sky(ccatalog)
    idx2, d2d2, d3d2 = match_coordinates_sky(cmatch, ccatalog)

    npt.assert_array_equal(idx1, idx2)
    npt.assert_array_equal(d2d1, d2d2)
    npt.assert_array_equal(d3d1, d3d2)


    assert len(idx1) == len(d2d1) == len(d3d1) == 20





