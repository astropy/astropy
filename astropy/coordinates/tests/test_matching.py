# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy import testing as npt

from ...tests.helper import assert_quantity_allclose as assert_allclose

from ... import units as u
from ...utils import minversion

from .. import matching

"""
These are the tests for coordinate matching.

Note that this requires scipy.
"""

try:
    import scipy
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

if HAS_SCIPY and minversion(scipy, '0.12.0', inclusive=False):
    OLDER_SCIPY = False
else:
    OLDER_SCIPY = True


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching_function():
    from .. import ICRS
    from ..matching import match_coordinates_3d
    # this only uses match_coordinates_3d because that's the actual implementation

    cmatch = ICRS([4, 2.1]*u.degree, [0, 0]*u.degree)
    ccatalog = ICRS([1, 2, 3, 4]*u.degree, [0, 0, 0, 0]*u.degree)

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

    cmatch = ICRS([4, 2.1]*u.degree, [0, 0]*u.degree, distance=[1, 5] * u.kpc)
    ccatalog = ICRS([1, 2, 3, 4]*u.degree, [0, 0, 0, 0]*u.degree, distance=[1, 1, 1, 5] * u.kpc)

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog)
    npt.assert_array_equal(idx, [2, 3])

    assert_allclose(d2d, [1, 1.9] * u.deg)
    assert np.abs(d3d[0].to_value(u.kpc) - np.radians(1)) < 1e-6
    assert np.abs(d3d[1].to_value(u.kpc) - 5*np.radians(1.9)) < 1e-5

    idx, d2d, d3d = match_coordinates_sky(cmatch, ccatalog)
    npt.assert_array_equal(idx, [3, 1])

    assert_allclose(d2d, [0, 0.1] * u.deg)
    assert_allclose(d3d, [4, 4.0000019] * u.kpc)


@pytest.mark.parametrize('functocheck, args, defaultkdtname, bothsaved',
                         [(matching.match_coordinates_3d, [], 'kdtree_3d', False),
                          (matching.match_coordinates_sky, [], 'kdtree_sky', False),
                          (matching.search_around_3d, [1*u.kpc], 'kdtree_3d', True),
                          (matching.search_around_sky, [1*u.deg], 'kdtree_sky', False)
                         ])
@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_kdtree_storage(functocheck, args, defaultkdtname, bothsaved):
    from .. import ICRS

    def make_scs():
        cmatch = ICRS([4, 2.1]*u.degree, [0, 0]*u.degree, distance=[1, 2]*u.kpc)
        ccatalog = ICRS([1, 2, 3, 4]*u.degree, [0, 0, 0, 0]*u.degree, distance=[1, 2, 3, 4]*u.kpc)
        return cmatch, ccatalog

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args, storekdtree=False)
    assert 'kdtree' not in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args)
    assert defaultkdtname in ccatalog.cache
    assert 'kdtree' not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args, storekdtree=True)
    assert 'kdtree' in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    assert 'tislit_cheese' not in ccatalog.cache
    functocheck(cmatch, ccatalog, *args, storekdtree='tislit_cheese')
    assert 'tislit_cheese' in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache
    assert 'kdtree' not in ccatalog.cache
    if bothsaved:
        assert 'tislit_cheese' in cmatch.cache
        assert defaultkdtname not in cmatch.cache
        assert 'kdtree' not in cmatch.cache
    else:
        assert 'tislit_cheese' not in cmatch.cache

    # now a bit of a hacky trick to make sure it at least tries to *use* it
    ccatalog.cache['tislit_cheese'] = 1
    cmatch.cache['tislit_cheese'] = 1
    with pytest.raises(TypeError) as e:
        functocheck(cmatch, ccatalog, *args, storekdtree='tislit_cheese')
    assert 'KD' in e.value.args[0]


@pytest.mark.skipif(str('not HAS_SCIPY'))
def test_matching_method():
    from .. import ICRS, SkyCoord
    from ...utils import NumpyRNGContext
    from ..matching import match_coordinates_3d, match_coordinates_sky

    with NumpyRNGContext(987654321):
        cmatch = ICRS(np.random.rand(20) * 360.*u.degree,
                      (np.random.rand(20) * 180. - 90.)*u.degree)
        ccatalog = ICRS(np.random.rand(100) * 360. * u.degree,
                       (np.random.rand(100) * 180. - 90.)*u.degree)

    idx1, d2d1, d3d1 = SkyCoord(cmatch).match_to_catalog_3d(ccatalog)
    idx2, d2d2, d3d2 = match_coordinates_3d(cmatch, ccatalog)

    npt.assert_array_equal(idx1, idx2)
    assert_allclose(d2d1, d2d2)
    assert_allclose(d3d1, d3d2)

    # should be the same as above because there's no distance, but just make sure this method works
    idx1, d2d1, d3d1 = SkyCoord(cmatch).match_to_catalog_sky(ccatalog)
    idx2, d2d2, d3d2 = match_coordinates_sky(cmatch, ccatalog)

    npt.assert_array_equal(idx1, idx2)
    assert_allclose(d2d1, d2d2)
    assert_allclose(d3d1, d3d2)

    assert len(idx1) == len(d2d1) == len(d3d1) == 20


@pytest.mark.skipif(str('not HAS_SCIPY'))
@pytest.mark.skipif(str('OLDER_SCIPY'))
def test_search_around():
    from .. import ICRS, SkyCoord
    from ..matching import search_around_sky, search_around_3d

    coo1 = ICRS([4, 2.1]*u.degree, [0, 0]*u.degree, distance=[1, 5] * u.kpc)
    coo2 = ICRS([1, 2, 3, 4]*u.degree, [0, 0, 0, 0]*u.degree, distance=[1, 1, 1, 5] * u.kpc)

    idx1_1deg, idx2_1deg, d2d_1deg, d3d_1deg = search_around_sky(coo1, coo2, 1.01*u.deg)
    idx1_0p05deg, idx2_0p05deg, d2d_0p05deg, d3d_0p05deg = search_around_sky(coo1, coo2, 0.05*u.deg)

    assert list(zip(idx1_1deg, idx2_1deg)) == [(0, 2), (0, 3), (1, 1), (1, 2)]
    assert d2d_1deg[0] == 1.0*u.deg
    assert_allclose(d2d_1deg, [1, 0, .1, .9]*u.deg)

    assert list(zip(idx1_0p05deg, idx2_0p05deg)) == [(0, 3)]

    idx1_1kpc, idx2_1kpc, d2d_1kpc, d3d_1kpc = search_around_3d(coo1, coo2, 1*u.kpc)
    idx1_sm, idx2_sm, d2d_sm, d3d_sm = search_around_3d(coo1, coo2, 0.05*u.kpc)

    assert list(zip(idx1_1kpc, idx2_1kpc)) == [(0, 0), (0, 1), (0, 2), (1, 3)]
    assert list(zip(idx1_sm, idx2_sm)) == [(0, 1), (0, 2)]
    assert_allclose(d2d_sm, [2, 1]*u.deg)

    # Test for the non-matches, #4877
    coo1 = ICRS([4.1, 2.1]*u.degree, [0, 0]*u.degree, distance=[1, 5] * u.kpc)
    idx1, idx2, d2d, d3d = search_around_sky(coo1, coo2, 1*u.arcsec)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    idx1, idx2, d2d, d3d = search_around_3d(coo1, coo2, 1*u.m)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc

    # Test when one or both of the coordinate arrays is empty, #4875
    empty = ICRS(ra=[] * u.degree, dec=[] * u.degree, distance=[] * u.kpc)
    idx1, idx2, d2d, d3d = search_around_sky(empty, coo2, 1*u.arcsec)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    idx1, idx2, d2d, d3d = search_around_sky(coo1, empty, 1*u.arcsec)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    empty = ICRS(ra=[] * u.degree, dec=[] * u.degree, distance=[] * u.kpc)
    idx1, idx2, d2d, d3d = search_around_sky(empty, empty[:], 1*u.arcsec)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    idx1, idx2, d2d, d3d = search_around_3d(empty, coo2, 1*u.m)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    idx1, idx2, d2d, d3d = search_around_3d(coo1, empty, 1*u.m)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc
    idx1, idx2, d2d, d3d = search_around_3d(empty, empty[:], 1*u.m)
    assert idx1.size == idx2.size == d2d.size == d3d.size == 0
    assert idx1.dtype == idx2.dtype == np.int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc

    # Test that input without distance units results in a
    # 'dimensionless_unscaled' unit
    cempty = SkyCoord(ra=[], dec=[], unit=u.deg)
    idx1, idx2, d2d, d3d = search_around_3d(cempty, cempty[:], 1*u.m)
    assert d2d.unit == u.deg
    assert d3d.unit == u.dimensionless_unscaled
    idx1, idx2, d2d, d3d = search_around_sky(cempty, cempty[:], 1*u.m)
    assert d2d.unit == u.deg
    assert d3d.unit == u.dimensionless_unscaled


@pytest.mark.skipif(str('not HAS_SCIPY'))
@pytest.mark.skipif(str('OLDER_SCIPY'))
def test_search_around_scalar():
    from astropy.coordinates import SkyCoord, Angle

    cat = SkyCoord([1, 2, 3], [-30, 45, 8], unit="deg")
    target = SkyCoord('1.1 -30.1', unit="deg")

    with pytest.raises(ValueError) as excinfo:
        cat.search_around_sky(target, Angle('2d'))

    # make sure the error message is *specific* to search_around_sky rather than
    # generic as reported in #3359
    assert 'search_around_sky' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        cat.search_around_3d(target, Angle('2d'))
    assert 'search_around_3d' in str(excinfo.value)


@pytest.mark.skipif(str('not HAS_SCIPY'))
@pytest.mark.skipif(str('OLDER_SCIPY'))
def test_match_catalog_empty():
    from astropy.coordinates import SkyCoord

    sc1 = SkyCoord(1, 2, unit="deg")
    cat0 = SkyCoord([], [], unit="deg")
    cat1 = SkyCoord([1.1], [2.1], unit="deg")
    cat2 = SkyCoord([1.1, 3], [2.1, 5], unit="deg")

    sc1.match_to_catalog_sky(cat2)
    sc1.match_to_catalog_3d(cat2)

    sc1.match_to_catalog_sky(cat1)
    sc1.match_to_catalog_3d(cat1)

    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_sky(cat1[0])
    assert 'catalog' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(cat1[0])
    assert 'catalog' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_sky(cat0)
    assert 'catalog' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(cat0)
    assert 'catalog' in str(excinfo.value)
