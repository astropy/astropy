# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
import pytest
from numpy import testing as npt

from astropy import units as u
from astropy.coordinates import (
    ICRS,
    Angle,
    CartesianRepresentation,
    Galactic,
    SkyCoord,
    match_coordinates_3d,
    match_coordinates_sky,
    search_around_3d,
    search_around_sky,
)
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.utils import NumpyRNGContext
from astropy.utils.compat.optional_deps import HAS_SCIPY

"""
These are the tests for coordinate matching.

Coordinate matching can involve caching, so it is best to recreate the
coordinate objects in every test instead of trying to reuse module-level
variables.
"""

if not HAS_SCIPY:
    pytest.skip("Coordinate matching requires scipy", allow_module_level=True)


def test_matching_function():
    # this only uses match_coordinates_3d because that's the actual implementation

    cmatch = ICRS([4, 2.1] * u.degree, [0, 0] * u.degree)
    ccatalog = ICRS([1, 2, 3, 4] * u.degree, [0, 0, 0, 0] * u.degree)

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog)
    npt.assert_array_equal(idx, [3, 1])
    npt.assert_array_almost_equal(d2d.degree, [0, 0.1])
    assert d3d.value[0] == 0

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog, nthneighbor=2)
    assert np.all(idx == 2)
    npt.assert_array_almost_equal(d2d.degree, [1, 0.9])
    npt.assert_array_less(d3d.value, 0.02)


def test_matching_function_3d_and_sky():
    cmatch = ICRS([4, 2.1] * u.degree, [0, 0] * u.degree, distance=[1, 5] * u.kpc)
    ccatalog = ICRS(
        [1, 2, 3, 4] * u.degree, [0, 0, 0, 0] * u.degree, distance=[1, 1, 1, 5] * u.kpc
    )

    idx, d2d, d3d = match_coordinates_3d(cmatch, ccatalog)
    npt.assert_array_equal(idx, [2, 3])

    assert_allclose(d2d, [1, 1.9] * u.deg)
    assert np.abs(d3d[0].to_value(u.kpc) - np.radians(1)) < 1e-6
    assert np.abs(d3d[1].to_value(u.kpc) - 5 * np.radians(1.9)) < 1e-5

    idx, d2d, d3d = match_coordinates_sky(cmatch, ccatalog)
    npt.assert_array_equal(idx, [3, 1])

    assert_allclose(d2d, [0, 0.1] * u.deg)
    assert_allclose(d3d, [4, 4.0000019] * u.kpc)


@pytest.mark.parametrize(
    "functocheck, args, defaultkdtname, bothsaved",
    [
        (match_coordinates_3d, [], "kdtree_3d", False),
        (match_coordinates_sky, [], "kdtree_sky", False),
        (search_around_3d, [1 * u.kpc], "kdtree_3d", True),
        (search_around_sky, [1 * u.deg], "kdtree_sky", False),
    ],
)
def test_kdtree_storage(functocheck, args, defaultkdtname, bothsaved):
    def make_scs():
        cmatch = ICRS([4, 2.1] * u.degree, [0, 0] * u.degree, distance=[1, 2] * u.kpc)
        ccatalog = ICRS(
            [1, 2, 3, 4] * u.degree,
            [0, 0, 0, 0] * u.degree,
            distance=[1, 2, 3, 4] * u.kpc,
        )
        return cmatch, ccatalog

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args, storekdtree=False)
    assert "kdtree" not in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args)
    assert defaultkdtname in ccatalog.cache
    assert "kdtree" not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    functocheck(cmatch, ccatalog, *args, storekdtree=True)
    assert "kdtree" in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache

    cmatch, ccatalog = make_scs()
    assert "tislit_cheese" not in ccatalog.cache
    functocheck(cmatch, ccatalog, *args, storekdtree="tislit_cheese")
    assert "tislit_cheese" in ccatalog.cache
    assert defaultkdtname not in ccatalog.cache
    assert "kdtree" not in ccatalog.cache
    if bothsaved:
        assert "tislit_cheese" in cmatch.cache
        assert defaultkdtname not in cmatch.cache
        assert "kdtree" not in cmatch.cache
    else:
        assert "tislit_cheese" not in cmatch.cache

    # now a bit of a hacky trick to make sure it at least tries to *use* it
    ccatalog.cache["tislit_cheese"] = 1
    cmatch.cache["tislit_cheese"] = 1
    with pytest.raises(TypeError) as e:
        functocheck(cmatch, ccatalog, *args, storekdtree="tislit_cheese")
    assert "KD" in e.value.args[0]


def test_matching_method():
    with NumpyRNGContext(987654321):
        cmatch = ICRS(
            np.random.rand(20) * 360.0 * u.degree,
            (np.random.rand(20) * 180.0 - 90.0) * u.degree,
        )
        ccatalog = ICRS(
            np.random.rand(100) * 360.0 * u.degree,
            (np.random.rand(100) * 180.0 - 90.0) * u.degree,
        )

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


@pytest.mark.parametrize(
    "search_limit,expected_idx1,expected_idx2,expected_d2d,expected_d3d",
    [
        pytest.param(
            1.01 * u.deg,
            [0, 0, 1, 1],
            [2, 3, 1, 2],
            [1, 0, 0.1, 0.9] * u.deg,
            [0.01745307, 4.0, 4.0000019, 4.00015421] * u.kpc,
            id="1.01_deg",
        ),
        pytest.param(
            [1.01, 1.01] * u.deg,
            [0, 0, 1, 1],
            [2, 3, 1, 2],
            [1, 0, 0.1, 0.9] * u.deg,
            [0.01745307, 4.0, 4.0000019, 4.00015421] * u.kpc,
            id="1.01_deg_as_multiple",
        ),
        pytest.param(
            [1.01, 0.05] * u.deg,
            [0, 0],
            [2, 3],
            [1, 0] * u.deg,
            [0.01745307, 4.0] * u.kpc,
            id="multiple",
        ),
        pytest.param(0.05 * u.deg, [0], [3], [0] * u.deg, [4] * u.kpc, id="0.05_deg"),
    ],
)
def test_search_around_sky(
    search_limit, expected_idx1, expected_idx2, expected_d2d, expected_d3d
):
    idx1, idx2, d2d, d3d = search_around_sky(
        ICRS([4, 2.1] * u.deg, [0, 0] * u.deg, distance=[1, 5] * u.kpc),
        ICRS([1, 2, 3, 4] * u.deg, [0, 0, 0, 0] * u.deg, distance=[1, 1, 1, 5] * u.kpc),
        search_limit,
    )
    npt.assert_array_equal(idx1, expected_idx1)
    npt.assert_array_equal(idx2, expected_idx2)
    assert_allclose(d2d, expected_d2d)
    assert_allclose(d3d, expected_d3d)


@pytest.mark.parametrize(
    "search_limit,expected_idx1,expected_idx2,expected_d2d,expected_d3d",
    [
        pytest.param(
            1 * u.kpc,
            [0, 0, 0, 1],
            [0, 1, 2, 3],
            [3, 2, 1, 1.9] * u.deg,
            [0.0523539, 0.03490481, 0.01745307, 0.16579868] * u.kpc,
            id="1_kpc",
        ),
        pytest.param(
            [1, 1] * u.kpc,
            [0, 0, 0, 1],
            [0, 1, 2, 3],
            [3, 2, 1, 1.9] * u.deg,
            [0.0523539, 0.03490481, 0.01745307, 0.16579868] * u.kpc,
            id="1_kpc",
        ),
        pytest.param(
            [1, 0.05] * u.kpc,
            [0, 0, 0],
            [0, 1, 2],
            [3, 2, 1] * u.deg,
            [0.0523539, 0.03490481, 0.01745307] * u.kpc,
            id="1_kpc",
        ),
        pytest.param(
            0.05 * u.kpc,
            [0, 0],
            [1, 2],
            [2, 1] * u.deg,
            [0.03490481, 0.01745307] * u.kpc,
            id="0.05_kpc",
        ),
    ],
)
def test_search_around_3d(
    search_limit, expected_idx1, expected_idx2, expected_d2d, expected_d3d
):
    idx1, idx2, d2d, d3d = search_around_3d(
        ICRS([4, 2.1] * u.deg, [0, 0] * u.deg, distance=[1, 5] * u.kpc),
        ICRS([1, 2, 3, 4] * u.deg, [0, 0, 0, 0] * u.deg, distance=[1, 1, 1, 5] * u.kpc),
        search_limit,
    )
    npt.assert_array_equal(idx1, expected_idx1)
    npt.assert_array_equal(idx2, expected_idx2)
    assert_allclose(d2d, expected_d2d)
    assert_allclose(d3d, expected_d3d)


@pytest.mark.parametrize(
    "function,search_limit",
    [
        pytest.param(func, limit, id=func.__name__)
        for func, limit in ([search_around_3d, 1 * u.m], [search_around_sky, 1 * u.deg])
    ],
)
def test_search_around_no_matches(function, search_limit):
    # Test for the non-matches, #4877
    idx1, idx2, d2d, d3d = function(
        ICRS([41, 21] * u.deg, [0, 0] * u.deg, distance=[1, 5] * u.kpc),
        ICRS([1, 2] * u.deg, [0, 0] * u.deg, distance=[1, 1] * u.kpc),
        search_limit,
    )
    assert idx1.size == 0
    assert idx2.size == 0
    assert d2d.size == 0
    assert d3d.size == 0
    assert idx1.dtype == int
    assert idx2.dtype == int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc


@pytest.mark.parametrize(
    "function,search_limit",
    [
        pytest.param(func, limit, id=func.__name__)
        for func, limit in ([search_around_3d, 1 * u.m], [search_around_sky, 1 * u.deg])
    ],
)
@pytest.mark.parametrize(
    "sources,catalog",
    [
        pytest.param(
            ICRS(ra=[] * u.deg, dec=[] * u.deg, distance=[] * u.kpc),
            ICRS([1] * u.deg, [0] * u.deg, distance=[1] * u.kpc),
            id="empty_sources",
        ),
        pytest.param(
            ICRS([1] * u.deg, [0] * u.deg, distance=[1] * u.kpc),
            ICRS(ra=[] * u.deg, dec=[] * u.deg, distance=[] * u.kpc),
            id="empty_catalog",
        ),
        pytest.param(
            ICRS(ra=[] * u.deg, dec=[] * u.deg, distance=[] * u.kpc),
            ICRS(ra=[] * u.deg, dec=[] * u.deg, distance=[] * u.kpc),
            id="empty_both",
        ),
    ],
)
def test_search_around_empty_input(sources, catalog, function, search_limit):
    # Test when one or both of the coordinate arrays is empty, #4875
    idx1, idx2, d2d, d3d = function(sources, catalog, search_limit)
    assert idx1.size == 0
    assert idx2.size == 0
    assert d2d.size == 0
    assert d3d.size == 0
    assert idx1.dtype == int
    assert idx2.dtype == int
    assert d2d.unit == u.deg
    assert d3d.unit == u.kpc


@pytest.mark.parametrize(
    "sources,catalog",
    [
        pytest.param(
            ICRS([1] * u.deg, [0] * u.deg),
            ICRS([1] * u.deg, [0] * u.deg),
            id="both_with_data",
        ),
        pytest.param(
            ICRS(ra=[] * u.deg, dec=[] * u.deg),
            ICRS([1] * u.deg, [0] * u.deg),
            id="empty_sources",
        ),
        pytest.param(
            ICRS([1] * u.deg, [0] * u.deg),
            ICRS(ra=[] * u.deg, dec=[] * u.deg),
            id="empty_catalog",
        ),
        pytest.param(
            ICRS(ra=[] * u.deg, dec=[] * u.deg),
            ICRS(ra=[] * u.deg, dec=[] * u.deg),
            id="both_empty",
        ),
    ],
)
def test_search_around_3d_no_dist_input(sources, catalog):
    # Regression test for #16280: UnitConversionError was not raised if at
    # least one of the coords was empty.
    with pytest.raises(
        u.UnitConversionError,
        match=r"^'pc' \(length\) and '' \(dimensionless\) are not convertible$",
    ):
        search_around_3d(sources, catalog, 1 * u.pc)


def test_search_around_sky_no_dist_input():
    # Test that input without distance units results in a
    # 'dimensionless_unscaled' unit
    empty_sc = SkyCoord([], [], unit=u.deg)
    idx1, idx2, d2d, d3d = search_around_sky(empty_sc, empty_sc[:], 1 * u.deg)
    assert d2d.unit == u.deg
    assert d3d.unit == u.dimensionless_unscaled


def test_search_around_scalar():
    cat = SkyCoord([1, 2, 3], [-30, 45, 8], unit="deg")
    target = SkyCoord("1.1 -30.1", unit="deg")

    with pytest.raises(ValueError) as excinfo:
        cat.search_around_sky(target, Angle("2d"))

    # make sure the error message is *specific* to search_around_sky rather than
    # generic as reported in #3359
    assert "search_around_sky" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        cat.search_around_3d(target, Angle("2d"))
    assert "search_around_3d" in str(excinfo.value)


def test_search_around_multidimensional():
    # search around methods only accept 1-dimensional coordinates, see #17824
    coo = SkyCoord([[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]], unit="deg")
    with pytest.raises(
        ValueError,
        match=("search_around_3d only supports 1-dimensional coordinate arrays.*"),
    ):
        coo.search_around_3d(coo, Angle("1d"))
    with pytest.raises(
        ValueError,
        match=("search_around_sky only supports 1-dimensional coordinate arrays.*"),
    ):
        coo.search_around_sky(coo, Angle("1d"))


def test_match_catalog_empty():
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
    assert "catalog" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(cat1[0])
    assert "catalog" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_sky(cat0)
    assert "catalog" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(cat0)
    assert "catalog" in str(excinfo.value)


@pytest.mark.filterwarnings(r"ignore:invalid value encountered in.*:RuntimeWarning")
def test_match_catalog_nan():
    sc1 = SkyCoord(1, 2, unit="deg")
    sc_with_nans = SkyCoord(1, np.nan, unit="deg")

    cat = SkyCoord([1.1, 3], [2.1, 5], unit="deg")
    cat_with_nans = SkyCoord([1.1, np.nan], [2.1, 5], unit="deg")
    galcat_with_nans = Galactic([1.2, np.nan] * u.deg, [5.6, 7.8] * u.deg)

    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_sky(cat_with_nans)
    assert "Catalog coordinates cannot contain" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(cat_with_nans)
    assert "Catalog coordinates cannot contain" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_sky(galcat_with_nans)
    assert "Catalog coordinates cannot contain" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc1.match_to_catalog_3d(galcat_with_nans)
    assert "Catalog coordinates cannot contain" in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        sc_with_nans.match_to_catalog_sky(cat)
    assert "Matching coordinates cannot contain" in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        sc_with_nans.match_to_catalog_3d(cat)
    assert "Matching coordinates cannot contain" in str(excinfo.value)


def test_match_catalog_nounit():
    i1 = ICRS([[1], [2], [3]], representation_type=CartesianRepresentation)
    i2 = ICRS([[1], [2], [4, 5]], representation_type=CartesianRepresentation)
    i, sep, sep3d = match_coordinates_sky(i1, i2)
    assert_allclose(sep3d, [1] * u.dimensionless_unscaled)
