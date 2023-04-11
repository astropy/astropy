# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
import pytest

from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates import transformations as t
from astropy.coordinates.attributes import Attribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, frame_transform_graph
from astropy.coordinates.builtin_frames import (
    FK4,
    FK5,
    HCRS,
    ICRS,
    AltAz,
    FK4NoETerms,
    Galactic,
)
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.tests.helper import assert_quantity_allclose as assert_allclose
from astropy.time import Time
from astropy.units import allclose as quantity_allclose
from astropy.utils.exceptions import AstropyWarning

CARTESIAN_POS = r.CartesianRepresentation([1, 2, 3] * u.kpc)
CARTESIAN_VEL = r.CartesianDifferential([8, 9, 10] * u.km / u.s)
CARTESIAN_POS_AND_VEL = CARTESIAN_POS.with_differentials(CARTESIAN_VEL)

RADIAL_VEL = r.RadialDifferential(1 * u.km / u.s)
SPHERICAL_COS_LAT_VEL = r.SphericalCosLatDifferential(
    1 * u.mas / u.yr, 2 * u.mas / u.yr, 3 * u.km / u.s
)
SPHERICAL_POS = r.SphericalRepresentation(
    lon=1 * u.deg, lat=2.0 * u.deg, distance=10 * u.pc
)
UNIT_SPHERICAL_POS = r.UnitSphericalRepresentation(lon=1 * u.deg, lat=2.0 * u.deg)

ROT_30 = rotation_matrix(30 * u.deg)
ROT_45 = rotation_matrix(45 * u.deg)
ROT_75 = rotation_matrix(75 * u.deg)
OFFSET_X = r.CartesianRepresentation([1, 0, 0])
OFFSET_Z = r.CartesianRepresentation([0, 0, 1])
OFFSET_123 = r.CartesianRepresentation([1, 2, 3])
OFFSET_456 = r.CartesianRepresentation([4, 5, 6])
OFFSET_579 = r.CartesianRepresentation([5, 7, 9])
SQRT_2 = np.sqrt(2)


# Coordinates just for these tests.
class TCoo1(ICRS):
    pass


class TCoo2(ICRS):
    pass


class TCoo3(ICRS):
    pass


def test_transform_classes():
    """
    Tests the class-based/OO syntax for creating transforms
    """

    def tfun(c, f):
        return f.__class__(ra=c.ra, dec=c.dec)

    _ = t.FunctionTransform(tfun, TCoo1, TCoo2, register_graph=frame_transform_graph)

    c1 = TCoo1(ra=1 * u.radian, dec=0.5 * u.radian)
    c2 = c1.transform_to(TCoo2())
    assert_allclose(c2.ra.radian, 1)
    assert_allclose(c2.dec.radian, 0.5)

    def matfunc(coo, fr):
        return [[1, 0, 0], [0, coo.ra.degree, 0], [0, 0, 1]]

    trans2 = t.DynamicMatrixTransform(matfunc, TCoo1, TCoo2)
    trans2.register(frame_transform_graph)

    c3 = TCoo1(ra=1 * u.deg, dec=2 * u.deg)
    c4 = c3.transform_to(TCoo2())

    assert_allclose(c4.ra.degree, 1)
    assert_allclose(c4.ra.degree, 1)

    # be sure to unregister the second one - no need for trans1 because it
    # already got unregistered when trans2 was created.
    trans2.unregister(frame_transform_graph)


def test_transform_decos():
    """
    Tests the decorator syntax for creating transforms
    """
    c1 = TCoo1(ra=1 * u.deg, dec=2 * u.deg)

    @frame_transform_graph.transform(t.FunctionTransform, TCoo1, TCoo2)
    def trans(coo1, f):
        return TCoo2(ra=coo1.ra, dec=coo1.dec * 2)

    c2 = c1.transform_to(TCoo2())
    assert_allclose(c2.ra.degree, 1)
    assert_allclose(c2.dec.degree, 4)

    c3 = TCoo1(r.CartesianRepresentation(x=1 * u.pc, y=1 * u.pc, z=2 * u.pc))

    @frame_transform_graph.transform(t.StaticMatrixTransform, TCoo1, TCoo2)
    def matrix():
        return [[2, 0, 0], [0, 1, 0], [0, 0, 1]]

    c4 = c3.transform_to(TCoo2())

    assert_allclose(c4.cartesian.x, 2 * u.pc)
    assert_allclose(c4.cartesian.y, 1 * u.pc)
    assert_allclose(c4.cartesian.z, 2 * u.pc)


def test_shortest_path():
    class FakeTransform:
        def __init__(self, pri):
            self.priority = pri

    g = t.TransformGraph()

    # cheating by adding graph elements directly that are not classes - the
    # graphing algorithm still works fine with integers - it just isn't a valid
    # TransformGraph

    # the graph looks is a down-going diamond graph with the lower-right slightly
    # heavier and a cycle from the bottom to the top
    # also, a pair of nodes isolated from 1

    g._graph[1][2] = FakeTransform(1)
    g._graph[1][3] = FakeTransform(1)
    g._graph[2][4] = FakeTransform(1)
    g._graph[3][4] = FakeTransform(2)
    g._graph[4][1] = FakeTransform(5)

    g._graph[5][6] = FakeTransform(1)

    path, d = g.find_shortest_path(1, 2)
    assert path == [1, 2]
    assert d == 1
    path, d = g.find_shortest_path(1, 3)
    assert path == [1, 3]
    assert d == 1
    path, d = g.find_shortest_path(1, 4)
    print("Cached paths:", g._shortestpaths)
    assert path == [1, 2, 4]
    assert d == 2

    # unreachable
    path, d = g.find_shortest_path(1, 5)
    assert path is None
    assert d == float("inf")

    path, d = g.find_shortest_path(5, 6)
    assert path == [5, 6]
    assert d == 1


def test_sphere_cart():
    """
    Tests the spherical <-> cartesian transform functions
    """
    from astropy.coordinates import cartesian_to_spherical, spherical_to_cartesian
    from astropy.utils import NumpyRNGContext

    x, y, z = spherical_to_cartesian(1, 0, 0)
    assert_allclose(x, 1)
    assert_allclose(y, 0)
    assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(0, 1, 1)
    assert_allclose(x, 0)
    assert_allclose(y, 0)
    assert_allclose(z, 0)

    x, y, z = spherical_to_cartesian(5, 0, np.arcsin(4.0 / 5.0))
    assert_allclose(x, 3)
    assert_allclose(y, 4)
    assert_allclose(z, 0)

    r, lat, lon = cartesian_to_spherical(0, 1, 0)
    assert_allclose(r, 1)
    assert_allclose(lat, 0 * u.deg)
    assert_allclose(lon, np.pi / 2 * u.rad)

    # test round-tripping
    with NumpyRNGContext(13579):
        x, y, z = np.random.randn(3, 5)

    x2, y2, z2 = spherical_to_cartesian(*cartesian_to_spherical(x, y, z))

    assert_allclose(x, x2)
    assert_allclose(y, y2)
    assert_allclose(z, z2)


def test_transform_path_pri():
    """
    This checks that the transformation path prioritization works by
    making sure the ICRS -> Gal transformation always goes through FK5
    and not FK4.
    """
    frame_transform_graph.invalidate_cache()
    tpath, td = frame_transform_graph.find_shortest_path(ICRS, Galactic)
    assert tpath == [ICRS, FK5, Galactic]
    assert td == 2

    # but direct from FK4 to Galactic should still be possible
    tpath, td = frame_transform_graph.find_shortest_path(FK4, Galactic)
    assert tpath == [FK4, FK4NoETerms, Galactic]
    assert td == 2


def test_obstime():
    """
    Checks to make sure observation time is
    accounted for at least in FK4 <-> ICRS transformations
    """
    b1950 = Time("B1950")
    j1975 = Time("J1975")

    fk4_50 = FK4(ra=1 * u.deg, dec=2 * u.deg, obstime=b1950)
    fk4_75 = FK4(ra=1 * u.deg, dec=2 * u.deg, obstime=j1975)

    icrs_50 = fk4_50.transform_to(ICRS())
    icrs_75 = fk4_75.transform_to(ICRS())

    # now check that the resulting coordinates are *different* - they should be,
    # because the obstime is different
    assert icrs_50.ra.degree != icrs_75.ra.degree
    assert icrs_50.dec.degree != icrs_75.dec.degree


# ------------------------------------------------------------------------------
# Affine transform tests and helpers:

# just acting as a namespace


class transfunc:
    rep = r.CartesianRepresentation(np.arange(3) * u.pc)
    dif = r.CartesianDifferential(*np.arange(3, 6) * u.pc / u.Myr)
    rep0 = r.CartesianRepresentation(np.zeros(3) * u.pc)

    @classmethod
    def both(cls, coo, fr):
        # exchange x <-> z and offset
        M = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        return M, cls.rep.with_differentials(cls.dif)

    @classmethod
    def just_matrix(cls, coo, fr):
        # exchange x <-> z and offset
        M = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
        return M, None

    @classmethod
    def no_matrix(cls, coo, fr):
        return None, cls.rep.with_differentials(cls.dif)

    @classmethod
    def no_pos(cls, coo, fr):
        return None, cls.rep0.with_differentials(cls.dif)

    @classmethod
    def no_vel(cls, coo, fr):
        return None, cls.rep


@pytest.mark.parametrize(
    "transfunc",
    [
        transfunc.both,
        transfunc.no_matrix,
        transfunc.no_pos,
        transfunc.no_vel,
        transfunc.just_matrix,
    ],
)
@pytest.mark.parametrize(
    "rep",
    (
        CARTESIAN_POS,
        CARTESIAN_POS_AND_VEL,
        CARTESIAN_POS_AND_VEL.represent_as(
            r.CylindricalRepresentation, r.CylindricalDifferential
        ),
    ),
)
def test_affine_transform_succeed(transfunc, rep):
    c = TCoo1(rep)

    # compute expected output
    M, offset = transfunc(c, TCoo2)

    expected_rep = rep.to_cartesian().with_differentials(
        {
            k: diff.represent_as(r.CartesianDifferential, rep)
            for k, diff in rep.differentials.items()
        }
    )

    if M is not None:
        expected_rep = expected_rep.transform(M)

    expected_pos = expected_rep.without_differentials()
    if offset is not None:
        expected_pos += offset.without_differentials()

    expected_vel = None
    if c.data.differentials:
        expected_vel = expected_rep.differentials["s"]

        if offset and offset.differentials:
            expected_vel += offset.differentials["s"]

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    c2 = c.transform_to(TCoo2())

    assert quantity_allclose(
        c2.data.to_cartesian().xyz, expected_pos.to_cartesian().xyz
    )

    if expected_vel is not None:
        diff = c2.data.differentials["s"].to_cartesian(base=c2.data)
        assert quantity_allclose(diff.xyz, expected_vel.d_xyz)

    trans.unregister(frame_transform_graph)


# these should fail
def transfunc_invalid_matrix(coo, fr):
    return np.eye(4), None


# Leaving this open in case we want to add more functions to check for failures


@pytest.mark.parametrize("transfunc", [transfunc_invalid_matrix])
def test_affine_transform_fail(transfunc):
    c = TCoo1(CARTESIAN_POS_AND_VEL)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    with pytest.raises(ValueError):
        c.transform_to(TCoo2())

    trans.unregister(frame_transform_graph)


def test_too_many_differentials():
    dif2 = r.CartesianDifferential(*np.arange(3, 6) * u.pc / u.Myr**2)
    rep = CARTESIAN_POS_AND_VEL.with_differentials(dif2)

    with pytest.raises(ValueError):
        c = TCoo1(rep)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc.both, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    # Check that if frame somehow gets through to transformation, multiple
    # differentials are caught
    c = TCoo1(rep.without_differentials())
    c._data = c._data.with_differentials({"s": CARTESIAN_VEL, "s2": dif2})
    with pytest.raises(ValueError):
        c.transform_to(TCoo2())

    trans.unregister(frame_transform_graph)


# A matrix transform of a unit spherical with differentials should work


@pytest.mark.parametrize(
    "rep",
    (
        UNIT_SPHERICAL_POS.with_differentials(SPHERICAL_COS_LAT_VEL),
        r.UnitSphericalRepresentation(
            UNIT_SPHERICAL_POS, differentials={"s": RADIAL_VEL}
        ),
        SPHERICAL_POS.with_differentials(RADIAL_VEL),
    ),
)
def test_unit_spherical_with_differentials(rep):
    c = TCoo1(rep)

    # register and do the transformation and check against expected
    trans = t.AffineTransform(transfunc.just_matrix, TCoo1, TCoo2)
    trans.register(frame_transform_graph)
    c2 = c.transform_to(TCoo2())

    assert "s" in rep.differentials
    assert isinstance(c2.data.differentials["s"], type(rep.differentials["s"]))

    if isinstance(rep.differentials["s"], r.RadialDifferential):
        assert c2.data.differentials["s"] is rep.differentials["s"]

    trans.unregister(frame_transform_graph)

    # should fail if we have to do offsets
    trans = t.AffineTransform(transfunc.both, TCoo1, TCoo2)
    trans.register(frame_transform_graph)

    with pytest.raises(TypeError):
        c.transform_to(TCoo2())

    trans.unregister(frame_transform_graph)


def test_vel_transformation_obstime_err():
    # TODO: replace after a final decision on PR #6280
    from astropy.coordinates.sites import get_builtin_sites

    diff = r.CartesianDifferential([0.1, 0.2, 0.3] * u.km / u.s)
    rep = r.CartesianRepresentation([1, 2, 3] * u.au, differentials=diff)

    loc = get_builtin_sites()["example_site"]

    aaf = AltAz(obstime="J2010", location=loc)
    aaf2 = AltAz(obstime=aaf.obstime + 3 * u.day, location=loc)
    aaf3 = AltAz(obstime=aaf.obstime + np.arange(3) * u.day, location=loc)
    aaf4 = AltAz(obstime=aaf.obstime, location=loc)

    aa = aaf.realize_frame(rep)

    with pytest.raises(NotImplementedError) as exc:
        aa.transform_to(aaf2)
    assert "cannot transform" in exc.value.args[0]

    with pytest.raises(NotImplementedError) as exc:
        aa.transform_to(aaf3)
    assert "cannot transform" in exc.value.args[0]

    aa.transform_to(aaf4)

    aa.transform_to(ICRS())


def test_function_transform_with_differentials():
    def tfun(c, f):
        return f.__class__(ra=c.ra, dec=c.dec)

    _ = t.FunctionTransform(tfun, TCoo3, TCoo2, register_graph=frame_transform_graph)

    t3 = TCoo3(
        ra=1 * u.deg,
        dec=2 * u.deg,
        pm_ra_cosdec=1 * u.marcsec / u.yr,
        pm_dec=1 * u.marcsec / u.yr,
    )

    with pytest.warns(AstropyWarning, match=r".*they have been dropped.*") as w:
        t3.transform_to(TCoo2())
    assert len(w) == 1


def test_frame_override_component_with_attribute():
    """
    It was previously possible to define a frame with an attribute with the
    same name as a component. We don't want to allow this!
    """

    class BorkedFrame(BaseCoordinateFrame):
        ra = Attribute(default=150)
        dec = Attribute(default=150)

    def trans_func(coo1, f):
        pass

    trans = t.FunctionTransform(trans_func, BorkedFrame, ICRS)
    with pytest.raises(ValueError) as exc:
        trans.register(frame_transform_graph)

    assert (
        "BorkedFrame" in exc.value.args[0]
        and "'ra'" in exc.value.args[0]
        and "'dec'" in exc.value.args[0]
    )


def test_static_matrix_combine_paths():
    """
    Check that combined staticmatrixtransform matrices provide the same
    transformation as using an intermediate transformation.

    This is somewhat of a regression test for #7706
    """

    class AFrame(BaseCoordinateFrame):
        default_representation = r.SphericalRepresentation
        default_differential = r.SphericalCosLatDifferential

    t1 = t.StaticMatrixTransform(rotation_matrix(30.0 * u.deg, "z"), ICRS, AFrame)
    t1.register(frame_transform_graph)
    t2 = t.StaticMatrixTransform(rotation_matrix(30.0 * u.deg, "z").T, AFrame, ICRS)
    t2.register(frame_transform_graph)

    class BFrame(BaseCoordinateFrame):
        default_representation = r.SphericalRepresentation
        default_differential = r.SphericalCosLatDifferential

    t3 = t.StaticMatrixTransform(rotation_matrix(30.0 * u.deg, "x"), ICRS, BFrame)
    t3.register(frame_transform_graph)
    t4 = t.StaticMatrixTransform(rotation_matrix(30.0 * u.deg, "x").T, BFrame, ICRS)
    t4.register(frame_transform_graph)

    c = Galactic(123 * u.deg, 45 * u.deg)
    c_direct = c.transform_to(BFrame())
    c_through_A = c.transform_to(AFrame()).transform_to(BFrame())
    c_through_ICRS = c.transform_to(ICRS()).transform_to(BFrame())

    assert quantity_allclose(c_direct.lon, c_through_A.lon)
    assert quantity_allclose(c_direct.lat, c_through_A.lat)

    assert quantity_allclose(c_direct.lon, c_through_ICRS.lon)
    assert quantity_allclose(c_direct.lat, c_through_ICRS.lat)

    for t_ in [t1, t2, t3, t4]:
        t_.unregister(frame_transform_graph)


def test_multiple_aliases():
    # Define a frame with multiple aliases
    class MultipleAliasesFrame(BaseCoordinateFrame):
        name = ["alias_1", "alias_2"]
        default_representation = r.SphericalRepresentation

    def tfun(c, f):
        return f.__class__(lon=c.lon, lat=c.lat)

    # Register a transform
    graph = t.TransformGraph()
    _ = t.FunctionTransform(
        tfun, MultipleAliasesFrame, MultipleAliasesFrame, register_graph=graph
    )

    # Test that both aliases have been added to the transform graph
    assert graph.lookup_name("alias_1") == MultipleAliasesFrame
    assert graph.lookup_name("alias_2") == MultipleAliasesFrame

    # Test that both aliases appear in the graphviz DOT format output
    dotstr = graph.to_dot_graph()
    assert "`alias_1`\\n`alias_2`" in dotstr


def test_remove_transform_and_unregister():
    def tfun(c, f):
        f.__class__(ra=c.ra, dec=c.dec)

    # Register transforms
    graph = t.TransformGraph()
    ftrans1 = t.FunctionTransform(tfun, TCoo1, TCoo1, register_graph=graph)
    ftrans2 = t.FunctionTransform(tfun, TCoo2, TCoo2, register_graph=graph)
    _ = t.FunctionTransform(tfun, TCoo1, TCoo2, register_graph=graph)

    # Confirm that the frames are part of the graph
    assert TCoo1 in graph.frame_set
    assert TCoo2 in graph.frame_set

    # Use all three ways to remove a transform

    # Remove the only transform with TCoo2 as the "from" frame
    ftrans2.unregister(graph)
    # TCoo2 should still be part of the graph because it is the "to" frame of a transform
    assert TCoo2 in graph.frame_set

    # Remove the remaining transform that involves TCoo2
    graph.remove_transform(TCoo1, TCoo2, None)
    # Now TCoo2 should not be part of the graph
    assert TCoo2 not in graph.frame_set

    # Remove the remaining  transform that involves TCoo1
    graph.remove_transform(None, None, ftrans1)
    # Now TCoo1 should not be part of the graph
    assert TCoo1 not in graph.frame_set


def test_remove_transform_errors():
    def tfun(c, f):
        return f.__class__(ra=c.ra, dec=c.dec)

    graph = t.TransformGraph()
    _ = t.FunctionTransform(tfun, TCoo1, TCoo1, register_graph=graph)

    # Test bad calls to remove_transform

    with pytest.raises(ValueError):
        graph.remove_transform(None, TCoo1, None)

    with pytest.raises(ValueError):
        graph.remove_transform(TCoo1, None, None)

    with pytest.raises(ValueError):
        graph.remove_transform(None, None, None)

    with pytest.raises(ValueError):
        graph.remove_transform(None, None, 1)

    with pytest.raises(ValueError):
        graph.remove_transform(TCoo1, TCoo1, 1)


def test_impose_finite_difference_dt():
    class H1(HCRS):
        pass

    class H2(HCRS):
        pass

    class H3(HCRS):
        pass

    graph = t.TransformGraph()
    tfun = lambda c, f: type(f)(ra=c.ra, dec=c.dec)

    # Set up a number of transforms with different time steps
    old_dt = 1 * u.min
    transform1 = t.FunctionTransformWithFiniteDifference(
        tfun, H1, H1, register_graph=graph, finite_difference_dt=old_dt
    )
    transform2 = t.FunctionTransformWithFiniteDifference(
        tfun, H2, H2, register_graph=graph, finite_difference_dt=old_dt * 2
    )
    transform3 = t.FunctionTransformWithFiniteDifference(
        tfun, H2, H3, register_graph=graph, finite_difference_dt=old_dt * 3
    )

    # Check that all of the transforms have the same new time step
    new_dt = 1 * u.yr
    with graph.impose_finite_difference_dt(new_dt):
        assert transform1.finite_difference_dt == new_dt
        assert transform2.finite_difference_dt == new_dt
        assert transform3.finite_difference_dt == new_dt

    # Check that all of the original time steps have been restored
    assert transform1.finite_difference_dt == old_dt
    assert transform2.finite_difference_dt == old_dt * 2
    assert transform3.finite_difference_dt == old_dt * 3


@pytest.mark.parametrize(
    "first,second,check",
    (
        ([ROT_30, None], [ROT_45, None], [ROT_75, None]),
        ([ROT_30, None], [ROT_45, OFFSET_Z], [ROT_75, OFFSET_Z]),
        ([ROT_30, OFFSET_123], [None, OFFSET_456], [ROT_30, OFFSET_579]),
        ([None, OFFSET_123], [None, OFFSET_456], [None, OFFSET_579]),
        ([ROT_30, OFFSET_X], [None, None], [ROT_30, OFFSET_X]),
        ([None, None], [ROT_45, OFFSET_Z], [ROT_45, OFFSET_Z]),
        ([None, None], [None, None], [None, None]),
        (
            [ROT_30, OFFSET_X],
            [ROT_45, None],
            [ROT_75, r.CartesianRepresentation([1 / SQRT_2, -1 / SQRT_2, 0])],
        ),
        (
            [ROT_30, OFFSET_X],
            [ROT_45, OFFSET_Z],
            [ROT_75, r.CartesianRepresentation([1 / SQRT_2, -1 / SQRT_2, 1])],
        ),
        (
            [None, OFFSET_123],
            [ROT_45, OFFSET_456],
            [ROT_45, r.CartesianRepresentation([3 / SQRT_2 + 4, 1 / SQRT_2 + 5, 9])],
        ),
    ),
)
def test_combine_affine_params(first, second, check):
    result = t._combine_affine_params(first, second)
    if check[0] is None:
        assert result[0] is None
    else:
        assert_allclose(result[0], check[0])
    if check[1] is None:
        assert result[1] is None
    else:
        assert_allclose(result[1].xyz, check[1].xyz)
