# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS
import functools
import numpy as np

from ... import units as u
from .. import (PhysicsSphericalRepresentation, CartesianRepresentation,
                CylindricalRepresentation, SphericalRepresentation,
                UnitSphericalRepresentation, Longitude, Latitude)
from ..angle_utilities import angular_separation
from ...tests.helper import pytest, assert_quantity_allclose


def assert_representation_allclose(actual, desired, rtol=1.e-7, atol=None,
                                   **kwargs):
    assert_quantity_allclose(actual.to_cartesian().xyz,
                             desired.to_cartesian().xyz,
                             rtol, atol, **kwargs)


def representation_equal(first, second):
    return functools.reduce(np.logical_and,
                            (getattr(first, component) ==
                             getattr(second, component)
                             for component in first.components))

class TestArithmetic():

    def setup(self):
        # Choose some specific coordinates, for which ``sum`` and ``dot``
        # works out nicely.
        self.lon = Longitude(np.arange(0, 12.1, 2), u.hourangle)
        self.lat = Latitude(np.arange(-90, 91, 30), u.deg)
        self.distance = [5., 12., 4., 2., 4., 12., 5.] * u.kpc
        self.spherical = SphericalRepresentation(self.lon, self.lat,
                                                 self.distance)
        self.unit_spherical = self.spherical.represent_as(
            UnitSphericalRepresentation)
        self.cartesian = self.spherical.to_cartesian()

    def test_norm_spherical(self):
        norm_s = self.spherical.norm()
        assert isinstance(norm_s, u.Quantity)
        # Just to be sure, test against getting object arrays.
        assert norm_s.dtype.kind == 'f'
        assert np.all(norm_s == self.distance)

    @pytest.mark.parametrize('representation',
                             (PhysicsSphericalRepresentation,
                              CartesianRepresentation,
                              CylindricalRepresentation))
    def test_norm(self, representation):
        in_rep = self.spherical.represent_as(representation)
        norm_rep = in_rep.norm()
        assert isinstance(norm_rep, u.Quantity)
        assert_quantity_allclose(norm_rep, self.distance)

    def test_norm_unitspherical(self):
        norm_rep = self.unit_spherical.norm()
        assert norm_rep.unit == u.dimensionless_unscaled
        assert np.all(norm_rep == 1. * u.dimensionless_unscaled)

    @pytest.mark.parametrize('representation',
                             (SphericalRepresentation,
                              PhysicsSphericalRepresentation,
                              CartesianRepresentation,
                              CylindricalRepresentation,
                              UnitSphericalRepresentation))
    def test_neg_pos(self, representation):
        in_rep = self.cartesian.represent_as(representation)
        pos_rep = +in_rep
        assert type(pos_rep) is type(in_rep)
        assert pos_rep is not in_rep
        assert np.all(representation_equal(pos_rep, in_rep))
        neg_rep = -in_rep
        assert type(neg_rep) is type(in_rep)
        assert np.all(neg_rep.norm() == in_rep.norm())
        in_rep_xyz = in_rep.to_cartesian().xyz
        assert_quantity_allclose(neg_rep.to_cartesian().xyz,
                                 -in_rep_xyz, atol=1.e-10*in_rep_xyz.unit)

    def test_mul_div_spherical(self):
        s0 = self.spherical / (1. * u.Myr)
        assert isinstance(s0, SphericalRepresentation)
        assert s0.distance.dtype.kind == 'f'
        assert np.all(s0.lon == self.spherical.lon)
        assert np.all(s0.lat == self.spherical.lat)
        assert np.all(s0.distance == self.distance / (1. * u.Myr))
        s1 = (1./u.Myr) * self.spherical
        assert isinstance(s1, SphericalRepresentation)
        assert np.all(representation_equal(s1, s0))
        s2 = self.spherical * np.array([[1.], [2.]])
        assert isinstance(s2, SphericalRepresentation)
        assert s2.shape == (2, self.spherical.shape[0])
        assert np.all(s2.lon == self.spherical.lon)
        assert np.all(s2.lat == self.spherical.lat)
        assert np.all(s2.distance ==
                      self.spherical.distance * np.array([[1.], [2.]]))
        s3 = np.array([[1.], [2.]]) * self.spherical
        assert isinstance(s3, SphericalRepresentation)
        assert np.all(representation_equal(s3, s2))
        s4 = -self.spherical
        assert isinstance(s4, SphericalRepresentation)
        assert np.all(s4.lon == self.spherical.lon)
        assert np.all(s4.lat == self.spherical.lat)
        assert np.all(s4.distance == -self.spherical.distance)
        s5 = +self.spherical
        assert s5 is not self.spherical
        assert np.all(representation_equal(s5, self.spherical))

    @pytest.mark.parametrize('representation',
                             (PhysicsSphericalRepresentation,
                              CartesianRepresentation,
                              CylindricalRepresentation))
    def test_mul_div(self, representation):
        in_rep = self.spherical.represent_as(representation)
        r1 = in_rep / (1. * u.Myr)
        assert isinstance(r1, representation)
        for component in in_rep.components:
            in_rep_comp = getattr(in_rep, component)
            r1_comp = getattr(r1, component)
            if in_rep_comp.unit == self.distance.unit:
                assert np.all(r1_comp == in_rep_comp / (1.*u.Myr))
            else:
                assert np.all(r1_comp == in_rep_comp)

        r2 = np.array([[1.], [2.]]) * in_rep
        assert isinstance(r2, representation)
        assert r2.shape == (2, in_rep.shape[0])
        assert_quantity_allclose(r2.norm(),
                                 self.distance * np.array([[1.], [2.]]))
        r3 = -in_rep
        assert np.all(representation_equal(r3, in_rep * -1.))
        with pytest.raises(TypeError):
            in_rep * in_rep
        with pytest.raises(TypeError):
            dict() * in_rep

    def test_mul_div_unit_spherical(self):
        s1 = self.unit_spherical * self.distance
        assert isinstance(s1, SphericalRepresentation)
        assert np.all(s1.lon == self.unit_spherical.lon)
        assert np.all(s1.lat == self.unit_spherical.lat)
        assert np.all(s1.distance == self.spherical.distance)
        s2 = self.unit_spherical / u.s
        assert isinstance(s2, SphericalRepresentation)
        assert np.all(s2.lon == self.unit_spherical.lon)
        assert np.all(s2.lat == self.unit_spherical.lat)
        assert np.all(s2.distance == 1./u.s)
        u3 = -self.unit_spherical
        assert isinstance(u3, UnitSphericalRepresentation)
        assert_quantity_allclose(u3.lon, self.unit_spherical.lon + 180.*u.deg)
        assert np.all(u3.lat == -self.unit_spherical.lat)
        assert_quantity_allclose(u3.to_cartesian().xyz,
                                 -self.unit_spherical.to_cartesian().xyz,
                                 atol=1.e-10*u.dimensionless_unscaled)
        u4 = +self.unit_spherical
        assert isinstance(u4, UnitSphericalRepresentation)
        assert u4 is not self.unit_spherical
        assert np.all(representation_equal(u4, self.unit_spherical))

    def test_add_sub_cartesian(self):
        c1 = self.cartesian + self.cartesian
        assert isinstance(c1, CartesianRepresentation)
        assert c1.x.dtype.kind == 'f'
        assert np.all(representation_equal(c1, 2. * self.cartesian))
        with pytest.raises(TypeError):
            self.cartesian + 10.*u.m
        with pytest.raises(u.UnitsError):
            self.cartesian + (self.cartesian / u.s)
        c2 = self.cartesian - self.cartesian
        assert isinstance(c2, CartesianRepresentation)
        assert np.all(representation_equal(
            c2, CartesianRepresentation(0.*u.m, 0.*u.m, 0.*u.m)))
        c3 = self.cartesian - self.cartesian / 2.
        assert isinstance(c3, CartesianRepresentation)
        assert np.all(representation_equal(c3, self.cartesian / 2.))

    @pytest.mark.parametrize('representation',
                             (PhysicsSphericalRepresentation,
                              SphericalRepresentation,
                              CylindricalRepresentation))
    def test_add_sub(self, representation):
        in_rep = self.cartesian.represent_as(representation)
        r1 = in_rep + in_rep
        assert isinstance(r1, representation)
        expected = 2. * in_rep
        for component in in_rep.components:
            assert_quantity_allclose(getattr(r1, component),
                                     getattr(expected, component))
        with pytest.raises(TypeError):
            10.*u.m + in_rep
        with pytest.raises(u.UnitsError):
            in_rep + (in_rep / u.s)
        r2 = in_rep - in_rep
        assert isinstance(r2, representation)
        assert np.all(representation_equal(
            r2.to_cartesian(), CartesianRepresentation(0.*u.m, 0.*u.m, 0.*u.m)))
        r3 = in_rep - in_rep / 2.
        assert isinstance(r3, representation)
        expected = in_rep / 2.
        assert_representation_allclose(r3, expected)

    def test_add_sub_unit_spherical(self):
        s1 = self.unit_spherical + self.unit_spherical
        assert isinstance(s1, SphericalRepresentation)
        expected = 2. * self.unit_spherical
        for component in s1.components:
            assert_quantity_allclose(getattr(s1, component),
                                     getattr(expected, component))
        with pytest.raises(TypeError):
            10.*u.m - self.unit_spherical
        with pytest.raises(u.UnitsError):
            self.unit_spherical + (self.unit_spherical / u.s)
        s2 = self.unit_spherical - self.unit_spherical / 2.
        assert isinstance(s2, SphericalRepresentation)
        expected = self.unit_spherical / 2.
        for component in s2.components:
            assert_quantity_allclose(getattr(s2, component),
                                     getattr(expected, component))

    @pytest.mark.parametrize('representation',
                             (CartesianRepresentation,
                              PhysicsSphericalRepresentation,
                              SphericalRepresentation,
                              CylindricalRepresentation))
    def test_sum_mean(self, representation):
        in_rep = self.spherical.represent_as(representation)
        r_sum = in_rep.sum()
        assert isinstance(r_sum, representation)
        expected = SphericalRepresentation(
            90. * u.deg, 0. * u.deg, 14. * u.kpc).represent_as(representation)
        for component in expected.components:
            exp_component = getattr(expected, component)
            assert_quantity_allclose(getattr(r_sum, component),
                                     exp_component,
                                     atol=1e-10*exp_component.unit)

        r_mean = in_rep.mean()
        assert isinstance(r_mean, representation)
        expected = expected / len(in_rep)
        for component in expected.components:
            exp_component = getattr(expected, component)
            assert_quantity_allclose(getattr(r_mean, component),
                                     exp_component,
                                     atol=1e-10*exp_component.unit)

    def test_sum_mean_unit_spherical(self):
        s_sum = self.unit_spherical.sum()
        assert isinstance(s_sum, SphericalRepresentation)
        expected = SphericalRepresentation(
            90. * u.deg, 0. * u.deg, 3. * u.dimensionless_unscaled)
        for component in expected.components:
            exp_component = getattr(expected, component)
            assert_quantity_allclose(getattr(s_sum, component),
                                     exp_component,
                                     atol=1e-10*exp_component.unit)

        s_mean = self.unit_spherical.mean()
        assert isinstance(s_mean, SphericalRepresentation)
        expected = expected / len(self.unit_spherical)
        for component in expected.components:
            exp_component = getattr(expected, component)
            assert_quantity_allclose(getattr(s_mean, component),
                                     exp_component,
                                     atol=1e-10*exp_component.unit)

    @pytest.mark.parametrize('representation',
                             (CartesianRepresentation,
                              PhysicsSphericalRepresentation,
                              SphericalRepresentation,
                              CylindricalRepresentation))
    def test_dot(self, representation):
        in_rep = self.cartesian.represent_as(representation)
        r_dot_r = in_rep.dot(in_rep)
        assert isinstance(r_dot_r, u.Quantity)
        assert r_dot_r.shape == in_rep.shape
        assert_quantity_allclose(np.sqrt(r_dot_r), self.distance)
        r_dot_r_rev = in_rep.dot(in_rep[::-1])
        assert isinstance(r_dot_r_rev, u.Quantity)
        assert r_dot_r_rev.shape == in_rep.shape
        expected = [-25., -126., 2., 4., 2., -126., -25.] * u.kpc**2
        assert_quantity_allclose(r_dot_r_rev, expected)
        for axis in 'xyz':
            project = CartesianRepresentation(*(
                (1. if axis == _axis else 0.) * u.dimensionless_unscaled
                for _axis in 'xyz'))
            assert_quantity_allclose(in_rep.dot(project),
                                     getattr(self.cartesian, axis),
                                     atol=1.*u.upc)
        with pytest.raises(TypeError):
            in_rep.dot(self.cartesian.xyz)

    def test_dot_unit_spherical(self):
        u_dot_u = self.unit_spherical.dot(self.unit_spherical)
        assert isinstance(u_dot_u, u.Quantity)
        assert u_dot_u.shape == self.unit_spherical.shape
        assert_quantity_allclose(u_dot_u, 1.*u.dimensionless_unscaled)
        cartesian = self.unit_spherical.to_cartesian()
        for axis in 'xyz':
            project = CartesianRepresentation(*(
                (1. if axis == _axis else 0.) * u.dimensionless_unscaled
                for _axis in 'xyz'))
            assert_quantity_allclose(self.unit_spherical.dot(project),
                                     getattr(cartesian, axis), atol=1.e-10)

    @pytest.mark.parametrize('representation',
                             (CartesianRepresentation,
                              PhysicsSphericalRepresentation,
                              SphericalRepresentation,
                              CylindricalRepresentation))
    def test_cross(self, representation):
        in_rep = self.cartesian.represent_as(representation)
        r_cross_r = in_rep.cross(in_rep)
        assert isinstance(r_cross_r, representation)
        assert_quantity_allclose(r_cross_r.norm(), 0.*u.kpc**2,
                                 atol=1.*u.mpc**2)
        r_cross_r_rev = in_rep.cross(in_rep[::-1])
        sep = angular_separation(self.lon, self.lat,
                                 self.lon[::-1], self.lat[::-1])
        expected = self.distance * self.distance[::-1] * np.sin(sep)
        assert_quantity_allclose(r_cross_r_rev.norm(), expected,
                                 atol=1.*u.mpc**2)
        unit_vectors = CartesianRepresentation(
            [1., 0., 0.]*u.one,
            [0., 1., 0.]*u.one,
            [0., 0., 1.]*u.one)[:, np.newaxis]
        r_cross_uv = in_rep.cross(unit_vectors)
        assert r_cross_uv.shape == (3, 7)
        assert_quantity_allclose(r_cross_uv.dot(unit_vectors), 0.*u.kpc,
                                 atol=1.*u.upc)
        assert_quantity_allclose(r_cross_uv.dot(in_rep), 0.*u.kpc**2,
                                 atol=1.*u.mpc**2)
        zeros = np.zeros(len(in_rep)) * u.kpc
        expected = CartesianRepresentation(
            u.Quantity((zeros, -self.cartesian.z, self.cartesian.y)),
            u.Quantity((self.cartesian.z, zeros, -self.cartesian.x)),
            u.Quantity((-self.cartesian.y, self.cartesian.x, zeros)))
        # Comparison with spherical is hard since some distances are zero,
        # implying the angles are undefined.
        r_cross_uv_cartesian = r_cross_uv.to_cartesian()
        assert_representation_allclose(r_cross_uv_cartesian,
                                       expected, atol=1.*u.upc)
        # A final check, with the side benefit of ensuring __div__ and norm
        # work on multi-D representations.
        r_cross_uv_by_distance = r_cross_uv / self.distance
        uv_sph = unit_vectors.represent_as(UnitSphericalRepresentation)
        sep = angular_separation(self.lon, self.lat, uv_sph.lon, uv_sph.lat)
        assert_quantity_allclose(r_cross_uv_by_distance.norm(), np.sin(sep),
                                 atol=1e-9)

        with pytest.raises(TypeError):
            in_rep.cross(self.cartesian.xyz)

    def test_cross_unit_spherical(self):
        u_cross_u = self.unit_spherical.cross(self.unit_spherical)
        assert isinstance(u_cross_u, SphericalRepresentation)
        assert_quantity_allclose(u_cross_u.norm(), 0.*u.one, atol=1.e-10*u.one)
        u_cross_u_rev = self.unit_spherical.cross(self.unit_spherical[::-1])
        assert isinstance(u_cross_u_rev, SphericalRepresentation)
        sep = angular_separation(self.lon, self.lat,
                                 self.lon[::-1], self.lat[::-1])
        expected = np.sin(sep)
        assert_quantity_allclose(u_cross_u_rev.norm(), expected,
                                 atol=1.e-10*u.one)
