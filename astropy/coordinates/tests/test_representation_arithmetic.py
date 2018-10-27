# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools

import pytest
import numpy as np

from ... import units as u
from .. import (PhysicsSphericalRepresentation, CartesianRepresentation,
                CylindricalRepresentation, SphericalRepresentation,
                UnitSphericalRepresentation, SphericalDifferential,
                CartesianDifferential, UnitSphericalDifferential,
                SphericalCosLatDifferential, UnitSphericalCosLatDifferential,
                PhysicsSphericalDifferential, CylindricalDifferential,
                RadialRepresentation, RadialDifferential, Longitude, Latitude)
from ..representation import DIFFERENTIAL_CLASSES
from ..angle_utilities import angular_separation
from ...tests.helper import assert_quantity_allclose, quantity_allclose


def assert_representation_allclose(actual, desired, rtol=1.e-7, atol=None,
                                   **kwargs):
    actual_xyz = actual.to_cartesian().get_xyz(xyz_axis=-1)
    desired_xyz = desired.to_cartesian().get_xyz(xyz_axis=-1)
    actual_xyz, desired_xyz = np.broadcast_arrays(actual_xyz, desired_xyz,
                                                  subok=True)
    assert_quantity_allclose(actual_xyz, desired_xyz, rtol, atol, **kwargs)


def assert_differential_allclose(actual, desired, rtol=1.e-7, **kwargs):
    assert actual.components == desired.components
    for component in actual.components:
        actual_c = getattr(actual, component)
        atol = 1.e-10 * actual_c.unit
        assert_quantity_allclose(actual_c, getattr(desired, component),
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
        assert quantity_allclose(s4.to_cartesian().xyz,
                                 -self.spherical.to_cartesian().xyz,
                                 atol=1e-15*self.spherical.distance.unit)
        assert np.all(s4.distance == self.spherical.distance)
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


class TestUnitVectorsAndScales():

    @staticmethod
    def check_unit_vectors(e):
        for v in e.values():
            assert type(v) is CartesianRepresentation
            assert_quantity_allclose(v.norm(), 1. * u.one)
        return e

    @staticmethod
    def check_scale_factors(sf, rep):
        unit = rep.norm().unit
        for c, f in sf.items():
            assert type(f) is u.Quantity
            assert (f.unit * getattr(rep, c).unit).is_equivalent(unit)

    def test_spherical(self):
        s = SphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                    lat=[0., -30., 85.] * u.deg,
                                    distance=[1, 2, 3] * u.kpc)
        e = s.unit_vectors()
        self.check_unit_vectors(e)
        sf = s.scale_factors()
        self.check_scale_factors(sf, s)

        s_lon = s + s.distance * 1e-5 * np.cos(s.lat) * e['lon']
        assert_quantity_allclose(s_lon.lon, s.lon + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        assert_quantity_allclose(s_lon.lat, s.lat, atol=1e-10*u.rad)
        assert_quantity_allclose(s_lon.distance, s.distance)
        s_lon2 = s + 1e-5 * u.radian * sf['lon'] * e['lon']
        assert_representation_allclose(s_lon2, s_lon)

        s_lat = s + s.distance * 1e-5 * e['lat']
        assert_quantity_allclose(s_lat.lon, s.lon)
        assert_quantity_allclose(s_lat.lat, s.lat + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        assert_quantity_allclose(s_lon.distance, s.distance)
        s_lat2 = s + 1.e-5 * u.radian * sf['lat'] * e['lat']
        assert_representation_allclose(s_lat2, s_lat)

        s_distance = s + 1. * u.pc * e['distance']
        assert_quantity_allclose(s_distance.lon, s.lon, atol=1e-10*u.rad)
        assert_quantity_allclose(s_distance.lat, s.lat, atol=1e-10*u.rad)
        assert_quantity_allclose(s_distance.distance, s.distance + 1.*u.pc)
        s_distance2 = s + 1. * u.pc * sf['distance'] * e['distance']
        assert_representation_allclose(s_distance2, s_distance)

    def test_unit_spherical(self):
        s = UnitSphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                        lat=[0., -30., 85.] * u.deg)

        e = s.unit_vectors()
        self.check_unit_vectors(e)
        sf = s.scale_factors()
        self.check_scale_factors(sf, s)

        s_lon = s + 1e-5 * np.cos(s.lat) * e['lon']
        assert_quantity_allclose(s_lon.lon, s.lon + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        assert_quantity_allclose(s_lon.lat, s.lat, atol=1e-10*u.rad)
        s_lon2 = s + 1e-5 * u.radian * sf['lon'] * e['lon']
        assert_representation_allclose(s_lon2, s_lon)

        s_lat = s + 1e-5 * e['lat']
        assert_quantity_allclose(s_lat.lon, s.lon)
        assert_quantity_allclose(s_lat.lat, s.lat + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        s_lat2 = s + 1.e-5 * u.radian * sf['lat'] * e['lat']
        assert_representation_allclose(s_lat2, s_lat)

    def test_radial(self):
        r = RadialRepresentation(10.*u.kpc)
        with pytest.raises(NotImplementedError):
            r.unit_vectors()
        sf = r.scale_factors()
        assert np.all(sf['distance'] == 1.*u.one)
        assert np.all(r.norm() == r.distance)
        with pytest.raises(TypeError):
            r + r

    def test_physical_spherical(self):

        s = PhysicsSphericalRepresentation(phi=[0., 6., 21.] * u.hourangle,
                                           theta=[90., 120., 5.] * u.deg,
                                           r=[1, 2, 3] * u.kpc)

        e = s.unit_vectors()
        self.check_unit_vectors(e)
        sf = s.scale_factors()
        self.check_scale_factors(sf, s)

        s_phi = s + s.r * 1e-5 * np.sin(s.theta) * e['phi']
        assert_quantity_allclose(s_phi.phi, s.phi + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        assert_quantity_allclose(s_phi.theta, s.theta, atol=1e-10*u.rad)
        assert_quantity_allclose(s_phi.r, s.r)
        s_phi2 = s + 1e-5 * u.radian * sf['phi'] * e['phi']
        assert_representation_allclose(s_phi2, s_phi)

        s_theta = s + s.r * 1e-5 * e['theta']
        assert_quantity_allclose(s_theta.phi, s.phi)
        assert_quantity_allclose(s_theta.theta, s.theta + 1e-5*u.rad,
                                 atol=1e-10*u.rad)
        assert_quantity_allclose(s_theta.r, s.r)
        s_theta2 = s + 1.e-5 * u.radian * sf['theta'] * e['theta']
        assert_representation_allclose(s_theta2, s_theta)

        s_r = s + 1. * u.pc * e['r']
        assert_quantity_allclose(s_r.phi, s.phi, atol=1e-10*u.rad)
        assert_quantity_allclose(s_r.theta, s.theta, atol=1e-10*u.rad)
        assert_quantity_allclose(s_r.r, s.r + 1.*u.pc)
        s_r2 = s + 1. * u.pc * sf['r'] * e['r']
        assert_representation_allclose(s_r2, s_r)

    def test_cartesian(self):

        s = CartesianRepresentation(x=[1, 2, 3] * u.pc,
                                    y=[2, 3, 4] * u.Mpc,
                                    z=[3, 4, 5] * u.kpc)

        e = s.unit_vectors()
        sf = s.scale_factors()
        for v, expected in zip(e.values(), ([1., 0., 0.] * u.one,
                                            [0., 1., 0.] * u.one,
                                            [0., 0., 1.] * u.one)):
            assert np.all(v.get_xyz(xyz_axis=-1) == expected)
        for f in sf.values():
            assert np.all(f == 1.*u.one)

    def test_cylindrical(self):

        s = CylindricalRepresentation(rho=[1, 2, 3] * u.pc,
                                      phi=[0., 90., -45.] * u.deg,
                                      z=[3, 4, 5] * u.kpc)
        e = s.unit_vectors()
        self.check_unit_vectors(e)
        sf = s.scale_factors()
        self.check_scale_factors(sf, s)

        s_rho = s + 1. * u.pc * e['rho']
        assert_quantity_allclose(s_rho.rho, s.rho + 1.*u.pc)
        assert_quantity_allclose(s_rho.phi, s.phi)
        assert_quantity_allclose(s_rho.z, s.z)
        s_rho2 = s + 1. * u.pc * sf['rho'] * e['rho']
        assert_representation_allclose(s_rho2, s_rho)

        s_phi = s + s.rho * 1e-5 * e['phi']
        assert_quantity_allclose(s_phi.rho, s.rho)
        assert_quantity_allclose(s_phi.phi, s.phi + 1e-5*u.rad)
        assert_quantity_allclose(s_phi.z, s.z)
        s_phi2 = s + 1e-5 * u.radian * sf['phi'] * e['phi']
        assert_representation_allclose(s_phi2, s_phi)

        s_z = s + 1. * u.pc * e['z']
        assert_quantity_allclose(s_z.rho, s.rho)
        assert_quantity_allclose(s_z.phi, s.phi, atol=1e-10*u.rad)
        assert_quantity_allclose(s_z.z, s.z + 1.*u.pc)
        s_z2 = s + 1. * u.pc * sf['z'] * e['z']
        assert_representation_allclose(s_z2, s_z)


@pytest.mark.parametrize('omit_coslat', [False, True], scope='class')
class TestSphericalDifferential():
    # these test cases are subclassed for SphericalCosLatDifferential,
    # hence some tests depend on omit_coslat.

    def _setup(self, omit_coslat):
        if omit_coslat:
            self.SD_cls = SphericalCosLatDifferential
        else:
            self.SD_cls = SphericalDifferential

        s = SphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                    lat=[0., -30., 85.] * u.deg,
                                    distance=[1, 2, 3] * u.kpc)
        self.s = s
        self.e = s.unit_vectors()
        self.sf = s.scale_factors(omit_coslat=omit_coslat)

    def test_name_coslat(self, omit_coslat):
        self._setup(omit_coslat)
        if omit_coslat:
            assert self.SD_cls is SphericalCosLatDifferential
            assert self.SD_cls.get_name() == 'sphericalcoslat'
        else:
            assert self.SD_cls is SphericalDifferential
            assert self.SD_cls.get_name() == 'spherical'
        assert self.SD_cls.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self, omit_coslat):
        self._setup(omit_coslat)
        s, e, sf = self.s, self.e, self.sf

        o_lon = self.SD_cls(1.*u.arcsec, 0.*u.arcsec, 0.*u.kpc)
        o_lonc = o_lon.to_cartesian(base=s)
        o_lon2 = self.SD_cls.from_cartesian(o_lonc, base=s)
        assert_differential_allclose(o_lon, o_lon2)
        # simple check by hand for first element.
        # lat[0] is 0, so cos(lat) term doesn't matter.
        assert_quantity_allclose(o_lonc[0].xyz,
                                 [0., np.pi/180./3600., 0.]*u.kpc)
        # check all using unit vectors and scale factors.
        s_lon = s + 1.*u.arcsec * sf['lon'] * e['lon']
        assert_representation_allclose(o_lonc, s_lon - s, atol=1*u.npc)
        s_lon2 = s + o_lon
        assert_representation_allclose(s_lon2, s_lon, atol=1*u.npc)

        o_lat = self.SD_cls(0.*u.arcsec, 1.*u.arcsec, 0.*u.kpc)
        o_latc = o_lat.to_cartesian(base=s)
        assert_quantity_allclose(o_latc[0].xyz,
                                 [0., 0., np.pi/180./3600.]*u.kpc,
                                 atol=1.*u.npc)
        s_lat = s + 1.*u.arcsec * sf['lat'] * e['lat']
        assert_representation_allclose(o_latc, s_lat - s, atol=1*u.npc)
        s_lat2 = s + o_lat
        assert_representation_allclose(s_lat2, s_lat, atol=1*u.npc)

        o_distance = self.SD_cls(0.*u.arcsec, 0.*u.arcsec, 1.*u.mpc)
        o_distancec = o_distance.to_cartesian(base=s)
        assert_quantity_allclose(o_distancec[0].xyz,
                                 [1e-6, 0., 0.]*u.kpc, atol=1.*u.npc)
        s_distance = s + 1.*u.mpc * sf['distance'] * e['distance']
        assert_representation_allclose(o_distancec, s_distance - s,
                                       atol=1*u.npc)
        s_distance2 = s + o_distance
        assert_representation_allclose(s_distance2, s_distance)

    def test_differential_arithmetic(self, omit_coslat):
        self._setup(omit_coslat)
        s = self.s

        o_lon = self.SD_cls(1.*u.arcsec, 0.*u.arcsec, 0.*u.kpc)
        o_lon_by_2 = o_lon / 2.
        assert_representation_allclose(o_lon_by_2.to_cartesian(s) * 2.,
                                       o_lon.to_cartesian(s), atol=1e-10*u.kpc)
        assert_representation_allclose(s + o_lon, s + 2 * o_lon_by_2,
                                       atol=1e-10*u.kpc)
        o_lon_rec = o_lon_by_2 + o_lon_by_2
        assert_representation_allclose(s + o_lon, s + o_lon_rec,
                                       atol=1e-10*u.kpc)
        o_lon_0 = o_lon - o_lon
        for c in o_lon_0.components:
            assert np.all(getattr(o_lon_0, c) == 0.)
        o_lon2 = self.SD_cls(1*u.mas/u.yr, 0*u.mas/u.yr, 0*u.km/u.s)
        assert_quantity_allclose(o_lon2.norm(s)[0], 4.74*u.km/u.s,
                                 atol=0.01*u.km/u.s)
        assert_representation_allclose(o_lon2.to_cartesian(s) * 1000.*u.yr,
                                       o_lon.to_cartesian(s), atol=1e-10*u.kpc)
        s_off = s + o_lon
        s_off2 = s + o_lon2 * 1000.*u.yr
        assert_representation_allclose(s_off, s_off2, atol=1e-10*u.kpc)

        factor = 1e5 * u.radian/u.arcsec
        if not omit_coslat:
            factor = factor / np.cos(s.lat)
        s_off_big = s + o_lon * factor

        assert_representation_allclose(
            s_off_big, SphericalRepresentation(s.lon + 90.*u.deg, 0.*u.deg,
                                               1e5*s.distance),
            atol=5.*u.kpc)

        o_lon3c = CartesianRepresentation(0., 4.74047, 0., unit=u.km/u.s)
        o_lon3 = self.SD_cls.from_cartesian(o_lon3c, base=s)
        expected0 = self.SD_cls(1.*u.mas/u.yr, 0.*u.mas/u.yr, 0.*u.km/u.s)
        assert_differential_allclose(o_lon3[0], expected0)
        s_off_big2 = s + o_lon3 * 1e5 * u.yr * u.radian/u.mas
        assert_representation_allclose(
            s_off_big2, SphericalRepresentation(90.*u.deg, 0.*u.deg,
                                                1e5*u.kpc), atol=5.*u.kpc)

        with pytest.raises(TypeError):
            o_lon - s
        with pytest.raises(TypeError):
            s.to_cartesian() + o_lon

    def test_differential_init_errors(self, omit_coslat):
        self._setup(omit_coslat)
        s = self.s
        with pytest.raises(u.UnitsError):
            self.SD_cls(1.*u.arcsec, 0., 0.)
        with pytest.raises(TypeError):
            self.SD_cls(1.*u.arcsec, 0.*u.arcsec, 0.*u.kpc,
                                  False, False)
        with pytest.raises(TypeError):
            self.SD_cls(1.*u.arcsec, 0.*u.arcsec, 0.*u.kpc,
                            copy=False, d_lat=0.*u.arcsec)
        with pytest.raises(TypeError):
            self.SD_cls(1.*u.arcsec, 0.*u.arcsec, 0.*u.kpc,
                            copy=False, flying='circus')
        with pytest.raises(ValueError):
            self.SD_cls(np.ones(2)*u.arcsec,
                            np.zeros(3)*u.arcsec, np.zeros(2)*u.kpc)
        with pytest.raises(u.UnitsError):
            self.SD_cls(1.*u.arcsec, 1.*u.s, 0.*u.kpc)
        with pytest.raises(u.UnitsError):
            self.SD_cls(1.*u.kpc, 1.*u.arcsec, 0.*u.kpc)
        o = self.SD_cls(1.*u.arcsec, 1.*u.arcsec, 0.*u.km/u.s)
        with pytest.raises(u.UnitsError):
            o.to_cartesian(s)
        with pytest.raises(AttributeError):
            o.d_lat = 0.*u.arcsec
        with pytest.raises(AttributeError):
            del o.d_lat

        o = self.SD_cls(1.*u.arcsec, 1.*u.arcsec, 0.*u.km)
        with pytest.raises(TypeError):
            o.to_cartesian()
        c = CartesianRepresentation(10., 0., 0., unit=u.km)
        with pytest.raises(TypeError):
            self.SD_cls.to_cartesian(c)
        with pytest.raises(TypeError):
            self.SD_cls.from_cartesian(c)
        with pytest.raises(TypeError):
            self.SD_cls.from_cartesian(c, SphericalRepresentation)
        with pytest.raises(TypeError):
            self.SD_cls.from_cartesian(c, c)


@pytest.mark.parametrize('omit_coslat', [False, True], scope='class')
class TestUnitSphericalDifferential():
    def _setup(self, omit_coslat):
        if omit_coslat:
            self.USD_cls = UnitSphericalCosLatDifferential
        else:
            self.USD_cls = UnitSphericalDifferential

        s = UnitSphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                        lat=[0., -30., 85.] * u.deg)
        self.s = s
        self.e = s.unit_vectors()
        self.sf = s.scale_factors(omit_coslat=omit_coslat)

    def test_name_coslat(self, omit_coslat):
        self._setup(omit_coslat)
        if omit_coslat:
            assert self.USD_cls is UnitSphericalCosLatDifferential
            assert self.USD_cls.get_name() == 'unitsphericalcoslat'
        else:
            assert self.USD_cls is UnitSphericalDifferential
            assert self.USD_cls.get_name() == 'unitspherical'
        assert self.USD_cls.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self, omit_coslat):
        self._setup(omit_coslat)
        s, e, sf = self.s, self.e, self.sf

        o_lon = self.USD_cls(1.*u.arcsec, 0.*u.arcsec)
        o_lonc = o_lon.to_cartesian(base=s)
        o_lon2 = self.USD_cls.from_cartesian(o_lonc, base=s)
        assert_differential_allclose(o_lon, o_lon2)
        # simple check by hand for first element
        # (lat[0]=0, so works for both normal and CosLat differential)
        assert_quantity_allclose(o_lonc[0].xyz,
                                 [0., np.pi/180./3600., 0.]*u.one)
        # check all using unit vectors and scale factors.
        s_lon = s + 1.*u.arcsec * sf['lon'] * e['lon']
        assert type(s_lon) is SphericalRepresentation
        assert_representation_allclose(o_lonc, s_lon - s, atol=1e-10*u.one)
        s_lon2 = s + o_lon
        assert_representation_allclose(s_lon2, s_lon, atol=1e-10*u.one)

        o_lat = self.USD_cls(0.*u.arcsec, 1.*u.arcsec)
        o_latc = o_lat.to_cartesian(base=s)
        assert_quantity_allclose(o_latc[0].xyz,
                                 [0., 0., np.pi/180./3600.]*u.one,
                                 atol=1e-10*u.one)
        s_lat = s + 1.*u.arcsec * sf['lat'] * e['lat']
        assert type(s_lat) is SphericalRepresentation
        assert_representation_allclose(o_latc, s_lat - s, atol=1e-10*u.one)
        s_lat2 = s + o_lat
        assert_representation_allclose(s_lat2, s_lat, atol=1e-10*u.one)

    def test_differential_arithmetic(self, omit_coslat):
        self._setup(omit_coslat)
        s = self.s

        o_lon = self.USD_cls(1.*u.arcsec, 0.*u.arcsec)
        o_lon_by_2 = o_lon / 2.
        assert type(o_lon_by_2) is self.USD_cls
        assert_representation_allclose(o_lon_by_2.to_cartesian(s) * 2.,
                                       o_lon.to_cartesian(s), atol=1e-10*u.one)
        s_lon = s + o_lon
        s_lon2 = s + 2 * o_lon_by_2
        assert type(s_lon) is SphericalRepresentation
        assert_representation_allclose(s_lon, s_lon2, atol=1e-10*u.one)
        o_lon_rec = o_lon_by_2 + o_lon_by_2
        assert type(o_lon_rec) is self.USD_cls
        assert representation_equal(o_lon, o_lon_rec)
        assert_representation_allclose(s + o_lon, s + o_lon_rec,
                                       atol=1e-10*u.one)
        o_lon_0 = o_lon - o_lon
        assert type(o_lon_0) is self.USD_cls
        for c in o_lon_0.components:
            assert np.all(getattr(o_lon_0, c) == 0.)

        o_lon2 = self.USD_cls(1.*u.mas/u.yr, 0.*u.mas/u.yr)
        kks = u.km/u.kpc/u.s
        assert_quantity_allclose(o_lon2.norm(s)[0], 4.74047*kks, atol=1e-4*kks)
        assert_representation_allclose(o_lon2.to_cartesian(s) * 1000.*u.yr,
                                       o_lon.to_cartesian(s), atol=1e-10*u.one)
        s_off = s + o_lon
        s_off2 = s + o_lon2 * 1000.*u.yr
        assert_representation_allclose(s_off, s_off2, atol=1e-10*u.one)

        factor = 1e5 * u.radian/u.arcsec
        if not omit_coslat:
            factor = factor / np.cos(s.lat)
        s_off_big = s + o_lon * factor

        assert_representation_allclose(
            s_off_big, SphericalRepresentation(s.lon + 90.*u.deg,
                                               0.*u.deg, 1e5),
            atol=5.*u.one)

        o_lon3c = CartesianRepresentation(0., 4.74047, 0., unit=kks)
        # This looses information!!
        o_lon3 = self.USD_cls.from_cartesian(o_lon3c, base=s)
        expected0 = self.USD_cls(1.*u.mas/u.yr, 0.*u.mas/u.yr)
        assert_differential_allclose(o_lon3[0], expected0)
        # Part of motion kept.
        part_kept = s.cross(CartesianRepresentation(0, 1, 0, unit=u.one)).norm()
        assert_quantity_allclose(o_lon3.norm(s), 4.74047*part_kept*kks,
                                 atol=1e-10*kks)
        # (lat[0]=0, so works for both normal and CosLat differential)
        s_off_big2 = s + o_lon3 * 1e5 * u.yr * u.radian/u.mas
        expected0 = SphericalRepresentation(90.*u.deg, 0.*u.deg,
                                            1e5*u.one)
        assert_representation_allclose(s_off_big2[0], expected0, atol=5.*u.one)

    def test_differential_init_errors(self, omit_coslat):
        self._setup(omit_coslat)
        with pytest.raises(u.UnitsError):
            self.USD_cls(0.*u.deg, 10.*u.deg/u.yr)


class TestRadialDifferential():
    def setup(self):
        s = SphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                    lat=[0., -30., 85.] * u.deg,
                                    distance=[1, 2, 3] * u.kpc)
        self.s = s
        self.r = s.represent_as(RadialRepresentation)
        self.e = s.unit_vectors()
        self.sf = s.scale_factors()

    def test_name(self):
        assert RadialDifferential.get_name() == 'radial'
        assert RadialDifferential.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self):
        r, s, e, sf = self.r, self.s, self.e, self.sf

        o_distance = RadialDifferential(1.*u.mpc)
        # Can be applied to RadialRepresentation, though not most useful.
        r_distance = r + o_distance
        assert_quantity_allclose(r_distance.distance,
                                 r.distance + o_distance.d_distance)
        r_distance2 = o_distance + r
        assert_quantity_allclose(r_distance2.distance,
                                 r.distance + o_distance.d_distance)
        # More sense to apply it relative to spherical representation.
        o_distancec = o_distance.to_cartesian(base=s)
        assert_quantity_allclose(o_distancec[0].xyz,
                                 [1e-6, 0., 0.]*u.kpc, atol=1.*u.npc)
        o_recover = RadialDifferential.from_cartesian(o_distancec, base=s)
        assert_quantity_allclose(o_recover.d_distance, o_distance.d_distance)

        s_distance = s + 1.*u.mpc * sf['distance'] * e['distance']
        assert_representation_allclose(o_distancec, s_distance - s,
                                       atol=1*u.npc)
        s_distance2 = s + o_distance
        assert_representation_allclose(s_distance2, s_distance)


class TestPhysicsSphericalDifferential():
    """Test copied from SphericalDifferential, so less extensive."""

    def setup(self):
        s = PhysicsSphericalRepresentation(phi=[0., 90., 315.] * u.deg,
                                           theta=[90., 120., 5.] * u.deg,
                                           r=[1, 2, 3] * u.kpc)
        self.s = s
        self.e = s.unit_vectors()
        self.sf = s.scale_factors()

    def test_name(self):
        assert PhysicsSphericalDifferential.get_name() == 'physicsspherical'
        assert PhysicsSphericalDifferential.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self):
        s, e, sf = self.s, self.e, self.sf

        o_phi = PhysicsSphericalDifferential(1*u.arcsec, 0*u.arcsec, 0*u.kpc)
        o_phic = o_phi.to_cartesian(base=s)
        o_phi2 = PhysicsSphericalDifferential.from_cartesian(o_phic, base=s)
        assert_quantity_allclose(o_phi.d_phi, o_phi2.d_phi, atol=1.*u.narcsec)
        assert_quantity_allclose(o_phi.d_theta, o_phi2.d_theta,
                                 atol=1.*u.narcsec)
        assert_quantity_allclose(o_phi.d_r, o_phi2.d_r, atol=1.*u.npc)
        # simple check by hand for first element.
        assert_quantity_allclose(o_phic[0].xyz,
                                 [0., np.pi/180./3600., 0.]*u.kpc,
                                 atol=1.*u.npc)
        # check all using unit vectors and scale factors.
        s_phi = s + 1.*u.arcsec * sf['phi'] * e['phi']
        assert_representation_allclose(o_phic, s_phi - s, atol=1e-10*u.kpc)

        o_theta = PhysicsSphericalDifferential(0*u.arcsec, 1*u.arcsec, 0*u.kpc)
        o_thetac = o_theta.to_cartesian(base=s)
        assert_quantity_allclose(o_thetac[0].xyz,
                                 [0., 0., -np.pi/180./3600.]*u.kpc,
                                 atol=1.*u.npc)
        s_theta = s + 1.*u.arcsec * sf['theta'] * e['theta']
        assert_representation_allclose(o_thetac, s_theta - s, atol=1e-10*u.kpc)
        s_theta2 = s + o_theta
        assert_representation_allclose(s_theta2, s_theta, atol=1e-10*u.kpc)

        o_r = PhysicsSphericalDifferential(0*u.arcsec, 0*u.arcsec, 1*u.mpc)
        o_rc = o_r.to_cartesian(base=s)
        assert_quantity_allclose(o_rc[0].xyz, [1e-6, 0., 0.]*u.kpc,
                                 atol=1.*u.npc)
        s_r = s + 1.*u.mpc * sf['r'] * e['r']
        assert_representation_allclose(o_rc, s_r - s, atol=1e-10*u.kpc)
        s_r2 = s + o_r
        assert_representation_allclose(s_r2, s_r)

    def test_differential_init_errors(self):
        with pytest.raises(u.UnitsError):
            PhysicsSphericalDifferential(1.*u.arcsec, 0., 0.)


class TestCylindricalDifferential():
    """Test copied from SphericalDifferential, so less extensive."""

    def setup(self):
        s = CylindricalRepresentation(rho=[1, 2, 3] * u.kpc,
                                      phi=[0., 90., 315.] * u.deg,
                                      z=[3, 2, 1] * u.kpc)
        self.s = s
        self.e = s.unit_vectors()
        self.sf = s.scale_factors()

    def test_name(self):
        assert CylindricalDifferential.get_name() == 'cylindrical'
        assert CylindricalDifferential.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self):
        s, e, sf = self.s, self.e, self.sf

        o_rho = CylindricalDifferential(1.*u.mpc, 0.*u.arcsec, 0.*u.kpc)
        o_rhoc = o_rho.to_cartesian(base=s)
        assert_quantity_allclose(o_rhoc[0].xyz, [1.e-6, 0., 0.]*u.kpc)
        s_rho = s + 1.*u.mpc * sf['rho'] * e['rho']
        assert_representation_allclose(o_rhoc, s_rho - s, atol=1e-10*u.kpc)
        s_rho2 = s + o_rho
        assert_representation_allclose(s_rho2, s_rho)

        o_phi = CylindricalDifferential(0.*u.kpc, 1.*u.arcsec, 0.*u.kpc)
        o_phic = o_phi.to_cartesian(base=s)
        o_phi2 = CylindricalDifferential.from_cartesian(o_phic, base=s)
        assert_quantity_allclose(o_phi.d_rho, o_phi2.d_rho, atol=1.*u.npc)
        assert_quantity_allclose(o_phi.d_phi, o_phi2.d_phi, atol=1.*u.narcsec)
        assert_quantity_allclose(o_phi.d_z, o_phi2.d_z, atol=1.*u.npc)
        # simple check by hand for first element.
        assert_quantity_allclose(o_phic[0].xyz,
                                 [0., np.pi/180./3600., 0.]*u.kpc)
        # check all using unit vectors and scale factors.
        s_phi = s + 1.*u.arcsec * sf['phi'] * e['phi']
        assert_representation_allclose(o_phic, s_phi - s, atol=1e-10*u.kpc)

        o_z = CylindricalDifferential(0.*u.kpc, 0.*u.arcsec, 1.*u.mpc)
        o_zc = o_z.to_cartesian(base=s)
        assert_quantity_allclose(o_zc[0].xyz, [0., 0., 1.e-6]*u.kpc)
        s_z = s + 1.*u.mpc * sf['z'] * e['z']
        assert_representation_allclose(o_zc, s_z - s, atol=1e-10*u.kpc)
        s_z2 = s + o_z
        assert_representation_allclose(s_z2, s_z)

    def test_differential_init_errors(self):
        with pytest.raises(u.UnitsError):
            CylindricalDifferential(1.*u.pc, 1.*u.arcsec, 3.*u.km/u.s)


class TestCartesianDifferential():
    """Test copied from SphericalDifferential, so less extensive."""

    def setup(self):
        s = CartesianRepresentation(x=[1, 2, 3] * u.kpc,
                                    y=[2, 3, 1] * u.kpc,
                                    z=[3, 1, 2] * u.kpc)
        self.s = s
        self.e = s.unit_vectors()
        self.sf = s.scale_factors()

    def test_name(self):
        assert CartesianDifferential.get_name() == 'cartesian'
        assert CartesianDifferential.get_name() in DIFFERENTIAL_CLASSES

    def test_simple_differentials(self):
        s, e, sf = self.s, self.e, self.sf

        for d, differential in (  # test different inits while we're at it.
                ('x', CartesianDifferential(1.*u.pc, 0.*u.pc, 0.*u.pc)),
                ('y', CartesianDifferential([0., 1., 0.], unit=u.pc)),
                ('z', CartesianDifferential(np.array([[0., 0., 1.]]) * u.pc,
                                            xyz_axis=1))):
            o_c = differential.to_cartesian(base=s)
            o_c2 = differential.to_cartesian()
            assert np.all(representation_equal(o_c, o_c2))
            assert all(np.all(getattr(differential, 'd_'+c) == getattr(o_c, c))
                       for c in ('x', 'y', 'z'))
            differential2 = CartesianDifferential.from_cartesian(o_c)
            assert np.all(representation_equal(differential2, differential))
            differential3 = CartesianDifferential.from_cartesian(o_c, base=o_c)
            assert np.all(representation_equal(differential3, differential))

            s_off = s + 1.*u.pc * sf[d] * e[d]
            assert_representation_allclose(o_c, s_off - s, atol=1e-10*u.kpc)
            s_off2 = s + differential
            assert_representation_allclose(s_off2, s_off)

    def test_init_failures(self):
        with pytest.raises(ValueError):
            CartesianDifferential(1.*u.kpc/u.s, 2.*u.kpc)
        with pytest.raises(u.UnitsError):
            CartesianDifferential(1.*u.kpc/u.s, 2.*u.kpc, 3.*u.kpc)
        with pytest.raises(ValueError):
            CartesianDifferential(1.*u.kpc, 2.*u.kpc, 3.*u.kpc, xyz_axis=1)


class TestDifferentialConversion():
    def setup(self):
        self.s = SphericalRepresentation(lon=[0., 6., 21.] * u.hourangle,
                                         lat=[0., -30., 85.] * u.deg,
                                         distance=[1, 2, 3] * u.kpc)

    @pytest.mark.parametrize('sd_cls', [SphericalDifferential,
                                        SphericalCosLatDifferential])
    def test_represent_as_own_class(self, sd_cls):
        so = sd_cls(1.*u.deg, 2.*u.deg, 0.1*u.kpc)
        so2 = so.represent_as(sd_cls)
        assert so2 is so

    def test_represent_other_coslat(self):
        s = self.s
        coslat = np.cos(s.lat)
        so = SphericalDifferential(1.*u.deg, 2.*u.deg, 0.1*u.kpc)
        so_coslat = so.represent_as(SphericalCosLatDifferential, base=s)
        assert_quantity_allclose(so.d_lon * coslat,
                                 so_coslat.d_lon_coslat)
        so2 = so_coslat.represent_as(SphericalDifferential, base=s)
        assert np.all(representation_equal(so2, so))
        so3 = SphericalDifferential.from_representation(so_coslat, base=s)
        assert np.all(representation_equal(so3, so))
        so_coslat2 = SphericalCosLatDifferential.from_representation(so, base=s)
        assert np.all(representation_equal(so_coslat2, so_coslat))
        # Also test UnitSpherical
        us = s.represent_as(UnitSphericalRepresentation)
        uo = so.represent_as(UnitSphericalDifferential)
        uo_coslat = so.represent_as(UnitSphericalCosLatDifferential, base=s)
        assert_quantity_allclose(uo.d_lon * coslat,
                                 uo_coslat.d_lon_coslat)
        uo2 = uo_coslat.represent_as(UnitSphericalDifferential, base=us)
        assert np.all(representation_equal(uo2, uo))
        uo3 = UnitSphericalDifferential.from_representation(uo_coslat, base=us)
        assert np.all(representation_equal(uo3, uo))
        uo_coslat2 = UnitSphericalCosLatDifferential.from_representation(
            uo, base=us)
        assert np.all(representation_equal(uo_coslat2, uo_coslat))
        uo_coslat3 = uo.represent_as(UnitSphericalCosLatDifferential, base=us)
        assert np.all(representation_equal(uo_coslat3, uo_coslat))

    @pytest.mark.parametrize('sd_cls', [SphericalDifferential,
                                        SphericalCosLatDifferential])
    @pytest.mark.parametrize('r_cls', (SphericalRepresentation,
                                       UnitSphericalRepresentation,
                                       PhysicsSphericalRepresentation,
                                       CylindricalRepresentation))
    def test_represent_regular_class(self, sd_cls, r_cls):
        so = sd_cls(1.*u.deg, 2.*u.deg, 0.1*u.kpc)
        r = so.represent_as(r_cls, base=self.s)
        c = so.to_cartesian(self.s)
        r_check = c.represent_as(r_cls)
        assert np.all(representation_equal(r, r_check))
        so2 = sd_cls.from_representation(r, base=self.s)
        so3 = sd_cls.from_cartesian(r.to_cartesian(), self.s)
        assert np.all(representation_equal(so2, so3))

    @pytest.mark.parametrize('sd_cls', [SphericalDifferential,
                                        SphericalCosLatDifferential])
    def test_convert_physics(self, sd_cls):
        # Conversion needs no base for SphericalDifferential, but does
        # need one (to get the latitude) for SphericalCosLatDifferential.
        if sd_cls is SphericalDifferential:
            usd_cls = UnitSphericalDifferential
            base_s = base_u = base_p = None
        else:
            usd_cls = UnitSphericalCosLatDifferential
            base_s = self.s[1]
            base_u = base_s.represent_as(UnitSphericalRepresentation)
            base_p = base_s.represent_as(PhysicsSphericalRepresentation)

        so = sd_cls(1.*u.deg, 2.*u.deg, 0.1*u.kpc)
        po = so.represent_as(PhysicsSphericalDifferential, base=base_s)
        so2 = sd_cls.from_representation(po, base=base_s)
        assert_differential_allclose(so, so2)
        po2 = PhysicsSphericalDifferential.from_representation(so, base=base_p)
        assert_differential_allclose(po, po2)
        so3 = po.represent_as(sd_cls, base=base_p)
        assert_differential_allclose(so, so3)

        s = self.s
        p = s.represent_as(PhysicsSphericalRepresentation)
        cso = so.to_cartesian(s[1])
        cpo = po.to_cartesian(p[1])
        assert_representation_allclose(cso, cpo)
        assert_representation_allclose(s[1] + so, p[1] + po)
        po2 = so.represent_as(PhysicsSphericalDifferential,
                              base=None if base_s is None else s)
        assert_representation_allclose(s + so, p + po2)

        suo = usd_cls.from_representation(so)
        puo = usd_cls.from_representation(po, base=base_u)
        assert_differential_allclose(suo, puo)
        suo2 = so.represent_as(usd_cls)
        puo2 = po.represent_as(usd_cls, base=base_p)
        assert_differential_allclose(suo2, puo2)
        assert_differential_allclose(puo, puo2)

        sro = RadialDifferential.from_representation(so)
        pro = RadialDifferential.from_representation(po)
        assert representation_equal(sro, pro)
        sro2 = so.represent_as(RadialDifferential)
        pro2 = po.represent_as(RadialDifferential)
        assert representation_equal(sro2, pro2)
        assert representation_equal(pro, pro2)

    @pytest.mark.parametrize(
        ('sd_cls', 'usd_cls'),
        [(SphericalDifferential, UnitSphericalDifferential),
         (SphericalCosLatDifferential, UnitSphericalCosLatDifferential)])
    def test_convert_unit_spherical_radial(self, sd_cls, usd_cls):
        s = self.s
        us = s.represent_as(UnitSphericalRepresentation)
        rs = s.represent_as(RadialRepresentation)
        assert_representation_allclose(rs * us, s)

        uo = usd_cls(2.*u.deg, 1.*u.deg)
        so = uo.represent_as(sd_cls, base=s)
        assert_quantity_allclose(so.d_distance, 0.*u.kpc, atol=1.*u.npc)
        uo2 = so.represent_as(usd_cls)
        assert_representation_allclose(uo.to_cartesian(us),
                                       uo2.to_cartesian(us))
        so1 = sd_cls(2.*u.deg, 1.*u.deg, 5.*u.pc)
        uo_r = so1.represent_as(usd_cls)
        ro_r = so1.represent_as(RadialDifferential)
        assert np.all(representation_equal(uo_r, uo))
        assert np.all(representation_equal(ro_r, RadialDifferential(5.*u.pc)))

    @pytest.mark.parametrize('sd_cls', [SphericalDifferential,
                                        SphericalCosLatDifferential])
    def test_convert_cylindrial(self, sd_cls):
        s = self.s
        so = sd_cls(1.*u.deg, 2.*u.deg, 0.1*u.kpc)
        cyo = so.represent_as(CylindricalDifferential, base=s)
        cy = s.represent_as(CylindricalRepresentation)
        so1 = cyo.represent_as(sd_cls, base=cy)
        assert_representation_allclose(so.to_cartesian(s),
                                       so1.to_cartesian(s))
        cyo2 = CylindricalDifferential.from_representation(so, base=cy)
        assert_representation_allclose(cyo2.to_cartesian(base=cy),
                                       cyo.to_cartesian(base=cy))
        so2 = sd_cls.from_representation(cyo2, base=s)
        assert_representation_allclose(so.to_cartesian(s),
                                       so2.to_cartesian(s))

    @pytest.mark.parametrize('sd_cls', [SphericalDifferential,
                                        SphericalCosLatDifferential])
    def test_combinations(self, sd_cls):
        if sd_cls is SphericalDifferential:
            uo = UnitSphericalDifferential(2.*u.deg, 1.*u.deg)
            uo_d_lon = uo.d_lon
        else:
            uo = UnitSphericalCosLatDifferential(2.*u.deg, 1.*u.deg)
            uo_d_lon = uo.d_lon_coslat
        ro = RadialDifferential(1.*u.mpc)
        so1 = uo + ro
        so1c = sd_cls(uo_d_lon, uo.d_lat, ro.d_distance)
        assert np.all(representation_equal(so1, so1c))

        so2 = uo - ro
        so2c = sd_cls(uo_d_lon, uo.d_lat, -ro.d_distance)
        assert np.all(representation_equal(so2, so2c))
        so3 = so2 + ro
        so3c = sd_cls(uo_d_lon, uo.d_lat, 0.*u.kpc)
        assert np.all(representation_equal(so3, so3c))
        so4 = so1 + ro
        so4c = sd_cls(uo_d_lon, uo.d_lat, 2*ro.d_distance)
        assert np.all(representation_equal(so4, so4c))
        so5 = so1 - uo
        so5c = sd_cls(0*u.deg, 0.*u.deg, ro.d_distance)
        assert np.all(representation_equal(so5, so5c))
        assert_representation_allclose(self.s + (uo+ro), self.s+so1)


@pytest.mark.parametrize('rep,dif', [
    [CartesianRepresentation([1, 2, 3]*u.kpc),
     CartesianDifferential([.1, .2, .3]*u.km/u.s)],
    [SphericalRepresentation(90*u.deg, 0.*u.deg, 14.*u.kpc),
     SphericalDifferential(1.*u.deg, 2.*u.deg, 0.1*u.kpc)]
])
def test_arithmetic_with_differentials_fail(rep, dif):

    rep = rep.with_differentials(dif)

    with pytest.raises(TypeError):
        rep + rep

    with pytest.raises(TypeError):
        rep - rep

    with pytest.raises(TypeError):
        rep * rep

    with pytest.raises(TypeError):
        rep / rep

    with pytest.raises(TypeError):
        10. * rep

    with pytest.raises(TypeError):
        rep / 10.

    with pytest.raises(TypeError):
        -rep
