# Licensed under a 3-clause BSD style license - see LICENSE.rst


import pytest
import numpy as np
from numpy.testing import assert_allclose

from ...wcs import wcs
from .. import models
from ... import units as u
from ...tests.helper import assert_quantity_allclose


@pytest.mark.parametrize(('inp'), [(0, 0), (4000, -20.56), (-2001.5, 45.9),
                                   (0, 90), (0, -90), (np.mgrid[:4, :6])])
def test_against_wcslib(inp):
    w = wcs.WCS()
    crval = [202.4823228, 47.17511893]
    w.wcs.crval = crval
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    lonpole = 180
    tan = models.Pix2Sky_TAN()
    n2c = models.RotateNative2Celestial(crval[0] * u.deg, crval[1] * u.deg, lonpole * u.deg)
    c2n = models.RotateCelestial2Native(crval[0] * u.deg, crval[1] * u.deg, lonpole * u.deg)
    m = tan | n2c
    minv = c2n | tan.inverse

    radec = w.wcs_pix2world(inp[0], inp[1], 1)
    xy = w.wcs_world2pix(radec[0], radec[1], 1)

    assert_allclose(m(*inp), radec, atol=1e-12)
    assert_allclose(minv(*radec), xy, atol=1e-12)


@pytest.mark.parametrize(('inp'), [(40 * u.deg, -0.057 * u.rad), (21.5 * u.arcsec, 45.9 * u.deg)])
def test_roundtrip_sky_rotaion(inp):
    lon, lat, lon_pole = 42 * u.deg, (43 * u.deg).to(u.arcsec), (44 * u.deg).to(u.rad)
    n2c = models.RotateNative2Celestial(lon, lat, lon_pole)
    c2n = models.RotateCelestial2Native(lon, lat, lon_pole)
    assert_quantity_allclose(n2c.inverse(*n2c(*inp)), inp, atol=1e-13 * u.deg)
    assert_quantity_allclose(c2n.inverse(*c2n(*inp)), inp, atol=1e-13 * u.deg)


def test_Rotation2D():
    model = models.Rotation2D(angle=90 * u.deg)
    a, b = 1 * u.deg, 0 * u.deg
    x, y = model(a, b)
    assert_quantity_allclose([x, y], [0 * u.deg, 1 * u.deg], atol=1e-10 * u.deg)


def test_Rotation2D_inverse():
    model = models.Rotation2D(angle=234.23494 * u.deg)
    x, y = model.inverse(*model(1 * u.deg, 0 * u.deg))
    assert_quantity_allclose([x, y], [1 * u.deg, 0 * u.deg], atol=1e-10 * u.deg)


def test_euler_angle_rotations():
    ydeg = (90 * u.deg, 0 * u.deg)
    y = (90, 0)
    z = (0, 90)

    # rotate y into minus z
    model = models.EulerAngleRotation(0 * u.rad, np.pi / 2 * u.rad, 0 * u.rad, 'zxz')
    assert_allclose(model(*z), y, atol=10**-12)
    model = models.EulerAngleRotation(0 * u.deg, 90 * u.deg, 0 * u.deg, 'zxz')
    assert_quantity_allclose(model(*(z * u.deg)), ydeg, atol=10**-12 * u.deg)


@pytest.mark.parametrize(('params'), [(60, 10, 25),
                                      (60 * u.deg, 10 * u.deg, 25 * u.deg),
                                      ((60 * u.deg).to(u.rad),
                                       (10 * u.deg).to(u.rad),
                                       (25 * u.deg).to(u.rad))])
def test_euler_rotations_with_units(params):
    x = 1 * u.deg
    y = 1 * u.deg
    phi, theta, psi = params

    urot = models.EulerAngleRotation(phi, theta, psi, axes_order='xyz')
    a, b = urot(x.value, y.value)
    assert_allclose((a, b), (-23.614457631192547, 9.631254579686113))
    a, b = urot(x, y)
    assert_quantity_allclose((a, b), (-23.614457631192547 * u.deg, 9.631254579686113 * u.deg))
    a, b = urot(x.to(u.rad), y.to(u.rad))
    assert_quantity_allclose((a, b), (-23.614457631192547 * u.deg, 9.631254579686113 * u.deg))


def test_attributes():
    n2c = models.RotateNative2Celestial(20016 * u.arcsec, -72.3 * u.deg, np.pi * u.rad)
    assert_allclose(n2c.lat.value, -72.3)
    assert_allclose(n2c.lat._raw_value, -1.2618730491919001)
    assert_allclose(n2c.lon.value, 20016)
    assert_allclose(n2c.lon._raw_value, 0.09704030641088472)
    assert_allclose(n2c.lon_pole.value, np.pi)
    assert_allclose(n2c.lon_pole._raw_value, np.pi)
    assert(n2c.lon.unit is u.Unit("arcsec"))
    assert(n2c._param_metrics['lon']['raw_unit'] is u.Unit("rad"))
    assert(n2c.lat.unit is u.Unit("deg"))
    assert(n2c._param_metrics['lat']['raw_unit'] is u.Unit("rad"))
    assert(n2c.lon_pole.unit is u.Unit("rad"))
    assert(n2c._param_metrics['lon_pole']['raw_unit'] is u.Unit("rad"))
