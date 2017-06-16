# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from math import cos, sin

import numpy as np
from numpy.testing import utils

from .. import models
from ... import units as u
from ...tests.helper import assert_quantity_allclose


def test_Rotation2D():
    model = models.Rotation2D(angle=90*u.deg)
    a, b = 1*u.deg, 0*u.deg
    x, y = model(a, b)
    assert_quantity_allclose([x, y], [0*u.deg, 1*u.deg], atol=1e-10*u.deg)


def test_Rotation2D_inverse():
    model = models.Rotation2D(angle=234.23494*u.deg)
    x, y = model.inverse(*model(1*u.deg, 0*u.deg))
    assert_quantity_allclose([x, y], [1*u.deg, 0*u.deg], atol=1e-10*u.deg)


def test_euler_angle_rotations():
    x = (0 * u.deg, 0 * u.deg)
    ydeg = (90 * u.deg, 0 * u.deg)
    y = (90, 0)
    z = (0, 90)

    # rotate y into minus z
    model = models.EulerAngleRotation(0*u.rad, (90*u.deg).to(u.rad), 0*u.rad, 'zxz')
    utils.assert_allclose(model(*z), y, atol=10**-12)
    model = models.EulerAngleRotation(0*u.deg, 90*u.deg, 0*u.deg, 'zxz')
    assert_quantity_allclose(model(*(z*u.deg)), ydeg, atol=10**-12*u.deg)


def test_euler_rotations_with_units():
    x = 1 * u.deg
    y = 1 * u.deg
    phi = 60 * u.deg
    theta = 10 * u.deg
    psi = 25 * u.deg

    rot = models.EulerAngleRotation(phi.value, theta.value, psi.value, axes_order='xyz')
    utils.assert_allclose(rot(x.value, y.value), (-23.614457631192547, 9.631254579686113))
    assert_quantity_allclose(rot(x.to(u.rad), y.to( u.rad)),
                             (-0.4121500367370107*u.rad, 0.16809710351330523*u.rad))
    assert_quantity_allclose(rot(x, y),
                             (-23.614457631192547*u.deg, 9.631254579686113*u.deg))

    urot = models.EulerAngleRotation(phi, theta, psi, axes_order='xyz')
    a, b = urot(x.value, y.value)
    utils.assert_allclose((a, b), (-23.614457631192547, 9.631254579686113))
    a, b = urot(x, y)
    assert_quantity_allclose((a, b), (-23.614457631192547*u.deg, 9.631254579686113*u.deg))
    a, b = urot(x.to(u.rad), y.to(u.rad))
    assert_quantity_allclose((a, b), (-0.4121500367370107*u.rad, 0.16809710351330523*u.rad))

    rurot = models.EulerAngleRotation(phi.to(u.rad), theta.to(u.rad), psi.to(u.rad), axes_order='xyz')
    a, b = rurot(x.value, y.value)
    utils.assert_allclose((a, b), (-23.614457631192547, 9.631254579686113))
    a, b = rurot(x, y)
    assert_quantity_allclose((a, b), (-23.614457631192547*u.deg, 9.631254579686113*u.deg))
    a, b = rurot(x.to(u.rad), y.to(u.rad))
    assert_quantity_allclose((a, b), (-0.4121500367370107*u.rad, 0.16809710351330523*u.rad))

