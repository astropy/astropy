# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import pytest
import numpy as np

from astropy import __minimum_asdf_version__

asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)
from asdf import util
from asdf.tests import helpers

import astropy.units as u
from astropy.modeling import models as astmodels

test_models = [
    astmodels.Identity(2), astmodels.Polynomial1D(2, c0=1, c1=2, c2=3),
    astmodels.Polynomial2D(1, c0_0=1, c0_1=2, c1_0=3), astmodels.Shift(2.),
    astmodels.Scale(3.4), astmodels.RotateNative2Celestial(5.63, -72.5, 180),
    astmodels.Multiply(3), astmodels.Multiply(10*u.m),
    astmodels.RotateCelestial2Native(5.63, -72.5, 180),
    astmodels.EulerAngleRotation(23, 14, 2.3, axes_order='xzx'),
    astmodels.Mapping((0, 1), n_inputs=3),
    astmodels.Shift(2.*u.deg),
    astmodels.Scale(3.4*u.deg),
    astmodels.RotateNative2Celestial(5.63*u.deg, -72.5*u.deg, 180*u.deg),
    astmodels.RotateCelestial2Native(5.63*u.deg, -72.5*u.deg, 180*u.deg),
]


def test_transforms_compound(tmpdir):
    tree = {
        'compound':
            astmodels.Shift(1) & astmodels.Shift(2) |
            astmodels.Sky2Pix_TAN() |
            astmodels.Rotation2D() |
            astmodels.AffineTransformation2D([[2, 0], [0, 2]], [42, 32]) +
            astmodels.Rotation2D(32)
    }

    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_inverse_transforms(tmpdir):
    rotation = astmodels.Rotation2D(32)
    rotation.inverse = astmodels.Rotation2D(45)

    real_rotation = astmodels.Rotation2D(32)

    tree = {
        'rotation': rotation,
        'real_rotation': real_rotation
    }

    def check(ff):
        assert ff.tree['rotation'].inverse.angle == 45

    helpers.assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)


@pytest.mark.parametrize(('model'), test_models)
def test_single_model(tmpdir, model):
    tree = {'single_model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_name(tmpdir):
    def check(ff):
        assert ff.tree['rot'].name == 'foo'

    tree = {'rot': astmodels.Rotation2D(23, name='foo')}
    helpers.assert_roundtrip_tree(tree, tmpdir, asdf_check_func=check)


def test_zenithal_with_arguments(tmpdir):
    tree = {
        'azp': astmodels.Sky2Pix_AZP(0.5, 0.3)
    }

    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_naming_of_compound_model(tmpdir):
    """Issue #87"""
    def asdf_check(ff):
        assert ff.tree['model'].name == 'compound_model'

    offx = astmodels.Shift(1)
    scl = astmodels.Scale(2)
    model = (offx | scl).rename('compound_model')
    tree = {
        'model': model
    }
    helpers.assert_roundtrip_tree(tree, tmpdir, asdf_check_func=asdf_check)


def test_generic_projections(tmpdir):
    from astropy.io.misc.asdf.tags.transform import projections

    for tag_name, (name, params, version) in projections._generic_projections.items():
        tree = {
            'forward': util.resolve_name(
                'astropy.modeling.projections.Sky2Pix_{0}'.format(name))(),
            'backward': util.resolve_name(
                'astropy.modeling.projections.Pix2Sky_{0}'.format(name))()
        }

        helpers.assert_roundtrip_tree(tree, tmpdir)


def test_tabular_model(tmpdir):
    points = np.arange(0, 5)
    values = [1., 10, 2, 45, -3]
    model = astmodels.Tabular1D(points=points, lookup_table=values)
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)
    table = np.array([[ 3.,  0.,  0.],
                      [ 0.,  2.,  0.],
                      [ 0.,  0.,  0.]])
    points = ([1, 2, 3], [1, 2, 3])
    model2 = astmodels.Tabular2D(points, lookup_table=table, bounds_error=False,
                                 fill_value=None, method='nearest')
    tree = {'model': model2}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_bounding_box(tmpdir):
    model = astmodels.Shift(1) & astmodels.Shift(2)
    model.bounding_box = ((1, 3), (2, 4))
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_linear1d(tmpdir):
    model = astmodels.Linear1D()
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_linear1d_quantity(tmpdir):
    model = astmodels.Linear1D(1*u.nm, 1*(u.nm/u.pixel))
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)


def test_tabular_model_units(tmpdir):
    points = np.arange(0, 5) * u.pix
    values = [1., 10, 2, 45, -3] * u.nm
    model = astmodels.Tabular1D(points=points, lookup_table=values)
    tree = {'model': model}
    helpers.assert_roundtrip_tree(tree, tmpdir)
    table = np.array([[3.,  0.,  0.],
                      [0.,  2.,  0.],
                      [0.,  0.,  0.]]) * u.nm
    points = ([1, 2, 3], [1, 2, 3]) * u.pix
    model2 = astmodels.Tabular2D(points, lookup_table=table,
                                 bounds_error=False, fill_value=None,
                                 method='nearest')
    tree = {'model': model2}
    helpers.assert_roundtrip_tree(tree, tmpdir)
