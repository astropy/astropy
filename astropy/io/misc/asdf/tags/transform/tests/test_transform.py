# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import warnings
from distutils.version import LooseVersion

import pytest
import numpy as np

from astropy import __minimum_asdf_version__

asdf = pytest.importorskip('asdf', minversion=__minimum_asdf_version__)
from asdf import util
from asdf.tests import helpers
from asdf import AsdfFile
import asdf

import astropy.units as u
from astropy.modeling.core import fix_inputs
from astropy.modeling import models as astmodels


def custom_and_analytical_inverse():
    p1 = astmodels.Polynomial1D(1)
    p2 = astmodels.Polynomial1D(1)
    p3 = astmodels.Polynomial1D(1)
    p4 = astmodels.Polynomial1D(1)
    m1 = p1 & p2
    m2 = p3 & p4
    m1.inverse = m2
    return m1


test_models = [
    astmodels.Identity(2), astmodels.Polynomial1D(2, c0=1, c1=2, c2=3),
    astmodels.Polynomial2D(1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Shift(2.),
    astmodels.Hermite1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Legendre1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Chebyshev1D(2, c0=2, c1=3, c2=0.5),
    astmodels.Chebyshev2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Legendre2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Hermite2D(1, 1, c0_0=1, c0_1=2, c1_0=3),
    astmodels.Scale(3.4), astmodels.RotateNative2Celestial(5.63, -72.5, 180),
    astmodels.Multiply(3), astmodels.Multiply(10*u.m),
    astmodels.RotateCelestial2Native(5.63, -72.5, 180),
    astmodels.EulerAngleRotation(23, 14, 2.3, axes_order='xzx'),
    astmodels.Mapping((0, 1), n_inputs=3),
    astmodels.Shift(2.*u.deg),
    astmodels.Scale(3.4*u.deg),
    astmodels.RotateNative2Celestial(5.63*u.deg, -72.5*u.deg, 180*u.deg),
    astmodels.RotateCelestial2Native(5.63*u.deg, -72.5*u.deg, 180*u.deg),
    astmodels.RotationSequence3D([1.2, 2.3, 3.4, .3], 'xyzx'),
    astmodels.SphericalRotationSequence([1.2, 2.3, 3.4, .3], 'xyzy'),
    custom_and_analytical_inverse(),
]


math_models = []

for kl in astmodels.math.__all__:
    klass = getattr(astmodels.math, kl)
    math_models.append(klass())

test_models.extend(math_models)


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
    with warnings.catch_warnings():
        # Some schema files are missing from asdf<=2.4.2 which causes warnings
        if LooseVersion(asdf.__version__) <= '2.4.2':
            warnings.filterwarnings('ignore', 'Unable to locate schema file')
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
                f'astropy.modeling.projections.Sky2Pix_{name}')(),
            'backward': util.resolve_name(
                f'astropy.modeling.projections.Pix2Sky_{name}')()
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


@pytest.mark.parametrize("standard_version", asdf.versioning.supported_versions)
@pytest.mark.parametrize("model", [
    astmodels.Polynomial1D(1, c0=5, c1=17),
    astmodels.Polynomial1D(1, c0=5, c1=17, domain=[-5, 4], window=[-2, 3]),
    astmodels.Polynomial2D(2, c0_0=3, c1_0=5, c0_1=7),
    astmodels.Polynomial2D(
        2, c0_0=3, c1_0=5, c0_1=7, x_domain=[-2, 2], y_domain=[-4, 4],
        x_window=[-6, 6], y_window=[-8, 8]
    ),
])
def test_polynomial(tmpdir, standard_version, model):
    helpers.assert_roundtrip_tree({"model": model}, tmpdir, init_options={"version": standard_version})


def test_domain_orthopoly(tmpdir):
    model1d = astmodels.Chebyshev1D(2, c0=2, c1=3, c2=0.5, domain=[-2, 2])
    model2d = astmodels.Chebyshev2D(1, 1, c0_0=1, c0_1=2, c1_0=3,
                                    x_domain=[-2, 2], y_domain=[-2, 2])
    fa = AsdfFile()
    fa.tree['model1d'] = model1d
    fa.tree['model2d'] = model2d

    file_path = str(tmpdir.join('orthopoly_domain.asdf'))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree['model1d'](1.8) == model1d(1.8)
        assert f.tree['model2d'](1.8, -1.5) == model2d(1.8, -1.5)


def test_window_orthopoly(tmpdir):
    model1d = astmodels.Chebyshev1D(2, c0=2, c1=3, c2=0.5,
                                    domain=[-2, 2], window=[-0.5, 0.5])
    model2d = astmodels.Chebyshev2D(1, 1, c0_0=1, c0_1=2, c1_0=3,
                                    x_domain=[-2, 2], y_domain=[-2, 2],
                                    x_window=[-0.5, 0.5], y_window=[-0.1, 0.5])
    fa = AsdfFile()
    fa.tree['model1d'] = model1d
    fa.tree['model2d'] = model2d

    file_path = str(tmpdir.join('orthopoly_window.asdf'))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree['model1d'](1.8) == model1d(1.8)
        assert f.tree['model2d'](1.8, -1.5) == model2d(1.8, -1.5)


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


def test_fix_inputs(tmpdir):

    with warnings.catch_warnings():
        # Some schema files are missing from asdf<=2.4.2 which causes warnings
        if LooseVersion(asdf.__version__) <= '2.5.1':
            warnings.filterwarnings('ignore', 'Unable to locate schema file')

        model = astmodels.Pix2Sky_TAN() | astmodels.Rotation2D()
        tree = {
            'compound': fix_inputs(model, {'x': 45}),
            'compound1': fix_inputs(model, {0: 45})
        }

        helpers.assert_roundtrip_tree(tree, tmpdir)


def test_fix_inputs_type():
    with pytest.raises(TypeError):
        tree = {
        'compound': fix_inputs(3, {'x': 45})
        }
        helpers.assert_roundtrip_tree(tree, tmpdir)

    with pytest.raises(AttributeError):
        tree = {
        'compound': astmodels.Pix2Sky_TAN() & {'x': 45}
        }
        helpers.assert_roundtrip_tree(tree, tmpdir)


comp_model = custom_and_analytical_inverse()


@pytest.mark.parametrize(('model'), [astmodels.Shift(1) & astmodels.Shift(2) | comp_model,
                                     comp_model | astmodels.Shift(1) & astmodels.Shift(2),
                                     astmodels.Shift(1) & comp_model,
                                     comp_model & astmodels.Shift(1)
                                     ])
def test_custom_and_analytical(model, tmpdir):
    fa = AsdfFile()
    fa.tree['model'] = model
    file_path = str(tmpdir.join('custom_and_analytical_inverse.asdf'))
    fa.write_to(file_path)
    with asdf.open(file_path) as f:
        assert f.tree['model'].inverse is not None
