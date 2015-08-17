# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import collections
import inspect

import pytest
import numpy as np
from numpy.testing.utils import assert_allclose
from ..core import Model, InputParameterError, custom_model, render_model
from ..parameters import Parameter
from .. import models


class NonFittableModel(Model):
    """An example class directly subclassing Model for testing."""

    a = Parameter()

    def __init__(self, a, model_set_axis=None):
        super(NonFittableModel, self).__init__(
            a, model_set_axis=model_set_axis)

    @staticmethod
    def evaluate():
        pass


def test_Model_instance_repr_and_str():
    m = NonFittableModel(42)
    assert repr(m) == "<NonFittableModel(a=42.0)>"
    assert (str(m) ==
        "Model: NonFittableModel\n"
        "Inputs: ()\n"
        "Outputs: ()\n"
        "Model set size: 1\n"
        "Parameters:\n"
        "     a  \n"
        "    ----\n"
        "    42.0")

    assert len(m) == 1


def test_Model_array_parameter():
    m = NonFittableModel([[42, 43], [1,2]])
    m = NonFittableModel([42, 43, 44, 45])

    phi, theta, psi = 42, 43, 44
    model = models.RotateNative2Celestial(phi, theta, psi)
    assert_allclose(model.param_sets, [[42], [43], [44]])


def test_inputless_model():
    """
    Regression test for
    https://github.com/astropy/astropy/pull/3772#issuecomment-101821641
    """

    class TestModel(Model):
        inputs = ()
        outputs = ('y',)
        a = Parameter()

        @staticmethod
        def evaluate(a):
            return a

    m = TestModel(1)
    assert m.a == 1
    assert m() == 1

    # Test array-like output
    m = TestModel([1, 2, 3], model_set_axis=False)
    assert len(m) == 1
    assert np.all(m() == [1, 2, 3])

    # Test a model set
    m = TestModel(a=[1, 2, 3], model_set_axis=0)
    assert len(m) == 3
    assert np.all(m() == [1, 2, 3])

    # Test a model set
    m = TestModel(a=[[1, 2, 3], [4, 5, 6]], model_set_axis=0)
    assert len(m) == 2
    assert np.all(m() == [[1, 2, 3], [4, 5, 6]])


def test_Model_add_model():
    m = models.Gaussian1D(1,2,3)
    m.add_model(m, 'p')
    m.add_model(m, 's')
    with pytest.raises(InputParameterError):
        m.add_model(m, 'q')
        m.add_model(m, 42)


def test_ParametericModel():
    with pytest.raises(TypeError):
        models.Gaussian1D(1, 2, 3, wrong=4)


def test_custom_model_signature():
    """
    Tests that the signatures for the __init__ and __call__
    methods of custom models are useful.
    """

    @custom_model
    def model_a(x):
        return x

    assert model_a.param_names == ()
    assert model_a.n_inputs == 1
    argspec = inspect.getargspec(model_a.__init__)
    assert argspec.args == ['self']
    argspec = inspect.getargspec(model_a.__call__)
    assert argspec.args == ['self', 'x', 'model_set_axis']

    @custom_model
    def model_b(x, a=1, b=2):
        return x + a + b

    assert model_b.param_names == ('a', 'b')
    assert model_b.n_inputs == 1
    argspec = inspect.getargspec(model_b.__init__)
    assert argspec.args == ['self', 'a', 'b']
    assert argspec.defaults == (1, 2)
    argspec = inspect.getargspec(model_b.__call__)
    assert argspec.args == ['self', 'x', 'model_set_axis']

    @custom_model
    def model_c(x, y, a=1, b=2):
        return x + y + a + b

    assert model_c.param_names == ('a', 'b')
    assert model_c.n_inputs == 2
    argspec = inspect.getargspec(model_c.__init__)
    assert argspec.args == ['self', 'a', 'b']
    assert argspec.defaults == (1, 2)
    argspec = inspect.getargspec(model_c.__call__)
    assert argspec.args == ['self', 'x', 'y', 'model_set_axis']


def test_custom_model_subclass():
    """Test that custom models can be subclassed."""

    @custom_model
    def model_a(x, a=1):
        return x * a

    class model_b(model_a):
        # Override the evaluate from model_a
        @classmethod
        def evaluate(cls, x, a):
            return -super(model_b, cls).evaluate(x, a)

    b = model_b()
    assert b.param_names == ('a',)
    assert b.a == 1
    assert b(1) == -1

    argspec = inspect.getargspec(model_b.__init__)
    assert argspec.args == ['self', 'a']
    argspec = inspect.getargspec(model_b.__call__)
    assert argspec.args == ['self', 'x', 'model_set_axis']


def test_custom_model_parametrized_decorator():
    """Tests using custom_model as a decorator with parameters."""

    def cosine(x, amplitude=1):
        return [amplitude * np.cos(x)]

    @custom_model(fit_deriv=cosine)
    def sine(x, amplitude=1):
        return amplitude * np.sin(x)

    assert issubclass(sine, Model)
    s = sine(2)
    assert_allclose(s(np.pi / 2), 2)
    assert_allclose(s.fit_deriv(0, 2), 2)


def test_custom_inverse():
    """Test setting a custom inverse on a model."""

    p = models.Polynomial1D(1, c0=-2, c1=3)
    # A trivial inverse for a trivial polynomial
    inv = models.Polynomial1D(1, c0=(2./3.), c1=(1./3.))

    with pytest.raises(NotImplementedError):
        p.inverse

    p.inverse = inv

    x = np.arange(100)

    assert_allclose(x, p(p.inverse(x)))
    assert_allclose(x, p.inverse(p(x)))

    p.inverse = None

    with pytest.raises(NotImplementedError):
        p.inverse


def test_render_model_2d():

    imshape = (71, 141)
    image = np.zeros(imshape)
    coords = y, x = np.indices(imshape)

    model = models.Gaussian2D(x_stddev=6.1, y_stddev=3.9, theta=np.pi / 4)

    # test points for edges
    ye, xe = [0, 35, 70], [0, 70, 140]
    # test points for floating point positions
    yf, xf = [35.1, 35.5, 35.9], [70.1, 70.5, 70.9]

    test_pts = [(a, b) for a in xe for b in ye] + [(a, b) for a in xf for b in yf]

    for x0, y0 in test_pts:
        model.x_mean = x0
        model.y_mean = y0
        expected = model(x, y)
        for im in [image, None]:
            for xy in [coords, None]:
                if (im is None) & (xy is None):
                    # this case is tested in Fittable2DModelTester
                    continue
                actual = render_model(model, arr=image, coords=xy)
                # assert images match
                assert_allclose(expected, actual, atol=2e-7)
                # assert flux conserved
                assert ((np.sum(expected) - np.sum(actual)) / np.sum(expected)) < 1e-7


def test_render_model_1d():

    npix = 101
    image = np.zeros(npix)
    coords = np.arange(npix)

    model = models.Gaussian1D(stddev=49.5)

    # test points
    test_pts = [0, 49.1, 49.5, 49.9, 100]

    # test widths
    test_stdv = np.arange(5.5, 6.7, .2)

    for x0, stdv in zip(test_pts, test_stdv):
        model.mean = x0
        model.stddev = stdv
        expected = model(coords)
        for im in [image, None]:
            for x in [coords, None]:
                if (im is None) & (x is None):
                    # this case is tested in Fittable1DModelTester
                    continue
                actual = render_model(model, arr=image, coords=x)
                # assert images match
                assert_allclose(expected, actual, atol=2e-7)
                # assert flux conserved
                assert ((np.sum(expected) - np.sum(actual)) / np.sum(expected)) < 1e-7


def test_render_model_3d():
    imshape = (17, 21, 27)
    image = np.zeros(imshape)
    coords = np.indices(imshape)

    def ellipsoid(x, y, z, x0=13., y0=10., z0=8., a=4., b=3., c=2., amp=1.):
        rsq = ((x - x0) / a) ** 2 + ((y - y0) / b) ** 2 + ((z - z0) / c) ** 2
        val = (rsq < 1) * amp
        return val

    Ellipsoid3D = models.custom_model(ellipsoid)

    model = Ellipsoid3D()

    # test points for edges
    ze, ye, xe = [0, 8, 16], [0, 10, 20], [0, 13, 26]
    # test points for floating point positions
    zf, yf, xf = [8.1, 8.5, 8.9], [10.1, 10.5, 10.9], [13.1, 13.5, 13.9]

    test_pts = [(x, y, z) for x in xe for y in ye for z in ze]
    test_pts += [(x, y, z) for x in xf for y in yf for z in zf]

    for x0, y0, z0 in [(8,10,13)]:#test_pts:
        model.x0 = x0
        model.y0 = y0
        model.z0 = z0
        expected = model(*coords[::-1])
        for im in [image, None]:
            for c in [coords, None]:
                if (im is None) & (c is None):
                    continue
                actual = render_model(model, arr=image, coords=c)
                # assert images match
                assert_allclose(expected, actual)
                # assert flux conserved
                assert ((np.sum(expected) - np.sum(actual)) / np.sum(expected)) == 0
