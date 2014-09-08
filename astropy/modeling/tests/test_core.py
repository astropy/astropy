# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import collections
import inspect

import pytest
import numpy as np
from numpy.testing.utils import assert_allclose
from ..core import Model, InputParameterError, custom_model
from ..parameters import Parameter
from .. import models


class NonFittableModel(Model):
    """An example class directly subclassing Model for testing"""

    a = Parameter()

    def __init__(self, a, model_set_axis=None):
        super(NonFittableModel, self).__init__(
            a, model_set_axis=model_set_axis)

    @staticmethod
    def evaluate():
        pass


def test_Model():
    """Some silly tests just to have all lines in the Model code covered by unit tests"""
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
