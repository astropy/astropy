# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import division
import collections
import pytest
import numpy as np
from numpy.testing.utils import assert_allclose
from ..core import Model
from ..parameters import Parameter
from .. import models
from ..utils import InputParameterError
 
 
class NonFittableModel(Model):
    """An example class directly subclassing Model for testing"""
    param_names = ['a']
 
    def __init__(self, a):
        if not isinstance(a, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(a)
        self._a = Parameter('a', a, self, param_dim)
        super(NonFittableModel, self).__init__(self.param_names, n_inputs=1, n_outputs=1,
                                         param_dim=param_dim)
 
    def __call__(self):
        pass


def test_Model():
    """Some silly tests just to have all lines in the Model code covered by unit tests"""
    m = NonFittableModel(42)
    assert repr(m) == 'NonFittableModel(\n            a=42.0,\n            )'
    assert str(m) == '\n        Model: NonFittableModel\n        Parameter sets: 1\n        Parameters:\n                   a: 42.0\n        '

    assert m.param_dim == 1
    m.param_dim = 42
    assert m.param_dim == 42
    m.param_dim = 1
    m.param_sets


def test_Model_array_parameter():
    m = NonFittableModel([[42, 43],[1,2]])
    m = NonFittableModel([42, 43, 44, 45])
    #assert m.param_dim == 2
    m.param_sets
    
    phi, theta, psi = 42, 43, 44
    model = models.RotateNative2Celestial(phi, theta, psi)
    assert_allclose(model.param_sets, [[0.73303829], [0.75049158], [0.76794487]])



def test_Model_add_model():
    m = models.Gaussian1DModel(1,2,3)
    m.add_model(m, 'p')
    m.add_model(m, 's')
    with pytest.raises(InputParameterError):
        m.add_model(m, 'q')
        m.add_model(m, 42)

def test_ParametericModel():
    with pytest.raises(TypeError):
        models.Gaussian1DModel(1, 2, 3, wrong=4)


def test_a():
    pass
