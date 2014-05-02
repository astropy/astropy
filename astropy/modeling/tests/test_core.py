# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import collections
import pytest
import numpy as np
from numpy.testing.utils import assert_allclose
from ..core import Model, InputParameterError
from ..parameters import Parameter
from .. import models


class NonFittableModel(Model):
    """An example class directly subclassing Model for testing"""

    a = Parameter()

    def __init__(self, a):
        if not isinstance(a, collections.Sequence):
            param_dim = 1
        else:
            param_dim = len(a)

        super(NonFittableModel, self).__init__(a, param_dim=param_dim)

    def __call__(self):
        pass


def test_Model():
    """Some silly tests just to have all lines in the Model code covered by unit tests"""
    m = NonFittableModel(42)
    assert repr(m) == "<NonFittableModel(a=42.0)>"
    assert (str(m) ==
        "Model: NonFittableModel\n"
        "Inputs: 1\n"
        "Outputs: 1\n"
        "Parameter sets: 1\n"
        "Parameters:\n"
        "     a  \n"
        "    ----\n"
        "    42.0")

    assert m.param_dim == 1

    with pytest.raises(AttributeError):
        m.param_dim = 42


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


def test_a():
    pass
