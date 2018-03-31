# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests models.parameters
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import itertools

import pytest
import numpy as np
from numpy.testing import (assert_allclose, assert_equal, assert_array_equal,
                           assert_almost_equal)

from . import irafutil
from .. import models, fitting
from ..core import Model, FittableModel
from ..parameters import Parameter, InputParameterError
from ...utils.data import get_pkg_data_filename


def setter1(val):
    return val


def setter2(val, model):
    model.do_something(val)
    return val * model.p


class SetterModel(FittableModel):

    inputs = ('x', 'y')
    outputs = ('z',)

    xc = Parameter(default=1, setter=setter1)
    yc = Parameter(default=1, setter=setter2)

    def __init__(self, xc, yc, p):
        self.p = p  # p is a value intended to be used by the setter
        super(SetterModel, self).__init__()
        self.xc = xc
        self.yc = yc

    def evaluate(self, x, y, xc, yc):
        return ((x - xc)**2 + (y - yc)**2)

    def do_something(self, v):
        pass


class TParModel(Model):
    """
    A toy model to test parameters machinery
    """

    coeff = Parameter()
    e = Parameter()

    def __init__(self, coeff, e, **kwargs):
        super(TParModel, self).__init__(coeff=coeff, e=e, **kwargs)

    @staticmethod
    def evaluate(coeff, e):
        pass


class MockModel(FittableModel):
    alpha = Parameter(name='alpha', default=42)

    @staticmethod
    def evaluate(*args):
        pass


def test_parameter_properties():
    """Test if getting / setting of Parameter properties works."""

    m = MockModel()
    p = m.alpha

    assert p.name == 'alpha'

    # Parameter names are immutable
    with pytest.raises(AttributeError):
        p.name = 'beta'

    assert p.fixed is False
    p.fixed = True
    assert p.fixed is True

    assert p.tied is False
    p.tied = lambda _: 0

    p.tied = False
    assert p.tied is False

    assert p.min is None
    p.min = 42
    assert p.min == 42
    p.min = None
    assert p.min is None

    assert p.max is None
    # TODO: shouldn't setting a max < min give an error?
    p.max = 41
    assert p.max == 41


def test_parameter_operators():
    """Test if the parameter arithmetic operators work."""

    m = MockModel()
    par = m.alpha
    num = 42.
    val = 3

    assert par - val == num - val
    assert val - par == val - num
    assert par / val == num / val
    assert val / par == val / num
    assert par ** val == num ** val
    assert val ** par == val ** num
    assert par < 45
    assert par > 41
    assert par <= par
    assert par >= par
    assert par == par
    assert -par == -num
    assert abs(par) == abs(num)


class TestParameters(object):

    def setup_class(self):
        """
        Unit tests for parameters

        Read an iraf database file created by onedspec.identify.  Use the
        information to create a 1D Chebyshev model and perform the same fit.

        Create also a gausian model.
        """
        test_file = get_pkg_data_filename('data/idcompspec.fits')
        f = open(test_file)
        lines = f.read()
        reclist = lines.split("begin")
        f.close()
        record = irafutil.IdentifyRecord(reclist[1])
        self.icoeff = record.coeff
        order = int(record.fields['order'])
        self.model = models.Chebyshev1D(order - 1)
        self.gmodel = models.Gaussian1D(2, mean=3, stddev=4)
        self.linear_fitter = fitting.LinearLSQFitter()
        self.x = record.x
        self.y = record.z
        self.yy = np.array([record.z, record.z])

    def test_set_slice(self):
        """
        Tests updating the parameters attribute with a slice.

        This is what fitters internally do.
        """

        self.model.parameters[:] = np.array([3, 4, 5, 6, 7])
        assert (self.model.parameters == [3., 4., 5., 6., 7.]).all()

    def test_set_parameters_as_list(self):
        """Tests updating parameters using a list."""

        self.model.parameters = [30, 40, 50, 60, 70]
        assert (self.model.parameters == [30., 40., 50., 60, 70]).all()

    def test_set_parameters_as_array(self):
        """Tests updating parameters using an array."""

        self.model.parameters = np.array([3, 4, 5, 6, 7])
        assert (self.model.parameters == [3., 4., 5., 6., 7.]).all()

    def test_set_as_tuple(self):
        """Tests updating parameters using a tuple."""

        self.model.parameters = (1, 2, 3, 4, 5)
        assert (self.model.parameters == [1, 2, 3, 4, 5]).all()

    def test_set_model_attr_seq(self):
        """
        Tests updating the parameters attribute when a model's
        parameter (in this case coeff) is updated.
        """

        self.model.parameters = [0, 0., 0., 0, 0]
        self.model.c0 = 7
        assert (self.model.parameters == [7, 0., 0., 0, 0]).all()

    def test_set_model_attr_num(self):
        """Update the parameter list when a model's parameter is updated."""

        self.gmodel.amplitude = 7
        assert (self.gmodel.parameters == [7, 3, 4]).all()

    def test_set_item(self):
        """Update the parameters using indexing."""

        self.model.parameters = [1, 2, 3, 4, 5]
        self.model.parameters[0] = 10.
        assert (self.model.parameters == [10, 2, 3, 4, 5]).all()
        assert self.model.c0 == 10

    def test_wrong_size1(self):
        """
        Tests raising an error when attempting to reset the parameters
        using a list of a different size.
        """

        with pytest.raises(InputParameterError):
            self.model.parameters = [1, 2, 3]

    def test_wrong_size2(self):
        """
        Tests raising an exception when attempting to update a model's
        parameter (in this case coeff) with a sequence of the wrong size.
        """

        with pytest.raises(InputParameterError):
            self.model.c0 = [1, 2, 3]

    def test_wrong_shape(self):
        """
        Tests raising an exception when attempting to update a model's
        parameter and the new value has the wrong shape.
        """

        with pytest.raises(InputParameterError):
            self.gmodel.amplitude = [1, 2]

    def test_par_against_iraf(self):
        """
        Test the fitter modifies model.parameters.

        Uses an iraf example.
        """

        new_model = self.linear_fitter(self.model, self.x, self.y)
        print(self.y, self.x)
        assert_allclose(new_model.parameters,
                        np.array(
                            [4826.1066602783685, 952.8943813407858,
                             12.641236013982386,
                             -1.7910672553339604,
                             0.90252884366711317]),
                        rtol=10 ** (-2))

    def testPolynomial1D(self):
        d = {'c0': 11, 'c1': 12, 'c2': 13, 'c3': 14}
        p1 = models.Polynomial1D(3, **d)
        assert_equal(p1.parameters, [11, 12, 13, 14])

    def test_poly1d_multiple_sets(self):
        p1 = models.Polynomial1D(3, n_models=3)
        assert_equal(p1.parameters, [0.0, 0.0, 0.0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0])
        assert_array_equal(p1.c0, [0, 0, 0])
        p1.c0 = [10, 10, 10]
        assert_equal(p1.parameters, [10.0, 10.0, 10.0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0])

    def test_par_slicing(self):
        """
        Test assigning to a parameter slice
        """
        p1 = models.Polynomial1D(3, n_models=3)
        p1.c0[:2] = [10, 10]
        assert_equal(p1.parameters, [10.0, 10.0, 0.0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0])

    def test_poly2d(self):
        p2 = models.Polynomial2D(degree=3)
        p2.c0_0 = 5
        assert_equal(p2.parameters, [5, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_poly2d_multiple_sets(self):
        kw = {'c0_0': [2, 3], 'c1_0': [1, 2], 'c2_0': [4, 5],
              'c0_1': [1, 1], 'c0_2': [2, 2], 'c1_1': [5, 5]}
        p2 = models.Polynomial2D(2, **kw)
        assert_equal(p2.parameters, [2, 3, 1, 2, 4, 5,
                                     1, 1, 2, 2, 5, 5])

    def test_shift_model_parameters1d(self):
        sh1 = models.Shift(2)
        sh1.offset = 3
        assert sh1.offset == 3
        assert sh1.offset.value == 3

    def test_scale_model_parametersnd(self):
        sc1 = models.Scale([2, 2])
        sc1.factor = [3, 3]
        assert np.all(sc1.factor == [3, 3])
        assert_array_equal(sc1.factor.value, [3, 3])

    def test_parameters_wrong_shape(self):
        sh1 = models.Shift(2)
        with pytest.raises(InputParameterError):
            sh1.offset = [3, 3]


class TestMultipleParameterSets(object):

    def setup_class(self):
        self.x1 = np.arange(1, 10, .1)
        self.y, self.x = np.mgrid[:10, :7]
        self.x11 = np.array([self.x1, self.x1]).T
        self.gmodel = models.Gaussian1D([12, 10], [3.5, 5.2], stddev=[.4, .7],
                                        n_models=2)

    def test_change_par(self):
        """
        Test that a change to one parameter as a set propagates to param_sets.
        """
        self.gmodel.amplitude = [1, 10]
        assert_almost_equal(
            self.gmodel.param_sets,
            np.array([[1.,
                       10],
                      [3.5,
                       5.2],
                      [0.4,
                       0.7]]))
        np.all(self.gmodel.parameters == [1.0, 10.0, 3.5, 5.2, 0.4, 0.7])

    def test_change_par2(self):
        """
        Test that a change to one single parameter in a set propagates to
        param_sets.
        """
        self.gmodel.amplitude[0] = 11
        assert_almost_equal(
            self.gmodel.param_sets,
            np.array([[11.,
                       10],
                      [3.5,
                       5.2],
                      [0.4,
                       0.7]]))
        np.all(self.gmodel.parameters == [11.0, 10.0, 3.5, 5.2, 0.4, 0.7])

    def test_change_parameters(self):
        self.gmodel.parameters = [13, 10, 9, 5.2, 0.4, 0.7]
        assert_almost_equal(self.gmodel.amplitude.value, [13., 10.])
        assert_almost_equal(self.gmodel.mean.value, [9., 5.2])


class TestParameterInitialization(object):
    """
    This suite of tests checks most if not all cases if instantiating a model
    with parameters of different shapes/sizes and with different numbers of
    parameter sets.
    """

    def test_single_model_scalar_parameters(self):
        t = TParModel(10, 1)
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(t.param_sets == [[10], [1]])
        assert np.all(t.parameters == [10, 1])
        assert t.coeff.shape == ()
        assert t.e.shape == ()

    def test_single_model_scalar_and_array_parameters(self):
        t = TParModel(10, [1, 2])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert len(t.param_sets) == 2
        assert np.all(t.param_sets[0] == [10])
        assert np.all(t.param_sets[1] == [[1, 2]])
        assert np.all(t.parameters == [10, 1, 2])
        assert t.coeff.shape == ()
        assert t.e.shape == (2,)

    def test_single_model_1d_array_parameters(self):
        t = TParModel([10, 20], [1, 2])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(t.param_sets == [[[10, 20]], [[1, 2]]])
        assert np.all(t.parameters == [10, 20, 1, 2])
        assert t.coeff.shape == (2,)
        assert t.e.shape == (2,)

    def test_single_model_1d_array_different_length_parameters(self):
        with pytest.raises(InputParameterError):
            # Not broadcastable
            t = TParModel([1, 2], [3, 4, 5])

    def test_single_model_2d_array_parameters(self):
        t = TParModel([[10, 20], [30, 40]], [[1, 2], [3, 4]])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(t.param_sets == [[[[10, 20], [30, 40]]],
                                       [[[1, 2], [3, 4]]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 1, 2, 3, 4])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2, 2)

    def test_single_model_2d_non_square_parameters(self):
        coeff = np.array([[10, 20], [30, 40], [50, 60]])
        e = np.array([[1, 2], [3, 4], [5, 6]])

        t = TParModel(coeff, e)
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(t.param_sets == [[[[10, 20], [30, 40], [50, 60]]],
                                       [[[1, 2], [3, 4], [5, 6]]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60,
                                       1, 2, 3, 4, 5, 6])
        assert t.coeff.shape == (3, 2)
        assert t.e.shape == (3, 2)

        t2 = TParModel(coeff.T, e.T)
        assert len(t2) == 1
        assert t2.model_set_axis is False
        assert np.all(t2.param_sets == [[[[10, 30, 50], [20, 40, 60]]],
                                        [[[1, 3, 5], [2, 4, 6]]]])
        assert np.all(t2.parameters == [10, 30, 50, 20, 40, 60,
                                        1, 3, 5, 2, 4, 6])
        assert t2.coeff.shape == (2, 3)
        assert t2.e.shape == (2, 3)

        # Not broadcastable
        with pytest.raises(InputParameterError):
            TParModel(coeff, e.T)

        with pytest.raises(InputParameterError):
            TParModel(coeff.T, e)

    def test_single_model_2d_broadcastable_parameters(self):
        t = TParModel([[10, 20, 30], [40, 50, 60]], [1, 2, 3])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[[10, 20, 30], [40, 50, 60]]])
        assert np.all(t.param_sets[1] == [[1, 2, 3]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 1, 2, 3])

    @pytest.mark.parametrize(('p1', 'p2'), [
        (1, 2), (1, [2, 3]), ([1, 2], 3), ([1, 2, 3], [4, 5]),
        ([1, 2], [3, 4, 5])])
    def test_two_model_incorrect_scalar_parameters(self, p1, p2):
        with pytest.raises(InputParameterError):
            TParModel(p1, p2, n_models=2)

    @pytest.mark.parametrize('kwargs', [
        {'n_models': 2}, {'model_set_axis': 0},
        {'n_models': 2, 'model_set_axis': 0}])
    def test_two_model_scalar_parameters(self, kwargs):
        t = TParModel([10, 20], [1, 2], **kwargs)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(t.param_sets == [[10, 20], [1, 2]])
        assert np.all(t.parameters == [10, 20, 1, 2])
        assert t.coeff.shape == ()
        assert t.e.shape == ()

    @pytest.mark.parametrize('kwargs', [
        {'n_models': 2}, {'model_set_axis': 0},
        {'n_models': 2, 'model_set_axis': 0}])
    def test_two_model_scalar_and_array_parameters(self, kwargs):
        t = TParModel([10, 20], [[1, 2], [3, 4]], **kwargs)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[10], [20]])
        assert np.all(t.param_sets[1] == [[1, 2], [3, 4]])
        assert np.all(t.parameters == [10, 20, 1, 2, 3, 4])
        assert t.coeff.shape == ()
        assert t.e.shape == (2,)

    def test_two_model_1d_array_parameters(self):
        t = TParModel([[10, 20], [30, 40]], [[1, 2], [3, 4]], n_models=2)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(t.param_sets == [[[10, 20], [30, 40]],
                                       [[1, 2], [3, 4]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 1, 2, 3, 4])
        assert t.coeff.shape == (2,)
        assert t.e.shape == (2,)

        t2 = TParModel([[10, 20, 30], [40, 50, 60]],
                       [[1, 2, 3], [4, 5, 6]], n_models=2)
        assert len(t2) == 2
        assert t2.model_set_axis == 0
        assert np.all(t2.param_sets == [[[10, 20, 30], [40, 50, 60]],
                                        [[1, 2, 3], [4, 5, 6]]])
        assert np.all(t2.parameters == [10, 20, 30, 40, 50, 60,
                                        1, 2, 3, 4, 5, 6])
        assert t2.coeff.shape == (3,)
        assert t2.e.shape == (3,)

    def test_two_model_mixed_dimension_array_parameters(self):
        with pytest.raises(InputParameterError):
            # Can't broadcast different array shapes
            TParModel([[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
                      [[9, 10, 11], [12, 13, 14]], n_models=2)

        t = TParModel([[[10, 20], [30, 40]], [[50, 60], [70, 80]]],
                      [[1, 2], [3, 4]], n_models=2)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[[10, 20], [30, 40]],
                                          [[50, 60], [70, 80]]])
        assert np.all(t.param_sets[1] == [[[1, 2]], [[3, 4]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 70, 80,
                                       1, 2, 3, 4])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2,)

    def test_two_model_2d_array_parameters(self):
        t = TParModel([[[10, 20], [30, 40]], [[50, 60], [70, 80]]],
                      [[[1, 2], [3, 4]], [[5, 6], [7, 8]]], n_models=2)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(t.param_sets == [[[[10, 20], [30, 40]],
                                        [[50, 60], [70, 80]]],
                                       [[[1, 2], [3, 4]],
                                        [[5, 6], [7, 8]]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 70, 80,
                                       1, 2, 3, 4, 5, 6, 7, 8])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2, 2)

    def test_two_model_nonzero_model_set_axis(self):
        # An example where the model set axis is the *last* axis of the
        # parameter arrays
        coeff = np.array([[[10, 20], [30, 40]], [[50, 60], [70, 80]]])
        coeff = np.rollaxis(coeff, 0, 3)
        e = np.array([[1, 2], [3, 4]])
        e = np.rollaxis(e, 0, 2)
        t = TParModel(coeff, e, model_set_axis=-1)
        assert len(t) == 2
        assert t.model_set_axis == -1
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[[10, 50], [20, 60]],
                                          [[30, 70], [40, 80]]])
        assert np.all(t.param_sets[1] == [[[1, 3], [2, 4]]])
        assert np.all(t.parameters == [10, 50, 20, 60, 30, 70, 40, 80,
                                       1, 3, 2, 4])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2,)

    def test_wrong_number_of_params(self):
        with pytest.raises(InputParameterError):
            TParModel(coeff=[[1, 2], [3, 4]], e=(2, 3, 4), n_models=2)
        with pytest.raises(InputParameterError):
            TParModel(coeff=[[1, 2], [3, 4]], e=(2, 3, 4), model_set_axis=0)

    def test_wrong_number_of_params2(self):
        with pytest.raises(InputParameterError):
            m = TParModel(coeff=[[1, 2], [3, 4]], e=4, n_models=2)
        with pytest.raises(InputParameterError):
            m = TParModel(coeff=[[1, 2], [3, 4]], e=4, model_set_axis=0)

    def test_array_parameter1(self):
        with pytest.raises(InputParameterError):
            t = TParModel(np.array([[1, 2], [3, 4]]), 1, model_set_axis=0)

    def test_array_parameter2(self):
        with pytest.raises(InputParameterError):
            m = TParModel(np.array([[1, 2], [3, 4]]), (1, 1, 11),
                          model_set_axis=0)

    def test_array_parameter4(self):
        """
        Test multiple parameter model with array-valued parameters of the same
        size as the number of parameter sets.
        """

        t4 = TParModel([[1, 2], [3, 4]], [5, 6], model_set_axis=False)
        assert len(t4) == 1
        assert t4.coeff.shape == (2, 2)
        assert t4.e.shape == (2,)
        assert np.issubdtype(t4.param_sets.dtype, np.object_)
        assert np.all(t4.param_sets[0] == [[1, 2], [3, 4]])
        assert np.all(t4.param_sets[1] == [5, 6])


def test_non_broadcasting_parameters():
    """
    Tests that in a model with 3 parameters that do not all mutually broadcast,
    this is determined correctly regardless of what order the parameters are
    in.
    """

    a = 3
    b = np.array([[1, 2, 3], [4, 5, 6]])
    c = np.array([[1, 2, 3, 4], [1, 2, 3, 4]])

    class TestModel(Model):
        p1 = Parameter()
        p2 = Parameter()
        p3 = Parameter()

        def evaluate(self, *args):
            return

    # a broadcasts with both b and c, but b does not broadcast with c
    for args in itertools.permutations((a, b, c)):
        with pytest.raises(InputParameterError):
            TestModel(*args)


def test_setter():
    pars = np.random.rand(20).reshape((10, 2))

    model = SetterModel(-1, 3, np.pi)

    for x, y in pars:
        model.x = x
        model.y = y
        assert_almost_equal(model(x, y), (x + 1)**2 + (y - np.pi * 3)**2)
