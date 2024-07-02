# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests models.parameters
"""
# pylint: disable=invalid-name

import functools
import itertools
import unittest.mock as mk

import numpy as np
import pytest

from astropy import units as u
from astropy.modeling import fitting, models
from astropy.modeling.core import FittableModel, Model
from astropy.modeling.parameters import (
    InputParameterError,
    Parameter,
    _tofloat,
    param_repr_oneline,
)
from astropy.utils.data import get_pkg_data_filename

from . import irafutil


def setter1(val):
    return val


def setter2(val, model):
    model.do_something(val)
    return val * model.p


def getter1(val):
    return val


class SetterModel(FittableModel):
    n_inputs = 2
    n_outputs = 1

    xc = Parameter(default=1, setter=setter1, getter=getter1)
    yc = Parameter(default=1, setter=setter2, getter=getter1)

    def do_something(self, v):
        pass

    def __init__(self, xc, yc, p):
        self.p = p  # p is a value intended to be used by the setter
        super().__init__()
        self.xc = xc
        self.yc = yc

    def evaluate(self, x, y, xc, yc):
        return (x - xc) ** 2 + (y - yc) ** 2


class TParModel(Model):
    """
    A toy model to test parameters machinery
    """

    coeff = Parameter()
    e = Parameter()

    def __init__(self, coeff, e, **kwargs):
        super().__init__(coeff=coeff, e=e, **kwargs)

    @staticmethod
    def evaluate(coeff, e):
        pass


class MockModel(FittableModel):
    alpha = Parameter(name="alpha", default=42)

    @staticmethod
    def evaluate(*args):
        pass


def test__tofloat():
    # iterable
    value = _tofloat([1, 2, 3])
    assert isinstance(value, np.ndarray)
    assert (value == np.array([1, 2, 3])).all()
    assert np.all([isinstance(val, float) for val in value])
    value = _tofloat(np.array([1, 2, 3]))
    assert isinstance(value, np.ndarray)
    assert (value == np.array([1, 2, 3])).all()
    assert np.all([isinstance(val, float) for val in value])
    MESSAGE = r"Parameter of .* could not be converted to float"
    with pytest.raises(InputParameterError, match=MESSAGE):
        _tofloat("test")

    # quantity
    assert _tofloat(1 * u.m) == 1 * u.m

    # dimensions/scalar array
    value = _tofloat(np.asanyarray(3))
    assert isinstance(value, float)
    assert value == 3

    # A regular number
    value = _tofloat(3)
    assert isinstance(value, float)
    assert value == 3
    value = _tofloat(3.0)
    assert isinstance(value, float)
    assert value == 3
    value = _tofloat(np.float32(3))
    assert isinstance(value, float)
    assert value == 3
    value = _tofloat(np.float64(3))
    assert isinstance(value, float)
    assert value == 3
    value = _tofloat(np.int32(3))
    assert isinstance(value, float)
    assert value == 3
    value = _tofloat(np.int64(3))
    assert isinstance(value, float)
    assert value == 3

    # boolean
    MESSAGE = r"Expected parameter to be of numerical type, not boolean"
    with pytest.raises(InputParameterError, match=MESSAGE):
        _tofloat(True)
    with pytest.raises(InputParameterError, match=MESSAGE):
        _tofloat(False)

    # other
    class Value:
        pass

    MESSAGE = r"Don't know how to convert parameter of .* to float"
    with pytest.raises(InputParameterError, match=MESSAGE):
        _tofloat(Value)


def test_parameter_properties():
    """Test if getting / setting of Parameter properties works."""

    p = Parameter("alpha", default=1)

    assert p.name == "alpha"

    # Parameter names are immutable
    with pytest.raises(AttributeError):
        p.name = "beta"

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
    p.max = 41
    assert p.max == 41


def test_parameter_operators():
    """Test if the parameter arithmetic operators work."""

    par = Parameter("alpha", default=42)
    num = 42.0
    val = 3

    assert par - val == num - val
    assert val - par == val - num
    assert par / val == num / val
    assert val / par == val / num
    assert par**val == num**val
    assert val**par == val**num
    assert par < 45
    assert par > 41
    assert par <= par
    assert par >= par
    assert par == par
    assert -par == -num
    assert abs(par) == abs(num)


# Test inherited models


class M1(Model):
    m1a = Parameter(default=1.0)
    m1b = Parameter(default=5.0)

    def evaluate():
        pass


class M2(M1):
    m2c = Parameter(default=11.0)


class M3(M2):
    m3d = Parameter(default=20.0)


def test_parameter_inheritance():
    mod = M3()
    assert mod.m1a == 1.0
    assert mod.m1b == 5.0
    assert mod.m2c == 11.0
    assert mod.m3d == 20.0
    for key in ["m1a", "m1b", "m2c", "m3d"]:
        assert key in mod.__dict__
    assert mod.param_names == ("m1a", "m1b", "m2c", "m3d")


def test_param_metric():
    mod = M3()
    assert mod._param_metrics["m1a"]["slice"] == slice(0, 1)
    assert mod._param_metrics["m1b"]["slice"] == slice(1, 2)
    assert mod._param_metrics["m2c"]["slice"] == slice(2, 3)
    assert mod._param_metrics["m3d"]["slice"] == slice(3, 4)
    mod._parameters_to_array()
    assert (mod._parameters == np.array([1.0, 5.0, 11.0, 20], dtype=np.float64)).all()


class TestParameters:
    def setup_class(self):
        """
        Unit tests for parameters

        Read an iraf database file created by onedspec.identify.  Use the
        information to create a 1D Chebyshev model and perform the same fit.

        Create also a gaussian model.
        """
        test_file = get_pkg_data_filename("data/idcompspec.fits")
        f = open(test_file)
        lines = f.read()
        reclist = lines.split("begin")
        f.close()
        record = irafutil.IdentifyRecord(reclist[1])
        self.icoeff = record.coeff
        order = int(record.fields["order"])
        self.model = models.Chebyshev1D(order - 1)
        self.gmodel = models.Gaussian1D(2, mean=3, stddev=4)
        self.linear_fitter = fitting.LinearLSQFitter()
        self.x = record.x
        self.y = record.z
        self.yy = np.array([record.z, record.z])

    def test_set_parameters_as_list(self):
        """Tests updating parameters using a list."""

        self.model.parameters = [30, 40, 50, 60, 70]
        assert (self.model.parameters == [30.0, 40.0, 50.0, 60, 70]).all()

    def test_set_parameters_as_array(self):
        """Tests updating parameters using an array."""

        self.model.parameters = np.array([3, 4, 5, 6, 7])
        assert (self.model.parameters == [3.0, 4.0, 5.0, 6.0, 7.0]).all()

    def test_set_as_tuple(self):
        """Tests updating parameters using a tuple."""

        self.model.parameters = (1, 2, 3, 4, 5)
        assert (self.model.parameters == [1, 2, 3, 4, 5]).all()

    def test_set_model_attr_seq(self):
        """
        Tests updating the parameters attribute when a model's
        parameter (in this case coeff) is updated.
        """

        self.model.parameters = [0, 0.0, 0.0, 0, 0]
        self.model.c0 = 7
        assert (self.model.parameters == [7, 0.0, 0.0, 0, 0]).all()

    def test_set_model_attr_num(self):
        """Update the parameter list when a model's parameter is updated."""

        self.gmodel.amplitude = 7
        assert (self.gmodel.parameters == [7, 3, 4]).all()

    def test_set_item(self):
        """Update the parameters using indexing."""

        self.model.parameters = [1, 2, 3, 4, 5]
        tpar = self.model.parameters
        tpar[0] = 10.0
        self.model.parameters = tpar
        assert (self.model.parameters == [10, 2, 3, 4, 5]).all()
        assert self.model.c0 == 10

    def test_wrong_size1(self):
        """
        Tests raising an error when attempting to reset the parameters
        using a list of a different size.
        """

        MESSAGE = (
            r"Input parameter values not compatible with the model parameters array: .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            self.model.parameters = [1, 2, 3]

    def test_wrong_size2(self):
        """
        Tests raising an exception when attempting to update a model's
        parameter (in this case coeff) with a sequence of the wrong size.
        """

        MESSAGE = (
            r"Value for parameter c0 does not match shape or size\nexpected by model .*"
            r" vs .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            self.model.c0 = [1, 2, 3]

    def test_wrong_shape(self):
        """
        Tests raising an exception when attempting to update a model's
        parameter and the new value has the wrong shape.
        """

        MESSAGE = (
            r"Value for parameter amplitude does not match shape or size\nexpected by"
            r" model .* vs .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            self.gmodel.amplitude = [1, 2]

    def test_par_against_iraf(self):
        """
        Test the fitter modifies model.parameters.

        Uses an iraf example.
        """

        new_model = self.linear_fitter(self.model, self.x, self.y)
        np.testing.assert_allclose(
            new_model.parameters,
            np.array(
                [
                    4826.1066602783685,
                    952.8943813407858,
                    12.641236013982386,
                    -1.7910672553339604,
                    0.90252884366711317,
                ]
            ),
            rtol=10 ** (-2),
        )

    def testPolynomial1D(self):
        d = {"c0": 11, "c1": 12, "c2": 13, "c3": 14}
        p1 = models.Polynomial1D(3, **d)
        np.testing.assert_equal(p1.parameters, [11, 12, 13, 14])

    def test_poly1d_multiple_sets(self):
        p1 = models.Polynomial1D(3, n_models=3)
        np.testing.assert_equal(
            p1.parameters, [0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )
        np.testing.assert_array_equal(p1.c0, [0, 0, 0])
        p1.c0 = [10, 10, 10]
        np.testing.assert_equal(
            p1.parameters, [10.0, 10.0, 10.0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )

    def test_par_slicing(self):
        """
        Test assigning to a parameter slice
        """
        p1 = models.Polynomial1D(3, n_models=3)
        p1.c0[:2] = [10, 10]
        np.testing.assert_equal(
            p1.parameters, [10.0, 10.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        )

    def test_poly2d(self):
        p2 = models.Polynomial2D(degree=3)
        p2.c0_0 = 5
        np.testing.assert_equal(p2.parameters, [5, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_poly2d_multiple_sets(self):
        kw = {
            "c0_0": [2, 3],
            "c1_0": [1, 2],
            "c2_0": [4, 5],
            "c0_1": [1, 1],
            "c0_2": [2, 2],
            "c1_1": [5, 5],
        }
        p2 = models.Polynomial2D(2, **kw)
        np.testing.assert_equal(p2.parameters, [2, 3, 1, 2, 4, 5, 1, 1, 2, 2, 5, 5])

    def test_shift_model_parameters1d(self):
        sh1 = models.Shift(2)
        sh1.offset = 3
        assert sh1.offset == 3
        assert sh1.offset.value == 3

    def test_scale_model_parametersnd(self):
        sc1 = models.Scale([2, 2])
        sc1.factor = [3, 3]
        assert np.all(sc1.factor == [3, 3])
        np.testing.assert_array_equal(sc1.factor.value, [3, 3])

    def test_bounds(self):
        # Valid __init__
        param = Parameter(bounds=(1, 2))
        assert param.bounds == (1, 2)
        param = Parameter(min=1, max=2)
        assert param.bounds == (1, 2)

        # Errors __init__
        MESSAGE = r"bounds may not be specified simultaneously with min or max .*"
        with pytest.raises(ValueError, match=MESSAGE):
            Parameter(bounds=(1, 2), min=1, name="test")
        with pytest.raises(ValueError, match=MESSAGE):
            Parameter(bounds=(1, 2), max=2, name="test")
        with pytest.raises(ValueError, match=MESSAGE):
            Parameter(bounds=(1, 2), min=1, max=2, name="test")

        # Setters
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.bounds == (None, None) == param._bounds

        # Set errors
        MESSAGE = "{} value must be a number or a Quantity"
        with pytest.raises(TypeError, match=MESSAGE.format("Min")):
            param.bounds = ("test", None)
        with pytest.raises(TypeError, match=MESSAGE.format("Max")):
            param.bounds = (None, "test")

        # Set number
        param.bounds = (1, 2)
        assert param.bounds == (1, 2) == param._bounds

        # Set Quantity
        param.bounds = (1 * u.m, 2 * u.m)
        assert param.bounds == (1, 2) == param._bounds

    def test_modify_value(self):
        param = Parameter(name="test", default=[1, 2, 3])
        assert (param.value == [1, 2, 3]).all()

        # Errors
        MESSAGE = r"Slice assignment outside the parameter dimensions for 'test'"
        with pytest.raises(InputParameterError, match=MESSAGE):
            param[slice(0, 0)] = 2

        MESSAGE = r"Input dimension 3 invalid for 'test' parameter with dimension 1"
        with pytest.raises(InputParameterError, match=MESSAGE):
            param[3] = np.array([5])

        # assignment of a slice
        param[slice(0, 2)] = [4, 5]
        assert (param.value == [4, 5, 3]).all()

        # assignment of a value
        param[2] = 6
        assert (param.value == [4, 5, 6]).all()

    def test__set_unit(self):
        param = Parameter(name="test", default=[1, 2, 3])
        assert param.unit is None

        # No force Error (no existing unit)
        MESSAGE = r"Cannot attach units to parameters that were .*"
        with pytest.raises(ValueError, match=MESSAGE):
            param._set_unit(u.m)

        # Force
        param._set_unit(u.m, True)
        assert param.unit == u.m

        # Force magnitude unit (mag=False)
        MESSAGE = r"This parameter does not support the magnitude units such as .*"
        with pytest.raises(ValueError, match=MESSAGE):
            param._set_unit(u.ABmag, True)

        # Force magnitude unit (mag=True)
        param._mag = True
        param._set_unit(u.ABmag, True)
        assert param._unit == u.ABmag

        # No force Error (existing unit)
        MESSAGE = r"Cannot change the unit attribute directly, instead change the .*"
        with pytest.raises(ValueError, match=MESSAGE):
            param._set_unit(u.K)

    def test_quantity(self):
        param = Parameter(name="test", default=[1, 2, 3])
        assert param.unit is None
        assert param.quantity is None

        param = Parameter(name="test", default=[1, 2, 3], unit=u.m)
        assert param.unit == u.m
        assert (param.quantity == np.array([1, 2, 3]) * u.m).all()

    def test_shape(self):
        # Array like
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.shape == (4,)
        # Reshape error
        MESSAGE = r"cannot reshape array of size 4 into shape .*"
        with pytest.raises(ValueError, match=MESSAGE):
            param.shape = (5,)
        # Reshape success
        param.shape = (2, 2)
        assert param.shape == (2, 2)
        assert (param.value == [[1, 2], [3, 4]]).all()

        # Scalar
        param = Parameter(name="test", default=1)
        assert param.shape == ()
        # Reshape error
        MESSAGE = r"Cannot assign this shape to a scalar quantity"
        with pytest.raises(ValueError, match=MESSAGE):
            param.shape = (5,)
        param.shape = (1,)

        # single value
        param = Parameter(name="test", default=np.array([1]))
        assert param.shape == (1,)
        # Reshape error
        with pytest.raises(ValueError, match=MESSAGE):
            param.shape = (5,)
        param.shape = ()

    def test_size(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.size == 4

        param = Parameter(name="test", default=[1])
        assert param.size == 1

        param = Parameter(name="test", default=1)
        assert param.size == 1

    def test_std(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.std is None
        assert param._std is None

        param.std = 5
        assert param.std == 5 == param._std

    def test_fixed(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.fixed is False
        assert param._fixed is False

        # Set error
        MESSAGE = r"Value must be boolean"
        with pytest.raises(ValueError, match=MESSAGE):
            param.fixed = 3

        # Set
        param.fixed = True
        assert param.fixed is True
        assert param._fixed is True

    def test_tied(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.tied is False
        assert param._tied is False

        # Set error
        MESSAGE = r"Tied must be a callable or set to False or None"
        with pytest.raises(TypeError, match=MESSAGE):
            param.tied = mk.NonCallableMagicMock()

        # Set None
        param.tied = None
        assert param.tied is None
        assert param._tied is None

        # Set False
        param.tied = False
        assert param.tied is False
        assert param._tied is False

        # Set other
        tied = mk.MagicMock()
        param.tied = tied
        assert param.tied == tied == param._tied

    def test_validator(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param._validator is None

        valid = mk.MagicMock()
        param.validator(valid)
        assert param._validator == valid

        MESSAGE = r"This decorator method expects a callable.*"
        with pytest.raises(ValueError, match=MESSAGE):
            param.validator(mk.NonCallableMagicMock())

    def test_validate(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param._validator is None
        assert param.model is None

        # Run without validator
        param.validate(mk.MagicMock())

        # Run with validator but no Model
        validator = mk.MagicMock()
        param.validator(validator)
        assert param._validator == validator
        param.validate(mk.MagicMock())
        assert validator.call_args_list == []

        # Full validate
        param._model = mk.MagicMock()
        value = mk.MagicMock()
        param.validate(value)
        assert validator.call_args_list == [mk.call(param._model, value)]

    def test_copy(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        copy_param = param.copy()

        assert (param == copy_param).all()
        assert id(param) != id(copy_param)

    def test_model(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param.model is None
        assert param._model is None
        assert param._model_required is False
        assert (param._value == [1, 2, 3, 4]).all()

        setter = mk.MagicMock()
        getter = mk.MagicMock()
        param._setter = setter
        param._getter = getter

        # No Model Required
        param._value = [5, 6, 7, 8]
        model0 = mk.MagicMock()
        setter0 = mk.MagicMock()
        getter0 = mk.MagicMock()
        with mk.patch.object(
            Parameter, "_create_value_wrapper", side_effect=[setter0, getter0]
        ) as mkCreate:
            param.model = model0
            assert param.model == model0 == param._model
            assert param._setter == setter0
            assert param._getter == getter0
            assert mkCreate.call_args_list == [
                mk.call(setter, model0),
                mk.call(getter, model0),
            ]
            assert param._value == [5, 6, 7, 8]

        param._setter = setter
        param._getter = getter

        # Model required
        param._model_required = True
        model1 = mk.MagicMock()
        setter1 = mk.MagicMock()
        getter1 = mk.MagicMock()
        setter1.return_value = np.array([9, 10, 11, 12])
        getter1.return_value = np.array([9, 10, 11, 12])
        with mk.patch.object(
            Parameter, "_create_value_wrapper", side_effect=[setter1, getter1]
        ) as mkCreate:
            param.model = model1
            assert param.model == model1 == param._model
            assert param._setter == setter1
            assert param._getter == getter1
            assert mkCreate.call_args_list == [
                mk.call(setter, model1),
                mk.call(getter, model1),
            ]
            assert (param.value == [9, 10, 11, 12]).all()

        param._setter = setter
        param._getter = getter
        param._default = None
        with mk.patch.object(
            Parameter, "_create_value_wrapper", side_effect=[setter1, getter1]
        ) as mkCreate:
            param.model = model1
            assert param.model == model1 == param._model
            assert param._setter == setter1
            assert param._getter == getter1
            assert mkCreate.call_args_list == [
                mk.call(setter, model1),
                mk.call(getter, model1),
            ]
            assert param._value is None

    def test_value(self):
        param = Parameter(name="test", default=1)
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(name="test", default=[1])
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(name="test", default=[[1]])
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(name="test", default=np.array([1]))
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(name="test", default=[1, 2, 3])
        assert isinstance(param.value, np.ndarray)
        assert (param.value == [1, 2, 3]).all()

        param = Parameter(name="test", default=[1], setter=setter1, getter=getter1)
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(name="test", default=[[1]], setter=setter1, getter=getter1)
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

        param = Parameter(
            name="test", default=np.array([1]), setter=setter1, getter=getter1
        )
        assert not isinstance(param.value, np.ndarray)
        assert param.value == 1

    def test_raw_value(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])

        # Normal case
        assert (param._raw_value == param.value).all()

        # Bad setter
        param._setter = True
        param._internal_value = 4
        assert param._raw_value == 4

    def test__create_value_wrapper(self):
        param = Parameter(name="test", default=[1, 2, 3, 4])

        # Bad ufunc
        MESSAGE = r"A numpy.ufunc used for Parameter getter/setter .*"
        with pytest.raises(TypeError, match=MESSAGE):
            param._create_value_wrapper(np.add, mk.MagicMock())
        # Good ufunc
        with mk.patch(
            "astropy.modeling.parameters._wrap_ufunc", autospec=True
        ) as mkWrap:
            assert (
                param._create_value_wrapper(np.negative, mk.MagicMock())
                == mkWrap.return_value
            )
            assert mkWrap.call_args_list == [mk.call(np.negative)]

        # None
        assert param._create_value_wrapper(None, mk.MagicMock()) is None

        # wrapper with one argument
        def wrapper1(a):
            pass

        assert param._create_value_wrapper(wrapper1, mk.MagicMock()) == wrapper1

        # wrapper with two argument2
        def wrapper2(a, b):
            pass

        # model is None
        assert param._model_required is False
        assert param._create_value_wrapper(wrapper2, None) == wrapper2
        assert param._model_required is True
        # model is not None
        param._model_required = False
        model = mk.MagicMock()
        partial_wrapper = param._create_value_wrapper(wrapper2, model)
        assert isinstance(partial_wrapper, functools.partial)
        assert partial_wrapper.func is wrapper2
        assert partial_wrapper.args == ()
        assert list(partial_wrapper.keywords.keys()) == ["b"]
        assert partial_wrapper.keywords["b"] is model

        # wrapper with more than 2 arguments
        def wrapper3(a, b, c):
            pass

        MESSAGE = r"Parameter getter/setter must be a function .*"
        with pytest.raises(TypeError, match=MESSAGE):
            param._create_value_wrapper(wrapper3, mk.MagicMock())

    def test_bool(self):
        # single value is true
        param = Parameter(name="test", default=1)
        assert param.value == 1
        assert np.all(param)
        assert param

        # single value is false
        param = Parameter(name="test", default=0)
        assert param.value == 0
        assert not np.all(param)
        assert not param

        # vector value all true
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert np.all(param.value == [1, 2, 3, 4])
        assert np.all(param)
        assert param

        # vector value at least one false
        param = Parameter(name="test", default=[1, 2, 0, 3, 4])
        assert np.all(param.value == [1, 2, 0, 3, 4])
        assert not np.all(param)
        assert not param

    def test_param_repr_oneline(self):
        # Single value no units
        param = Parameter(name="test", default=1)
        assert param_repr_oneline(param) == "1."

        # Vector value no units
        param = Parameter(name="test", default=[1, 2, 3, 4])
        assert param_repr_oneline(param) == "[1., 2., 3., 4.]"

        # Single value units
        param = Parameter(name="test", default=1 * u.m)
        assert param_repr_oneline(param) == "1. m"

        # Vector value units
        param = Parameter(name="test", default=[1, 2, 3, 4] * u.m)
        assert param_repr_oneline(param) == "[1., 2., 3., 4.] m"

    def test_getter_setter(self):
        msg = "setter and getter must both be input"
        with pytest.raises(ValueError, match=msg):
            Parameter(name="test", default=1, getter=getter1)
        with pytest.raises(ValueError, match=msg):
            Parameter(name="test", default=1, setter=setter1)


class TestMultipleParameterSets:
    def setup_class(self):
        self.x1 = np.arange(1, 10, 0.1)
        self.y, self.x = np.mgrid[:10, :7]
        self.x11 = np.array([self.x1, self.x1]).T
        self.gmodel = models.Gaussian1D(
            [12, 10], [3.5, 5.2], stddev=[0.4, 0.7], n_models=2
        )

    def test_change_par(self):
        """
        Test that a change to one parameter as a set propagates to param_sets.
        """
        self.gmodel.amplitude = [1, 10]
        np.testing.assert_almost_equal(
            self.gmodel.param_sets,
            np.array(
                [
                    [1.0, 10],
                    [3.5, 5.2],
                    [0.4, 0.7],
                ]
            ),
        )
        np.all(self.gmodel.parameters == [1.0, 10.0, 3.5, 5.2, 0.4, 0.7])

    def test_change_par2(self):
        """
        Test that a change to one single parameter in a set propagates to
        param_sets.
        """
        self.gmodel.amplitude[0] = 11
        np.testing.assert_almost_equal(
            self.gmodel.param_sets,
            np.array(
                [
                    [11.0, 10],
                    [3.5, 5.2],
                    [0.4, 0.7],
                ]
            ),
        )
        np.all(self.gmodel.parameters == [11.0, 10.0, 3.5, 5.2, 0.4, 0.7])

    def test_change_parameters(self):
        self.gmodel.parameters = [13, 10, 9, 5.2, 0.4, 0.7]
        np.testing.assert_almost_equal(self.gmodel.amplitude.value, [13.0, 10.0])
        np.testing.assert_almost_equal(self.gmodel.mean.value, [9.0, 5.2])


class TestParameterInitialization:
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
        MESSAGE = (
            r"Parameter .* of shape .* cannot be broadcast with parameter .* of"
            r" shape .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            # Not broadcastable
            TParModel([1, 2], [3, 4, 5])

    def test_single_model_2d_array_parameters(self):
        t = TParModel([[10, 20], [30, 40]], [[1, 2], [3, 4]])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(
            t.param_sets
            == [
                [[[10, 20], [30, 40]]],
                [[[1, 2], [3, 4]]],
            ]
        )
        assert np.all(t.parameters == [10, 20, 30, 40, 1, 2, 3, 4])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2, 2)

    def test_single_model_2d_non_square_parameters(self):
        coeff = np.array(
            [
                [10, 20],
                [30, 40],
                [50, 60],
            ]
        )
        e = np.array([[1, 2], [3, 4], [5, 6]])

        t = TParModel(coeff, e)
        assert len(t) == 1
        assert t.model_set_axis is False
        assert np.all(
            t.param_sets
            == [
                [[[10, 20], [30, 40], [50, 60]]],
                [[[1, 2], [3, 4], [5, 6]]],
            ]
        )
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 1, 2, 3, 4, 5, 6])
        assert t.coeff.shape == (3, 2)
        assert t.e.shape == (3, 2)

        t2 = TParModel(coeff.T, e.T)
        assert len(t2) == 1
        assert t2.model_set_axis is False
        assert np.all(
            t2.param_sets
            == [
                [[[10, 30, 50], [20, 40, 60]]],
                [[[1, 3, 5], [2, 4, 6]]],
            ]
        )
        assert np.all(t2.parameters == [10, 30, 50, 20, 40, 60, 1, 3, 5, 2, 4, 6])
        assert t2.coeff.shape == (2, 3)
        assert t2.e.shape == (2, 3)

        # Not broadcastable
        MESSAGE = (
            r"Parameter .* of shape .* cannot be broadcast with parameter .* of"
            r" shape .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff, e.T)

        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff.T, e)

    def test_single_model_2d_broadcastable_parameters(self):
        t = TParModel([[10, 20, 30], [40, 50, 60]], [1, 2, 3])
        assert len(t) == 1
        assert t.model_set_axis is False
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(
            t.param_sets[0]
            == [
                [[10, 20, 30], [40, 50, 60]],
            ]
        )
        assert np.all(t.param_sets[1] == [[1, 2, 3]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 1, 2, 3])

    @pytest.mark.parametrize(
        ("p1", "p2"),
        [
            (1, 2),
            (1, [2, 3]),
            ([1, 2], 3),
            ([1, 2, 3], [4, 5]),
            ([1, 2], [3, 4, 5]),
        ],
    )
    def test_two_model_incorrect_scalar_parameters(self, p1, p2):
        with pytest.raises(InputParameterError, match=r".*"):
            TParModel(p1, p2, n_models=2)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"n_models": 2},
            {"model_set_axis": 0},
            {"n_models": 2, "model_set_axis": 0},
        ],
    )
    def test_two_model_scalar_parameters(self, kwargs):
        t = TParModel([10, 20], [1, 2], **kwargs)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(t.param_sets == [[10, 20], [1, 2]])
        assert np.all(t.parameters == [10, 20, 1, 2])
        assert t.coeff.shape == (2,)
        assert t.e.shape == (2,)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"n_models": 2},
            {"model_set_axis": 0},
            {"n_models": 2, "model_set_axis": 0},
        ],
    )
    def test_two_model_scalar_and_array_parameters(self, kwargs):
        t = TParModel([10, 20], [[1, 2], [3, 4]], **kwargs)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[10], [20]])
        assert np.all(t.param_sets[1] == [[1, 2], [3, 4]])
        assert np.all(t.parameters == [10, 20, 1, 2, 3, 4])
        assert t.coeff.shape == (2,)
        assert t.e.shape == (2, 2)

    def test_two_model_1d_array_parameters(self):
        t = TParModel([[10, 20], [30, 40]], [[1, 2], [3, 4]], n_models=2)
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(
            t.param_sets
            == [
                [[10, 20], [30, 40]],
                [[1, 2], [3, 4]],
            ]
        )
        assert np.all(t.parameters == [10, 20, 30, 40, 1, 2, 3, 4])
        assert t.coeff.shape == (2, 2)
        assert t.e.shape == (2, 2)

        t2 = TParModel([[10, 20, 30], [40, 50, 60]], [[1, 2, 3], [4, 5, 6]], n_models=2)
        assert len(t2) == 2
        assert t2.model_set_axis == 0
        assert np.all(
            t2.param_sets
            == [
                [[10, 20, 30], [40, 50, 60]],
                [[1, 2, 3], [4, 5, 6]],
            ]
        )
        assert np.all(t2.parameters == [10, 20, 30, 40, 50, 60, 1, 2, 3, 4, 5, 6])
        assert t2.coeff.shape == (2, 3)
        assert t2.e.shape == (2, 3)

    def test_two_model_mixed_dimension_array_parameters(self):
        MESSAGE = (
            r"Parameter .* of shape .* cannot be broadcast with parameter .* of"
            r" shape .*"
        )
        with pytest.raises(InputParameterError, match=MESSAGE):
            # Can't broadcast different array shapes
            TParModel(
                [[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
                [[9, 10, 11], [12, 13, 14]],
                n_models=2,
            )

        t = TParModel(
            [[[10, 20], [30, 40]], [[50, 60], [70, 80]]], [[1, 2], [3, 4]], n_models=2
        )
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(t.param_sets[0] == [[[10, 20], [30, 40]], [[50, 60], [70, 80]]])
        assert np.all(t.param_sets[1] == [[[1, 2]], [[3, 4]]])
        assert np.all(t.parameters == [10, 20, 30, 40, 50, 60, 70, 80, 1, 2, 3, 4])
        assert t.coeff.shape == (2, 2, 2)
        assert t.e.shape == (2, 2)

    def test_two_model_2d_array_parameters(self):
        t = TParModel(
            [[[10, 20], [30, 40]], [[50, 60], [70, 80]]],
            [[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
            n_models=2,
        )
        assert len(t) == 2
        assert t.model_set_axis == 0
        assert np.all(
            t.param_sets
            == [
                [[[10, 20], [30, 40]], [[50, 60], [70, 80]]],
                [[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
            ]
        )
        assert np.all(
            t.parameters == [10, 20, 30, 40, 50, 60, 70, 80, 1, 2, 3, 4, 5, 6, 7, 8]
        )
        assert t.coeff.shape == (2, 2, 2)
        assert t.e.shape == (2, 2, 2)

    def test_two_model_nonzero_model_set_axis(self):
        # An example where the model set axis is the *last* axis of the
        # parameter arrays
        coeff = np.array([[[10, 20, 30], [30, 40, 50]], [[50, 60, 70], [70, 80, 90]]])
        coeff = np.rollaxis(coeff, 0, 3)
        e = np.array([[1, 2, 3], [3, 4, 5]])
        e = np.rollaxis(e, 0, 2)
        t = TParModel(coeff, e, n_models=2, model_set_axis=-1)
        assert len(t) == 2
        assert t.model_set_axis == -1
        assert len(t.param_sets) == 2
        assert np.issubdtype(t.param_sets.dtype, np.object_)
        assert np.all(
            t.param_sets[0]
            == [
                [[10, 50], [20, 60], [30, 70]],
                [[30, 70], [40, 80], [50, 90]],
            ]
        )
        assert np.all(t.param_sets[1] == [[[1, 3], [2, 4], [3, 5]]])
        assert np.all(
            t.parameters
            == [10, 50, 20, 60, 30, 70, 30, 70, 40, 80, 50, 90, 1, 3, 2, 4, 3, 5]
        )
        assert t.coeff.shape == (2, 3, 2)  # note change in api
        assert t.e.shape == (3, 2)  # note change in api

    def test_wrong_number_of_params(self):
        MESSAGE = r"Inconsistent dimensions for parameter .* for 2 model sets.*"
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff=[[1, 2], [3, 4]], e=(2, 3, 4), n_models=2)
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff=[[1, 2], [3, 4]], e=(2, 3, 4), model_set_axis=0)

    def test_wrong_number_of_params2(self):
        MESSAGE = r"All parameter values must be arrays of dimension at .*"
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff=[[1, 2], [3, 4]], e=4, n_models=2)
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(coeff=[[1, 2], [3, 4]], e=4, model_set_axis=0)

    def test_array_parameter1(self):
        MESSAGE = r"All parameter values must be arrays of dimension at .*"
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(np.array([[1, 2], [3, 4]]), 1, model_set_axis=0)

    def test_array_parameter2(self):
        MESSAGE = r"Inconsistent dimensions for parameter .* for 2 model sets.*"
        with pytest.raises(InputParameterError, match=MESSAGE):
            TParModel(np.array([[1, 2], [3, 4]]), (1, 1, 11), model_set_axis=0)

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
    MESSAGE = (
        r"Parameter '.*' of shape .* cannot be broadcast with parameter '.*' of"
        r" shape .*"
    )
    for args in itertools.permutations((a, b, c)):
        with pytest.raises(InputParameterError, match=MESSAGE):
            TestModel(*args)


def test_setter():
    pars = np.random.rand(20).reshape((10, 2))

    model = SetterModel(xc=-1, yc=3, p=np.pi)

    for x, y in pars:
        np.testing.assert_almost_equal(model(x, y), (x + 1) ** 2 + (y - np.pi * 3) ** 2)
