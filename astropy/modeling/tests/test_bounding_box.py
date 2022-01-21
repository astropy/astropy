# Licensed under a 3-clause BSD style license - see LICENSE.rst

import unittest.mock as mk

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SpectralCoord
from astropy.modeling.bounding_box import (CompoundBoundingBox, ModelBoundingBox, _BaseInterval,
                                           _BaseSelectorArgument, _BoundingDomain,
                                           _ignored_interval, _Interval, _SelectorArgument,
                                           _SelectorArguments)
from astropy.modeling.core import Model, fix_inputs
from astropy.modeling.models import Gaussian1D, Gaussian2D, Identity, Scale, Shift


class Test_Interval:
    def test_create(self):
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = _Interval(lower, upper)
        assert isinstance(interval, _BaseInterval)
        assert interval.lower == lower
        assert interval.upper == upper
        assert interval == (lower, upper)

        assert interval.__repr__() == \
            f"Interval(lower={lower}, upper={upper})"

    def test_copy(self):
        interval = _Interval(0.5, 1.5)
        copy = interval.copy()

        assert interval == copy
        assert id(interval) != id(copy)

        # Same float values have will have same id
        assert interval.lower == copy.lower
        assert id(interval.lower) == id(copy.lower)

        # Same float values have will have same id
        assert interval.upper == copy.upper
        assert id(interval.upper) == id(copy.upper)

    def test__validate_shape(self):
        message = "An interval must be some sort of sequence of length 2"
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = _Interval(lower, upper)

        # Passes (2,)
        interval._validate_shape((1, 2))
        interval._validate_shape([1, 2])
        interval._validate_shape((1*u.m, 2*u.m))
        interval._validate_shape([1*u.m, 2*u.m])

        # Passes (1, 2)
        interval._validate_shape(((1, 2),))
        interval._validate_shape(([1, 2],))
        interval._validate_shape([(1, 2)])
        interval._validate_shape([[1, 2]])
        interval._validate_shape(((1*u.m, 2*u.m),))
        interval._validate_shape(([1*u.m, 2*u.m],))
        interval._validate_shape([(1*u.m, 2*u.m)])
        interval._validate_shape([[1*u.m, 2*u.m]])

        # Passes (2, 0)
        interval._validate_shape((mk.MagicMock(), mk.MagicMock()))
        interval._validate_shape([mk.MagicMock(), mk.MagicMock()])

        # Passes with array inputs:
        interval._validate_shape((np.array([-2.5, -3.5]), np.array([2.5, 3.5])))
        interval._validate_shape((np.array([-2.5, -3.5, -4.5]),
                                  np.array([2.5, 3.5, 4.5])))

        # Fails shape (no units)
        with pytest.raises(ValueError) as err:
            interval._validate_shape((1, 2, 3))
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape([1, 2, 3])
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape([[1, 2, 3], [4, 5, 6]])
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape(1)
        assert str(err.value) == message

        # Fails shape (units)
        message = "An interval must be some sort of sequence of length 2"
        with pytest.raises(ValueError) as err:
            interval._validate_shape((1*u.m, 2*u.m, 3*u.m))
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape([1*u.m, 2*u.m, 3*u.m])
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape([[1*u.m, 2*u.m, 3*u.m], [4*u.m, 5*u.m, 6*u.m]])
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape(1*u.m)
        assert str(err.value) == message

        # Fails shape (arrays):
        with pytest.raises(ValueError) as err:
            interval._validate_shape((np.array([-2.5, -3.5]),
                                      np.array([2.5, 3.5]),
                                      np.array([3, 4])))
        assert str(err.value) == message
        with pytest.raises(ValueError) as err:
            interval._validate_shape((np.array([-2.5, -3.5]), [2.5, 3.5]))
        assert str(err.value) == message

    def test__validate_bounds(self):
        # Passes
        assert _Interval._validate_bounds(1, 2) == (1, 2)
        assert _Interval._validate_bounds(1*u.m, 2*u.m) == (1*u.m, 2*u.m)

        interval = _Interval._validate_bounds(np.array([-2.5, -3.5]), np.array([2.5, 3.5]))
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()

        # Fails
        with pytest.warns(RuntimeWarning,
                          match="Invalid interval: upper bound 1 is strictly "
                          r"less than lower bound 2\."):
            _Interval._validate_bounds(2, 1)
        with pytest.warns(RuntimeWarning,
                          match=r"Invalid interval: upper bound 1\.0 m is strictly "
                          r"less than lower bound 2\.0 m\."):
            _Interval._validate_bounds(2*u.m, 1*u.m)

    def test_validate(self):
        # Passes
        assert _Interval.validate((1, 2)) == (1, 2)
        assert _Interval.validate([1, 2]) == (1, 2)
        assert _Interval.validate((1*u.m, 2*u.m)) == (1*u.m, 2*u.m)
        assert _Interval.validate([1*u.m, 2*u.m]) == (1*u.m, 2*u.m)

        assert _Interval.validate(((1, 2),)) == (1, 2)
        assert _Interval.validate(([1, 2],)) == (1, 2)
        assert _Interval.validate([(1, 2)]) == (1, 2)
        assert _Interval.validate([[1, 2]]) == (1, 2)
        assert _Interval.validate(((1*u.m, 2*u.m),)) == (1*u.m, 2*u.m)
        assert _Interval.validate(([1*u.m, 2*u.m],)) == (1*u.m, 2*u.m)
        assert _Interval.validate([(1*u.m, 2*u.m)]) == (1*u.m, 2*u.m)
        assert _Interval.validate([[1*u.m, 2*u.m]]) == (1*u.m, 2*u.m)

        interval = _Interval.validate((np.array([-2.5, -3.5]),
                                      np.array([2.5, 3.5])))
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()
        interval = _Interval.validate((np.array([-2.5, -3.5, -4.5]),
                                     np.array([2.5, 3.5, 4.5])))
        assert (interval.lower == np.array([-2.5, -3.5, -4.5])).all()
        assert (interval.upper == np.array([2.5, 3.5, 4.5])).all()

        # Fail shape
        with pytest.raises(ValueError):
            _Interval.validate((1, 2, 3))

        # Fail bounds
        with pytest.warns(RuntimeWarning):
            _Interval.validate((2, 1))

    def test_outside(self):
        interval = _Interval.validate((0, 1))

        assert (interval.outside(np.linspace(-1, 2, 13)) ==
                [True, True, True, True,
                 False, False, False, False, False,
                 True, True, True, True]).all()

    def test_domain(self):
        interval = _Interval.validate((0, 1))
        assert (interval.domain(0.25) == np.linspace(0, 1, 5)).all()

    def test__ignored_interval(self):
        assert _ignored_interval.lower == -np.inf
        assert _ignored_interval.upper == np.inf

        for num in [0, -1, -100, 3.14, 10**100, -10**100]:
            assert not num < _ignored_interval[0]
            assert num > _ignored_interval[0]

            assert not num > _ignored_interval[1]
            assert num < _ignored_interval[1]

            assert not (_ignored_interval.outside(np.array([num]))).all()

    def test_validate_with_SpectralCoord(self):
        """Regression test for issue #12439"""

        lower = SpectralCoord(1, u.um)
        upper = SpectralCoord(10, u.um)

        interval = _Interval.validate((lower, upper))
        assert interval.lower == lower
        assert interval.upper == upper


class Test_BoundingDomain:
    def setup(self):
        class BoundingDomain(_BoundingDomain):
            def _verify(self, _external_ignored = None):
                pass

            def fix_inputs(self, model, fix_inputs):
                pass

            def prepare_inputs(self, input_shape, inputs, ignored=[]):
                pass

        self.BoundingDomain = BoundingDomain

    def test_create(self):
        model = mk.MagicMock()
        model.inputs = ['x']

        # Defaults
        bounding_box = self.BoundingDomain()
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # make sure verify is called when using model
        with mk.patch.object(self.BoundingDomain, 'verify',
                             autospec=True) as mkVerify:
            bounding_box = self.BoundingDomain()
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

        # Only model
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # make sure verify is called when using model
        with mk.patch.object(self.BoundingDomain, 'verify',
                             autospec=True) as mkVerify:
            bounding_box = self.BoundingDomain(model)
            assert mkVerify.call_args_list == [mk.call(bounding_box, model)]

        # Model and ignored
        bounding_box = self.BoundingDomain(model, ['x'])
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'C'
        bounding_box = self.BoundingDomain(model, [0])
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'C'

        # Model and order
        bounding_box = self.BoundingDomain(model, order='F')
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'F'

        # Model, ignored, and order
        bounding_box = self.BoundingDomain(model, ['x'], 'F')
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'F'
        bounding_box = self.BoundingDomain(model, [0], 'F')
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'F'

        # Only ignored
        bounding_box = self.BoundingDomain(ignored=['x'])
        assert bounding_box._model is None
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'C'
        bounding_box = self.BoundingDomain(ignored=[0])
        assert bounding_box._model is None
        assert bounding_box._ignored == [0]
        assert bounding_box._order == 'C'

        # Ignored and order
        bounding_box = self.BoundingDomain(ignored=['x'], order='F')
        assert bounding_box._model is None
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'F'
        bounding_box = self.BoundingDomain(ignored=[0], order='F')
        assert bounding_box._model is None
        assert bounding_box._ignored == [0]
        assert bounding_box._order == 'F'

        # Only order
        bounding_box = self.BoundingDomain(order='F')
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'F'

        # Error from bad order
        with pytest.raises(ValueError, match=r"order must be*"):
            self.BoundingDomain(order=mk.MagicMock())

        # Error from bad ignored for model
        with pytest.raises(ValueError, match=r"Key value:*"):
            self.BoundingDomain(model, [mk.MagicMock()])
        with pytest.raises(ValueError, match=r"'*' is not *"):
            self.BoundingDomain(model, ['y'])
        with pytest.raises(IndexError, match=r"Integer key:*"):
            self.BoundingDomain(model, [1])

    def test_model(self):
        model = mk.MagicMock()

        # Test get with no error
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        assert bounding_box.model == model

        # Test get with error
        bounding_box = self.BoundingDomain()
        assert bounding_box._model is None
        with pytest.raises(RuntimeError) as err:
            bounding_box.model
        assert str(err.value) == \
            "Method requires a model to function, please attach to a model"

        # Test set
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        with mk.patch.object(self.BoundingDomain, 'verify',
                             autospec=True) as mkVerify:
            # None model
            bounding_box.model = None
            assert bounding_box._model == model
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

            mkVerify.reset_mock()

            # Model
            bounding_box.model = model
            assert bounding_box._model == model
            assert bounding_box.model == model
            assert mkVerify.call_args_list == [mk.call(bounding_box, model)]

    def test__has_model(self):
        bounding_box = self.BoundingDomain()
        assert bounding_box._model is None
        assert bounding_box._has_model is False

        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        assert bounding_box._has_model is True

    def test_order(self):
        bounding_box = self.BoundingDomain(mk.MagicMock(), order='C')
        assert bounding_box._order == 'C'
        assert bounding_box.order == 'C'

        bounding_box = self.BoundingDomain(mk.MagicMock(), order='F')
        assert bounding_box._order == 'F'
        assert bounding_box.order == 'F'

        bounding_box._order = 'test'
        assert bounding_box.order == 'test'

    def test_ignored(self):
        ignored = ['x']
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = self.BoundingDomain(model, ignored=ignored)

        assert bounding_box._ignored == ignored
        assert bounding_box.ignored == ignored

    def test__get_order(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # Success (default 'C')
        assert bounding_box._order == 'C'
        assert bounding_box._get_order() == 'C'
        assert bounding_box._get_order('C') == 'C'
        assert bounding_box._get_order('F') == 'F'

        # Success (default 'F')
        bounding_box._order = 'F'
        assert bounding_box._order == 'F'
        assert bounding_box._get_order() == 'F'
        assert bounding_box._get_order('C') == 'C'
        assert bounding_box._get_order('F') == 'F'

        # Error
        order = mk.MagicMock()
        with pytest.raises(ValueError) as err:
            bounding_box._get_order(order)
        assert str(err.value) ==\
            "order must be either 'C' (C/python order) or " +\
            f"'F' (Fortran/mathematical order), got: {order}."

    def test__get_index(self):
        bounding_box = self.BoundingDomain(Gaussian2D())

        # Pass input name
        assert bounding_box._get_index('x') == 0
        assert bounding_box._get_index('y') == 1

        # Pass invalid input name
        with pytest.raises(ValueError) as err:
            bounding_box._get_index('z')
        assert str(err.value) ==\
            "'z' is not one of the inputs: ('x', 'y')."

        # Pass valid index
        assert bounding_box._get_index(0) == 0
        assert bounding_box._get_index(1) == 1
        assert bounding_box._get_index(np.int32(0)) == 0
        assert bounding_box._get_index(np.int32(1)) == 1
        assert bounding_box._get_index(np.int64(0)) == 0
        assert bounding_box._get_index(np.int64(1)) == 1

        # Pass invalid index
        MESSAGE = "Integer key: 2 must be non-negative and < 2."
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(2)
        assert str(err.value) == MESSAGE
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(np.int32(2))
        assert str(err.value) == MESSAGE
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(np.int64(2))
        assert str(err.value) == MESSAGE
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(-1)
        assert str(err.value) ==\
            "Integer key: -1 must be non-negative and < 2."

        # Pass invalid key
        value = mk.MagicMock()
        with pytest.raises(ValueError) as err:
            bounding_box._get_index(value)
        assert str(err.value) ==\
            f"Key value: {value} must be string or integer."

    def test__get_name(self):
        import astropy.modeling.bounding_box as bbox
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = self.BoundingDomain(model)

        key = mk.MagicMock()
        index = mk.MagicMock()
        name = mk.MagicMock()
        model.inputs = mk.MagicMock()
        model.inputs.__getitem__.return_value = name

        with mk.patch('astropy.modeling.bounding_box.get_index', attribute=True,
                      return_value=index) as mkIndex:
            assert bounding_box._get_name(key) == name
            assert mkIndex.call_args_list == [mk.call(model, key)]
            assert model.inputs.__getitem__.call_args_list == [mk.call(index)]

    def test_ignored_inputs(self):
        model = mk.MagicMock()
        ignored = [f"x{index}" for index in range(4, 8)]
        model.n_inputs = 8
        model.inputs = [f"x{index}" for index in range(8)]
        bounding_box = self.BoundingDomain(model, ignored=ignored)

        ignored_inputs = bounding_box.ignored_inputs
        assert isinstance(ignored_inputs, list)

        # Check each ignored input is valid
        for index, _input in enumerate(ignored_inputs):
            assert 0 < _input <= len(model.inputs)
            assert ignored[index] == model.inputs[_input]

        # Check that each input is sorted correctly
        for index, _input in enumerate(model.inputs):
            if _input in ignored:
                assert index >= 4
                assert index in ignored_inputs
            else:
                assert index < 4
                assert index not in ignored_inputs

    def test_verify(self):
        model = mk.MagicMock()
        model.inputs = ['x', 'y']
        external_ignored = mk.MagicMock()

        # Pass, No ignored, default args
        bounding_box = self.BoundingDomain()
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        with mk.patch.object(self.BoundingDomain, '_verify',
                             autospec=True) as mkVerify:
            bounding_box.verify(None)
            assert bounding_box._model is None
            assert bounding_box._ignored == []
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]
            mkVerify.reset_mock()
            bounding_box.verify(model)
            assert bounding_box._model == model
            assert bounding_box._ignored == []
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

        # Pass, No ignored, non-default args
        bounding_box = self.BoundingDomain()
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        with mk.patch.object(self.BoundingDomain, '_verify',
                             autospec=True) as mkVerify:
            bounding_box.verify(None, external_ignored)
            assert bounding_box._model is None
            assert bounding_box._ignored == []
            assert mkVerify.call_args_list == [mk.call(bounding_box, external_ignored)]
            mkVerify.reset_mock()
            bounding_box.verify(model, external_ignored)
            assert bounding_box._model == model
            assert bounding_box._ignored == []
            assert mkVerify.call_args_list == [mk.call(bounding_box, external_ignored)]

        # Pass, str ignored, default args
        bounding_box = self.BoundingDomain(ignored=['x', 'y'])
        assert bounding_box._model is None
        assert bounding_box._ignored == ['x', 'y']
        with mk.patch.object(self.BoundingDomain, '_verify',
                             autospec=True) as mkVerify:
            bounding_box.verify(None)
            assert bounding_box._model is None
            assert bounding_box._ignored == ['x', 'y']
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]
            mkVerify.reset_mock()
            bounding_box.verify(model)
            assert bounding_box._model == model
            assert bounding_box._ignored == ['x', 'y']
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

        # Pass, integer ignored, default args
        bounding_box = self.BoundingDomain(ignored=[0, 1])
        assert bounding_box._model is None
        assert bounding_box._ignored == [0, 1]
        with mk.patch.object(self.BoundingDomain, '_verify',
                             autospec=True) as mkVerify:
            bounding_box.verify(None)
            assert bounding_box._model is None
            assert bounding_box._ignored == [0, 1]
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]
            mkVerify.reset_mock()
            bounding_box.verify(model)
            assert bounding_box._model == model
            assert bounding_box._ignored == ['x', 'y']
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

        # Pass, numpy integer ignored, default args
        bounding_box = self.BoundingDomain(ignored=[np.int32(0), np.int64(1)])
        assert bounding_box._model is None
        assert bounding_box._ignored == [np.int32(0), np.int64(1)]
        with mk.patch.object(self.BoundingDomain, '_verify',
                             autospec=True) as mkVerify:
            bounding_box.verify(None)
            assert bounding_box._model is None
            assert bounding_box._ignored == [np.int32(0), np.int64(1)]
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]
            mkVerify.reset_mock()
            bounding_box.verify(model)
            assert bounding_box._model == model
            assert bounding_box._ignored == ['x', 'y']
            assert mkVerify.call_args_list == [mk.call(bounding_box, None)]

        # Fail, arbitrary value
        bounding_box = self.BoundingDomain(ignored=[mk.MagicMock()])
        with pytest.raises(ValueError):
            bounding_box.verify(model)

        # Fail, non input str
        bounding_box = self.BoundingDomain(ignored=['z'])
        with pytest.raises(ValueError):
            bounding_box.verify(model)

        # Fail, non input integer index
        bounding_box = self.BoundingDomain(ignored=[3])
        with pytest.raises(IndexError):
            bounding_box.verify(model)
        bounding_box = self.BoundingDomain(ignored=[np.int32(3)])
        with pytest.raises(IndexError):
            bounding_box.verify(model)
        bounding_box = self.BoundingDomain(ignored=[np.int64(3)])
        with pytest.raises(IndexError):
            bounding_box.verify(model)

    def test___call__(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        args = tuple([mk.MagicMock() for _ in range(3)])
        kwargs = {f"test{idx}": mk.MagicMock() for idx in range(3)}

        with pytest.raises(RuntimeError) as err:
            bounding_box(*args, **kwargs)
        assert str(err.value) ==\
            "This bounding box is fixed by the model and does not have " +\
            "adjustable parameters."

    def test_fix_inputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())
        model = mk.MagicMock()
        fixed_inputs = mk.MagicMock()
        bounding_box.fix_inputs(model, fixed_inputs)

    def test_prepare_inputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        bounding_box.prepare_inputs(mk.MagicMock(), mk.MagicMock())
        bounding_box.prepare_inputs(mk.MagicMock(), mk.MagicMock(), mk.MagicMock())

    def test__base_ouput(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # Simple shape
        input_shape = (13,)
        output = bounding_box._base_output(input_shape, 0)
        assert (output == 0).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, np.nan)
        assert (np.isnan(output)).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, 14)
        assert (output == 14).all()
        assert output.shape == input_shape

        # Complex shape
        input_shape = (13, 7)
        output = bounding_box._base_output(input_shape, 0)
        assert (output == 0).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, np.nan)
        assert (np.isnan(output)).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, 14)
        assert (output == 14).all()
        assert output.shape == input_shape

    def test__all_out_output(self):
        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)

        # Simple shape
        model.n_outputs = 1
        input_shape = (13,)
        output, output_unit = bounding_box._all_out_output(input_shape, 0)
        assert (np.array(output) == 0).all()
        assert np.array(output).shape == (1, 13)
        assert output_unit is None

        # Complex shape
        model.n_outputs = 6
        input_shape = (13, 7)
        output, output_unit = bounding_box._all_out_output(input_shape, 0)
        assert (np.array(output) == 0).all()
        assert np.array(output).shape == (6, 13, 7)
        assert output_unit is None

    def test__modify_output(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        # Simple shape
        with mk.patch.object(_BoundingDomain, '_base_output', autospec=True,
                             return_value=np.asanyarray(0)) as mkBase:
            assert (np.array([1, 2, 3]) ==
                    bounding_box._modify_output([1, 2, 3], valid_index, input_shape, fill_value)).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

        # Replacement
        with mk.patch.object(_BoundingDomain, '_base_output', autospec=True,
                             return_value=np.array([1, 2, 3, 4, 5, 6])) as mkBase:
            assert (np.array([7, 2, 8, 4, 9, 6]) ==
                    bounding_box._modify_output([7, 8, 9], np.array([[0, 2, 4]]), input_shape, fill_value)).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

    def test__prepare_outputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        valid_outputs = [mk.MagicMock() for _ in range(3)]
        effects = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(_BoundingDomain, '_modify_output', autospec=True,
                             side_effect=effects) as mkModify:
            assert effects == bounding_box._prepare_outputs(valid_outputs, valid_index,
                                                            input_shape, fill_value)
            assert mkModify.call_args_list == \
                [mk.call(bounding_box, valid_outputs[idx], valid_index, input_shape, fill_value)
                 for idx in range(3)]

    def test_prepare_outputs(self):
        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)

        valid_outputs = mk.MagicMock()
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        with mk.patch.object(_BoundingDomain, '_prepare_outputs', autospec=True) as mkPrepare:
            # Reshape valid_outputs
            model.n_outputs = 1
            assert mkPrepare.return_value == \
                bounding_box.prepare_outputs(valid_outputs, valid_index, input_shape, fill_value)
            assert mkPrepare.call_args_list == \
                [mk.call(bounding_box, [valid_outputs], valid_index, input_shape, fill_value)]
            mkPrepare.reset_mock()

            # No reshape valid_outputs
            model.n_outputs = 2
            assert mkPrepare.return_value == \
                bounding_box.prepare_outputs(valid_outputs, valid_index, input_shape, fill_value)
            assert mkPrepare.call_args_list == \
                [mk.call(bounding_box, valid_outputs, valid_index, input_shape, fill_value)]

    def test__get_valid_outputs_unit(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # Don't get unit
        assert bounding_box._get_valid_outputs_unit(mk.MagicMock(), False) is None

        # Get unit from unitless
        assert bounding_box._get_valid_outputs_unit(7, True) is None

        # Get unit
        assert bounding_box._get_valid_outputs_unit(25 * u.m, True) == u.m

    def test__evaluate_model(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        evaluate = mk.MagicMock()
        valid_inputs = mk.MagicMock()
        input_shape = mk.MagicMock()
        valid_index = mk.MagicMock()
        fill_value = mk.MagicMock()
        with_units = mk.MagicMock()

        with mk.patch.object(_BoundingDomain, '_get_valid_outputs_unit',
                             autospec=True) as mkGet:
            with mk.patch.object(_BoundingDomain, 'prepare_outputs',
                                 autospec=True) as mkPrepare:
                assert bounding_box._evaluate_model(evaluate, valid_inputs,
                                                    valid_index, input_shape,
                                                    fill_value, with_units) == \
                    (mkPrepare.return_value, mkGet.return_value)
                assert mkPrepare.call_args_list == \
                    [mk.call(bounding_box, evaluate.return_value, valid_index,
                             input_shape, fill_value)]
                assert mkGet.call_args_list == \
                    [mk.call(evaluate.return_value, with_units)]
                assert evaluate.call_args_list == \
                    [mk.call(valid_inputs)]

    def test__evaluate(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        evaluate = mk.MagicMock()
        inputs = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()
        with_units = mk.MagicMock()

        valid_inputs = mk.MagicMock()
        valid_index = mk.MagicMock()

        effects = [(valid_inputs, valid_index, True), (valid_inputs, valid_index, False)]
        with mk.patch.object(self.BoundingDomain, 'prepare_inputs', autospec=True,
                             side_effect=effects) as mkPrepare:
            with mk.patch.object(_BoundingDomain, '_all_out_output',
                                 autospec=True) as mkAll:
                with mk.patch.object(_BoundingDomain, '_evaluate_model',
                                     autospec=True) as mkEvaluate:
                    # all_out
                    assert bounding_box._evaluate(evaluate, inputs, input_shape,
                                                  fill_value, with_units) == \
                        mkAll.return_value
                    assert mkAll.call_args_list == \
                        [mk.call(bounding_box, input_shape, fill_value)]
                    assert mkEvaluate.call_args_list == []
                    assert mkPrepare.call_args_list == \
                        [mk.call(bounding_box, input_shape, inputs)]

                    mkAll.reset_mock()
                    mkPrepare.reset_mock()

                    # not all_out
                    assert bounding_box._evaluate(evaluate, inputs, input_shape,
                                                  fill_value, with_units) == \
                        mkEvaluate.return_value
                    assert mkAll.call_args_list == []
                    assert mkEvaluate.call_args_list == \
                        [mk.call(bounding_box, evaluate, valid_inputs, valid_index,
                                 input_shape, fill_value, with_units)]
                    assert mkPrepare.call_args_list == \
                        [mk.call(bounding_box, input_shape, inputs)]

    def test__set_outputs_unit(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # set no unit
        assert 27 == bounding_box._set_outputs_unit(27, None)

        # set unit
        assert 27 * u.m == bounding_box._set_outputs_unit(27, u.m)

    def test_evaluate(self):
        bounding_box = self.BoundingDomain(Gaussian2D())

        evaluate = mk.MagicMock()
        inputs = mk.MagicMock()
        fill_value = mk.MagicMock()

        outputs = mk.MagicMock()
        valid_outputs_unit = mk.MagicMock()
        value = (outputs, valid_outputs_unit)
        with mk.patch.object(_BoundingDomain, '_evaluate',
                             autospec=True, return_value=value) as mkEvaluate:
            with mk.patch.object(_BoundingDomain, '_set_outputs_unit',
                                 autospec=True) as mkSet:
                with mk.patch.object(Model, 'input_shape', autospec=True) as mkShape:
                    with mk.patch.object(Model, 'bbox_with_units',
                                         new_callable=mk.PropertyMock) as mkUnits:
                        assert tuple(mkSet.return_value) == \
                            bounding_box.evaluate(evaluate, inputs, fill_value)
                        assert mkSet.call_args_list == \
                            [mk.call(outputs, valid_outputs_unit)]
                        assert mkEvaluate.call_args_list == \
                            [mk.call(bounding_box, evaluate, inputs, mkShape.return_value,
                                     fill_value, mkUnits.return_value)]
                        assert mkShape.call_args_list == \
                            [mk.call(bounding_box._model, inputs)]
                        assert mkUnits.call_args_list == [mk.call()]


class TestModelBoundingBox:
    def test_create(self):
        # Intervals only, 0D-tuple
        bounding_box = ModelBoundingBox(())
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 0D-dictionary
        bounding_box = ModelBoundingBox({})
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 1D-tuple
        bounding_box = ModelBoundingBox((0, 1))
        assert bounding_box._intervals == (0, 1)
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 1D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1)})
        assert bounding_box._intervals == {'x': (0, 1)}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 2D-tuple
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)))
        assert bounding_box._intervals == ((0, 1), (2, 3))
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 2D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1), 'y': (2, 3)})
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3)}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 3D-tuple
        bounding_box = ModelBoundingBox(((0, 1), (2, 3), (4, 5)))
        assert bounding_box._intervals == ((0, 1), (2, 3), (4, 5))
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, 3D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1), 'y': (2, 3), 'z': (4, 5)})
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3), 'z': (4, 5)}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals and ignored
        bounding_box = ModelBoundingBox({}, ignored=['a', 'b'])
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == ['a', 'b']
        assert bounding_box._order == 'C'
        bounding_box = ModelBoundingBox({}, ignored=[1, 2])
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == [1, 2]
        assert bounding_box._order == 'C'

        # Intervals and order
        bounding_box = ModelBoundingBox({}, order='F')
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'F'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), order='F')
        assert bounding_box._intervals == ((0, 1), (2, 3))
        assert bounding_box._model is None
        assert bounding_box._ignored == []
        assert bounding_box._order == 'F'

        # Intervals, ignored, and order
        bounding_box = ModelBoundingBox({}, ignored=['a', 'b'], order='F')
        assert bounding_box._intervals == {}
        assert bounding_box._model is None
        assert bounding_box._ignored == ['a', 'b']
        assert bounding_box._order == 'F'

        # Order error
        with pytest.raises(ValueError, match=r"order must be*"):
            ModelBoundingBox({}, ignored=['a', 'b'], order=mk.MagicMock())

        model = mk.MagicMock()

        # Intervals and model, 0D-dictionary
        bounding_box = ModelBoundingBox({}, model)
        assert bounding_box._intervals == {}
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        model.inputs = ['x']

        # Intervals and model, 1D-tuple
        bounding_box = ModelBoundingBox((0, 1), model)
        assert bounding_box._intervals == {'x': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals only, bad tuple
        with pytest.raises(ValueError, match=r"The intervals:*"):
            ModelBoundingBox((1,), model)

        # Intervals and ignored Error
        with pytest.raises(ValueError, match=r"At least one*"):
            ModelBoundingBox({'x': (0, 1)}, model, ignored=['x'])

        # Intervals and model, 1D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1)}, model)
        assert bounding_box._intervals == {'x': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        model.inputs = ['x', 'y']

        # Intervals and model, 2D-tuple
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model)
        assert bounding_box._intervals == {'x': (2, 3), 'y': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals and model, 2D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1), 'y': (2, 3)}, model)
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        model.inputs = ['x', 'y', 'z']

        # Intervals and model, 3D-tuple
        bounding_box = ModelBoundingBox(((0, 1), (2, 3), (4, 5)), model)
        assert bounding_box._intervals == {'x': (4, 5), 'y': (2, 3), 'z': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Intervals and model, 3D-dictionary
        bounding_box = ModelBoundingBox({'x': (0, 1), 'y': (2, 3), 'z': (4, 5)}, model)
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3), 'z': (4, 5)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == 'C'

        # Tuple interval with ignored 1D
        model = mk.MagicMock()
        model.inputs = ['x', 'y']
        bounding_box = ModelBoundingBox((0, 1), model, ignored=['x'])
        assert bounding_box._intervals == {'y': (0, 1)}
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'C'
        bounding_box = ModelBoundingBox((0, 1), model, ignored=['y'])
        assert bounding_box._intervals == {'x': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['y']
        assert bounding_box._order == 'C'

        # Tuple interval with ignored 2D
        model.inputs = ['x', 'y', 'z']
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['x'])
        assert bounding_box._intervals == {'y': (2, 3), 'z': (0, 1)}
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'C'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['x'], order='F')
        assert bounding_box._intervals == {'y': (0, 1), 'z': (2, 3)}
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['x']
        assert bounding_box._order == 'F'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['y'])
        assert bounding_box._intervals == {'x': (2, 3), 'z': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['y']
        assert bounding_box._order == 'C'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['y'], order='F')
        assert bounding_box._intervals == {'x': (0, 1), 'z': (2, 3)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['z'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['y']
        assert bounding_box._order == 'F'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['z'])
        assert bounding_box._intervals == {'x': (2, 3), 'y': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['z']
        assert bounding_box._order == 'C'
        bounding_box = ModelBoundingBox(((0, 1), (2, 3)), model, ignored=['z'], order='F')
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)
        assert bounding_box._model == model
        assert bounding_box._ignored == ['z']
        assert bounding_box._order == 'F'

        # Intervals, model, and ignored error (Incorrect number of intervals)
        with pytest.raises(ValueError, match=r"An interval must be*"):
            ModelBoundingBox((0, 1), model, ignored=['x'])
        model.inputs = ['x', 'y', 'z', 'a']
        with pytest.raises(ValueError, match=r"Given .*, need .* intervals!"):
            ModelBoundingBox(((0, 1), (2, 3)), model, ignored=[0])
        with pytest.raises(ValueError, match=r"Given .*, need .* intervals!"):
            ModelBoundingBox({'x': (0, 1)}, model, ignored=[0])
        with pytest.raises(ValueError, match=r"Given .*, need .* intervals!"):
            ModelBoundingBox({'x': (0, 1), 'y': (2, 3)}, model, ignored=['z'])

        # Intervals, model, and ignored error (shared input name)
        model.inputs = ['x']
        with pytest.raises(ValueError, match=r"At least one*"):
            ModelBoundingBox({'x': (0, 1)}, model, ignored=[0])
        with pytest.raises(ValueError, match=r"At least one*"):
            ModelBoundingBox({0: (0, 1)}, model, ignored=['x'])

    def test__pop_intervals(self):
        intervals = {0: (1, 2)}
        bounding_box = ModelBoundingBox(intervals)
        assert bounding_box._intervals == intervals
        assert bounding_box._pop_intervals() == intervals
        assert intervals != {}
        assert bounding_box._intervals == {}

    def test__verify_intervals_dict(self):
        bounding_box = ModelBoundingBox({}, Gaussian2D())
        assert bounding_box._intervals == {}

        # Success integer keys
        intervals = {0: (0, 1), 1: (2, 3)}
        bounding_box._intervals = intervals
        assert not isinstance(bounding_box._intervals[0], _Interval)
        assert not isinstance(bounding_box._intervals[1], _Interval)
        bounding_box._verify_intervals_dict()
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)

        # Success integer keys (ignored key)
        intervals = {0: (0, 1)}
        bounding_box._intervals = intervals
        assert not isinstance(bounding_box._intervals[0], _Interval)
        bounding_box._verify_intervals_dict(['y'])
        assert bounding_box._intervals == {'x': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)

        # Fail integer keys (ignored key)
        bounding_box._intervals = {0: (0, 1), 1: (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['x'])
        bounding_box._intervals = {0: (0, 1), 1: (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['y'])
        bounding_box._intervals = {0: (0, 1)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['x'])
        bounding_box._intervals = {1: (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['y'])

        # Success str keys
        intervals = {'x': (0, 1), 'y': (2, 3)}
        bounding_box._intervals = intervals
        assert not isinstance(bounding_box._intervals['x'], _Interval)
        assert not isinstance(bounding_box._intervals['y'], _Interval)
        bounding_box._verify_intervals_dict()
        assert bounding_box._intervals == intervals
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)

        # Fail str keys (ignored key)
        bounding_box._intervals = {'x': (0, 1), 'y': (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['x'])
        bounding_box._intervals = {'x': (0, 1), 'y': (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['y'])
        bounding_box._intervals = {'x': (0, 1)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['x'])
        bounding_box._intervals = {'y': (2, 3)}
        with pytest.raises(ValueError, match=r"Interval: .* is being externally ignored!"):
            bounding_box._verify_intervals_dict(['y'])

        # Fail, bad interval
        bounding_box._intervals = {0: (0, 1), 1: (2,)}
        with pytest.raises(ValueError, match=r"An interval must*"):
            bounding_box._verify_intervals_dict()

        # Fail, bad key
        bounding_box._intervals = {0: (0, 1), 2: (2, 3)}
        with pytest.raises(IndexError, match=r"Integer key:*"):
            bounding_box._verify_intervals_dict()
        bounding_box._intervals = {0: (0, 1), 'z': (2, 3)}
        with pytest.raises(ValueError, match=r"'*' is not one*"):
            bounding_box._verify_intervals_dict()

        model = mk.MagicMock()
        model.inputs = ['x', 'y', 'z', 'a']
        bounding_box._ignored = ['z']
        bounding_box._model = model
        bounding_box._intervals = {'x': (0, 1), 'y': (2, 3)}
        with pytest.raises(ValueError, match=r"Given 2, need 3 intervals!"):
            bounding_box._verify_intervals_dict()

    def test__verify_intervals_sequence(self):
        bounding_box = ModelBoundingBox({})
        assert bounding_box._intervals == {}

        # Success 1D with model
        bounding_box._model = Gaussian1D()
        bounding_box._intervals = (0, 1)
        bounding_box._verify_intervals_sequence()
        assert bounding_box._intervals == {'x': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)

        # Success 2D with model, order 'C'
        bounding_box._order = 'C'
        bounding_box._model = Gaussian2D()
        bounding_box._intervals = ((0, 1), (2, 3))
        bounding_box._verify_intervals_sequence()
        assert bounding_box._intervals == {'x': (2, 3), 'y': (0, 1)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)

        # Success 2D with model, order 'F'
        bounding_box._order = 'F'
        bounding_box._model = Gaussian2D()
        bounding_box._intervals = ((0, 1), (2, 3))
        bounding_box._verify_intervals_sequence()
        assert bounding_box._intervals == {'x': (0, 1), 'y': (2, 3)}
        assert isinstance(bounding_box._intervals['x'], _Interval)
        assert isinstance(bounding_box._intervals['y'], _Interval)

        # Fail, bad interval
        bounding_box._intervals = ((0, 1), (2,))
        with pytest.raises(ValueError, match=r"An interval must*"):
            bounding_box._verify_intervals_sequence()

        # Fail, ignoring all intervals but passing one.
        bounding_box._model = Gaussian1D()
        bounding_box._intervals = (0, 1)
        bounding_box._ignored = ['x']
        with pytest.raises(ValueError, match="All intervals have been ignored!"):
            bounding_box._verify_intervals_sequence()
        bounding_box._model = Gaussian2D()
        bounding_box._intervals = ((0, 1), (2, 3))
        bounding_box._ignored = ['x', 'y']
        with pytest.raises(ValueError, match="All intervals have been ignored!"):
            bounding_box._verify_intervals_sequence()

        # Intervals only, tuple is not long enough
        bounding_box._intervals = (1,)
        bounding_box._ignored = []
        with pytest.raises(ValueError) as err:
            bounding_box._verify_intervals_sequence()
        assert str(err.value) ==\
            "The intervals: (1,), do not contain enough information to construct a bounding_box!"

        model = mk.MagicMock()
        model.inputs = ['x', 'y', 'z', 'a']
        bounding_box._model = model
        # Difference in available intervals and provided intervals
        bounding_box._intervals = ((0, 1), (2, 3))
        bounding_box._ignored = ['a']
        with pytest.raises(ValueError, match=r"Given 2, need 3 intervals!"):
            bounding_box._verify_intervals_sequence()
        bounding_box._intervals = ((0, 1), (2, 3), (4, 5))
        bounding_box._ignored = ['a', 'z']
        with pytest.raises(ValueError, match=r"Given 3, need 2 intervals!"):
            bounding_box._verify_intervals_sequence()

    def test__verify_intervals(self):
        bounding_box = ModelBoundingBox({}, Gaussian2D())
        assert bounding_box._intervals == {}

        with mk.patch.object(ModelBoundingBox, '_verify_intervals_dict',
                             autospec=True) as mkDict:
            with mk.patch.object(ModelBoundingBox, '_verify_intervals_sequence',
                                 autospec=True) as mkSeq:
                # Verify dictionary intervals
                bounding_box._verify_intervals()
                assert mkDict.call_args_list == [mk.call(bounding_box, None)]
                assert mkSeq.call_args_list == []

                mkDict.reset_mock()

                # Verify tuple intervals
                bounding_box._intervals = ((0, 1), (2, 3))
                bounding_box._verify_intervals()
                assert mkDict.call_args_list == []
                assert mkSeq.call_args_list == [mk.call(bounding_box, None)]

        # Test error from ignored and intervals sharing keys
        bounding_box._intervals = {'x': (0, 1)}
        bounding_box._ignored = ['x']
        bounding_box._model = Gaussian1D()
        with pytest.raises(ValueError, match="At least one interval is being ignored"):
            bounding_box._verify_intervals()

    def test__verify(self):
        external_ignored = mk.MagicMock()

        bounding_box_0 = ModelBoundingBox({})
        bounding_box_1 = ModelBoundingBox({}, mk.MagicMock())
        with mk.patch.object(ModelBoundingBox, '_verify_intervals',
                             autospec=True) as mkIntervals:
            # No model, no _external_ignored
            bounding_box_0._verify()
            assert mkIntervals.call_args_list == []

            # No model, with _external_ignored
            bounding_box_0._verify(external_ignored)
            assert mkIntervals.call_args_list == []

            # With model, no _external_ignored
            bounding_box_1._verify()
            assert mkIntervals.call_args_list == \
                [mk.call(bounding_box_1, None)]

            mkIntervals.reset_mock()

            # With model, with _external_ignored
            bounding_box_1._verify(external_ignored)
            assert mkIntervals.call_args_list == \
                [mk.call(bounding_box_1, external_ignored)]

    def test_copy(self):
        bounding_box = ModelBoundingBox.validate(Gaussian2D(), ((-4.5, 4.5), (-1.4, 1.4)))
        copy = bounding_box.copy()

        assert bounding_box == copy
        assert id(bounding_box) != id(copy)

        assert bounding_box.ignored == copy.ignored
        assert id(bounding_box.ignored) != id(copy.ignored)

        # model is not copied to prevent infinite recursion
        assert bounding_box._model == copy._model
        assert id(bounding_box._model) == id(copy._model)

        # Same string values have will have same id
        assert bounding_box._order == copy._order
        assert id(bounding_box._order) == id(copy._order)

        # Check interval objects
        for index, interval in bounding_box.intervals.items():
            assert interval == copy.intervals[index]
            assert id(interval) != id(copy.intervals[index])

            # Same float values have will have same id
            assert interval.lower == copy.intervals[index].lower
            assert id(interval.lower) == id(copy.intervals[index].lower)

            # Same float values have will have same id
            assert interval.upper == copy.intervals[index].upper
            assert id(interval.upper) == id(copy.intervals[index].upper)
        assert len(bounding_box.intervals) == len(copy.intervals)
        assert bounding_box.intervals.keys() == copy.intervals.keys()

    def test_intervals(self):
        intervals = {0: _Interval(1, 2)}
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = ModelBoundingBox(intervals, model)

        truth_intervals = {'x': _Interval(1, 2)}
        assert bounding_box._intervals == truth_intervals
        assert bounding_box.intervals == truth_intervals

    def test_indexed_intervals(self):
        intervals = {f"x{idx}": _Interval(idx, idx + 1) for idx in range(4)}
        model = mk.MagicMock()
        model.n_inputs = 4
        model.inputs = [f"x{idx}" for idx in range(4)]
        bounding_box = ModelBoundingBox(intervals, model)

        indexed = bounding_box.indexed_intervals
        assert isinstance(indexed, dict)

        for index, interval in indexed.items():
            assert 0 <= index <= len(model.inputs)
            assert interval == intervals[model.inputs[index]]

        for index, name in enumerate(model.inputs):
            assert intervals[name] == indexed[index]

    def test___repr__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        assert bounding_box.__repr__() ==\
            "ModelBoundingBox(\n" +\
            "    intervals={\n" +\
            "        x: Interval(lower=-1, upper=1)\n" +\
            "        y: Interval(lower=-4, upper=4)\n" +\
            "    }\n" +\
            "    model=Gaussian2D(inputs=('x', 'y'))\n" +\
            "    order='C'\n" +\
            ")"

        intervals = {0: _Interval(-1, 1)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals, ignored=['y'])

        assert bounding_box.__repr__() ==\
            "ModelBoundingBox(\n" +\
            "    intervals={\n" +\
            "        x: Interval(lower=-1, upper=1)\n" +\
            "    }\n" +\
            "    ignored=['y']\n" +\
            "    model=Gaussian2D(inputs=('x', 'y'))\n" +\
            "    order='C'\n" +\
            ")"

    def test___len__(self):
        intervals = {0: _Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert len(bounding_box) == 1 == len(bounding_box._intervals)

        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert len(bounding_box) == 2 == len(bounding_box._intervals)

        bounding_box._intervals = {}
        assert len(bounding_box) == 0 == len(bounding_box._intervals)

    def test___contains__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Contains with keys
        assert 'x' in bounding_box
        assert 'y' in bounding_box
        assert 'z' not in bounding_box

        # Contains with index
        assert 0 in bounding_box
        assert 1 in bounding_box
        assert 2 not in bounding_box

        # General not in
        assert mk.MagicMock() not in bounding_box

        # Contains with ignored
        del bounding_box['y']

        # Contains with keys
        assert 'x' in bounding_box
        assert 'y' in bounding_box
        assert 'z' not in bounding_box

        # Contains with index
        assert 0 in bounding_box
        assert 1 in bounding_box
        assert 2 not in bounding_box

    def test___getitem__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Get using input key
        assert bounding_box['x'] == (-1, 1)
        assert bounding_box['y'] == (-4, 4)

        # Fail with input key
        with pytest.raises(ValueError):
            bounding_box['z']

        # Get using index
        assert bounding_box[0] == (-1, 1)
        assert bounding_box[1] == (-4, 4)
        assert bounding_box[np.int32(0)] == (-1, 1)
        assert bounding_box[np.int32(1)] == (-4, 4)
        assert bounding_box[np.int64(0)] == (-1, 1)
        assert bounding_box[np.int64(1)] == (-4, 4)

        # Fail with index
        with pytest.raises(IndexError):
            bounding_box[2]
        with pytest.raises(IndexError):
            bounding_box[np.int32(2)]
        with pytest.raises(IndexError):
            bounding_box[np.int64(2)]

        # get ignored interval
        del bounding_box[0]
        assert bounding_box[0] == _ignored_interval
        assert bounding_box[1] == (-4, 4)

        del bounding_box[1]
        assert bounding_box[0] == _ignored_interval
        assert bounding_box[1] == _ignored_interval

    def test_bounding_box(self):
        # 0D
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model, {}, ignored=['x'])
        assert bounding_box.bounding_box() == (-np.inf, np.inf)
        assert bounding_box.bounding_box('C') == (-np.inf, np.inf)
        assert bounding_box.bounding_box('F') == (-np.inf, np.inf)

        # 1D
        intervals = {0: _Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert bounding_box.bounding_box() == (-1, 1)
        assert bounding_box.bounding_box(mk.MagicMock()) == (-1, 1)

        # > 1D
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert bounding_box.bounding_box() == ((-4, 4), (-1, 1))
        assert bounding_box.bounding_box('C') == ((-4, 4), (-1, 1))
        assert bounding_box.bounding_box('F') == ((-1, 1), (-4, 4))

    def test___eq__(self):
        intervals = {0: _Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model.copy(), intervals.copy())

        assert bounding_box == bounding_box
        assert bounding_box == ModelBoundingBox.validate(model.copy(), intervals.copy())
        assert bounding_box == (-1, 1)

        assert not (bounding_box == mk.MagicMock())
        assert not (bounding_box == (-2, 2))
        assert not (bounding_box == ModelBoundingBox.validate(model, {0: _Interval(-2, 2)}))

        # Respect ordering
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box_1 = ModelBoundingBox.validate(model, intervals)
        bounding_box_2 = ModelBoundingBox.validate(model, intervals, order='F')
        assert bounding_box_1._order == 'C'
        assert bounding_box_1 == ((-4, 4), (-1, 1))
        assert not (bounding_box_1 == ((-1, 1), (-4, 4)))

        assert bounding_box_2._order == 'F'
        assert not (bounding_box_2 == ((-4, 4), (-1, 1)))
        assert bounding_box_2 == ((-1, 1), (-4, 4))

        assert bounding_box_1 == bounding_box_2

        # Respect ignored
        model = Gaussian2D()
        bounding_box_1._ignored = [mk.MagicMock()]
        bounding_box_2._ignored = [mk.MagicMock()]
        assert bounding_box_1._ignored != bounding_box_2._ignored
        assert not (bounding_box_1 == bounding_box_2)

        # Empty == (-np.inf, np.inf)
        bounding_box = ModelBoundingBox((), model)
        assert bounding_box == (-np.inf, np.inf)

    def test__setitem__(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, {}, ignored=[0, 1])
        assert bounding_box._ignored == ['x', 'y']

        # USING Intervals directly
        # Set interval using key
        assert 'x' not in bounding_box.intervals
        assert 'x' in bounding_box.ignored
        bounding_box['x'] = _Interval(-1, 1)
        assert 'x' in bounding_box.intervals
        assert 'x' not in bounding_box.ignored
        assert isinstance(bounding_box['x'], _Interval)
        assert bounding_box['x'] == (-1, 1)

        assert 'y' not in bounding_box.intervals
        assert 'y' in bounding_box.ignored
        bounding_box['y'] = _Interval(-4, 4)
        assert 'y' in bounding_box.intervals
        assert 'y' not in bounding_box.ignored
        assert isinstance(bounding_box['y'], _Interval)
        assert bounding_box['y'] == (-4, 4)

        assert bounding_box.ignored == []
        del bounding_box['x']
        del bounding_box['y']
        assert bounding_box.ignored == ['x', 'y']

        # Set interval using index
        assert 'x' not in bounding_box.intervals
        assert 'x' in bounding_box.ignored
        bounding_box[0] = _Interval(-1, 1)
        assert 'x' in bounding_box.intervals
        assert 'x' not in bounding_box.ignored
        assert isinstance(bounding_box[0], _Interval)
        assert bounding_box[0] == (-1, 1)

        assert 'y' not in bounding_box.intervals
        assert 'y' in bounding_box.ignored
        bounding_box[1] = _Interval(-4, 4)
        assert 'y' in bounding_box.intervals
        assert 'y' not in bounding_box.ignored
        assert isinstance(bounding_box[1], _Interval)
        assert bounding_box[1] == (-4, 4)

        assert bounding_box.ignored == []
        del bounding_box[0]
        del bounding_box[1]
        assert bounding_box.ignored == ['x', 'y']

        # USING tuples
        # Set interval using key
        assert 'x' not in bounding_box.intervals
        assert 'x' in bounding_box.ignored
        bounding_box['x'] = (-1, 1)
        assert 'x' in bounding_box.intervals
        assert 'x' not in bounding_box.ignored
        assert isinstance(bounding_box['x'], _Interval)
        assert bounding_box['x'] == (-1, 1)

        assert 'y' not in bounding_box.intervals
        assert 'y' in bounding_box.ignored
        bounding_box['y'] = (-4, 4)
        assert 'y' in bounding_box.intervals
        assert 'y' not in bounding_box.ignored
        assert isinstance(bounding_box['y'], _Interval)
        assert bounding_box['y'] == (-4, 4)

        assert bounding_box.ignored == []
        del bounding_box['x']
        del bounding_box['y']
        assert bounding_box.ignored == ['x', 'y']

        # Set interval using index
        assert 'x' not in bounding_box.intervals
        assert 'x' in bounding_box.ignored
        bounding_box[0] = (-1, 1)
        assert 'x' in bounding_box.intervals
        assert 'x' not in bounding_box.ignored
        assert isinstance(bounding_box[0], _Interval)
        assert bounding_box[0] == (-1, 1)

        assert 'y' not in bounding_box.intervals
        assert 'y' in bounding_box.ignored
        bounding_box[1] = (-4, 4)
        assert 'y' in bounding_box.intervals
        assert 'y' not in bounding_box.ignored
        assert isinstance(bounding_box[1], _Interval)
        assert bounding_box[1] == (-4, 4)

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = ModelBoundingBox({}, model)
        # USING Intervals directly
        # Set interval using key
        assert 'x' not in bounding_box
        bounding_box['x'] = _Interval(np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], _Interval)
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        bounding_box[0] = _Interval(np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box[0], _Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()
        # USING tuples
        # Set interval using key
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        bounding_box['x'] = (np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], _Interval)
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        bounding_box[0] = (np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box[0], _Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()

        # No model
        bounding_box = ModelBoundingBox({})
        bounding_box[0] = (0, 1)
        assert bounding_box._intervals[0] == (0, 1)
        assert isinstance(bounding_box._intervals[0], _Interval)

    def test___delitem__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Using index
        assert 'x' in bounding_box.intervals
        assert 'x' not in bounding_box.ignored
        assert 'x' in bounding_box
        assert 'x' in bounding_box
        del bounding_box[0]
        assert 'x' not in bounding_box.intervals
        assert 'x' in bounding_box.ignored
        assert 'x' in bounding_box
        assert 'x' in bounding_box

        # Delete an ignored item
        with pytest.raises(RuntimeError) as err:
            del bounding_box[0]
        assert str(err.value) ==\
            "Cannot delete ignored input: 0!"

        # Using key
        assert 'y' in bounding_box.intervals
        assert 'y' not in bounding_box.ignored
        assert 'x' in bounding_box
        assert 'y' in bounding_box
        del bounding_box['y']
        assert 'y' not in bounding_box.intervals
        assert 'y' in bounding_box.ignored
        assert 'x' in bounding_box
        assert 'y' in bounding_box

        # Delete an ignored item
        with pytest.raises(RuntimeError) as err:
            del bounding_box['y']
        assert str(err.value) ==\
            "Cannot delete ignored input: y!"

        # No model
        bounding_box = ModelBoundingBox({0: (0, 1)})
        assert bounding_box._intervals[0] == (0, 1)
        del bounding_box[0]
        assert 0 not in bounding_box._intervals

    def test_validate(self):
        model = Gaussian2D()
        kwargs = {'test': mk.MagicMock()}

        # Pass sequence Default order
        bounding_box = ModelBoundingBox.validate(model, ((-4, 4), (-1, 1)), **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box = ModelBoundingBox.validate(model, ((-4, 4), (-1, 1)), order='F', **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-4, 4)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        bounding_box = ModelBoundingBox.validate(model, intervals, order='F', **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2
        assert bounding_box.order == 'F'

        # Pass ModelBoundingBox
        bbox = bounding_box
        bounding_box = ModelBoundingBox.validate(model, bbox, **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2
        assert bounding_box.order == 'F'

        # Pass single ignored
        intervals = {0: _Interval(-1, 1)}
        bounding_box = ModelBoundingBox.validate(model, intervals, ignored=['y'], **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == _ignored_interval
        assert len(bounding_box.intervals) == 1

        # Pass single
        bounding_box = ModelBoundingBox.validate(Gaussian1D(), (-1, 1), **kwargs)
        assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert len(bounding_box.intervals) == 1

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        sequence = (np.array([-1, -2]), np.array([1, 2]))
        bounding_box = ModelBoundingBox.validate(model, sequence, **kwargs)
        assert 'x' in bounding_box
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()

    def test_fix_inputs(self):
        bounding_box = ModelBoundingBox.validate(Gaussian2D(), ((-4, 4), (-1, 1)))

        new_bounding_box = bounding_box.fix_inputs(Gaussian1D(), {1: mk.MagicMock()})
        assert not (bounding_box == new_bounding_box)

        assert (new_bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert 'x' in new_bounding_box
        assert new_bounding_box['x'] == (-1, 1)
        assert 'y' not in new_bounding_box
        assert len(new_bounding_box.intervals) == 1
        assert new_bounding_box.ignored == []

    def test_dimension(self):
        intervals = {0: _Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert bounding_box.dimension == 1 == len(bounding_box._intervals)

        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert bounding_box.dimension == 2 == len(bounding_box._intervals)

        bounding_box._intervals = {}
        assert bounding_box.dimension == 0 == len(bounding_box._intervals)

    def test_domain(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # test defaults
        assert (np.array(bounding_box.domain(0.25)) ==
                np.array([np.linspace(0, 2, 9), np.linspace(-1, 1, 9)])).all()

        # test C order
        assert (np.array(bounding_box.domain(0.25, 'C')) ==
                np.array([np.linspace(0, 2, 9), np.linspace(-1, 1, 9)])).all()

        # test Fortran order
        assert (np.array(bounding_box.domain(0.25, 'F')) ==
                np.array([np.linspace(-1, 1, 9), np.linspace(0, 2, 9)])).all()

        # test error order
        order = mk.MagicMock()
        with pytest.raises(ValueError):
            bounding_box.domain(0.25, order)

    def test__get_interval(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Get non-ignored
        for index in range(2):
            interval = bounding_box._get_interval(index, ['a'])
            assert isinstance(interval, _Interval)
            assert interval == intervals[index]

        # Get always ignored
        for index in range(2):
            interval = bounding_box._get_interval(index, ['x', 'y'])
            assert isinstance(interval, _Interval)
            assert interval == _ignored_interval

        # ignore only x
        interval = bounding_box._get_interval(0, ['x'])
        assert isinstance(interval, _Interval)
        assert interval == _ignored_interval
        interval = bounding_box._get_interval(1, ['x'])
        assert isinstance(interval, _Interval)
        assert interval == intervals[1]

        # ignore only y
        interval = bounding_box._get_interval(0, ['y'])
        assert isinstance(interval, _Interval)
        assert interval == intervals[0]
        interval = bounding_box._get_interval(1, ['y'])
        assert isinstance(interval, _Interval)
        assert interval == _ignored_interval

    def test__outside(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, [])
        assert (outside_index == [False for _ in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, [])
        assert (outside_index ==
                [True, True, True, True,
                 False, False, False, False, False,
                 True, True, True, True]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, [])
        assert (outside_index == [True for _ in range(13)]).all()
        assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, [])
        assert (outside_index == [False]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, [])
        assert (outside_index == [True]).all()
        assert all_out and isinstance(all_out, bool)

    def test__outside_with_ignored(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside, ignore x or y
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            outside_index, all_out = bounding_box._outside(input_shape, inputs, [ignore])
            assert (outside_index == [False for _ in range(13)]).all()
            assert not all_out and isinstance(all_out, bool)

        # normal array input, some inside and some outside, ignore x
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs, ['x'])
        assert (outside_index ==
                [False, False, False, False,
                 False, False, False, False, False,
                 True, True, True, True]).all()
        assert not all_out and isinstance(all_out, bool)
        # normal array input, some inside and some outside, ignore y
        outside_index, all_out = bounding_box._outside(input_shape, inputs, ['y'])
        assert (outside_index ==
                [True, True, True, True,
                 False, False, False, False, False,
                 False, False, False, False]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside, ignore x or y
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            outside_index, all_out = bounding_box._outside(input_shape, inputs, [ignore])
            assert (outside_index == [True for _ in range(13)]).all()
            assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box, ignore x or y
        inputs = (0.5, 0.5)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            outside_index, all_out = bounding_box._outside(input_shape, inputs, [ignore])
            assert (outside_index == [False]).all()
            assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box, ignore x or y
        inputs = (2, -1)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            outside_index, all_out = bounding_box._outside(input_shape, inputs, [ignore])
            assert (outside_index == [True]).all()
            assert all_out and isinstance(all_out, bool)

    def test__valid_index(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [])
        assert len(valid_index) == 1
        assert (valid_index[0] == [idx for idx in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [])
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [])
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [])
        assert len(valid_index) == 1
        assert (valid_index[0] == [0]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [])
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

    def test__valid_index_with_ignored(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside, ignore x or y
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [ignore])
            assert len(valid_index) == 1
            assert (valid_index[0] == [idx for idx in range(13)]).all()
            assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside, ignore x
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, ['x'])
        assert len(valid_index) == 1
        assert (valid_index[0] == [0, 1, 2, 3, 4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)
        # Normal array input, some inside and some outside, ignore y
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs, ['y'])
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8, 9, 10, 11, 12]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside, ignore x or y
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [ignore])
            assert len(valid_index) == 1
            assert (valid_index[0] == []).all()
            assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box, ignore x or y
        inputs = (0.5, 0.5)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [ignore])
            assert len(valid_index) == 1
            assert (valid_index[0] == [0]).all()
            assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box, ignore x or y
        inputs = (2, -1)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            valid_index, all_out = bounding_box._valid_index(input_shape, inputs, [ignore])
            assert len(valid_index) == 1
            assert (valid_index[0] == []).all()
            assert all_out and isinstance(all_out, bool)

    def test_prepare_inputs(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs)
        assert (np.array(new_inputs) == np.array(inputs)).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [idx for idx in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs)
        assert (np.array(new_inputs) ==
                np.array(
                    [
                        [x[4], x[5], x[6], x[7], x[8]],
                        [y[4], y[5], y[6], y[7], y[8]],
                    ]
                )).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs)
        assert new_inputs == ()
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs)
        assert (np.array(new_inputs) == np.array([[0.5], [0.5]])).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [0]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs)
        assert new_inputs == ()
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

    def test_prepare_inputs_with_ignored(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside, ignore x or y
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, [ignore])
            assert (np.array(new_inputs) == np.array(inputs)).all()
            assert len(valid_index) == 1
            assert (valid_index[0] == [idx for idx in range(13)]).all()
            assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside, ignore x
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, ['x'])
        assert (np.array(new_inputs) ==
                np.array(
                    [
                        [x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]],
                        [y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]],
                    ]
                )).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [0, 1, 2, 3, 4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)
        # Normal array input, some inside and some outside, ignore y
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, ['y'])
        assert (np.array(new_inputs) ==
                np.array(
                    [
                        [x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12]],
                        [y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11], y[12]],
                    ]
                )).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8, 9, 10, 11, 12]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside, ignore x or y
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        for ignore in ['x', 'y']:
            new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, [ignore])
            assert new_inputs == ()
            assert len(valid_index) == 1
            assert (valid_index[0] == []).all()
            assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box, ignore x or y
        inputs = (0.5, 0.5)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, [ignore])
            assert (np.array(new_inputs) == np.array([[0.5], [0.5]])).all()
            assert len(valid_index) == 1
            assert (valid_index[0] == [0]).all()
            assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box, ignore x or y
        inputs = (2, -1)
        input_shape = (1,)
        for ignore in ['x', 'y']:
            new_inputs, valid_index, all_out = bounding_box.prepare_inputs(input_shape, inputs, [ignore])
            assert new_inputs == ()
            assert len(valid_index) == 1
            assert (valid_index[0] == []).all()
            assert all_out and isinstance(all_out, bool)


class Test_SelectorArgument:
    def test_create(self):
        name = mk.MagicMock()
        ignore = mk.MagicMock()
        argument = _SelectorArgument(name, ignore)

        assert isinstance(argument, _BaseSelectorArgument)
        assert argument.name == name
        assert argument.ignore == ignore
        assert argument == (name, ignore)

    def test_validate(self):
        model = Gaussian2D()

        # default integer
        assert _SelectorArgument.validate(model, 0) == ('x', True)
        assert _SelectorArgument.validate(model, 1) == ('y', True)

        # default string
        assert _SelectorArgument.validate(model, 'x') == ('x', True)
        assert _SelectorArgument.validate(model, 'y') == ('y', True)

        ignore = mk.MagicMock()
        # non-default integer
        assert _SelectorArgument.validate(model, 0, ignore) == ('x', ignore)
        assert _SelectorArgument.validate(model, 1, ignore) == ('y', ignore)

        # non-default string
        assert _SelectorArgument.validate(model, 'x', ignore) == ('x', ignore)
        assert _SelectorArgument.validate(model, 'y', ignore) == ('y', ignore)

        # Fail
        with pytest.raises(ValueError):
            _SelectorArgument.validate(model, 'z')
        with pytest.raises(ValueError):
            _SelectorArgument.validate(model, mk.MagicMock())
        with pytest.raises(IndexError):
            _SelectorArgument.validate(model, 2)

    def test_index(self):
        model = Gaussian2D()
        for index, name in enumerate(model.inputs):
            assert _SelectorArgument(name, mk.MagicMock()).index(model) == index

    def test_get_selector(self):
        model = mk.MagicMock()
        model.inputs = ['x', 'y', 'z']

        # single inputs
        inputs = [idx + 17 for idx in range(len(model.inputs))]
        for index, name in enumerate(model.inputs):
            assert _SelectorArgument(name, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]
            assert _SelectorArgument(index, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]

        # numpy array of single inputs
        inputs = [np.array([idx + 11]) for idx in range(len(model.inputs))]
        for index, name in enumerate(model.inputs):
            assert _SelectorArgument(name, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]
            assert _SelectorArgument(index, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]
        inputs = [np.asanyarray(idx + 13) for idx in range(len(model.inputs))]
        for index, name in enumerate(model.inputs):
            assert _SelectorArgument(name, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]
            assert _SelectorArgument(index, mk.MagicMock()).get_selector(model, *inputs) == inputs[index]

        # multi entry numpy array
        inputs = [np.array([idx + 27, idx - 31]) for idx in range(len(model.inputs))]
        for index, name in enumerate(model.inputs):
            assert _SelectorArgument(name, mk.MagicMock()).get_selector(model, *inputs) == tuple(inputs[index])
            assert _SelectorArgument(index, mk.MagicMock()).get_selector(model, *inputs) == tuple(inputs[index])

    def test___repr__(self):
        model = Gaussian2D()

        assert _SelectorArgument('x', False).__repr__() ==\
            "Argument(name='x', ignore=False)"
        assert _SelectorArgument('x', True).__repr__() ==\
            "Argument(name='x', ignore=True)"

        assert _SelectorArgument('y', False).__repr__() ==\
            "Argument(name='y', ignore=False)"
        assert _SelectorArgument('y', True).__repr__() ==\
            "Argument(name='y', ignore=True)"

    def test_get_fixed_value(self):
        model = Gaussian2D()
        values = {0: 5, 'y': 7}

        # Get index value
        assert _SelectorArgument('x', mk.MagicMock()).get_fixed_value(model, values) == 5

        # Get name value
        assert _SelectorArgument('y', mk.MagicMock()).get_fixed_value(model, values) == 7

        # Fail
        values = {0: 5}
        with pytest.raises(RuntimeError) as err:
            _SelectorArgument('y', True).get_fixed_value(model, values)
        assert str(err.value) == \
            "Argument(name='y', ignore=True) was not found in {0: 5}"

    def test_is_argument(self):
        model = Gaussian2D()
        argument = _SelectorArgument.validate(model, 0)

        # Is true
        assert argument.is_argument(model, 0) == True
        assert argument.is_argument(model, 'x') == True

        # Is false
        assert argument.is_argument(model, 1) == False
        assert argument.is_argument(model, 'y') == False

        # Fail
        with pytest.raises(ValueError):
            argument.is_argument(model, 'z')
        with pytest.raises(ValueError):
            argument.is_argument(model, mk.MagicMock())
        with pytest.raises(IndexError):
            argument.is_argument(model, 2)

    def test_index_tuple(self):
        model = Gaussian2D()
        for index, name in enumerate(model.inputs):
            ignore = mk.MagicMock()
            assert _SelectorArgument(name, ignore).index_tuple(model) == \
                (index, ignore)


class Test_SelectorArguments:
    def test_create(self):
        arguments = _SelectorArguments((_SelectorArgument(0, True), _SelectorArgument(1, False)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, False))

    def test___repr__(self):
        arguments = _SelectorArguments((_SelectorArgument('x', True), _SelectorArgument('y', False)))

        assert arguments.__repr__() ==\
            "SelectorArguments(\n" +\
            "    Argument(name='x', ignore=True)\n" +\
            "    Argument(name='y', ignore=False)\n" +\
            ")"

    def test_ignored(self):
        assert _SelectorArguments((_SelectorArgument('x', True),
                                  _SelectorArgument('y', True))).ignored == ['x', 'y']
        assert _SelectorArguments((_SelectorArgument('x', True),
                                   _SelectorArgument('y', False))).ignored == ['x']
        assert _SelectorArguments((_SelectorArgument('x', False),
                                   _SelectorArgument('y', True))).ignored == ['y']
        assert _SelectorArguments((_SelectorArgument('x', False),
                                   _SelectorArgument('y', False))).ignored == []

    def test_validate(self):
        # Integer key and passed ignore
        arguments = _SelectorArguments.validate(Gaussian2D(), ((0, True), (1, False)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == (('x', True), ('y', False))

        # Default ignore
        arguments = _SelectorArguments.validate(Gaussian2D(), ((0,), (1,)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == (('x', True), ('y', True))

        # String key and passed ignore
        arguments = _SelectorArguments.validate(Gaussian2D(), (('x', True), ('y', False)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == (('x', True), ('y', False))

        # Invalid, bad argument
        with pytest.raises(ValueError):
            _SelectorArguments.validate(Gaussian2D(), ((0, True), ('z', False)))
        with pytest.raises(ValueError):
            _SelectorArguments.validate(Gaussian2D(), ((mk.MagicMock(), True), (1, False)))
        with pytest.raises(IndexError):
            _SelectorArguments.validate(Gaussian2D(), ((0, True), (2, False)))

        # Invalid, repeated argument
        with pytest.raises(ValueError) as err:
            _SelectorArguments.validate(Gaussian2D(), ((0, True), (0, False)))
        assert str(err.value) == \
            "Input: 'x' has been repeated."

        # Invalid, no arguments
        with pytest.raises(ValueError) as err:
            _SelectorArguments.validate(Gaussian2D(), ())
        assert str(err.value) == \
            "There must be at least one selector argument."

    def test_get_selector(self):
        model = mk.MagicMock()
        model.inputs = ['a', 'b', 'c', 'd']

        inputs = [idx + 19 for idx in range(len(model.inputs))]

        assert _SelectorArguments.validate(model,
                                       (('a', True), ('b', False))).get_selector(model, *inputs) ==\
            tuple(inputs[:2])
        assert _SelectorArguments.validate(model,
                                       (('b', True), ('a', False))).get_selector(model, *inputs) ==\
            tuple(inputs[:2][::-1])
        assert _SelectorArguments.validate(model,
                                       (('b', False),)).get_selector(model, *inputs) ==\
            (inputs[1],)
        assert _SelectorArguments.validate(model,
                                       (('a', True),)).get_selector(model, *inputs) ==\
            (inputs[0],)

    def test_is_selector(self):
        # Is Selector
        assert _SelectorArguments.validate(Gaussian2D(),
                                       ((0, True), (1, False))).is_selector((0.5, 2.5))
        assert _SelectorArguments.validate(Gaussian2D(),
                                       ((0, True),)).is_selector((0.5,))

        # Is not selector
        assert not _SelectorArguments.validate(Gaussian2D(),
                                           ((0, True), (1, False))).is_selector((0.5, 2.5, 3.5))
        assert not _SelectorArguments.validate(Gaussian2D(),
                                           ((0, True), (1, False))).is_selector((0.5,))
        assert not _SelectorArguments.validate(Gaussian2D(),
                                           ((0, True), (1, False))).is_selector(0.5)
        assert not _SelectorArguments.validate(Gaussian2D(),
                                           ((0, True),)).is_selector((0.5, 2.5))
        assert not _SelectorArguments.validate(Gaussian2D(),
                                           ((0, True),)).is_selector(2.5)

    def test_get_fixed_values(self):
        model = Gaussian2D()

        assert _SelectorArguments.validate(model,
                                       ((0, True), (1, False))).get_fixed_values(model, {0: 11, 1: 7}) \
            == (11, 7)
        assert _SelectorArguments.validate(model,
                                       ((0, True), (1, False))).get_fixed_values(model, {0: 5, 'y': 47}) \
            == (5, 47)
        assert _SelectorArguments.validate(model,
                                       ((0, True), (1, False))).get_fixed_values(model, {'x': 2, 'y': 9}) \
            == (2, 9)
        assert _SelectorArguments.validate(model,
                                       ((0, True), (1, False))).get_fixed_values(model, {'x': 12, 1: 19}) \
            == (12, 19)

    def test_is_argument(self):
        model = Gaussian2D()

        # Is true
        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.is_argument(model, 0) == True
        assert arguments.is_argument(model, 'x') == True
        assert arguments.is_argument(model, 1) == True
        assert arguments.is_argument(model, 'y') == True

        # Is true and false
        arguments = _SelectorArguments.validate(model, ((0, True),))
        assert arguments.is_argument(model, 0) == True
        assert arguments.is_argument(model, 'x') == True
        assert arguments.is_argument(model, 1) == False
        assert arguments.is_argument(model, 'y') == False

        arguments = _SelectorArguments.validate(model, ((1, False),))
        assert arguments.is_argument(model, 0) == False
        assert arguments.is_argument(model, 'x') == False
        assert arguments.is_argument(model, 1) == True
        assert arguments.is_argument(model, 'y') == True

    def test_selector_index(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.selector_index(model, 0) == 0
        assert arguments.selector_index(model, 'x') == 0
        assert arguments.selector_index(model, 1) == 1
        assert arguments.selector_index(model, 'y') == 1

        arguments = _SelectorArguments.validate(model, ((1, True), (0, False)))
        assert arguments.selector_index(model, 0) == 1
        assert arguments.selector_index(model, 'x') == 1
        assert arguments.selector_index(model, 1) == 0
        assert arguments.selector_index(model, 'y') == 0

        # Error
        arguments = _SelectorArguments.validate(model, ((0, True),))
        with pytest.raises(ValueError) as err:
            arguments.selector_index(model, 'y')
        assert str(err.value) ==\
            "y does not correspond to any selector argument."

    def test_reduce(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))

        new_arguments = arguments.reduce(model, 0)
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == (('y', False),)

        new_arguments = arguments.reduce(model, 'x')
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == (('y', False),)

        new_arguments = arguments.reduce(model, 1)
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == (('x', True),)

        new_arguments = arguments.reduce(model, 'y')
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == (('x', True),)

    def test_index_tuple(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.index_tuple(model) == ((0, True), (1, False))

        arguments = _SelectorArguments.validate(model, (('x', True), ('y', False)))
        assert arguments.index_tuple(model) == ((0, True), (1, False))


class TestCompoundBoundingBox:
    def test_create(self):
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        model = mk.MagicMock()
        model.inputs = ['x', 'y']
        selector_args = (('x', True),)
        ignored = ['x']

        # No bounding_boxes
        bounding_box = CompoundBoundingBox({})
        assert bounding_box._model is None
        assert bounding_box._selector_args is None
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == {}

        # bounding_boxes only:
        bounding_box = CompoundBoundingBox(bounding_boxes)
        assert bounding_box._model is None
        assert bounding_box._selector_args is None
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == bounding_boxes

        # model only:
        bounding_box = CompoundBoundingBox({}, model)
        assert bounding_box._model == model
        assert bounding_box._selector_args is None
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == {}

        # bounding_boxes and model only:
        bounding_box = CompoundBoundingBox(bounding_boxes, model)
        assert bounding_box._model == model
        assert bounding_box._selector_args is None
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == bounding_boxes

        # bounding_boxes and selector_args only:
        bounding_box = CompoundBoundingBox(bounding_boxes, selector_args=selector_args)
        assert bounding_box._model is None
        assert bounding_box._selector_args == selector_args
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == bounding_boxes

        # bounding_boxes and ignored only:
        bounding_box = CompoundBoundingBox(bounding_boxes, ignored=ignored)
        assert bounding_box._model is None
        assert bounding_box._selector_args is None
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == ignored
        assert bounding_box._bounding_boxes == bounding_boxes

        # bounding_box, model, and selector_args only:
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)
        assert bounding_box._model == model
        assert bounding_box._selector_args == selector_args
        assert not isinstance(selector_args, _SelectorArguments)
        assert isinstance(bounding_box._selector_args, _SelectorArguments)
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == []
        assert bounding_box._bounding_boxes == bounding_boxes
        for selector, bbox in bounding_box._bounding_boxes.items():
            assert not isinstance(bounding_boxes[selector], ModelBoundingBox)
            assert isinstance(bbox, ModelBoundingBox)
            assert bbox._intervals == {'y': bounding_boxes[selector]}
            assert bbox._ignored == []

        # bounding_box, model, selector_args, and ignored only:
        model.inputs = ['x', 'y', 'z']
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args, ignored=['z'])
        assert bounding_box._model == model
        assert bounding_box._selector_args == selector_args
        assert not isinstance(selector_args, _SelectorArguments)
        assert isinstance(bounding_box._selector_args, _SelectorArguments)
        assert bounding_box._create_selector is None
        assert bounding_box._order == 'C'
        assert bounding_box._ignored == ['z']
        assert bounding_box._bounding_boxes == bounding_boxes
        for selector, bbox in bounding_box._bounding_boxes.items():
            assert not isinstance(bounding_boxes[selector], ModelBoundingBox)
            assert isinstance(bbox, ModelBoundingBox)
            assert bbox._intervals == {'y': bounding_boxes[selector]}
            assert bbox._ignored == []

        # More practical example
        model = Gaussian2D()
        selector_args = (('x', True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        create_selector = mk.MagicMock()

        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args, create_selector, order='F')
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        for selector, bbox in bounding_boxes.items():
            assert selector in bounding_box._bounding_boxes
            assert bounding_box._bounding_boxes[selector] == bbox
        for selector, bbox in bounding_box._bounding_boxes.items():
            assert selector in bounding_boxes
            assert bounding_boxes[selector] == bbox
            assert isinstance(bbox, ModelBoundingBox)
            assert bbox._intervals == {'y': bounding_boxes[selector]}
            assert bbox._ignored == []
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == 'F'

    def test__verify_selector_args(self):
        bounding_box = CompoundBoundingBox({}, Gaussian2D())
        bounding_box._verify_selector_args()

        bounding_box._selector_args = (('x', True),)
        assert not isinstance(bounding_box._selector_args, _SelectorArguments)
        bounding_box._verify_selector_args()
        assert isinstance(bounding_box._selector_args, _SelectorArguments)
        assert bounding_box._selector_args == (('x', True),)

    def test__pop_bounding_boxes(self):
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes)
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._pop_bounding_boxes() == bounding_boxes
        assert bounding_boxes != {}
        assert bounding_box._bounding_boxes == {}

    def test__verify_bounding_boxes(self):
        bounding_box = CompoundBoundingBox({}, Gaussian2D(), (('x', True),))
        assert bounding_box._bounding_boxes == {}

        # With selector args
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box._bounding_boxes = bounding_boxes
        bounding_box._verify_bounding_boxes()
        assert bounding_box._bounding_boxes == bounding_boxes
        for selector, bbox in bounding_box._bounding_boxes.items():
            assert not isinstance(bounding_boxes[selector], ModelBoundingBox)
            assert isinstance(bbox, ModelBoundingBox)

        # Without selector args
        bounding_box._selector_args = None
        bounding_box._bounding_boxes = bounding_boxes
        bounding_box._verify_bounding_boxes()
        assert bounding_box._bounding_boxes == bounding_boxes
        for selector, bbox in bounding_box._bounding_boxes.items():
            assert not isinstance(bounding_boxes[selector], ModelBoundingBox)
            assert not isinstance(bbox, ModelBoundingBox)

    def test__verify(self):
        external_ignored = mk.MagicMock()

        bounding_box_0 = CompoundBoundingBox({})
        bounding_box_1 = CompoundBoundingBox({}, mk.MagicMock())
        with mk.patch.object(CompoundBoundingBox, '_verify_selector_args') as mkArgs:
            with mk.patch.object(CompoundBoundingBox, '_verify_bounding_boxes') as mkBbox:
                main = mk.MagicMock()
                main.attach_mock(mkArgs, 'select')
                main.attach_mock(mkBbox, 'bbox')

                # No model, no _external_ignored
                bounding_box_0._verify()
                assert main.mock_calls == []

                # No model, with _external_ignored
                bounding_box_0._verify(external_ignored)
                assert main.mock_calls == []

                # With model, no _external_ignored
                bounding_box_1._verify()
                assert main.mock_calls == [
                    mk.call.select(),
                    mk.call.bbox(None)
                ]

                main.reset_mock()

                # With model, with _external_ignored
                bounding_box_1._verify(external_ignored)
                assert main.mock_calls == [
                    mk.call.select(),
                    mk.call.bbox(external_ignored)
                ]

    def test_copy(self):
        bounding_box = CompoundBoundingBox.validate(Gaussian2D(), {(1,): (-1.5, 1.3), (2,): (-2.7, 2.4)},
                                                    ((0, True),), mk.MagicMock())
        copy = bounding_box.copy()

        assert bounding_box == copy
        assert id(bounding_box) != id(copy)

        # model is not copied to prevent infinite recursion
        assert bounding_box._model == copy._model
        assert id(bounding_box._model) == id(copy._model)

        # Same string values have will have same id
        assert bounding_box._order == copy._order
        assert id(bounding_box._order) == id(copy._order)

        assert bounding_box._create_selector == copy._create_selector
        assert id(bounding_box._create_selector) != id(copy._create_selector)

        # Check selector_args
        for index, argument in enumerate(bounding_box.selector_args):
            assert argument == copy.selector_args[index]
            assert id(argument) != id(copy.selector_args[index])

            # Same integer values have will have same id
            assert argument.name == copy.selector_args[index].name
            assert id(argument.name) == id(copy.selector_args[index].name)

            # Same boolean values have will have same id
            assert argument.ignore == copy.selector_args[index].ignore
            assert id(argument.ignore) == id(copy.selector_args[index].ignore)
        assert len(bounding_box.selector_args) == len(copy.selector_args)

        # Check bounding_boxes
        for selector, bbox in bounding_box.bounding_boxes.items():
            assert bbox == copy.bounding_boxes[selector]
            assert id(bbox) != id(copy.bounding_boxes[selector])

            assert bbox.ignored == copy.bounding_boxes[selector].ignored
            assert id(bbox.ignored) != id(copy.bounding_boxes[selector].ignored)

            # model is not copied to prevent infinite recursion
            assert bbox._model == copy.bounding_boxes[selector]._model
            assert id(bbox._model) == id(copy.bounding_boxes[selector]._model)

            # Same string values have will have same id
            assert bbox._order == copy.bounding_boxes[selector]._order
            assert id(bbox._order) == id(copy.bounding_boxes[selector]._order)

            # Check interval objects
            for index, interval in bbox.intervals.items():
                assert interval == copy.bounding_boxes[selector].intervals[index]
                assert id(interval) != id(copy.bounding_boxes[selector].intervals[index])

                # Same float values have will have same id
                assert interval.lower == copy.bounding_boxes[selector].intervals[index].lower
                assert id(interval.lower) == id(copy.bounding_boxes[selector].intervals[index].lower)

                # Same float values have will have same id
                assert interval.upper == copy.bounding_boxes[selector].intervals[index].upper
                assert id(interval.upper) == id(copy.bounding_boxes[selector].intervals[index].upper)
            assert len(bbox.intervals) == len(copy.bounding_boxes[selector].intervals)
            assert bbox.intervals.keys() == copy.bounding_boxes[selector].intervals.keys()
        assert len(bounding_box.bounding_boxes) == len(copy.bounding_boxes)
        assert bounding_box.bounding_boxes.keys() == copy.bounding_boxes.keys()

    def test___repr__(self):
        model = Gaussian2D()
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}

        # No global ignore
        selector_args = ((0, True),)
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        assert bounding_box.__repr__() ==\
            "CompoundBoundingBox(\n" + \
            "    bounding_boxes={\n" + \
            "        (1,) = ModelBoundingBox(\n" + \
            "                intervals={\n" + \
            "                    y: Interval(lower=-1, upper=1)\n" + \
            "                }\n" + \
            "                model=Gaussian2D(inputs=('x', 'y'))\n" + \
            "                order='C'\n" + \
            "            )\n" + \
            "        (2,) = ModelBoundingBox(\n" + \
            "                intervals={\n" + \
            "                    y: Interval(lower=-2, upper=2)\n" + \
            "                }\n" + \
            "                model=Gaussian2D(inputs=('x', 'y'))\n" + \
            "                order='C'\n" + \
            "            )\n" + \
            "    }\n" + \
            "    selector_args = SelectorArguments(\n" + \
            "            Argument(name='x', ignore=True)\n" + \
            "        )\n" + \
            ")"

        # Global ignore
        selector_args = ((0, False),)
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args, ignored=['x'])

        assert bounding_box.__repr__() ==\
            "CompoundBoundingBox(\n" + \
            "    bounding_boxes={\n" + \
            "        (1,) = ModelBoundingBox(\n" + \
            "                intervals={\n" + \
            "                    y: Interval(lower=-1, upper=1)\n" + \
            "                }\n" + \
            "                model=Gaussian2D(inputs=('x', 'y'))\n" + \
            "                order='C'\n" + \
            "            )\n" + \
            "        (2,) = ModelBoundingBox(\n" + \
            "                intervals={\n" + \
            "                    y: Interval(lower=-2, upper=2)\n" + \
            "                }\n" + \
            "                model=Gaussian2D(inputs=('x', 'y'))\n" + \
            "                order='C'\n" + \
            "            )\n" + \
            "    }\n" + \
            "    ignored=['x']\n" +\
            "    selector_args = SelectorArguments(\n" + \
            "            Argument(name='x', ignore=False)\n" + \
            "        )\n" + \
            ")"

    def test_bounding_boxes(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box.bounding_boxes == bounding_boxes

    def test_selector_args(self):
        model = Gaussian2D()
        selector_args = (('x', True),)
        bounding_box = CompoundBoundingBox({}, model)

        # Get error
        assert bounding_box._selector_args is None
        with pytest.raises(RuntimeError, match="selector_args must be specified for a fully functional*"):
            bounding_box.selector_args

        bounding_box._selector_args = selector_args

        # Get
        assert bounding_box._selector_args == selector_args
        assert bounding_box.selector_args == selector_args

        # Set without override
        bounding_box = CompoundBoundingBox({}, model)
        bounding_box.selector_args = selector_args
        assert bounding_box._selector_args == selector_args
        assert bounding_box.selector_args == selector_args

        # Set with override
        selector_args = (('y', False),)
        with pytest.warns(RuntimeWarning, match=r"Overriding selector_args.*"):
            bounding_box.selector_args = selector_args
        assert bounding_box._selector_args == selector_args
        assert bounding_box.selector_args == selector_args

    def test_create_selector(self):
        model = Gaussian2D()
        create_selector = mk.MagicMock()
        bounding_box = CompoundBoundingBox({}, model, ((1,),), create_selector)

        assert bounding_box._create_selector == create_selector
        assert bounding_box.create_selector == create_selector

    def test__get_selector_key(self):
        bounding_box = CompoundBoundingBox({}, Gaussian2D(), ((1, True),))
        assert len(bounding_box.bounding_boxes) == 0

        # Singlar
        assert bounding_box._get_selector_key(5) == (5,)
        assert bounding_box._get_selector_key((5,)) == (5,)
        assert bounding_box._get_selector_key([5]) == (5,)
        assert bounding_box._get_selector_key(np.asanyarray(5)) == (5,)
        assert bounding_box._get_selector_key(np.array([5])) == (5,)

        # multiple
        assert bounding_box._get_selector_key((5, 19)) == (5, 19)
        assert bounding_box._get_selector_key([5, 19]) == (5, 19)
        assert bounding_box._get_selector_key(np.array([5, 19])) == (5, 19)

    def test___setitem__(self):
        model = Gaussian2D()

        # Ignored argument
        bounding_box = CompoundBoundingBox({}, model, ((1, True),), order='F')
        assert len(bounding_box.bounding_boxes) == 0
        # Valid
        bounding_box[(15, )] = (-15, 15)
        assert len(bounding_box.bounding_boxes) == 1
        assert (15,) in bounding_box._bounding_boxes
        assert isinstance(bounding_box._bounding_boxes[(15,)], ModelBoundingBox)
        assert bounding_box._bounding_boxes[(15,)] == (-15, 15)
        assert bounding_box._bounding_boxes[(15,)].order == 'F'
        # Invalid key
        assert (7, 13) not in bounding_box._bounding_boxes
        with pytest.raises(ValueError) as err:
            bounding_box[(7, 13)] = (-7, 7)
        assert str(err.value) == \
            "(7, 13) is not a selector!"
        assert (7, 13) not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1
        # Invalid bounding box
        assert 13 not in bounding_box._bounding_boxes
        with pytest.raises(ValueError):
            bounding_box[(13,)] = ((-13, 13), (-3, 3))
        assert 13 not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1

        # No ignored argument
        bounding_box = CompoundBoundingBox({}, model, ((1, False),), order='F')
        assert len(bounding_box.bounding_boxes) == 0
        # Valid
        bounding_box[(15, )] = ((-15, 15), (-6, 6))
        assert len(bounding_box.bounding_boxes) == 1
        assert (15,) in bounding_box._bounding_boxes
        assert isinstance(bounding_box._bounding_boxes[(15,)], ModelBoundingBox)
        assert bounding_box._bounding_boxes[(15,)] == ((-15, 15), (-6, 6))
        assert bounding_box._bounding_boxes[(15,)].order == 'F'
        # Invalid key
        assert (14, 11) not in bounding_box._bounding_boxes
        with pytest.raises(ValueError) as err:
            bounding_box[(14, 11)] = ((-7, 7), (-12, 12))
        assert str(err.value) == \
            "(14, 11) is not a selector!"
        assert (14, 11) not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1
        # Invalid bounding box
        assert 13 not in bounding_box._bounding_boxes
        with pytest.raises(ValueError):
            bounding_box[(13,)] = (-13, 13)
        assert 13 not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1

        # _external_ignored
        model = mk.MagicMock()
        model.inputs = ['x', 'y', 'z']
        bounding_box = CompoundBoundingBox({}, model, ((1, True),), order='F')
        assert len(bounding_box.bounding_boxes) == 0
        bounding_box.__setitem__((12,), (-7, 4), _external_ignored=['z'])
        assert (12,) in bounding_box
        assert isinstance(bounding_box._bounding_boxes[(12,)], ModelBoundingBox)
        assert bounding_box._bounding_boxes[(12,)] == (-7, 4)
        assert bounding_box._bounding_boxes[(12,)].order == 'F'

    def test___eq__(self):
        bounding_box_1 = CompoundBoundingBox({(1,): (-1, 1), (2,): (-2, 2)}, Gaussian2D(), ((0, True),))
        bounding_box_2 = CompoundBoundingBox({(1,): (-1, 1), (2,): (-2, 2)}, Gaussian2D(), ((0, True),))

        # Equal
        assert bounding_box_1 == bounding_box_2

        # Not equal to non-compound bounding_box
        assert not bounding_box_1 == mk.MagicMock()
        assert not bounding_box_2 == mk.MagicMock()

        # Not equal bounding_boxes
        bounding_box_2[(15,)] = (-15, 15)
        assert not bounding_box_1 == bounding_box_2
        del bounding_box_2._bounding_boxes[(15,)]
        assert bounding_box_1 == bounding_box_2

        # Not equal selector_args
        bounding_box_2._selector_args = _SelectorArguments.validate(Gaussian2D(), ((0, False),))
        assert not bounding_box_1 == bounding_box_2
        bounding_box_2._selector_args = _SelectorArguments.validate(Gaussian2D(), ((0, True),))
        assert bounding_box_1 == bounding_box_2

        # Not equal create_selector
        bounding_box_2._create_selector = mk.MagicMock()
        assert not bounding_box_1 == bounding_box_2

    def test_validate(self):
        model = Gaussian2D()
        selector_args = (('x', True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        create_selector = mk.MagicMock()

        # Fail selector_args
        with pytest.raises(ValueError) as err:
            CompoundBoundingBox.validate(model, bounding_boxes)
        assert str(err.value) ==\
            "Selector arguments must be provided (can be passed as part of bounding_box argument)!"

        # Normal validate
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args,
                                                    create_selector, order='F')
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == 'F'

        # Re-validate
        new_bounding_box = CompoundBoundingBox.validate(model, bounding_box)
        assert bounding_box == new_bounding_box
        assert new_bounding_box._order == 'F'

        # Default order
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args,
                                                    create_selector)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == 'C'

        # ignored
        model = mk.MagicMock()
        model.inputs = ['x', 'y', 'z']
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args,
                                                    ignored=['z'])
        assert bounding_box._ignored == ['z']
        model.inputs = ['x', 'y', 'z', 'a']
        bounding_box = CompoundBoundingBox.validate(model, bounding_box, ignored=['a'])
        assert bounding_box._ignored == ['a', 'z']

    def test___contains__(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        assert (1,) in bounding_box
        assert (2,) in bounding_box

        assert (3,) not in bounding_box
        assert 1 not in bounding_box
        assert 2 not in bounding_box

    def test__create_bounding_box(self):
        model = Gaussian2D()
        create_selector = mk.MagicMock()
        bounding_box = CompoundBoundingBox({}, model, ((1, False),),
                                           create_selector)

        # Create is successful
        create_selector.return_value = ((-15, 15), (-23, 23))
        assert len(bounding_box._bounding_boxes) == 0
        bbox = bounding_box._create_bounding_box((7,))
        assert isinstance(bbox, ModelBoundingBox)
        assert bbox == ((-15, 15), (-23, 23))
        assert len(bounding_box._bounding_boxes) == 1
        assert (7,) in bounding_box
        assert isinstance(bounding_box[(7,)], ModelBoundingBox)
        assert bounding_box[(7,)] == bbox

        # Create is unsuccessful
        create_selector.return_value = (-42, 42)
        with pytest.raises(ValueError):
            bounding_box._create_bounding_box((27,))

    def test___getitem__(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        # already exists
        assert isinstance(bounding_box[1], ModelBoundingBox)
        assert bounding_box[1] == (-1, 1)
        assert isinstance(bounding_box[(2,)], ModelBoundingBox)
        assert bounding_box[2] == (-2, 2)
        assert isinstance(bounding_box[(1,)], ModelBoundingBox)
        assert bounding_box[(1,)] == (-1, 1)
        assert isinstance(bounding_box[(2,)], ModelBoundingBox)
        assert bounding_box[(2,)] == (-2, 2)

        # no selector
        with pytest.raises(RuntimeError) as err:
            bounding_box[(3,)]
        assert str(err.value) == \
            "No bounding box is defined for selector: (3,)."

        # Create a selector
        bounding_box._create_selector = mk.MagicMock()
        with mk.patch.object(CompoundBoundingBox, '_create_bounding_box',
                             autospec=True) as mkCreate:
            assert bounding_box[(3,)] == mkCreate.return_value
            assert mkCreate.call_args_list == \
                [mk.call(bounding_box, (3,))]

    def test__select_bounding_box(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        inputs = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(_SelectorArguments, 'get_selector',
                             autospec=True) as mkSelector:
            with mk.patch.object(CompoundBoundingBox, '__getitem__',
                                 autospec=True) as mkGet:
                assert bounding_box._select_bounding_box(inputs) == mkGet.return_value
                assert mkGet.call_args_list == \
                    [mk.call(bounding_box, mkSelector.return_value)]
                assert mkSelector.call_args_list == \
                    [mk.call(bounding_box.selector_args, model, *inputs)]

    def test_prepare_inputs(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        input_shape = mk.MagicMock()
        with mk.patch.object(ModelBoundingBox, 'prepare_inputs',
                             autospec=True) as mkPrepare:
            assert bounding_box.prepare_inputs(input_shape, [1, 2, 3]) == mkPrepare.return_value
            assert mkPrepare.call_args_list == \
                [mk.call(bounding_box[(1,)], input_shape, [1, 2, 3], ['x'])]
            mkPrepare.reset_mock()
            assert bounding_box.prepare_inputs(input_shape, [2, 2, 3], []) == mkPrepare.return_value
            assert mkPrepare.call_args_list == \
                [mk.call(bounding_box[(2,)], input_shape, [2, 2, 3], ['x'])]
            mkPrepare.reset_mock()

    def test__fix_inputs_matching_bounding_boxes(self):
        # Single selector index
        selector_args = ((0, False),)
        bounding_boxes = {(1,): ((-1, 1), (-2, 2)), (2,): ((-2, 2), (-3, 3)), (3,): ((-3, 3), (-4, 4))}
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            matching = bounding_box._fix_inputs_matching_bounding_boxes('x', value)
            assert isinstance(matching, dict)
            assert () in matching
            bbox = matching[()]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert 'y' in bbox
            assert bbox['y'] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == []

        # Multiple selector index
        selector_args = ((0, False), (1, False))
        bounding_boxes = {(1, 3): ((-1, 1), (-2, 2)), (2, 2): ((-2, 2), (-3, 3)), (3, 1): ((-3, 3), (-4, 4))}
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            matching = bounding_box._fix_inputs_matching_bounding_boxes('x', value)
            assert isinstance(matching, dict)
            assert (4 - value,) in matching
            bbox = matching[(4 - value,)]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert 'y' in bbox
            assert bbox['y'] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == []

            matching = bounding_box._fix_inputs_matching_bounding_boxes('y', value)
            assert isinstance(matching, dict)
            assert (4 - value,) in matching
            bbox = matching[(4 - value,)]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert 'x' in bbox
            assert bbox['x'] == (-(5 - value), (5 - value))
            assert len(bbox.intervals) == 1
            assert bbox.ignored == []

        # Real fix input of slicing input
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ('x', 'y', 'slit_id')
        bounding_boxes = {(0,): ((-0.5, 1047.5), (-0.5, 2047.5)), (1,): ((-0.5, 3047.5), (-0.5, 4047.5)), }
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args=[('slit_id', True)], order='F')

        matching = bounding_box._fix_inputs_matching_bounding_boxes('slit_id', 0)
        assert isinstance(matching, dict)
        assert () in matching
        bbox = matching[()]
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored == []
        assert bbox.intervals == {'x': (-0.5, 1047.5),
                                  'y': (-0.5, 2047.5)}
        assert bbox.order == 'F'

        matching = bounding_box._fix_inputs_matching_bounding_boxes('slit_id', 1)
        assert isinstance(matching, dict)
        assert () in matching
        bbox = matching[()]
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored == []
        assert bbox.intervals == {'x': (-0.5, 3047.5),
                                  'y': (-0.5, 4047.5)}
        assert bbox.order == 'F'

        # Errors
        with pytest.raises(ValueError) as err:
            bounding_box._fix_inputs_matching_bounding_boxes('slit_id', 2)
        assert str(err.value) ==\
            "Attempting to fix input slit_id, but there are no bounding boxes for argument value 2."

    def test__fix_input_selector_arg(self):
        # Single selector index
        selector_args = ((0, False),)
        bounding_boxes = {(1,): ((-1, 1), (-2, 2)), (2,): ((-2, 2), (-3, 3)), (3,): ((-3, 3), (-4, 4))}
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            bbox = bounding_box._fix_input_selector_arg('x', value)
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert 'y' in bbox
            assert bbox['y'] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == []

        # Multiple selector index
        selector_args = ((0, False), (1, False))
        bounding_boxes = {(1, 3): ((-1, 1), (-2, 2)), (2, 2): ((-2, 2), (-3, 3)), (3, 1): ((-3, 3), (-4, 4))}
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            bbox = bounding_box._fix_input_selector_arg('x', value)
            assert isinstance(bbox, CompoundBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert bbox.selector_args == (('y', False),)
            assert (4 - value,) in bbox
            bbox_selector = bbox[(4 - value,)]
            assert isinstance(bbox_selector, ModelBoundingBox)
            assert (bbox_selector._model.parameters == Gaussian2D().parameters).all()
            assert 'y' in bbox_selector
            assert bbox_selector['y'] == (-value, value)
            assert len(bbox_selector.intervals) == 1
            assert bbox_selector.ignored == []

            bbox = bounding_box._fix_input_selector_arg('y', value)
            assert isinstance(bbox, CompoundBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert bbox.selector_args == (('x', False),)
            assert (4 - value,) in bbox
            bbox_selector = bbox[(4 - value,)]
            assert isinstance(bbox_selector, ModelBoundingBox)
            assert (bbox_selector._model.parameters == Gaussian2D().parameters).all()
            assert 'x' in bbox_selector
            assert bbox_selector['x'] == (-(5 - value), (5 - value))
            assert len(bbox_selector.intervals) == 1
            assert bbox_selector.ignored == []

        # Real fix input of slicing input
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ('x', 'y', 'slit_id')
        bounding_boxes = {(0,): ((-0.5, 1047.5), (-0.5, 2047.5)), (1,): ((-0.5, 3047.5), (-0.5, 4047.5)), }
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args=[('slit_id', True)], order='F')

        bbox = bounding_box._fix_input_selector_arg('slit_id', 0)
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored == []
        assert bbox.intervals == {'x': (-0.5, 1047.5),
                                  'y': (-0.5, 2047.5)}
        assert bbox.order == 'F'

        bbox = bounding_box._fix_input_selector_arg('slit_id', 1)
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored == []
        assert bbox.intervals == {'x': (-0.5, 3047.5),
                                  'y': (-0.5, 4047.5)}
        assert bbox.order == 'F'

    def test__fix_input_bbox_arg(self):
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ('x', 'y', 'slit_id')
        bounding_boxes = {(0,): ((-0.5, 1047.5), (-0.5, 2047.5)), (1,): ((-0.5, 3047.5), (-0.5, 4047.5)), }
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args=[('slit_id', True)], order='F')

        bbox = bounding_box._fix_input_bbox_arg('x', 5)
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == (('slit_id', True),)
        assert bbox._bounding_boxes[(0,)] == (-0.5, 2047.5)
        assert bbox._bounding_boxes[(1,)] == (-0.5, 4047.5)
        assert len(bbox._bounding_boxes) == 2

        bbox = bounding_box._fix_input_bbox_arg('y', 5)
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == (('slit_id', True),)
        assert bbox._bounding_boxes[(0,)] == (-0.5, 1047.5)
        assert bbox._bounding_boxes[(1,)] == (-0.5, 3047.5)
        assert len(bbox._bounding_boxes) == 2

    def test_fix_inputs(self):
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ('x', 'y', 'slit_id')
        bounding_boxes = {(0,): ((-0.5, 1047.5), (-0.5, 2047.5)), (1,): ((-0.5, 3047.5), (-0.5, 4047.5)), }
        bounding_box = CompoundBoundingBox.validate(model, bounding_boxes, selector_args=[('slit_id', True)], order='F')
        model.bounding_box = bounding_box

        # Fix selector argument
        new_model = fix_inputs(model, {'slit_id': 0})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.intervals == {'x': (-0.5, 1047.5),
                                  'y': (-0.5, 2047.5)}
        assert bbox.order == 'F'

        # Fix a bounding_box field
        new_model = fix_inputs(model, {'x': 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == (('slit_id', True),)
        assert bbox._bounding_boxes[(0,)] == (-0.5, 2047.5)
        assert bbox._bounding_boxes[(0,)].order == 'F'
        assert bbox._bounding_boxes[(1,)] == (-0.5, 4047.5)
        assert bbox._bounding_boxes[(1,)].order == 'F'
        assert len(bbox._bounding_boxes) == 2
        new_model = fix_inputs(model, {'y': 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == (('slit_id', True),)
        assert bbox._bounding_boxes[(0,)] == (-0.5, 1047.5)
        assert bbox._bounding_boxes[(0,)].order == 'F'
        assert bbox._bounding_boxes[(1,)] == (-0.5, 3047.5)
        assert bbox._bounding_boxes[(1,)].order == 'F'
        assert len(bbox._bounding_boxes) == 2

        # Fix selector argument and a bounding_box field
        new_model = fix_inputs(model, {'slit_id': 0, 'x': 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.intervals == {'y': (-0.5, 2047.5)}
        assert bbox.order == 'F'
        new_model = fix_inputs(model, {'y': 5, 'slit_id': 1})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.intervals == {'x': (-0.5, 3047.5)}
        assert bbox.order == 'F'

        # Fix two bounding_box fields
        new_model = fix_inputs(model, {'x': 5, 'y': 7})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert bbox.selector_args == (('slit_id', True),)
        assert bbox._bounding_boxes[(0,)] == (-np.inf, np.inf)
        assert bbox._bounding_boxes[(0,)].order == 'F'
        assert bbox._bounding_boxes[(1,)] == (-np.inf, np.inf)
        assert bbox._bounding_boxes[(1,)].order == 'F'
        assert len(bbox._bounding_boxes) == 2
