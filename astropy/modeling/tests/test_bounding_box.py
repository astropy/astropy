# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
import unittest.mock as mk

from astropy.modeling.bounding_box import (_BaseInterval, Interval, BoundingBox)
from astropy.modeling.models import Gaussian1D, Gaussian2D
import astropy.units as u


class TestInterval:
    def test_create(self):
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = Interval(lower, upper)
        assert isinstance(interval, _BaseInterval)
        assert interval.lower == lower
        assert interval.upper == upper
        assert interval == (lower, upper)

        assert interval.__repr__() == \
            f"Interval(lower={lower}, upper={upper})"

    def test__validate_shape(self):
        message = "An interval must be some sort of sequence of length 2"
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = Interval(lower, upper)

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
        assert Interval._validate_bounds(1, 2) == (1, 2)
        assert Interval._validate_bounds(1*u.m, 2*u.m) == (1*u.m, 2*u.m)

        interval = Interval._validate_bounds(np.array([-2.5, -3.5]), np.array([2.5, 3.5]))
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()

        # Fails
        with pytest.warns(RuntimeWarning,
                          match="Invalid interval: upper bound 1 is strictly "
                          "less than lower bound 2."):
            Interval._validate_bounds(2, 1)
        with pytest.warns(RuntimeWarning,
                          match=r"Invalid interval: upper bound 1\.0 m is strictly "
                          r"less than lower bound 2\.0 m\."):
            Interval._validate_bounds(2*u.m, 1*u.m)

    def test_validate(self):
        # Passes
        assert Interval.validate((1, 2)) == (1, 2)
        assert Interval.validate([1, 2]) == (1, 2)
        assert Interval.validate((1*u.m, 2*u.m)) == (1*u.m, 2*u.m)
        assert Interval.validate([1*u.m, 2*u.m]) == (1*u.m, 2*u.m)

        assert Interval.validate(((1, 2),)) == (1, 2)
        assert Interval.validate(([1, 2],)) == (1, 2)
        assert Interval.validate([(1, 2)]) == (1, 2)
        assert Interval.validate([[1, 2]]) == (1, 2)
        assert Interval.validate(((1*u.m, 2*u.m),)) == (1*u.m, 2*u.m)
        assert Interval.validate(([1*u.m, 2*u.m],)) == (1*u.m, 2*u.m)
        assert Interval.validate([(1*u.m, 2*u.m)]) == (1*u.m, 2*u.m)
        assert Interval.validate([[1*u.m, 2*u.m]]) == (1*u.m, 2*u.m)

        interval = Interval.validate((np.array([-2.5, -3.5]),
                                      np.array([2.5, 3.5])))
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()
        interval = Interval.validate((np.array([-2.5, -3.5, -4.5]),
                                     np.array([2.5, 3.5, 4.5])))
        assert (interval.lower == np.array([-2.5, -3.5, -4.5])).all()
        assert (interval.upper == np.array([2.5, 3.5, 4.5])).all()

        # Fail shape
        with pytest.raises(ValueError):
            Interval.validate((1, 2, 3))

        # Fail bounds
        with pytest.warns(RuntimeWarning):
            Interval.validate((2, 1))

    def test_outside(self):
        interval = Interval.validate((0, 1))

        assert (interval.outside(np.linspace(-1, 2, 13)) ==
                [True, True, True, True,
                 False, False, False, False, False,
                 True, True, True, True]).all()

    def test_domain(self):
        interval = Interval.validate((0, 1))
        assert (interval.domain(0.25) == np.linspace(0, 1, 5)).all()


class TestBoundingBox:
    def test_create(self):
        intervals = ()
        model = mk.MagicMock()
        bounding_box = BoundingBox(intervals, model)

        assert bounding_box._intervals == {}
        assert bounding_box._model == model
        assert bounding_box._order == 'C'

        intervals = {}
        model = mk.MagicMock()
        bounding_box = BoundingBox(intervals, model, 'test')

        assert bounding_box._intervals == {}
        assert bounding_box._model == model
        assert bounding_box._order == 'test'

        intervals = (1, 2)
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = BoundingBox(intervals, model)

        assert bounding_box._intervals == {0: (1, 2)}
        assert bounding_box._model == model

    def test_intervals(self):
        intervals = {0: Interval(1, 2)}
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = BoundingBox(intervals, model)

        assert bounding_box._intervals == intervals
        assert bounding_box.intervals == intervals

    def test__get_name(self):
        intervals = {0: Interval(1, 2)}
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ['x']
        bounding_box = BoundingBox(intervals, model)

        index = mk.MagicMock()
        name = mk.MagicMock()
        model.inputs = mk.MagicMock()
        model.inputs.__getitem__.return_value = name
        assert bounding_box._get_name(index) == name
        assert model.inputs.__getitem__.call_args_list == [mk.call(index)]

    def test_named_intervals(self):
        intervals = {idx: Interval(idx, idx + 1) for idx in range(4)}
        model = mk.MagicMock()
        model.n_inputs = 4
        model.inputs = [mk.MagicMock() for _ in range(4)]
        bounding_box = BoundingBox(intervals, model)

        named = bounding_box.named_intervals
        assert isinstance(named, dict)
        for name, interval in named.items():
            assert name in model.inputs
            assert intervals[model.inputs.index(name)] == interval
        for index, name in enumerate(model.inputs):
            assert index in intervals
            assert name in named
            assert intervals[index] == named[name]

    def test___repr__(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        assert bounding_box.__repr__() ==\
            "BoundingBox(\n" +\
            "    intervals={\n" +\
            "        x: Interval(lower=-1, upper=1)\n" +\
            "        y: Interval(lower=-4, upper=4)\n" +\
            "    }\n" +\
            "    model=Gaussian2D(inputs=('x', 'y'))\n" +\
            "    order='C'\n" +\
            ")"

    def test___call__(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        args = tuple([mk.MagicMock() for _ in range(3)])
        kwargs = {f"test{idx}": mk.MagicMock() for idx in range(3)}

        with pytest.raises(RuntimeError) as err:
            bounding_box(*args, **kwargs)
        assert str(err.value) ==\
            "This bounding box is fixed by the model and does not have " +\
            "adjustable parameters."

    def test__get_index(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(2)
        assert str(err.value) ==\
            "Integer key: 2 must be < 2."
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(np.int32(2))
        assert str(err.value) ==\
            "Integer key: 2 must be < 2."
        with pytest.raises(IndexError) as err:
            bounding_box._get_index(np.int64(2))
        assert str(err.value) ==\
            "Integer key: 2 must be < 2."

        # Pass invalid key
        value = mk.MagicMock()
        with pytest.raises(ValueError) as err:
            bounding_box._get_index(value)
        assert str(err.value) ==\
            f"Key value: {value} must be string or integer."

    def test___len__(self):
        intervals = {0: Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert len(bounding_box) == 1 == len(bounding_box._intervals)

        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert len(bounding_box) == 2 == len(bounding_box._intervals)

        bounding_box._intervals = {}
        assert len(bounding_box) == 0 == len(bounding_box._intervals)

    def test___contains__(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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

    def test___getitem__(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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

    def test__get_order(self):
        intervals = {0: Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = BoundingBox.validate(model, intervals)

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

    def test_bounding_box(self):
        # 1D
        intervals = {0: Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert bounding_box.bounding_box() == (-1, 1)
        assert bounding_box.bounding_box(mk.MagicMock()) == (-1, 1)

        # > 1D
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert bounding_box.bounding_box() == ((-4, 4), (-1, 1))
        assert bounding_box.bounding_box('C') == ((-4, 4), (-1, 1))
        assert bounding_box.bounding_box('F') == ((-1, 1), (-4, 4))

    def test___eq__(self):
        intervals = {0: Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = BoundingBox.validate(model.copy(), intervals.copy())

        assert bounding_box == bounding_box
        assert bounding_box == BoundingBox.validate(model.copy(), intervals.copy())
        assert bounding_box == (-1, 1)

        assert not (bounding_box == mk.MagicMock())
        assert not (bounding_box == (-2, 2))
        assert not (bounding_box == BoundingBox.validate(model, {0: Interval(-2, 2)}))

        # Respect ordering
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box_1 = BoundingBox.validate(model, intervals)
        bounding_box_2 = BoundingBox.validate(model, intervals, order='F')
        assert bounding_box_1._order == 'C'
        assert bounding_box_1 == ((-4, 4), (-1, 1))
        assert not (bounding_box_1 == ((-1, 1), (-4, 4)))

        assert bounding_box_2._order == 'F'
        assert not (bounding_box_2 == ((-4, 4), (-1, 1)))
        assert bounding_box_2 == ((-1, 1), (-4, 4))

        assert bounding_box_1 == bounding_box_2

    def test__setitem__(self):
        model = Gaussian2D()
        bounding_box = BoundingBox({}, model)

        # USING Intervals directly
        # Set interval using key
        assert 'x' not in bounding_box
        bounding_box['x'] = Interval(-1, 1)
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], Interval)
        assert bounding_box['x'] == (-1, 1)

        assert 'y' not in bounding_box
        bounding_box['y'] = Interval(-4, 4)
        assert 'y' in bounding_box
        assert isinstance(bounding_box['y'], Interval)
        assert bounding_box['y'] == (-4, 4)

        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = Interval(-1, 1)
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], Interval)
        assert bounding_box[0] == (-1, 1)

        assert 1 not in bounding_box
        bounding_box[1] = Interval(-4, 4)
        assert 1 in bounding_box
        assert isinstance(bounding_box[1], Interval)
        assert bounding_box[1] == (-4, 4)

        # USING tuples
        # Set interval using key
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        bounding_box['x'] = (-1, 1)
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], Interval)
        assert bounding_box['x'] == (-1, 1)

        assert 'y' not in bounding_box
        bounding_box['y'] = (-4, 4)
        assert 'y' in bounding_box
        assert isinstance(bounding_box['y'], Interval)
        assert bounding_box['y'] == (-4, 4)

        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = (-1, 1)
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], Interval)
        assert bounding_box[0] == (-1, 1)

        assert 1 not in bounding_box
        bounding_box[1] = (-4, 4)
        assert 1 in bounding_box
        assert isinstance(bounding_box[1], Interval)
        assert bounding_box[1] == (-4, 4)

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = BoundingBox({}, model)
        # USING Intervals directly
        # Set interval using key
        assert 'x' not in bounding_box
        bounding_box['x'] = Interval(np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], Interval)
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = Interval(np.array([-1, -2]), np.array([1, 2]))
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()
        # USING tuples
        # Set interval using key
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        bounding_box['x'] = (np.array([-1, -2]), np.array([1, 2]))
        assert 'x' in bounding_box
        assert isinstance(bounding_box['x'], Interval)
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = (np.array([-1, -2]), np.array([1, 2]))
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()

    def test___delitem__(self):
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        # Using index
        assert 0 in bounding_box
        assert 'x' in bounding_box
        del bounding_box[0]
        assert 0 not in bounding_box
        assert 'x' not in bounding_box

        # Using key
        assert 1 in bounding_box
        assert 'y' in bounding_box
        del bounding_box['y']
        assert 1 not in bounding_box
        assert 'y' not in bounding_box

    def test__validate_dict(self):
        model = Gaussian2D()
        bounding_box = BoundingBox({}, model)

        # Input name keys
        intervals = {'x': Interval(-1, 1), 'y': Interval(-4, 4)}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_dict(intervals)
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Input index
        bounding_box._intervals = {}
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        assert 0 not in bounding_box
        assert 1 not in bounding_box
        bounding_box._validate_dict(intervals)
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = BoundingBox({}, model)
        # name keys
        intervals = {'x': Interval(np.array([-1, -2]), np.array([1, 2]))}
        assert 'x' not in bounding_box
        bounding_box._validate_dict(intervals)
        assert 'x' in bounding_box
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()
        # input index
        bounding_box._intervals = {}
        intervals = {0: Interval(np.array([-1, -2]), np.array([1, 2]))}
        assert 0 not in bounding_box
        bounding_box._validate_dict(intervals)
        assert 0 in bounding_box
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()

    def test__validate_sequence(self):
        model = Gaussian2D()
        bounding_box = BoundingBox({}, model)

        # Default order
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)))
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # C order
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)), order='C')
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Fortran order
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)), order='F')
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-4, 4)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Invalid order
        bounding_box._intervals = {}
        order = mk.MagicMock()
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        with pytest.raises(ValueError):
            bounding_box._validate_sequence(((-4, 4), (-1, 1)), order=order)
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        assert len(bounding_box.intervals) == 0

    def test__validate_iterable(self):
        model = Gaussian2D()
        bounding_box = BoundingBox({}, model)

        # Pass sequence Default order
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_iterable(((-4, 4), (-1, 1)))
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate_iterable(((-4, 4), (-1, 1)), order='F')
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-4, 4)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        bounding_box._intervals = {}
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        assert 0 not in bounding_box
        assert 1 not in bounding_box
        bounding_box._validate_iterable(intervals)
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Invalid iterable
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        with pytest.raises(ValueError) as err:
            bounding_box._validate_iterable(((-4, 4), (-1, 1), (-3, 3)))
        assert str(err.value) ==\
            "Found 3 intervals, but must have exactly 2."
        assert len(bounding_box.intervals) == 0
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        intervals = {0: Interval(-1, 1)}
        with pytest.raises(ValueError) as err:
            bounding_box._validate_iterable(intervals)
        assert str(err.value) ==\
            "Found 1 intervals, but must have exactly 2."
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        assert len(bounding_box.intervals) == 0

    def test__validate(self):
        model = Gaussian2D()
        bounding_box = BoundingBox({}, model)

        # Pass sequence Default order
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate(((-4, 4), (-1, 1)))
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box._intervals = {}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate(((-4, 4), (-1, 1)), order='F')
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-4, 4)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        bounding_box._intervals = {}
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        assert 'x' not in bounding_box
        assert 'y' not in bounding_box
        bounding_box._validate(intervals)
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass single
        model = Gaussian1D()
        bounding_box = BoundingBox({}, model)

        assert 'x' not in bounding_box
        bounding_box._validate((-1, 1))
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert len(bounding_box.intervals) == 1

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = BoundingBox({}, model)
        sequence = (np.array([-1, -2]), np.array([1, 2]))
        assert 'x' not in bounding_box
        bounding_box._validate(sequence)
        assert 'x' in bounding_box
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()

    def test_validate(self):
        model = Gaussian2D()

        # Pass sequence Default order
        bounding_box = BoundingBox.validate(model, ((-4, 4), (-1, 1)))
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box = BoundingBox.validate(model, ((-4, 4), (-1, 1)), order='F')
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-4, 4)
        assert 'y' in bounding_box
        assert bounding_box['y'] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        bounding_box = BoundingBox.validate(model, intervals)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass BoundingBox
        bbox = bounding_box
        bounding_box = BoundingBox.validate(model, bbox)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass single
        bounding_box = BoundingBox.validate(Gaussian1D(), (-1, 1))
        assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert 'x' in bounding_box
        assert bounding_box['x'] == (-1, 1)
        assert len(bounding_box.intervals) == 1

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        sequence = (np.array([-1, -2]), np.array([1, 2]))
        bounding_box = BoundingBox.validate(model, sequence)
        assert 'x' in bounding_box
        assert (bounding_box['x'].lower == np.array([-1, -2])).all()
        assert (bounding_box['x'].upper == np.array([1, 2])).all()

    def test_copy(self):
        bounding_box = BoundingBox.validate(Gaussian2D(), ((-4, 4), (-1, 1)))

        new_bounding_box = bounding_box.copy()
        assert bounding_box == new_bounding_box
        assert id(bounding_box) != id(new_bounding_box)

    def test_fix_inputs(self):
        bounding_box = BoundingBox.validate(Gaussian2D(), ((-4, 4), (-1, 1)))

        new_bounding_box = bounding_box.fix_inputs(Gaussian1D(), [1])
        assert not (bounding_box == new_bounding_box)

        assert (new_bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert 'x' in new_bounding_box
        assert new_bounding_box['x'] == (-1, 1)
        assert 'y' not in new_bounding_box
        assert len(new_bounding_box.intervals) == 1

    def test_dimension(self):
        intervals = {0: Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert bounding_box.dimension == 1 == len(bounding_box._intervals)

        intervals = {0: Interval(-1, 1), 1: Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)
        assert bounding_box.dimension == 2 == len(bounding_box._intervals)

        bounding_box._intervals = {}
        assert bounding_box.dimension == 0 == len(bounding_box._intervals)

    def test_domain(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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

    def test__outside(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
        assert (outside_index == [False for _ in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
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
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
        assert (outside_index == [True for _ in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
        assert (outside_index == [False]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
        assert (outside_index == [True]).all()
        assert all_out and isinstance(all_out, bool)

    def test__valid_index(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs)
        assert len(valid_index) == 1
        assert (valid_index[0] == [idx for idx in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs)
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs)
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs)
        assert len(valid_index) == 1
        assert (valid_index[0] == [0]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        valid_index, all_out = bounding_box._valid_index(input_shape, inputs)
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

    def test_prepare_inputs(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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
        assert new_inputs == []
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
        assert new_inputs == []
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

    def test__base_ouput(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

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
        input_shape = (13,7)
        output = bounding_box._base_output(input_shape, 0)
        assert (output == 0).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, np.nan)
        assert (np.isnan(output)).all()
        assert output.shape == input_shape
        output = bounding_box._base_output(input_shape, 14)
        assert (output == 14).all()
        assert output.shape == input_shape

    def test__modify_output(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        # Simple shape
        with mk.patch.object(BoundingBox, '_base_output', autospec=True,
                             return_value=np.asanyarray(0)) as mkBase:
            assert (np.array([1, 2, 3]) ==
                    bounding_box._modify_output([1, 2, 3], valid_index, input_shape, fill_value)).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

        # Replacement
        with mk.patch.object(BoundingBox, '_base_output', autospec=True,
                             return_value=np.array([1, 2, 3, 4, 5, 6])) as mkBase:
            assert (np.array([7, 2, 8, 4, 9, 6]) ==
                    bounding_box._modify_output([7, 8, 9], np.array([[0, 2, 4]]), input_shape, fill_value)).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

    def test__prepare_outputs(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        valid_outputs = [mk.MagicMock() for _ in range(3)]
        effects = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(BoundingBox, '_modify_output', autospec=True,
                             side_effect=effects) as mkModify:
            assert effects == bounding_box._prepare_outputs(valid_outputs, valid_index, input_shape, fill_value)
            assert mkModify.call_args_list == \
                [mk.call(bounding_box, valid_outputs[idx], valid_index, input_shape, fill_value)
                 for idx in range(3)]

    def test_prepare_outputs(self):
        intervals = {0: Interval(-1, 1), 1: Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = BoundingBox.validate(model, intervals)

        valid_outputs = mk.MagicMock()
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        with mk.patch.object(BoundingBox, '_prepare_outputs', autospec=True) as mkPrepare:
            # Reshape valid_outputs
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
