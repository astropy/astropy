# Licensed under a 3-clause BSD style license - see LICENSE.rst
import unittest.mock as mk

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SpectralCoord
from astropy.modeling.bounding_box import (
    CompoundBoundingBox,
    ModelBoundingBox,
    _BaseInterval,
    _BaseSelectorArgument,
    _BoundingDomain,
    _ignored_interval,
    _Interval,
    _SelectorArgument,
    _SelectorArguments,
)
from astropy.modeling.core import Model, fix_inputs
from astropy.modeling.models import (
    Gaussian1D,
    Gaussian2D,
    Identity,
    Polynomial2D,
    Scale,
    Shift,
)


class Test_Interval:
    def test_create(self):
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = _Interval(lower, upper)
        assert isinstance(interval, _BaseInterval)
        assert interval.lower == lower
        assert interval.upper == upper
        assert interval == (lower, upper)

        assert interval.__repr__() == f"Interval(lower={lower}, upper={upper})"

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
        MESSAGE = r"An interval must be some sort of sequence of length 2"
        lower = mk.MagicMock()
        upper = mk.MagicMock()
        interval = _Interval(lower, upper)

        # Passes (2,)
        interval._validate_shape((1, 2))
        interval._validate_shape([1, 2])
        interval._validate_shape((1 * u.m, 2 * u.m))
        interval._validate_shape([1 * u.m, 2 * u.m])

        # Passes (1, 2)
        interval._validate_shape(((1, 2),))
        interval._validate_shape(([1, 2],))
        interval._validate_shape([(1, 2)])
        interval._validate_shape([[1, 2]])
        interval._validate_shape(((1 * u.m, 2 * u.m),))
        interval._validate_shape(([1 * u.m, 2 * u.m],))
        interval._validate_shape([(1 * u.m, 2 * u.m)])
        interval._validate_shape([[1 * u.m, 2 * u.m]])

        # Passes (2, 0)
        interval._validate_shape((mk.MagicMock(), mk.MagicMock()))
        interval._validate_shape([mk.MagicMock(), mk.MagicMock()])

        # Passes with array inputs:
        interval._validate_shape((np.array([-2.5, -3.5]), np.array([2.5, 3.5])))
        interval._validate_shape(
            (np.array([-2.5, -3.5, -4.5]), np.array([2.5, 3.5, 4.5]))
        )

        # Fails shape (no units)
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape((1, 2, 3))
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape([1, 2, 3])
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape(1)

        # Fails shape (units)
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape((1 * u.m, 2 * u.m, 3 * u.m))
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape([1 * u.m, 2 * u.m, 3 * u.m])
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape(
                [[1 * u.m, 2 * u.m, 3 * u.m], [4 * u.m, 5 * u.m, 6 * u.m]]
            )
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape(1 * u.m)

        # Fails shape (arrays):
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape(
                (np.array([-2.5, -3.5]), np.array([2.5, 3.5]), np.array([3, 4]))
            )
        with pytest.raises(ValueError, match=MESSAGE):
            interval._validate_shape((np.array([-2.5, -3.5]), [2.5, 3.5]))

    def test__validate_bounds(self):
        # Passes
        assert _Interval._validate_bounds(1, 2) == (1, 2)
        assert _Interval._validate_bounds(1 * u.m, 2 * u.m) == (1 * u.m, 2 * u.m)

        interval = _Interval._validate_bounds(
            np.array([-2.5, -3.5]), np.array([2.5, 3.5])
        )
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()

        # Fails
        with pytest.warns(
            RuntimeWarning,
            match=r"Invalid interval: upper bound 1 is strictly "
            r"less than lower bound 2\.",
        ):
            _Interval._validate_bounds(2, 1)
        with pytest.warns(
            RuntimeWarning,
            match=r"Invalid interval: upper bound 1\.0 m is strictly "
            r"less than lower bound 2\.0 m\.",
        ):
            _Interval._validate_bounds(2 * u.m, 1 * u.m)

    def test_validate(self):
        # Passes
        assert _Interval.validate((1, 2)) == (1, 2)
        assert _Interval.validate([1, 2]) == (1, 2)
        assert _Interval.validate((1 * u.m, 2 * u.m)) == (1 * u.m, 2 * u.m)
        assert _Interval.validate([1 * u.m, 2 * u.m]) == (1 * u.m, 2 * u.m)

        assert _Interval.validate(((1, 2),)) == (1, 2)
        assert _Interval.validate(([1, 2],)) == (1, 2)
        assert _Interval.validate([(1, 2)]) == (1, 2)
        assert _Interval.validate([[1, 2]]) == (1, 2)
        assert _Interval.validate(((1 * u.m, 2 * u.m),)) == (1 * u.m, 2 * u.m)
        assert _Interval.validate(([1 * u.m, 2 * u.m],)) == (1 * u.m, 2 * u.m)
        assert _Interval.validate([(1 * u.m, 2 * u.m)]) == (1 * u.m, 2 * u.m)
        assert _Interval.validate([[1 * u.m, 2 * u.m]]) == (1 * u.m, 2 * u.m)

        interval = _Interval.validate((np.array([-2.5, -3.5]), np.array([2.5, 3.5])))
        assert (interval.lower == np.array([-2.5, -3.5])).all()
        assert (interval.upper == np.array([2.5, 3.5])).all()
        interval = _Interval.validate(
            (np.array([-2.5, -3.5, -4.5]), np.array([2.5, 3.5, 4.5]))
        )
        assert (interval.lower == np.array([-2.5, -3.5, -4.5])).all()
        assert (interval.upper == np.array([2.5, 3.5, 4.5])).all()

        # Fail shape
        MESSAGE = r"An interval must be some sort of sequence of length 2"
        with pytest.raises(ValueError, match=MESSAGE):
            _Interval.validate((1, 2, 3))

        # Fail bounds
        with pytest.warns(RuntimeWarning):
            _Interval.validate((2, 1))

    def test_outside(self):
        interval = _Interval.validate((0, 1))

        # fmt: off
        assert (
            interval.outside(np.linspace(-1, 2, 13))
            == [
                True, True, True, True,
                False, False, False, False, False,
                True, True, True, True
            ]
        ).all()
        # fmt: on

    def test_domain(self):
        interval = _Interval.validate((0, 1))
        assert (interval.domain(0.25) == np.linspace(0, 1, 5)).all()

    def test__ignored_interval(self):
        assert _ignored_interval.lower == -np.inf
        assert _ignored_interval.upper == np.inf

        for num in [0, -1, -100, 3.14, 10**100, -(10**100)]:
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
    def setup_method(self):
        class BoundingDomain(_BoundingDomain):
            def fix_inputs(self, model, fix_inputs):
                super().fix_inputs(model, fixed_inputs=fix_inputs)

            def prepare_inputs(self, input_shape, inputs):
                super().prepare_inputs(input_shape, inputs)

        self.BoundingDomain = BoundingDomain

    def test_create(self):
        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == "C"

        bounding_box = self.BoundingDomain(model, order="F")
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == "F"

        bounding_box = self.BoundingDomain(Gaussian2D(), ["x"])
        assert bounding_box._ignored == [0]
        assert bounding_box._order == "C"

        # Error
        MESSAGE = r"order must be either 'C' .* or 'F' .*, got: .*"
        with pytest.raises(ValueError, match=MESSAGE):
            self.BoundingDomain(model, order=mk.MagicMock())

    def test_model(self):
        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)
        assert bounding_box._model == model
        assert bounding_box.model == model

    def test_order(self):
        bounding_box = self.BoundingDomain(mk.MagicMock(), order="C")
        assert bounding_box._order == "C"
        assert bounding_box.order == "C"

        bounding_box = self.BoundingDomain(mk.MagicMock(), order="F")
        assert bounding_box._order == "F"
        assert bounding_box.order == "F"

        bounding_box._order = "test"
        assert bounding_box.order == "test"

    def test_ignored(self):
        ignored = [0]
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ["x"]
        bounding_box = self.BoundingDomain(model, ignored=ignored)

        assert bounding_box._ignored == ignored
        assert bounding_box.ignored == ignored

    def test__get_order(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # Success (default 'C')
        assert bounding_box._order == "C"
        assert bounding_box._get_order() == "C"
        assert bounding_box._get_order("C") == "C"
        assert bounding_box._get_order("F") == "F"

        # Success (default 'F')
        bounding_box._order = "F"
        assert bounding_box._order == "F"
        assert bounding_box._get_order() == "F"
        assert bounding_box._get_order("C") == "C"
        assert bounding_box._get_order("F") == "F"

        # Error
        MESSAGE = r"order must be either 'C' .* or 'F' .*, got: .*"
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box._get_order(mk.MagicMock())

    def test__get_index(self):
        bounding_box = self.BoundingDomain(Gaussian2D())

        # Pass input name
        assert bounding_box._get_index("x") == 0
        assert bounding_box._get_index("y") == 1

        # Pass invalid input name
        MESSAGE = r"'z' is not one of the inputs: .*"
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box._get_index("z")

        # Pass valid index
        assert bounding_box._get_index(0) == 0
        assert bounding_box._get_index(1) == 1
        assert bounding_box._get_index(np.int32(0)) == 0
        assert bounding_box._get_index(np.int32(1)) == 1
        assert bounding_box._get_index(np.int64(0)) == 0
        assert bounding_box._get_index(np.int64(1)) == 1

        # Pass invalid index
        MESSAGE = r"Integer key: .* must be non-negative and < 2"
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._get_index(2)
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._get_index(np.int32(2))
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._get_index(np.int64(2))
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._get_index(-1)

        # Pass invalid key
        MESSAGE = r"Key value: .* must be string or integer"
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box._get_index(mk.MagicMock())

    def test__get_name(self):
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ["x"]
        bounding_box = self.BoundingDomain(model)

        index = mk.MagicMock()
        name = mk.MagicMock()
        model.inputs = mk.MagicMock()
        model.inputs.__getitem__.return_value = name
        assert bounding_box._get_name(index) == name
        assert model.inputs.__getitem__.call_args_list == [mk.call(index)]

    def test_ignored_inputs(self):
        model = mk.MagicMock()
        ignored = list(range(4, 8))
        model.n_inputs = 8
        model.inputs = [mk.MagicMock() for _ in range(8)]
        bounding_box = self.BoundingDomain(model, ignored=ignored)

        inputs = bounding_box.ignored_inputs
        assert isinstance(inputs, list)
        for index, _input in enumerate(inputs):
            assert _input in model.inputs
            assert model.inputs[index + 4] == _input
        for index, _input in enumerate(model.inputs):
            if _input in inputs:
                assert inputs[index - 4] == _input
            else:
                assert index < 4

    def test__validate_ignored(self):
        bounding_box = self.BoundingDomain(Gaussian2D())

        # Pass
        assert bounding_box._validate_ignored(None) == []
        assert bounding_box._validate_ignored(["x", "y"]) == [0, 1]
        assert bounding_box._validate_ignored([0, 1]) == [0, 1]
        assert bounding_box._validate_ignored([np.int32(0), np.int64(1)]) == [0, 1]

        # Fail
        with pytest.raises(
            ValueError, match=r"Key value: .* must be string or integer"
        ):
            bounding_box._validate_ignored([mk.MagicMock()])
        with pytest.raises(ValueError, match=r"'.*' is not one of the inputs: .*"):
            bounding_box._validate_ignored(["z"])

        MESSAGE = r"Integer key: 3 must be non-negative and < 2"
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._validate_ignored([3])
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._validate_ignored([np.int32(3)])
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box._validate_ignored([np.int64(3)])

    def test___call__(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        args = tuple(mk.MagicMock() for _ in range(3))
        kwargs = {f"test{idx}": mk.MagicMock() for idx in range(3)}

        MESSAGE = (
            r"This bounding box is fixed by the model and does not have adjustable"
            r" parameters"
        )
        with pytest.raises(RuntimeError, match=MESSAGE):
            bounding_box(*args, **kwargs)

    def test_fix_inputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())
        model = mk.MagicMock()
        fixed_inputs = mk.MagicMock()

        with pytest.raises(
            NotImplementedError, match=r"This should be implemented by a child class"
        ):
            bounding_box.fix_inputs(model, fixed_inputs)

    def test__prepare_inputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        with pytest.raises(
            NotImplementedError,
            match=r"This has not been implemented for BoundingDomain",
        ):
            bounding_box.prepare_inputs(mk.MagicMock(), mk.MagicMock())

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
        with mk.patch.object(
            _BoundingDomain,
            "_base_output",
            autospec=True,
            return_value=np.asanyarray(0),
        ) as mkBase:
            assert (
                np.array([1, 2, 3])
                == bounding_box._modify_output(
                    [1, 2, 3], valid_index, input_shape, fill_value
                )
            ).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

        # Replacement
        with mk.patch.object(
            _BoundingDomain,
            "_base_output",
            autospec=True,
            return_value=np.array([1, 2, 3, 4, 5, 6]),
        ) as mkBase:
            assert (
                np.array([7, 2, 8, 4, 9, 6])
                == bounding_box._modify_output(
                    [7, 8, 9], np.array([[0, 2, 4]]), input_shape, fill_value
                )
            ).all()
            assert mkBase.call_args_list == [mk.call(input_shape, fill_value)]

    def test__prepare_outputs(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        valid_outputs = [mk.MagicMock() for _ in range(3)]
        effects = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(
            _BoundingDomain, "_modify_output", autospec=True, side_effect=effects
        ) as mkModify:
            assert effects == bounding_box._prepare_outputs(
                valid_outputs, valid_index, input_shape, fill_value
            )
            assert mkModify.call_args_list == [
                mk.call(
                    bounding_box,
                    valid_outputs[idx],
                    valid_index,
                    input_shape,
                    fill_value,
                )
                for idx in range(3)
            ]

    def test_prepare_outputs(self):
        model = mk.MagicMock()
        bounding_box = self.BoundingDomain(model)

        valid_outputs = mk.MagicMock()
        valid_index = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()

        with mk.patch.object(
            _BoundingDomain, "_prepare_outputs", autospec=True
        ) as mkPrepare:
            # Reshape valid_outputs
            model.n_outputs = 1
            assert mkPrepare.return_value == bounding_box.prepare_outputs(
                valid_outputs, valid_index, input_shape, fill_value
            )
            assert mkPrepare.call_args_list == [
                mk.call(
                    bounding_box, [valid_outputs], valid_index, input_shape, fill_value
                )
            ]
            mkPrepare.reset_mock()

            # No reshape valid_outputs
            model.n_outputs = 2
            assert mkPrepare.return_value == bounding_box.prepare_outputs(
                valid_outputs, valid_index, input_shape, fill_value
            )
            assert mkPrepare.call_args_list == [
                mk.call(
                    bounding_box, valid_outputs, valid_index, input_shape, fill_value
                )
            ]

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

        with mk.patch.object(
            _BoundingDomain, "_get_valid_outputs_unit", autospec=True
        ) as mkGet:
            with mk.patch.object(
                _BoundingDomain, "prepare_outputs", autospec=True
            ) as mkPrepare:
                assert bounding_box._evaluate_model(
                    evaluate,
                    valid_inputs,
                    valid_index,
                    input_shape,
                    fill_value,
                    with_units,
                ) == (mkPrepare.return_value, mkGet.return_value)
                assert mkPrepare.call_args_list == [
                    mk.call(
                        bounding_box,
                        evaluate.return_value,
                        valid_index,
                        input_shape,
                        fill_value,
                    )
                ]
                assert mkGet.call_args_list == [
                    mk.call(evaluate.return_value, with_units)
                ]
                assert evaluate.call_args_list == [mk.call(valid_inputs)]

    def test__evaluate(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        evaluate = mk.MagicMock()
        inputs = mk.MagicMock()
        input_shape = mk.MagicMock()
        fill_value = mk.MagicMock()
        with_units = mk.MagicMock()

        valid_inputs = mk.MagicMock()
        valid_index = mk.MagicMock()

        effects = [
            (valid_inputs, valid_index, True),
            (valid_inputs, valid_index, False),
        ]
        with mk.patch.object(
            self.BoundingDomain, "prepare_inputs", autospec=True, side_effect=effects
        ) as mkPrepare:
            with mk.patch.object(
                _BoundingDomain, "_all_out_output", autospec=True
            ) as mkAll:
                with mk.patch.object(
                    _BoundingDomain, "_evaluate_model", autospec=True
                ) as mkEvaluate:
                    # all_out
                    assert (
                        bounding_box._evaluate(
                            evaluate, inputs, input_shape, fill_value, with_units
                        )
                        == mkAll.return_value
                    )
                    assert mkAll.call_args_list == [
                        mk.call(bounding_box, input_shape, fill_value)
                    ]
                    assert mkEvaluate.call_args_list == []
                    assert mkPrepare.call_args_list == [
                        mk.call(bounding_box, input_shape, inputs)
                    ]

                    mkAll.reset_mock()
                    mkPrepare.reset_mock()

                    # not all_out
                    assert (
                        bounding_box._evaluate(
                            evaluate, inputs, input_shape, fill_value, with_units
                        )
                        == mkEvaluate.return_value
                    )
                    assert mkAll.call_args_list == []
                    assert mkEvaluate.call_args_list == [
                        mk.call(
                            bounding_box,
                            evaluate,
                            valid_inputs,
                            valid_index,
                            input_shape,
                            fill_value,
                            with_units,
                        )
                    ]
                    assert mkPrepare.call_args_list == [
                        mk.call(bounding_box, input_shape, inputs)
                    ]

    def test__set_outputs_unit(self):
        bounding_box = self.BoundingDomain(mk.MagicMock())

        # set no unit
        assert bounding_box._set_outputs_unit(27, None) == 27

        # set unit
        assert bounding_box._set_outputs_unit(27, u.m) == 27 * u.m

    def test_evaluate(self):
        bounding_box = self.BoundingDomain(Gaussian2D())

        evaluate = mk.MagicMock()
        inputs = mk.MagicMock()
        fill_value = mk.MagicMock()

        outputs = mk.MagicMock()
        valid_outputs_unit = mk.MagicMock()
        value = (outputs, valid_outputs_unit)
        with mk.patch.object(
            _BoundingDomain, "_evaluate", autospec=True, return_value=value
        ) as mkEvaluate:
            with mk.patch.object(
                _BoundingDomain, "_set_outputs_unit", autospec=True
            ) as mkSet:
                with mk.patch.object(Model, "input_shape", autospec=True) as mkShape:
                    with mk.patch.object(
                        Model, "bbox_with_units", new_callable=mk.PropertyMock
                    ) as mkUnits:
                        assert tuple(mkSet.return_value) == bounding_box.evaluate(
                            evaluate, inputs, fill_value
                        )
                        assert mkSet.call_args_list == [
                            mk.call(outputs, valid_outputs_unit)
                        ]
                        assert mkEvaluate.call_args_list == [
                            mk.call(
                                bounding_box,
                                evaluate,
                                inputs,
                                mkShape.return_value,
                                fill_value,
                                mkUnits.return_value,
                            )
                        ]
                        assert mkShape.call_args_list == [
                            mk.call(bounding_box._model, inputs)
                        ]
                        assert mkUnits.call_args_list == [mk.call()]


class TestModelBoundingBox:
    def test_create(self):
        intervals = ()
        model = mk.MagicMock()
        bounding_box = ModelBoundingBox(intervals, model)

        assert isinstance(bounding_box, _BoundingDomain)
        assert bounding_box._intervals == {}
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == "C"

        # Set optional
        intervals = {}
        model = mk.MagicMock()
        bounding_box = ModelBoundingBox(intervals, model, order="F")

        assert isinstance(bounding_box, _BoundingDomain)
        assert bounding_box._intervals == {}
        assert bounding_box._model == model
        assert bounding_box._ignored == []
        assert bounding_box._order == "F"

        # Set interval
        intervals = (1, 2)
        model = mk.MagicMock()
        model.n_inputs = 1
        model.inputs = ["x"]
        bounding_box = ModelBoundingBox(intervals, model)

        assert isinstance(bounding_box, _BoundingDomain)
        assert bounding_box._intervals == {0: (1, 2)}
        assert bounding_box._model == model

        # Set ignored
        intervals = (1, 2)
        model = mk.MagicMock()
        model.n_inputs = 2
        model.inputs = ["x", "y"]
        bounding_box = ModelBoundingBox(intervals, model, ignored=[1])

        assert isinstance(bounding_box, _BoundingDomain)
        assert bounding_box._intervals == {0: (1, 2)}
        assert bounding_box._model == model
        assert bounding_box._ignored == [1]

        intervals = ((1, 2), (3, 4))
        model = mk.MagicMock()
        model.n_inputs = 3
        model.inputs = ["x", "y", "z"]
        bounding_box = ModelBoundingBox(intervals, model, ignored=[2], order="F")

        assert isinstance(bounding_box, _BoundingDomain)
        assert bounding_box._intervals == {0: (1, 2), 1: (3, 4)}
        assert bounding_box._model == model
        assert bounding_box._ignored == [2]
        assert bounding_box._order == "F"

    def test_copy(self):
        bounding_box = ModelBoundingBox.validate(
            Gaussian2D(), ((-4.5, 4.5), (-1.4, 1.4))
        )
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
        model.inputs = ["x"]
        bounding_box = ModelBoundingBox(intervals, model)

        assert bounding_box._intervals == intervals
        assert bounding_box.intervals == intervals

    def test_named_intervals(self):
        intervals = {idx: _Interval(idx, idx + 1) for idx in range(4)}
        model = mk.MagicMock()
        model.n_inputs = 4
        model.inputs = [mk.MagicMock() for _ in range(4)]
        bounding_box = ModelBoundingBox(intervals, model)

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
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        assert (
            bounding_box.__repr__() == "ModelBoundingBox(\n"
            "    intervals={\n"
            "        x: Interval(lower=-1, upper=1)\n"
            "        y: Interval(lower=-4, upper=4)\n"
            "    }\n"
            "    model=Gaussian2D(inputs=('x', 'y'))\n"
            "    order='C'\n"
            ")"
        )

        intervals = {0: _Interval(-1, 1)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals, ignored=["y"])

        assert (
            bounding_box.__repr__() == "ModelBoundingBox(\n"
            "    intervals={\n"
            "        x: Interval(lower=-1, upper=1)\n"
            "    }\n"
            "    ignored=['y']\n"
            "    model=Gaussian2D(inputs=('x', 'y'))\n"
            "    order='C'\n"
            ")"
        )

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
        assert "x" in bounding_box
        assert "y" in bounding_box
        assert "z" not in bounding_box

        # Contains with index
        assert 0 in bounding_box
        assert 1 in bounding_box
        assert 2 not in bounding_box

        # General not in
        assert mk.MagicMock() not in bounding_box

        # Contains with ignored
        del bounding_box["y"]

        # Contains with keys
        assert "x" in bounding_box
        assert "y" in bounding_box
        assert "z" not in bounding_box

        # Contains with index
        assert 0 in bounding_box
        assert 1 in bounding_box
        assert 2 not in bounding_box

    def test___getitem__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Get using input key
        assert bounding_box["x"] == (-1, 1)
        assert bounding_box["y"] == (-4, 4)

        # Fail with input key
        with pytest.raises(ValueError, match=r"'.*' is not one of the inputs: .*"):
            bounding_box["z"]

        # Get using index
        assert bounding_box[0] == (-1, 1)
        assert bounding_box[1] == (-4, 4)
        assert bounding_box[np.int32(0)] == (-1, 1)
        assert bounding_box[np.int32(1)] == (-4, 4)
        assert bounding_box[np.int64(0)] == (-1, 1)
        assert bounding_box[np.int64(1)] == (-4, 4)

        # Fail with index
        MESSAGE = r"Integer key: 2 must be non-negative and < 2"
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box[2]
        with pytest.raises(IndexError, match=MESSAGE):
            bounding_box[np.int32(2)]
        with pytest.raises(IndexError, match=MESSAGE):
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
        bounding_box = ModelBoundingBox.validate(model, {}, ignored=["x"])
        assert bounding_box.bounding_box() == (-np.inf, np.inf)
        assert bounding_box.bounding_box("C") == (-np.inf, np.inf)
        assert bounding_box.bounding_box("F") == (-np.inf, np.inf)

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
        assert bounding_box.bounding_box("C") == ((-4, 4), (-1, 1))
        assert bounding_box.bounding_box("F") == ((-1, 1), (-4, 4))

    def test___eq__(self):
        intervals = {0: _Interval(-1, 1)}
        model = Gaussian1D()
        bounding_box = ModelBoundingBox.validate(model.copy(), intervals.copy())

        assert bounding_box == bounding_box
        assert bounding_box == ModelBoundingBox.validate(model.copy(), intervals.copy())
        assert bounding_box == (-1, 1)

        assert not (bounding_box == mk.MagicMock())
        assert not (bounding_box == (-2, 2))
        assert not (
            bounding_box == ModelBoundingBox.validate(model, {0: _Interval(-2, 2)})
        )

        # Respect ordering
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box_1 = ModelBoundingBox.validate(model, intervals)
        bounding_box_2 = ModelBoundingBox.validate(model, intervals, order="F")
        assert bounding_box_1._order == "C"
        assert bounding_box_1 == ((-4, 4), (-1, 1))
        assert not (bounding_box_1 == ((-1, 1), (-4, 4)))

        assert bounding_box_2._order == "F"
        assert not (bounding_box_2 == ((-4, 4), (-1, 1)))
        assert bounding_box_2 == ((-1, 1), (-4, 4))

        assert bounding_box_1 == bounding_box_2

        # Respect ignored
        model = Gaussian2D()
        bounding_box_1._ignored = [mk.MagicMock()]
        bounding_box_2._ignored = [mk.MagicMock()]
        assert bounding_box_1._ignored != bounding_box_2._ignored
        assert not (bounding_box_1 == bounding_box_2)

    def test__setitem__(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, {}, ignored=[0, 1])
        assert bounding_box._ignored == [0, 1]

        # USING Intervals directly
        # Set interval using key
        assert 0 not in bounding_box.intervals
        assert 0 in bounding_box.ignored
        bounding_box["x"] = _Interval(-1, 1)
        assert 0 in bounding_box.intervals
        assert 0 not in bounding_box.ignored
        assert isinstance(bounding_box["x"], _Interval)
        assert bounding_box["x"] == (-1, 1)

        assert 1 not in bounding_box.intervals
        assert 1 in bounding_box.ignored
        bounding_box["y"] = _Interval(-4, 4)
        assert 1 in bounding_box.intervals
        assert 1 not in bounding_box.ignored
        assert isinstance(bounding_box["y"], _Interval)
        assert bounding_box["y"] == (-4, 4)

        del bounding_box["x"]
        del bounding_box["y"]

        # Set interval using index
        assert 0 not in bounding_box.intervals
        assert 0 in bounding_box.ignored
        bounding_box[0] = _Interval(-1, 1)
        assert 0 in bounding_box.intervals
        assert 0 not in bounding_box.ignored
        assert isinstance(bounding_box[0], _Interval)
        assert bounding_box[0] == (-1, 1)

        assert 1 not in bounding_box.intervals
        assert 1 in bounding_box.ignored
        bounding_box[1] = _Interval(-4, 4)
        assert 1 in bounding_box.intervals
        assert 1 not in bounding_box.ignored
        assert isinstance(bounding_box[1], _Interval)
        assert bounding_box[1] == (-4, 4)

        del bounding_box[0]
        del bounding_box[1]

        # USING tuples
        # Set interval using key
        assert 0 not in bounding_box.intervals
        assert 0 in bounding_box.ignored
        bounding_box["x"] = (-1, 1)
        assert 0 in bounding_box.intervals
        assert 0 not in bounding_box.ignored
        assert isinstance(bounding_box["x"], _Interval)
        assert bounding_box["x"] == (-1, 1)

        assert 1 not in bounding_box.intervals
        assert 1 in bounding_box.ignored
        bounding_box["y"] = (-4, 4)
        assert 1 in bounding_box.intervals
        assert 1 not in bounding_box.ignored
        assert isinstance(bounding_box["y"], _Interval)
        assert bounding_box["y"] == (-4, 4)

        del bounding_box["x"]
        del bounding_box["y"]

        # Set interval using index
        assert 0 not in bounding_box.intervals
        assert 0 in bounding_box.ignored
        bounding_box[0] = (-1, 1)
        assert 0 in bounding_box.intervals
        assert 0 not in bounding_box.ignored
        assert isinstance(bounding_box[0], _Interval)
        assert bounding_box[0] == (-1, 1)

        assert 1 not in bounding_box.intervals
        assert 1 in bounding_box.ignored
        bounding_box[1] = (-4, 4)
        assert 1 in bounding_box.intervals
        assert 1 not in bounding_box.ignored
        assert isinstance(bounding_box[1], _Interval)
        assert bounding_box[1] == (-4, 4)

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = ModelBoundingBox({}, model)
        # USING Intervals directly
        # Set interval using key
        assert "x" not in bounding_box
        bounding_box["x"] = _Interval(np.array([-1, -2]), np.array([1, 2]))
        assert "x" in bounding_box
        assert isinstance(bounding_box["x"], _Interval)
        assert (bounding_box["x"].lower == np.array([-1, -2])).all()
        assert (bounding_box["x"].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = _Interval(np.array([-1, -2]), np.array([1, 2]))
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], _Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()
        # USING tuples
        # Set interval using key
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        bounding_box["x"] = (np.array([-1, -2]), np.array([1, 2]))
        assert "x" in bounding_box
        assert isinstance(bounding_box["x"], _Interval)
        assert (bounding_box["x"].lower == np.array([-1, -2])).all()
        assert (bounding_box["x"].upper == np.array([1, 2])).all()
        # Set interval using index
        bounding_box._intervals = {}
        assert 0 not in bounding_box
        bounding_box[0] = (np.array([-1, -2]), np.array([1, 2]))
        assert 0 in bounding_box
        assert isinstance(bounding_box[0], _Interval)
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()

    def test___delitem__(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Using index
        assert 0 in bounding_box.intervals
        assert 0 not in bounding_box.ignored
        assert 0 in bounding_box
        assert "x" in bounding_box
        del bounding_box[0]
        assert 0 not in bounding_box.intervals
        assert 0 in bounding_box.ignored
        assert 0 in bounding_box
        assert "x" in bounding_box

        # Delete an ignored item
        with pytest.raises(RuntimeError, match=r"Cannot delete ignored input: 0!"):
            del bounding_box[0]

        # Using key
        assert 1 in bounding_box.intervals
        assert 1 not in bounding_box.ignored
        assert 0 in bounding_box
        assert "y" in bounding_box
        del bounding_box["y"]
        assert 1 not in bounding_box.intervals
        assert 1 in bounding_box.ignored
        assert 0 in bounding_box
        assert "y" in bounding_box

        # Delete an ignored item
        with pytest.raises(RuntimeError, match=r"Cannot delete ignored input: y!"):
            del bounding_box["y"]

    def test__validate_dict(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox({}, model)

        # Input name keys
        intervals = {"x": _Interval(-1, 1), "y": _Interval(-4, 4)}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_dict(intervals)
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Input index
        bounding_box._intervals = {}
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
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
        bounding_box = ModelBoundingBox({}, model)
        # name keys
        intervals = {"x": _Interval(np.array([-1, -2]), np.array([1, 2]))}
        assert "x" not in bounding_box
        bounding_box._validate_dict(intervals)
        assert "x" in bounding_box
        assert (bounding_box["x"].lower == np.array([-1, -2])).all()
        assert (bounding_box["x"].upper == np.array([1, 2])).all()
        # input index
        bounding_box._intervals = {}
        intervals = {0: _Interval(np.array([-1, -2]), np.array([1, 2]))}
        assert 0 not in bounding_box
        bounding_box._validate_dict(intervals)
        assert 0 in bounding_box
        assert (bounding_box[0].lower == np.array([-1, -2])).all()
        assert (bounding_box[0].upper == np.array([1, 2])).all()

    def test__validate_sequence(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox({}, model)

        # Default order
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)))
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # C order
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)), order="C")
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Fortran order
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_sequence(((-4, 4), (-1, 1)), order="F")
        assert "x" in bounding_box
        assert bounding_box["x"] == (-4, 4)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Invalid order
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        assert "y" not in bounding_box

        MESSAGE = r"order must be either 'C' .* or 'F' .*, got: .*"
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box._validate_sequence(((-4, 4), (-1, 1)), order=mk.MagicMock())
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        assert len(bounding_box.intervals) == 0

    def test__n_inputs(self):
        model = Gaussian2D()

        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        bounding_box = ModelBoundingBox.validate(model, intervals)
        assert bounding_box._n_inputs == 2

        intervals = {0: _Interval(-1, 1)}
        bounding_box = ModelBoundingBox.validate(model, intervals, ignored=["y"])
        assert bounding_box._n_inputs == 1

        bounding_box = ModelBoundingBox.validate(model, {}, ignored=["x", "y"])
        assert bounding_box._n_inputs == 0

        bounding_box._ignored = ["x", "y", "z"]
        assert bounding_box._n_inputs == 0

    def test__validate_iterable(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox({}, model)

        # Pass sequence Default order
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_iterable(((-4, 4), (-1, 1)))
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate_iterable(((-4, 4), (-1, 1)), order="F")
        assert "x" in bounding_box
        assert bounding_box["x"] == (-4, 4)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        bounding_box._intervals = {}
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        assert 0 not in bounding_box
        assert 1 not in bounding_box
        bounding_box._validate_iterable(intervals)
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass with ignored
        bounding_box._intervals = {}
        bounding_box._ignored = [1]
        intervals = {0: _Interval(-1, 1)}
        assert 0 not in bounding_box.intervals
        bounding_box._validate_iterable(intervals)
        assert 0 in bounding_box.intervals
        assert bounding_box[0] == (-1, 1)

        # Invalid iterable
        MESSAGE = "Found {} intervals, but must have exactly {}"
        bounding_box._intervals = {}
        bounding_box._ignored = []
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        with pytest.raises(ValueError, match=MESSAGE.format(3, 2)):
            bounding_box._validate_iterable(((-4, 4), (-1, 1), (-3, 3)))
        assert len(bounding_box.intervals) == 0
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._ignored = [1]
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        with pytest.raises(ValueError, match=MESSAGE.format(2, 1)):
            bounding_box._validate_iterable(intervals)
        assert len(bounding_box.intervals) == 0
        bounding_box._ignored = []
        intervals = {0: _Interval(-1, 1)}
        with pytest.raises(ValueError, match=MESSAGE.format(1, 2)):
            bounding_box._validate_iterable(intervals)
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        assert len(bounding_box.intervals) == 0

    def test__validate(self):
        model = Gaussian2D()
        bounding_box = ModelBoundingBox({}, model)

        # Pass sequence Default order
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate(((-4, 4), (-1, 1)))
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box._intervals = {}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate(((-4, 4), (-1, 1)), order="F")
        assert "x" in bounding_box
        assert bounding_box["x"] == (-4, 4)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        bounding_box._intervals = {}
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        assert "x" not in bounding_box
        assert "y" not in bounding_box
        bounding_box._validate(intervals)
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass single with ignored
        intervals = {0: _Interval(-1, 1)}
        bounding_box = ModelBoundingBox({}, model, ignored=[1])

        assert 0 not in bounding_box.intervals
        assert 1 not in bounding_box.intervals
        bounding_box._validate(intervals)
        assert 0 in bounding_box.intervals
        assert bounding_box[0] == (-1, 1)
        assert 1 not in bounding_box.intervals
        assert len(bounding_box.intervals) == 1

        # Pass single
        model = Gaussian1D()
        bounding_box = ModelBoundingBox({}, model)

        assert "x" not in bounding_box
        bounding_box._validate((-1, 1))
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert len(bounding_box.intervals) == 1

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        bounding_box = ModelBoundingBox({}, model)
        sequence = (np.array([-1, -2]), np.array([1, 2]))
        assert "x" not in bounding_box
        bounding_box._validate(sequence)
        assert "x" in bounding_box
        assert (bounding_box["x"].lower == np.array([-1, -2])).all()
        assert (bounding_box["x"].upper == np.array([1, 2])).all()

    def test_validate(self):
        model = Gaussian2D()
        kwargs = {"test": mk.MagicMock()}

        # Pass sequence Default order
        bounding_box = ModelBoundingBox.validate(model, ((-4, 4), (-1, 1)), **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-4, 4)
        assert len(bounding_box.intervals) == 2

        # Pass sequence
        bounding_box = ModelBoundingBox.validate(
            model, ((-4, 4), (-1, 1)), order="F", **kwargs
        )
        assert (bounding_box._model.parameters == model.parameters).all()
        assert "x" in bounding_box
        assert bounding_box["x"] == (-4, 4)
        assert "y" in bounding_box
        assert bounding_box["y"] == (-1, 1)
        assert len(bounding_box.intervals) == 2

        # Pass Dict
        intervals = {0: _Interval(-1, 1), 1: _Interval(-4, 4)}
        bounding_box = ModelBoundingBox.validate(model, intervals, order="F", **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2
        assert bounding_box.order == "F"

        # Pass ModelBoundingBox
        bbox = bounding_box
        bounding_box = ModelBoundingBox.validate(model, bbox, **kwargs)
        assert (bounding_box._model.parameters == model.parameters).all()
        assert 0 in bounding_box
        assert bounding_box[0] == (-1, 1)
        assert 1 in bounding_box
        assert bounding_box[1] == (-4, 4)
        assert len(bounding_box.intervals) == 2
        assert bounding_box.order == "F"

        # Pass single ignored
        intervals = {0: _Interval(-1, 1)}
        bounding_box = ModelBoundingBox.validate(
            model, intervals, ignored=["y"], **kwargs
        )
        assert (bounding_box._model.parameters == model.parameters).all()
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert "y" in bounding_box
        assert bounding_box["y"] == _ignored_interval
        assert len(bounding_box.intervals) == 1

        # Pass single
        bounding_box = ModelBoundingBox.validate(Gaussian1D(), (-1, 1), **kwargs)
        assert (bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert "x" in bounding_box
        assert bounding_box["x"] == (-1, 1)
        assert len(bounding_box.intervals) == 1

        # Model set support
        model = Gaussian1D([0.1, 0.2], [0, 0], [5, 7], n_models=2)
        sequence = (np.array([-1, -2]), np.array([1, 2]))
        bounding_box = ModelBoundingBox.validate(model, sequence, **kwargs)
        assert "x" in bounding_box
        assert (bounding_box["x"].lower == np.array([-1, -2])).all()
        assert (bounding_box["x"].upper == np.array([1, 2])).all()

    def test_fix_inputs(self):
        bounding_box = ModelBoundingBox.validate(Gaussian2D(), ((-4, 4), (-1, 1)))

        # keep_ignored = False (default)
        new_bounding_box = bounding_box.fix_inputs(Gaussian1D(), {1: mk.MagicMock()})
        assert not (bounding_box == new_bounding_box)

        assert (new_bounding_box._model.parameters == Gaussian1D().parameters).all()
        assert "x" in new_bounding_box
        assert new_bounding_box["x"] == (-1, 1)
        assert "y" not in new_bounding_box
        assert len(new_bounding_box.intervals) == 1
        assert new_bounding_box.ignored == []

        # keep_ignored = True
        new_bounding_box = bounding_box.fix_inputs(
            Gaussian2D(), {1: mk.MagicMock()}, _keep_ignored=True
        )
        assert not (bounding_box == new_bounding_box)

        assert (new_bounding_box._model.parameters == Gaussian2D().parameters).all()
        assert "x" in new_bounding_box
        assert new_bounding_box["x"] == (-1, 1)
        assert "y" in new_bounding_box
        assert "y" in new_bounding_box.ignored_inputs
        assert len(new_bounding_box.intervals) == 1
        assert new_bounding_box.ignored == [1]

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
        assert (
            np.array(bounding_box.domain(0.25))
            == np.array([np.linspace(0, 2, 9), np.linspace(-1, 1, 9)])
        ).all()

        # test C order
        assert (
            np.array(bounding_box.domain(0.25, "C"))
            == np.array([np.linspace(0, 2, 9), np.linspace(-1, 1, 9)])
        ).all()

        # test Fortran order
        assert (
            np.array(bounding_box.domain(0.25, "F"))
            == np.array([np.linspace(-1, 1, 9), np.linspace(0, 2, 9)])
        ).all()

        # test error order
        MESSAGE = r"order must be either 'C' .* or 'F' .*, got: .*"
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box.domain(0.25, mk.MagicMock())

    def test__outside(self):
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

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
        # fmt: off
        assert (
            outside_index
            == [
                True, True, True, True,
                False, False, False, False, False,
                True, True, True, True
            ]
        ).all()
        # fmt: on
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        outside_index, all_out = bounding_box._outside(input_shape, inputs)
        assert (outside_index == [True for _ in range(13)]).all()
        assert all_out and isinstance(all_out, bool)

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
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

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
        intervals = {0: _Interval(-1, 1), 1: _Interval(0, 2)}
        model = Gaussian2D()
        bounding_box = ModelBoundingBox.validate(model, intervals)

        # Normal array input, all inside
        x = np.linspace(-1, 1, 13)
        y = np.linspace(0, 2, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(
            input_shape, inputs
        )
        assert (np.array(new_inputs) == np.array(inputs)).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [idx for idx in range(13)]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, some inside and some outside
        x = np.linspace(-2, 1, 13)
        y = np.linspace(0, 3, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(
            input_shape, inputs
        )
        assert (
            np.array(new_inputs)
            == np.array(
                [
                    [x[4], x[5], x[6], x[7], x[8]],
                    [y[4], y[5], y[6], y[7], y[8]],
                ]
            )
        ).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [4, 5, 6, 7, 8]).all()
        assert not all_out and isinstance(all_out, bool)

        # Normal array input, all outside
        x = np.linspace(2, 3, 13)
        y = np.linspace(-2, -1, 13)
        input_shape = x.shape
        inputs = (x, y)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(
            input_shape, inputs
        )
        assert new_inputs == ()
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

        # Scalar input inside bounding_box
        inputs = (0.5, 0.5)
        input_shape = (1,)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(
            input_shape, inputs
        )
        assert (np.array(new_inputs) == np.array([[0.5], [0.5]])).all()
        assert len(valid_index) == 1
        assert (valid_index[0] == [0]).all()
        assert not all_out and isinstance(all_out, bool)

        # Scalar input outside bounding_box
        inputs = (2, -1)
        input_shape = (1,)
        new_inputs, valid_index, all_out = bounding_box.prepare_inputs(
            input_shape, inputs
        )
        assert new_inputs == ()
        assert len(valid_index) == 1
        assert (valid_index[0] == []).all()
        assert all_out and isinstance(all_out, bool)

    def test_bounding_box_ignore(self):
        """Regression test for #13028"""

        bbox_x = ModelBoundingBox((9, 10), Polynomial2D(1), ignored=["x"])
        assert bbox_x.ignored_inputs == ["x"]

        bbox_y = ModelBoundingBox((11, 12), Polynomial2D(1), ignored=["y"])
        assert bbox_y.ignored_inputs == ["y"]


class Test_SelectorArgument:
    def test_create(self):
        index = mk.MagicMock()
        ignore = mk.MagicMock()
        argument = _SelectorArgument(index, ignore)

        assert isinstance(argument, _BaseSelectorArgument)
        assert argument.index == index
        assert argument.ignore == ignore
        assert argument == (index, ignore)

    def test_validate(self):
        model = Gaussian2D()

        # default integer
        assert _SelectorArgument.validate(model, 0) == (0, True)
        assert _SelectorArgument.validate(model, 1) == (1, True)

        # default string
        assert _SelectorArgument.validate(model, "x") == (0, True)
        assert _SelectorArgument.validate(model, "y") == (1, True)

        ignore = mk.MagicMock()
        # non-default integer
        assert _SelectorArgument.validate(model, 0, ignore) == (0, ignore)
        assert _SelectorArgument.validate(model, 1, ignore) == (1, ignore)

        # non-default string
        assert _SelectorArgument.validate(model, "x", ignore) == (0, ignore)
        assert _SelectorArgument.validate(model, "y", ignore) == (1, ignore)

        # Fail
        with pytest.raises(ValueError, match=r"'.*' is not one of the inputs: .*"):
            _SelectorArgument.validate(model, "z")
        with pytest.raises(
            ValueError, match=r"Key value: .* must be string or integer."
        ):
            _SelectorArgument.validate(model, mk.MagicMock())
        with pytest.raises(
            IndexError, match=r"Integer key: .* must be non-negative and < .*"
        ):
            _SelectorArgument.validate(model, 2)

    def test_get_selector(self):
        # single inputs
        inputs = [idx + 17 for idx in range(3)]
        for index in range(3):
            assert (
                _SelectorArgument(index, mk.MagicMock()).get_selector(*inputs)
                == inputs[index]
            )

        # numpy array of single inputs
        inputs = [np.array([idx + 11]) for idx in range(3)]
        for index in range(3):
            assert (
                _SelectorArgument(index, mk.MagicMock()).get_selector(*inputs)
                == inputs[index]
            )
        inputs = [np.asanyarray(idx + 13) for idx in range(3)]
        for index in range(3):
            assert (
                _SelectorArgument(index, mk.MagicMock()).get_selector(*inputs)
                == inputs[index]
            )

        # multi entry numpy array
        inputs = [np.array([idx + 27, idx - 31]) for idx in range(3)]
        for index in range(3):
            assert _SelectorArgument(index, mk.MagicMock()).get_selector(
                *inputs
            ) == tuple(inputs[index])

    def test_name(self):
        model = Gaussian2D()
        for index in range(model.n_inputs):
            assert (
                _SelectorArgument(index, mk.MagicMock()).name(model)
                == model.inputs[index]
            )

    def test_pretty_repr(self):
        model = Gaussian2D()

        assert (
            _SelectorArgument(0, False).pretty_repr(model)
            == "Argument(name='x', ignore=False)"
        )
        assert (
            _SelectorArgument(0, True).pretty_repr(model)
            == "Argument(name='x', ignore=True)"
        )
        assert (
            _SelectorArgument(1, False).pretty_repr(model)
            == "Argument(name='y', ignore=False)"
        )
        assert (
            _SelectorArgument(1, True).pretty_repr(model)
            == "Argument(name='y', ignore=True)"
        )

    def test_get_fixed_value(self):
        model = Gaussian2D()
        values = {0: 5, "y": 7}

        # Get index value
        assert _SelectorArgument(0, mk.MagicMock()).get_fixed_value(model, values) == 5

        # Get name value
        assert _SelectorArgument(1, mk.MagicMock()).get_fixed_value(model, values) == 7

        # Fail
        MESSAGE = r".* was not found in .*"
        with pytest.raises(RuntimeError, match=MESSAGE) as err:
            _SelectorArgument(1, True).get_fixed_value(model, {0: 5})

    def test_is_argument(self):
        model = Gaussian2D()
        argument = _SelectorArgument.validate(model, 0)

        # Is true
        assert argument.is_argument(model, 0) is True
        assert argument.is_argument(model, "x") is True

        # Is false
        assert argument.is_argument(model, 1) is False
        assert argument.is_argument(model, "y") is False

        # Fail
        with pytest.raises(ValueError, match=r"'.*' is not one of the inputs: .*"):
            argument.is_argument(model, "z")
        with pytest.raises(
            ValueError, match=r"Key value: .* must be string or integer"
        ):
            argument.is_argument(model, mk.MagicMock())
        with pytest.raises(
            IndexError, match=r"Integer key: .* must be non-negative and < .*"
        ):
            argument.is_argument(model, 2)

    def test_named_tuple(self):
        model = Gaussian2D()
        for index in range(model.n_inputs):
            ignore = mk.MagicMock()
            assert _SelectorArgument(index, ignore).named_tuple(model) == (
                model.inputs[index],
                ignore,
            )


class Test_SelectorArguments:
    def test_create(self):
        arguments = _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, False))
        )
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, False))
        assert arguments._kept_ignore == []

        kept_ignore = mk.MagicMock()
        arguments = _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, False)), kept_ignore
        )
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, False))
        assert arguments._kept_ignore == kept_ignore

    def test_pretty_repr(self):
        model = Gaussian2D()
        arguments = _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, False))
        )

        assert (
            arguments.pretty_repr(model) == "SelectorArguments(\n"
            "    Argument(name='x', ignore=True)\n"
            "    Argument(name='y', ignore=False)\n"
            ")"
        )

    def test_ignore(self):
        assert _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, True))
        ).ignore == [0, 1]
        assert _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, True)), [13, 4]
        ).ignore == [0, 1, 13, 4]
        assert _SelectorArguments(
            (_SelectorArgument(0, True), _SelectorArgument(1, False))
        ).ignore == [0]
        assert _SelectorArguments(
            (_SelectorArgument(0, False), _SelectorArgument(1, True))
        ).ignore == [1]
        assert (
            _SelectorArguments(
                (_SelectorArgument(0, False), _SelectorArgument(1, False))
            ).ignore
            == []
        )
        assert _SelectorArguments(
            (_SelectorArgument(0, False), _SelectorArgument(1, False)), [17, 14]
        ).ignore == [17, 14]

    def test_validate(self):
        # Integer key and passed ignore
        arguments = _SelectorArguments.validate(Gaussian2D(), ((0, True), (1, False)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, False))
        assert arguments.kept_ignore == []

        # Default ignore
        arguments = _SelectorArguments.validate(Gaussian2D(), ((0,), (1,)))
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, True))
        assert arguments.kept_ignore == []

        # String key and passed ignore
        arguments = _SelectorArguments.validate(
            Gaussian2D(), (("x", True), ("y", False))
        )
        assert isinstance(arguments, _SelectorArguments)
        assert arguments == ((0, True), (1, False))
        assert arguments.kept_ignore == []

        # Test kept_ignore option
        new_arguments = _SelectorArguments.validate(Gaussian2D(), arguments, [11, 5, 8])
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((0, True), (1, False))
        assert new_arguments.kept_ignore == [11, 5, 8]

        arguments._kept_ignore = [13, 17, 14]
        new_arguments = _SelectorArguments.validate(Gaussian2D(), arguments)
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((0, True), (1, False))
        assert new_arguments.kept_ignore == [13, 17, 14]

        # Invalid, bad argument
        with pytest.raises(ValueError, match=r"'.*' is not one of the inputs: .*"):
            _SelectorArguments.validate(Gaussian2D(), ((0, True), ("z", False)))
        with pytest.raises(
            ValueError, match=r"Key value: .* must be string or integer"
        ):
            _SelectorArguments.validate(
                Gaussian2D(), ((mk.MagicMock(), True), (1, False))
            )
        with pytest.raises(
            IndexError, match=r"Integer key: .* must be non-negative and < .*"
        ):
            _SelectorArguments.validate(Gaussian2D(), ((0, True), (2, False)))

        # Invalid, repeated argument
        with pytest.raises(ValueError, match=r"Input: 'x' has been repeated"):
            _SelectorArguments.validate(Gaussian2D(), ((0, True), (0, False)))

        # Invalid, no arguments
        with pytest.raises(
            ValueError, match=r"There must be at least one selector argument"
        ):
            _SelectorArguments.validate(Gaussian2D(), ())

    def test_get_selector(self):
        inputs = [idx + 19 for idx in range(4)]

        assert _SelectorArguments.validate(
            Gaussian2D(), ((0, True), (1, False))
        ).get_selector(*inputs) == tuple(inputs[:2])
        assert _SelectorArguments.validate(
            Gaussian2D(), ((1, True), (0, False))
        ).get_selector(*inputs) == tuple(inputs[:2][::-1])
        assert _SelectorArguments.validate(Gaussian2D(), ((1, False),)).get_selector(
            *inputs
        ) == (inputs[1],)
        assert _SelectorArguments.validate(Gaussian2D(), ((0, True),)).get_selector(
            *inputs
        ) == (inputs[0],)

    def test_is_selector(self):
        # Is Selector
        assert _SelectorArguments.validate(
            Gaussian2D(), ((0, True), (1, False))
        ).is_selector((0.5, 2.5))
        assert _SelectorArguments.validate(Gaussian2D(), ((0, True),)).is_selector(
            (0.5,)
        )

        # Is not selector
        assert not _SelectorArguments.validate(
            Gaussian2D(), ((0, True), (1, False))
        ).is_selector((0.5, 2.5, 3.5))
        assert not _SelectorArguments.validate(
            Gaussian2D(), ((0, True), (1, False))
        ).is_selector((0.5,))
        assert not _SelectorArguments.validate(
            Gaussian2D(), ((0, True), (1, False))
        ).is_selector(0.5)
        assert not _SelectorArguments.validate(Gaussian2D(), ((0, True),)).is_selector(
            (0.5, 2.5)
        )
        assert not _SelectorArguments.validate(Gaussian2D(), ((0, True),)).is_selector(
            2.5
        )

    def test_get_fixed_values(self):
        model = Gaussian2D()

        assert _SelectorArguments.validate(
            model, ((0, True), (1, False))
        ).get_fixed_values(model, {0: 11, 1: 7}) == (11, 7)
        assert _SelectorArguments.validate(
            model, ((0, True), (1, False))
        ).get_fixed_values(model, {0: 5, "y": 47}) == (5, 47)
        assert _SelectorArguments.validate(
            model, ((0, True), (1, False))
        ).get_fixed_values(model, {"x": 2, "y": 9}) == (2, 9)
        assert _SelectorArguments.validate(
            model, ((0, True), (1, False))
        ).get_fixed_values(model, {"x": 12, 1: 19}) == (12, 19)

    def test_is_argument(self):
        model = Gaussian2D()

        # Is true
        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.is_argument(model, 0) is True
        assert arguments.is_argument(model, "x") is True
        assert arguments.is_argument(model, 1) is True
        assert arguments.is_argument(model, "y") is True

        # Is true and false
        arguments = _SelectorArguments.validate(model, ((0, True),))
        assert arguments.is_argument(model, 0) is True
        assert arguments.is_argument(model, "x") is True
        assert arguments.is_argument(model, 1) is False
        assert arguments.is_argument(model, "y") is False

        arguments = _SelectorArguments.validate(model, ((1, False),))
        assert arguments.is_argument(model, 0) is False
        assert arguments.is_argument(model, "x") is False
        assert arguments.is_argument(model, 1) is True
        assert arguments.is_argument(model, "y") is True

    def test_selector_index(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.selector_index(model, 0) == 0
        assert arguments.selector_index(model, "x") == 0
        assert arguments.selector_index(model, 1) == 1
        assert arguments.selector_index(model, "y") == 1

        arguments = _SelectorArguments.validate(model, ((1, True), (0, False)))
        assert arguments.selector_index(model, 0) == 1
        assert arguments.selector_index(model, "x") == 1
        assert arguments.selector_index(model, 1) == 0
        assert arguments.selector_index(model, "y") == 0

        # Error
        arguments = _SelectorArguments.validate(model, ((0, True),))
        with pytest.raises(
            ValueError, match=r"y does not correspond to any selector argument"
        ):
            arguments.selector_index(model, "y")

    def test_add_ignore(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True),))
        assert arguments == ((0, True),)
        assert arguments._kept_ignore == []

        new_arguments0 = arguments.add_ignore(model, 1)
        assert new_arguments0 == arguments
        assert new_arguments0._kept_ignore == [1]
        assert arguments._kept_ignore == []

        assert arguments._kept_ignore == []
        new_arguments1 = new_arguments0.add_ignore(model, "y")
        assert new_arguments1 == arguments == new_arguments0
        assert new_arguments0._kept_ignore == [1]
        assert new_arguments1._kept_ignore == [1, 1]
        assert arguments._kept_ignore == []

        # Error
        with pytest.raises(
            ValueError, match=r"0: is a selector argument and cannot be ignored"
        ):
            arguments.add_ignore(model, 0)

    def test_reduce(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))

        new_arguments = arguments.reduce(model, 0)
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((1, False),)
        assert new_arguments._kept_ignore == [0]
        assert arguments._kept_ignore == []

        new_arguments = arguments.reduce(model, "x")
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((1, False),)
        assert new_arguments._kept_ignore == [0]
        assert arguments._kept_ignore == []

        new_arguments = arguments.reduce(model, 1)
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((0, True),)
        assert new_arguments._kept_ignore == [1]
        assert arguments._kept_ignore == []

        new_arguments = arguments.reduce(model, "y")
        assert isinstance(new_arguments, _SelectorArguments)
        assert new_arguments == ((0, True),)
        assert new_arguments._kept_ignore == [1]
        assert arguments._kept_ignore == []

    def test_named_tuple(self):
        model = Gaussian2D()

        arguments = _SelectorArguments.validate(model, ((0, True), (1, False)))
        assert arguments.named_tuple(model) == (("x", True), ("y", False))


class TestCompoundBoundingBox:
    def test_create(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        create_selector = mk.MagicMock()

        bounding_box = CompoundBoundingBox(
            bounding_boxes, model, selector_args, create_selector, order="F"
        )
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        for _selector, bbox in bounding_boxes.items():
            assert _selector in bounding_box._bounding_boxes
            assert bounding_box._bounding_boxes[_selector] == bbox
        for _selector, bbox in bounding_box._bounding_boxes.items():
            assert _selector in bounding_boxes
            assert bounding_boxes[_selector] == bbox
            assert isinstance(bbox, ModelBoundingBox)
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == "F"

    def test_copy(self):
        bounding_box = CompoundBoundingBox.validate(
            Gaussian2D(),
            {(1,): (-1.5, 1.3), (2,): (-2.7, 2.4)},
            ((0, True),),
            mk.MagicMock(),
        )
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
            assert argument.index == copy.selector_args[index].index
            assert id(argument.index) == id(copy.selector_args[index].index)

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
                assert id(interval) != id(
                    copy.bounding_boxes[selector].intervals[index]
                )

                # Same float values have will have same id
                assert (
                    interval.lower
                    == copy.bounding_boxes[selector].intervals[index].lower
                )
                assert id(interval.lower) == id(
                    copy.bounding_boxes[selector].intervals[index].lower
                )

                # Same float values have will have same id
                assert (
                    interval.upper
                    == copy.bounding_boxes[selector].intervals[index].upper
                )
                assert id(interval.upper) == id(
                    copy.bounding_boxes[selector].intervals[index].upper
                )
            assert len(bbox.intervals) == len(copy.bounding_boxes[selector].intervals)
            assert (
                bbox.intervals.keys() == copy.bounding_boxes[selector].intervals.keys()
            )
        assert len(bounding_box.bounding_boxes) == len(copy.bounding_boxes)
        assert bounding_box.bounding_boxes.keys() == copy.bounding_boxes.keys()

    def test___repr__(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        assert (
            bounding_box.__repr__() == "CompoundBoundingBox(\n"
            "    bounding_boxes={\n"
            "        (1,) = ModelBoundingBox(\n"
            "                intervals={\n"
            "                    y: Interval(lower=-1, upper=1)\n"
            "                }\n"
            "                ignored=['x']\n"
            "                model=Gaussian2D(inputs=('x', 'y'))\n"
            "                order='C'\n"
            "            )\n"
            "        (2,) = ModelBoundingBox(\n"
            "                intervals={\n"
            "                    y: Interval(lower=-2, upper=2)\n"
            "                }\n"
            "                ignored=['x']\n"
            "                model=Gaussian2D(inputs=('x', 'y'))\n"
            "                order='C'\n"
            "            )\n"
            "    }\n"
            "    selector_args = SelectorArguments(\n"
            "            Argument(name='x', ignore=True)\n"
            "        )\n"
            ")"
        )

    def test_bounding_boxes(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box.bounding_boxes == bounding_boxes

    def test_selector_args(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_box = CompoundBoundingBox({}, model, selector_args)

        # Get
        assert bounding_box._selector_args == selector_args
        assert bounding_box.selector_args == selector_args

        # Set
        selector_args = ((1, False),)
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

        # Singular
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
        bounding_box = CompoundBoundingBox({}, model, ((1, True),), order="F")
        assert len(bounding_box.bounding_boxes) == 0
        # Valid
        bounding_box[(15,)] = (-15, 15)
        assert len(bounding_box.bounding_boxes) == 1
        assert (15,) in bounding_box._bounding_boxes
        assert isinstance(bounding_box._bounding_boxes[(15,)], ModelBoundingBox)
        assert bounding_box._bounding_boxes[(15,)] == (-15, 15)
        assert bounding_box._bounding_boxes[(15,)].order == "F"
        # Invalid key
        assert (7, 13) not in bounding_box._bounding_boxes
        with pytest.raises(ValueError, match=".* is not a selector!"):
            bounding_box[(7, 13)] = (-7, 7)
        assert (7, 13) not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1
        # Invalid bounding box
        assert 13 not in bounding_box._bounding_boxes
        with pytest.raises(
            ValueError, match="An interval must be some sort of sequence of length 2"
        ):
            bounding_box[(13,)] = ((-13, 13), (-3, 3))
        assert 13 not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1

        # No ignored argument
        bounding_box = CompoundBoundingBox({}, model, ((1, False),), order="F")
        assert len(bounding_box.bounding_boxes) == 0
        # Valid
        bounding_box[(15,)] = ((-15, 15), (-6, 6))
        assert len(bounding_box.bounding_boxes) == 1
        assert (15,) in bounding_box._bounding_boxes
        assert isinstance(bounding_box._bounding_boxes[(15,)], ModelBoundingBox)
        assert bounding_box._bounding_boxes[(15,)] == ((-15, 15), (-6, 6))
        assert bounding_box._bounding_boxes[(15,)].order == "F"
        # Invalid key
        assert (14, 11) not in bounding_box._bounding_boxes
        with pytest.raises(ValueError, match=".* is not a selector!"):
            bounding_box[(14, 11)] = ((-7, 7), (-12, 12))
        assert (14, 11) not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1
        # Invalid bounding box
        assert 13 not in bounding_box._bounding_boxes
        with pytest.raises(
            ValueError, match="An interval must be some sort of sequence of length 2"
        ):
            bounding_box[(13,)] = (-13, 13)
        assert 13 not in bounding_box._bounding_boxes
        assert len(bounding_box.bounding_boxes) == 1

    def test__validate(self):
        model = Gaussian2D()
        selector_args = ((0, True),)

        # Tuple selector_args
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox({}, model, selector_args)
        bounding_box._validate(bounding_boxes)
        for _selector, bbox in bounding_boxes.items():
            assert _selector in bounding_box._bounding_boxes
            assert bounding_box._bounding_boxes[_selector] == bbox
        for _selector, bbox in bounding_box._bounding_boxes.items():
            assert _selector in bounding_boxes
            assert bounding_boxes[_selector] == bbox
            assert isinstance(bbox, ModelBoundingBox)
        assert bounding_box._bounding_boxes == bounding_boxes

    def test___eq__(self):
        bounding_box_1 = CompoundBoundingBox(
            {(1,): (-1, 1), (2,): (-2, 2)}, Gaussian2D(), ((0, True),)
        )
        bounding_box_2 = CompoundBoundingBox(
            {(1,): (-1, 1), (2,): (-2, 2)}, Gaussian2D(), ((0, True),)
        )

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
        bounding_box_2._selector_args = _SelectorArguments.validate(
            Gaussian2D(), ((0, False),)
        )
        assert not bounding_box_1 == bounding_box_2
        bounding_box_2._selector_args = _SelectorArguments.validate(
            Gaussian2D(), ((0, True),)
        )
        assert bounding_box_1 == bounding_box_2

        # Not equal create_selector
        bounding_box_2._create_selector = mk.MagicMock()
        assert not bounding_box_1 == bounding_box_2

    def test_validate(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        create_selector = mk.MagicMock()

        # Fail selector_args
        MESSAGE = r"Selector arguments must be provided .*"
        with pytest.raises(ValueError, match=MESSAGE):
            CompoundBoundingBox.validate(model, bounding_boxes)

        # Normal validate
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args, create_selector, order="F"
        )
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == "F"

        # Re-validate
        new_bounding_box = CompoundBoundingBox.validate(model, bounding_box)
        assert bounding_box == new_bounding_box
        assert new_bounding_box._order == "F"

        # Default order
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args, create_selector
        )
        assert (bounding_box._model.parameters == model.parameters).all()
        assert bounding_box._selector_args == selector_args
        assert bounding_box._bounding_boxes == bounding_boxes
        assert bounding_box._create_selector == create_selector
        assert bounding_box._order == "C"

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
        bounding_box = CompoundBoundingBox({}, model, ((1, False),), create_selector)

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
        with pytest.raises(
            ValueError, match="An interval must be some sort of sequence of length 2"
        ):
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
        with pytest.raises(
            RuntimeError, match="No bounding box is defined for selector: .*"
        ):
            bounding_box[(3,)]

        # Create a selector
        bounding_box._create_selector = mk.MagicMock()
        with mk.patch.object(
            CompoundBoundingBox, "_create_bounding_box", autospec=True
        ) as mkCreate:
            assert bounding_box[(3,)] == mkCreate.return_value
            assert mkCreate.call_args_list == [mk.call(bounding_box, (3,))]

    def test__select_bounding_box(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        inputs = [mk.MagicMock() for _ in range(3)]
        with mk.patch.object(
            _SelectorArguments, "get_selector", autospec=True
        ) as mkSelector:
            with mk.patch.object(
                CompoundBoundingBox, "__getitem__", autospec=True
            ) as mkGet:
                assert bounding_box._select_bounding_box(inputs) == mkGet.return_value
                assert mkGet.call_args_list == [
                    mk.call(bounding_box, mkSelector.return_value)
                ]
                assert mkSelector.call_args_list == [
                    mk.call(bounding_box.selector_args, *inputs)
                ]

    def test_prepare_inputs(self):
        model = Gaussian2D()
        selector_args = ((0, True),)
        bounding_boxes = {(1,): (-1, 1), (2,): (-2, 2)}
        bounding_box = CompoundBoundingBox(bounding_boxes, model, selector_args)

        input_shape = mk.MagicMock()
        with mk.patch.object(
            ModelBoundingBox, "prepare_inputs", autospec=True
        ) as mkPrepare:
            assert (
                bounding_box.prepare_inputs(input_shape, [1, 2, 3])
                == mkPrepare.return_value
            )
            assert mkPrepare.call_args_list == [
                mk.call(bounding_box[(1,)], input_shape, [1, 2, 3])
            ]
            mkPrepare.reset_mock()
            assert (
                bounding_box.prepare_inputs(input_shape, [2, 2, 3])
                == mkPrepare.return_value
            )
            assert mkPrepare.call_args_list == [
                mk.call(bounding_box[(2,)], input_shape, [2, 2, 3])
            ]
            mkPrepare.reset_mock()

    def test__matching_bounding_boxes(self):
        # Single selector index
        selector_args = ((0, False),)
        bounding_boxes = {
            (1,): ((-1, 1), (-2, 2)),
            (2,): ((-2, 2), (-3, 3)),
            (3,): ((-3, 3), (-4, 4)),
        }
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            matching = bounding_box._matching_bounding_boxes("x", value)
            assert isinstance(matching, dict)
            assert () in matching
            bbox = matching[()]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert "x" in bbox
            assert "x" in bbox.ignored_inputs
            assert "y" in bbox
            assert bbox["y"] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == [0]

        # Multiple selector index
        selector_args = ((0, False), (1, False))
        bounding_boxes = {
            (1, 3): ((-1, 1), (-2, 2)),
            (2, 2): ((-2, 2), (-3, 3)),
            (3, 1): ((-3, 3), (-4, 4)),
        }
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            matching = bounding_box._matching_bounding_boxes("x", value)
            assert isinstance(matching, dict)
            assert (4 - value,) in matching
            bbox = matching[(4 - value,)]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert "x" in bbox
            assert "x" in bbox.ignored_inputs
            assert "y" in bbox
            assert bbox["y"] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == [0]

            matching = bounding_box._matching_bounding_boxes("y", value)
            assert isinstance(matching, dict)
            assert (4 - value,) in matching
            bbox = matching[(4 - value,)]
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert "y" in bbox
            assert "y" in bbox.ignored_inputs
            assert "x" in bbox
            assert bbox["x"] == (-(5 - value), (5 - value))
            assert len(bbox.intervals) == 1
            assert bbox.ignored == [1]

        # Real fix input of slicing input
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ("x", "y", "slit_id")
        bounding_boxes = {
            (0,): ((-0.5, 1047.5), (-0.5, 2047.5)),
            (1,): ((-0.5, 3047.5), (-0.5, 4047.5)),
        }
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args=[("slit_id", True)], order="F"
        )

        matching = bounding_box._matching_bounding_boxes("slit_id", 0)
        assert isinstance(matching, dict)
        assert () in matching
        bbox = matching[()]
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored_inputs == ["slit_id"]
        assert bbox.named_intervals == {"x": (-0.5, 1047.5), "y": (-0.5, 2047.5)}
        assert bbox.order == "F"

        matching = bounding_box._matching_bounding_boxes("slit_id", 1)
        assert isinstance(matching, dict)
        assert () in matching
        bbox = matching[()]
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored_inputs == ["slit_id"]
        assert bbox.named_intervals == {"x": (-0.5, 3047.5), "y": (-0.5, 4047.5)}
        assert bbox.order == "F"

        # Errors
        MESSAGE = (
            r"Attempting to fix input .*, but there are no bounding boxes for argument"
            r" value .*"
        )
        with pytest.raises(ValueError, match=MESSAGE):
            bounding_box._matching_bounding_boxes("slit_id", 2)

    def test__fix_input_selector_arg(self):
        # Single selector index
        selector_args = ((0, False),)
        bounding_boxes = {
            (1,): ((-1, 1), (-2, 2)),
            (2,): ((-2, 2), (-3, 3)),
            (3,): ((-3, 3), (-4, 4)),
        }
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            bbox = bounding_box._fix_input_selector_arg("x", value)
            assert isinstance(bbox, ModelBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert "x" in bbox
            assert "x" in bbox.ignored_inputs
            assert "y" in bbox
            assert bbox["y"] == (-value, value)
            assert len(bbox.intervals) == 1
            assert bbox.ignored == [0]

        # Multiple selector index
        selector_args = ((0, False), (1, False))
        bounding_boxes = {
            (1, 3): ((-1, 1), (-2, 2)),
            (2, 2): ((-2, 2), (-3, 3)),
            (3, 1): ((-3, 3), (-4, 4)),
        }
        bounding_box = CompoundBoundingBox(bounding_boxes, Gaussian2D(), selector_args)

        for value in [1, 2, 3]:
            bbox = bounding_box._fix_input_selector_arg("x", value)
            assert isinstance(bbox, CompoundBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert bbox.selector_args == ((1, False),)
            assert (4 - value,) in bbox
            bbox_selector = bbox[(4 - value,)]
            assert isinstance(bbox_selector, ModelBoundingBox)
            assert (bbox_selector._model.parameters == Gaussian2D().parameters).all()
            assert "x" in bbox_selector
            assert "x" in bbox_selector.ignored_inputs
            assert "y" in bbox_selector
            assert bbox_selector["y"] == (-value, value)
            assert len(bbox_selector.intervals) == 1
            assert bbox_selector.ignored == [0]

            bbox = bounding_box._fix_input_selector_arg("y", value)
            assert isinstance(bbox, CompoundBoundingBox)
            assert (bbox._model.parameters == Gaussian2D().parameters).all()
            assert bbox.selector_args == ((0, False),)
            assert (4 - value,) in bbox
            bbox_selector = bbox[(4 - value,)]
            assert isinstance(bbox_selector, ModelBoundingBox)
            assert (bbox_selector._model.parameters == Gaussian2D().parameters).all()
            assert "y" in bbox_selector
            assert "y" in bbox_selector.ignored_inputs
            assert "x" in bbox_selector
            assert bbox_selector["x"] == (-(5 - value), (5 - value))
            assert len(bbox_selector.intervals) == 1
            assert bbox_selector.ignored == [1]

        # Real fix input of slicing input
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ("x", "y", "slit_id")
        bounding_boxes = {
            (0,): ((-0.5, 1047.5), (-0.5, 2047.5)),
            (1,): ((-0.5, 3047.5), (-0.5, 4047.5)),
        }
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args=[("slit_id", True)], order="F"
        )

        bbox = bounding_box._fix_input_selector_arg("slit_id", 0)
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored_inputs == ["slit_id"]
        assert bbox.named_intervals == {"x": (-0.5, 1047.5), "y": (-0.5, 2047.5)}
        assert bbox.order == "F"

        bbox = bounding_box._fix_input_selector_arg("slit_id", 1)
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.ignored_inputs == ["slit_id"]
        assert bbox.named_intervals == {"x": (-0.5, 3047.5), "y": (-0.5, 4047.5)}
        assert bbox.order == "F"

    def test__fix_input_bbox_arg(self):
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ("x", "y", "slit_id")
        bounding_boxes = {
            (0,): ((-0.5, 1047.5), (-0.5, 2047.5)),
            (1,): ((-0.5, 3047.5), (-0.5, 4047.5)),
        }
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args=[("slit_id", True)], order="F"
        )

        bbox = bounding_box._fix_input_bbox_arg("x", 5)
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == ((2, True),)
        assert bbox.selector_args._kept_ignore == [0]
        assert bbox._bounding_boxes[(0,)] == (-0.5, 2047.5)
        assert bbox._bounding_boxes[(1,)] == (-0.5, 4047.5)
        assert len(bbox._bounding_boxes) == 2

        bbox = bounding_box._fix_input_bbox_arg("y", 5)
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == ((2, True),)
        assert bbox.selector_args._kept_ignore == [1]
        assert bbox._bounding_boxes[(0,)] == (-0.5, 1047.5)
        assert bbox._bounding_boxes[(1,)] == (-0.5, 3047.5)
        assert len(bbox._bounding_boxes) == 2

    def test_fix_inputs(self):
        model = Shift(1) & Scale(2) & Identity(1)
        model.inputs = ("x", "y", "slit_id")
        bounding_boxes = {
            (0,): ((-0.5, 1047.5), (-0.5, 2047.5)),
            (1,): ((-0.5, 3047.5), (-0.5, 4047.5)),
        }
        bounding_box = CompoundBoundingBox.validate(
            model, bounding_boxes, selector_args=[("slit_id", True)], order="F"
        )
        model.bounding_box = bounding_box

        # Fix selector argument
        new_model = fix_inputs(model, {"slit_id": 0})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.named_intervals == {"x": (-0.5, 1047.5), "y": (-0.5, 2047.5)}
        assert bbox.order == "F"

        # Fix a bounding_box field
        new_model = fix_inputs(model, {"x": 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == ((1, True),)
        assert bbox.selector_args._kept_ignore == []
        assert bbox._bounding_boxes[(0,)] == (-0.5, 2047.5)
        assert bbox._bounding_boxes[(0,)].order == "F"
        assert bbox._bounding_boxes[(1,)] == (-0.5, 4047.5)
        assert bbox._bounding_boxes[(1,)].order == "F"
        assert len(bbox._bounding_boxes) == 2
        new_model = fix_inputs(model, {"y": 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert (bbox._model.parameters == model.parameters).all()
        assert bbox.selector_args == ((1, True),)
        assert bbox.selector_args._kept_ignore == []
        assert bbox._bounding_boxes[(0,)] == (-0.5, 1047.5)
        assert bbox._bounding_boxes[(0,)].order == "F"
        assert bbox._bounding_boxes[(1,)] == (-0.5, 3047.5)
        assert bbox._bounding_boxes[(1,)].order == "F"
        assert len(bbox._bounding_boxes) == 2

        # Fix selector argument and a bounding_box field
        new_model = fix_inputs(model, {"slit_id": 0, "x": 5})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.named_intervals == {"y": (-0.5, 2047.5)}
        assert bbox.order == "F"
        new_model = fix_inputs(model, {"y": 5, "slit_id": 1})
        bbox = new_model.bounding_box
        assert isinstance(bbox, ModelBoundingBox)
        assert (bbox._model.parameters == new_model.parameters).all()
        assert bbox.ignored_inputs == []
        assert bbox.named_intervals == {"x": (-0.5, 3047.5)}
        assert bbox.order == "F"

        # Fix two bounding_box fields
        new_model = fix_inputs(model, {"x": 5, "y": 7})
        bbox = new_model.bounding_box
        assert isinstance(bbox, CompoundBoundingBox)
        assert bbox.selector_args == ((0, True),)
        assert bbox.selector_args._kept_ignore == []
        assert bbox._bounding_boxes[(0,)] == (-np.inf, np.inf)
        assert bbox._bounding_boxes[(0,)].order == "F"
        assert bbox._bounding_boxes[(1,)] == (-np.inf, np.inf)
        assert bbox._bounding_boxes[(1,)].order == "F"
        assert len(bbox._bounding_boxes) == 2

    def test_complex_compound_bounding_box(self):
        model = Identity(4)
        bounding_boxes = {
            (2.5, 1.3): ((-1, 1), (-3, 3)),
            (2.5, 2.71): ((-3, 3), (-1, 1)),
        }
        selector_args = (("x0", True), ("x1", True))

        bbox = CompoundBoundingBox.validate(model, bounding_boxes, selector_args)
        assert bbox[(2.5, 1.3)] == ModelBoundingBox(
            ((-1, 1), (-3, 3)), model, ignored=["x0", "x1"]
        )
        assert bbox[(2.5, 2.71)] == ModelBoundingBox(
            ((-3, 3), (-1, 1)), model, ignored=["x0", "x1"]
        )
