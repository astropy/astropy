# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module is to contain an improved bounding box.
"""

from __future__ import annotations

import abc
import copy
import warnings
from typing import Any, Callable, NamedTuple

import numpy as np

from astropy.units import Quantity
from astropy.utils import isiterable

__all__ = ["ModelBoundingBox", "CompoundBoundingBox"]


class _BaseInterval(NamedTuple):
    lower: float
    upper: float


class _Interval(_BaseInterval):
    """
    A single input's bounding box interval.

    Parameters
    ----------
    lower : float
        The lower bound of the interval

    upper : float
        The upper bound of the interval

    Methods
    -------
    validate :
        Constructs a valid interval

    outside :
        Determine which parts of an input array are outside the interval.

    domain :
        Constructs a discretization of the points inside the interval.
    """

    def __repr__(self):
        return f"Interval(lower={self.lower}, upper={self.upper})"

    def copy(self):
        return copy.deepcopy(self)

    @staticmethod
    def _validate_shape(interval):
        """Validate the shape of an interval representation."""
        MESSAGE = """An interval must be some sort of sequence of length 2"""

        try:
            shape = np.shape(interval)
        except TypeError:
            try:
                # np.shape does not work with lists of Quantities
                if len(interval) == 1:
                    interval = interval[0]
                shape = np.shape([b.to_value() for b in interval])
            except (ValueError, TypeError, AttributeError):
                raise ValueError(MESSAGE)

        valid_shape = shape in ((2,), (1, 2), (2, 0))
        if not valid_shape:
            valid_shape = (
                len(shape) > 0
                and shape[0] == 2
                and all(isinstance(b, np.ndarray) for b in interval)
            )

        if not isiterable(interval) or not valid_shape:
            raise ValueError(MESSAGE)

    @classmethod
    def _validate_bounds(cls, lower, upper):
        """Validate the bounds are reasonable and construct an interval from them."""
        if (np.asanyarray(lower) > np.asanyarray(upper)).all():
            warnings.warn(
                f"Invalid interval: upper bound {upper} "
                f"is strictly less than lower bound {lower}.",
                RuntimeWarning,
            )

        return cls(lower, upper)

    @classmethod
    def validate(cls, interval):
        """
        Construct and validate an interval.

        Parameters
        ----------
        interval : iterable
            A representation of the interval.

        Returns
        -------
        A validated interval.
        """
        cls._validate_shape(interval)

        if len(interval) == 1:
            interval = tuple(interval[0])
        else:
            interval = tuple(interval)

        return cls._validate_bounds(interval[0], interval[1])

    def outside(self, _input: np.ndarray):
        """
        Parameters
        ----------
        _input : np.ndarray
            The evaluation input in the form of an array.

        Returns
        -------
        Boolean array indicating which parts of _input are outside the interval:
            True  -> position outside interval
            False -> position inside  interval
        """
        return np.logical_or(_input < self.lower, _input > self.upper)

    def domain(self, resolution):
        return np.arange(self.lower, self.upper + resolution, resolution)


# The interval where all ignored inputs can be found.
_ignored_interval = _Interval.validate((-np.inf, np.inf))


def get_index(model, key) -> int:
    """
    Get the input index corresponding to the given key.
        Can pass in either:
            the string name of the input or
            the input index itself.
    """
    if isinstance(key, str):
        if key in model.inputs:
            index = model.inputs.index(key)
        else:
            raise ValueError(f"'{key}' is not one of the inputs: {model.inputs}.")
    elif np.issubdtype(type(key), np.integer):
        if 0 <= key < len(model.inputs):
            index = key
        else:
            raise IndexError(
                f"Integer key: {key} must be non-negative and < {len(model.inputs)}."
            )
    else:
        raise ValueError(f"Key value: {key} must be string or integer.")

    return index


def get_name(model, index: int):
    """Get the input name corresponding to the input index."""
    return model.inputs[index]


class _BoundingDomain(abc.ABC):
    """
    Base class for ModelBoundingBox and CompoundBoundingBox.
        This is where all the `~astropy.modeling.core.Model` evaluation
        code for evaluating with a bounding box is because it is common
        to both types of bounding box.

    Parameters
    ----------
    model : `~astropy.modeling.Model`
        The Model this bounding domain is for.

    prepare_inputs :
        Generates the necessary input information so that model can
        be evaluated only for input points entirely inside bounding_box.
        This needs to be implemented by a subclass. Note that most of
        the implementation is in ModelBoundingBox.

    prepare_outputs :
        Fills the output values in for any input points outside the
        bounding_box.

    evaluate :
        Performs a complete model evaluation while enforcing the bounds
        on the inputs and returns a complete output.
    """

    def __init__(self, model, ignored: list[int] | None = None, order: str = "C"):
        self._model = model
        self._ignored = self._validate_ignored(ignored)
        self._order = self._get_order(order)

    @property
    def model(self):
        return self._model

    @property
    def order(self) -> str:
        return self._order

    @property
    def ignored(self) -> list[int]:
        return self._ignored

    def _get_order(self, order: str | None = None) -> str:
        """
        Get if bounding_box is C/python ordered or Fortran/mathematically
        ordered.
        """
        if order is None:
            order = self._order

        if order not in ("C", "F"):
            raise ValueError(
                "order must be either 'C' (C/python order) or "
                f"'F' (Fortran/mathematical order), got: {order}."
            )

        return order

    def _get_index(self, key) -> int:
        """
        Get the input index corresponding to the given key.
            Can pass in either:
                the string name of the input or
                the input index itself.
        """
        return get_index(self._model, key)

    def _get_name(self, index: int):
        """Get the input name corresponding to the input index."""
        return get_name(self._model, index)

    @property
    def ignored_inputs(self) -> list[str]:
        return [self._get_name(index) for index in self._ignored]

    def _validate_ignored(self, ignored: list) -> list[int]:
        if ignored is None:
            return []
        else:
            return [self._get_index(key) for key in ignored]

    def __call__(self, *args, **kwargs):
        raise NotImplementedError(
            "This bounding box is fixed by the model and does not have "
            "adjustable parameters."
        )

    @abc.abstractmethod
    def fix_inputs(self, model, fixed_inputs: dict):
        """
        Fix the bounding_box for a `fix_inputs` compound model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The new model for which this will be a bounding_box
        fixed_inputs : dict
            Dictionary of inputs which have been fixed by this bounding box.
        """
        raise NotImplementedError("This should be implemented by a child class.")

    @abc.abstractmethod
    def prepare_inputs(self, input_shape, inputs) -> tuple[Any, Any, Any]:
        """
        Get prepare the inputs with respect to the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs

        Returns
        -------
        valid_inputs : list
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : array_like
            array of all indices inside the bounding box
        all_out: bool
            if all of the inputs are outside the bounding_box
        """
        raise NotImplementedError("This has not been implemented for BoundingDomain.")

    @staticmethod
    def _base_output(input_shape, fill_value):
        """
        Create a baseline output, assuming that the entire input is outside
        the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box

        Returns
        -------
        An array of the correct shape containing all fill_value
        """
        return np.zeros(input_shape) + fill_value

    def _all_out_output(self, input_shape, fill_value):
        """
        Create output if all inputs are outside the domain.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box

        Returns
        -------
        A full set of outputs for case that all inputs are outside domain.
        """
        return [
            self._base_output(input_shape, fill_value)
            for _ in range(self._model.n_outputs)
        ], None

    def _modify_output(self, valid_output, valid_index, input_shape, fill_value):
        """
        For a single output fill in all the parts corresponding to inputs
        outside the bounding box.

        Parameters
        ----------
        valid_output : numpy array
            The output from the model corresponding to inputs inside the
            bounding box
        valid_index : numpy array
            array of all indices of inputs inside the bounding box
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box

        Returns
        -------
        An output array with all the indices corresponding to inputs
        outside the bounding box filled in by fill_value
        """
        output = self._base_output(input_shape, fill_value)
        if not output.shape:
            output = np.array(valid_output)
        else:
            output[valid_index] = valid_output

        if np.isscalar(valid_output):
            output = output.item(0)

        return output

    def _prepare_outputs(self, valid_outputs, valid_index, input_shape, fill_value):
        """
        Fill in all the outputs of the model corresponding to inputs
        outside the bounding_box.

        Parameters
        ----------
        valid_outputs : list of numpy array
            The list of outputs from the model corresponding to inputs
            inside the bounding box
        valid_index : numpy array
            array of all indices of inputs inside the bounding box
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box

        Returns
        -------
        List of filled in output arrays.
        """
        outputs = []
        for valid_output in valid_outputs:
            outputs.append(
                self._modify_output(valid_output, valid_index, input_shape, fill_value)
            )

        return outputs

    def prepare_outputs(self, valid_outputs, valid_index, input_shape, fill_value):
        """
        Fill in all the outputs of the model corresponding to inputs
        outside the bounding_box, adjusting any single output model so that
        its output becomes a list of containing that output.

        Parameters
        ----------
        valid_outputs : list
            The list of outputs from the model corresponding to inputs
            inside the bounding box
        valid_index : array_like
            array of all indices of inputs inside the bounding box
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box
        """
        if self._model.n_outputs == 1:
            valid_outputs = [valid_outputs]

        return self._prepare_outputs(
            valid_outputs, valid_index, input_shape, fill_value
        )

    @staticmethod
    def _get_valid_outputs_unit(valid_outputs, with_units: bool):
        """
        Get the unit for outputs if one is required.

        Parameters
        ----------
        valid_outputs : list of numpy array
            The list of outputs from the model corresponding to inputs
            inside the bounding box
        with_units : bool
            whether or not a unit is required
        """
        if with_units:
            return getattr(valid_outputs, "unit", None)

    def _evaluate_model(
        self,
        evaluate: Callable,
        valid_inputs,
        valid_index,
        input_shape,
        fill_value,
        with_units: bool,
    ):
        """
        Evaluate the model using the given evaluate routine.

        Parameters
        ----------
        evaluate : Callable
            callable which takes in the valid inputs to evaluate model
        valid_inputs : list of numpy arrays
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : numpy array
            array of all indices inside the bounding box
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box
        with_units : bool
            whether or not a unit is required

        Returns
        -------
        outputs :
            list containing filled in output values
        valid_outputs_unit :
            the unit that will be attached to the outputs
        """
        valid_outputs = evaluate(valid_inputs)
        valid_outputs_unit = self._get_valid_outputs_unit(valid_outputs, with_units)

        return (
            self.prepare_outputs(valid_outputs, valid_index, input_shape, fill_value),
            valid_outputs_unit,
        )

    def _evaluate(
        self, evaluate: Callable, inputs, input_shape, fill_value, with_units: bool
    ):
        """Evaluate model with steps: prepare_inputs -> evaluate -> prepare_outputs.

        Parameters
        ----------
        evaluate : Callable
            callable which takes in the valid inputs to evaluate model
        valid_inputs : list of numpy arrays
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : numpy array
            array of all indices inside the bounding box
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box
        with_units : bool
            whether or not a unit is required

        Returns
        -------
        outputs :
            list containing filled in output values
        valid_outputs_unit :
            the unit that will be attached to the outputs
        """
        valid_inputs, valid_index, all_out = self.prepare_inputs(input_shape, inputs)

        if all_out:
            return self._all_out_output(input_shape, fill_value)
        else:
            return self._evaluate_model(
                evaluate, valid_inputs, valid_index, input_shape, fill_value, with_units
            )

    @staticmethod
    def _set_outputs_unit(outputs, valid_outputs_unit):
        """
        Set the units on the outputs
            prepare_inputs -> evaluate -> prepare_outputs -> set output units.

        Parameters
        ----------
        outputs :
            list containing filled in output values
        valid_outputs_unit :
            the unit that will be attached to the outputs

        Returns
        -------
        List containing filled in output values and units
        """
        if valid_outputs_unit is not None:
            return Quantity(outputs, valid_outputs_unit, copy=False, subok=True)

        return outputs

    def evaluate(self, evaluate: Callable, inputs, fill_value):
        """
        Perform full model evaluation steps:
            prepare_inputs -> evaluate -> prepare_outputs -> set output units.

        Parameters
        ----------
        evaluate : callable
            callable which takes in the valid inputs to evaluate model
        valid_inputs : list
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : array_like
            array of all indices inside the bounding box
        fill_value : float
            The value which will be assigned to inputs which are outside
            the bounding box
        """
        input_shape = self._model.input_shape(inputs)

        # NOTE: CompoundModel does not currently support units during
        #   evaluation for bounding_box so this feature is turned off
        #   for CompoundModel(s).
        outputs, valid_outputs_unit = self._evaluate(
            evaluate, inputs, input_shape, fill_value, self._model.bbox_with_units
        )
        return tuple(self._set_outputs_unit(outputs, valid_outputs_unit))


class ModelBoundingBox(_BoundingDomain):
    """
    A model's bounding box.

    Parameters
    ----------
    intervals : dict
        A dictionary containing all the intervals for each model input
            keys   -> input index
            values -> interval for that index

    model : `~astropy.modeling.Model`
        The Model this bounding_box is for.

    ignored : list
        A list containing all the inputs (index) which will not be
        checked for whether or not their elements are in/out of an interval.

    order : optional, str
        The ordering that is assumed for the tuple representation of this
        bounding_box. Options: 'C': C/Python order, e.g. z, y, x.
        (default), 'F': Fortran/mathematical notation order, e.g. x, y, z.
    """

    def __init__(
        self,
        intervals: dict[int, _Interval],
        model,
        ignored: list[int] | None = None,
        order: str = "C",
    ):
        super().__init__(model, ignored, order)

        self._intervals = {}
        if intervals != () and intervals != {}:
            self._validate(intervals, order=order)

    def copy(self, ignored=None):
        intervals = {
            index: interval.copy() for index, interval in self._intervals.items()
        }

        if ignored is None:
            ignored = self._ignored.copy()

        return ModelBoundingBox(
            intervals, self._model, ignored=ignored, order=self._order
        )

    @property
    def intervals(self) -> dict[int, _Interval]:
        """Return bounding_box labeled using input positions."""
        return self._intervals

    @property
    def named_intervals(self) -> dict[str, _Interval]:
        """Return bounding_box labeled using input names."""
        return {self._get_name(index): bbox for index, bbox in self._intervals.items()}

    def __repr__(self):
        parts = ["ModelBoundingBox(", "    intervals={"]

        for name, interval in self.named_intervals.items():
            parts.append(f"        {name}: {interval}")

        parts.append("    }")
        if len(self._ignored) > 0:
            parts.append(f"    ignored={self.ignored_inputs}")

        parts.append(
            f"    model={self._model.__class__.__name__}(inputs={self._model.inputs})"
        )
        parts.append(f"    order='{self._order}'")
        parts.append(")")

        return "\n".join(parts)

    def __len__(self):
        return len(self._intervals)

    def __contains__(self, key):
        try:
            return self._get_index(key) in self._intervals or self._ignored
        except (IndexError, ValueError):
            return False

    def has_interval(self, key):
        return self._get_index(key) in self._intervals

    def __getitem__(self, key):
        """Get bounding_box entries by either input name or input index."""
        index = self._get_index(key)
        if index in self._ignored:
            return _ignored_interval
        else:
            return self._intervals[self._get_index(key)]

    def bounding_box(self, order: str | None = None):
        """
        Return the old tuple of tuples representation of the bounding_box
            order='C' corresponds to the old bounding_box ordering
            order='F' corresponds to the gwcs bounding_box ordering.
        """
        if len(self._intervals) == 1:
            return tuple(next(iter(self._intervals.values())))
        else:
            order = self._get_order(order)
            inputs = self._model.inputs
            if order == "C":
                inputs = inputs[::-1]

            bbox = tuple(tuple(self[input_name]) for input_name in inputs)
            if len(bbox) == 1:
                bbox = bbox[0]

            return bbox

    def __eq__(self, value):
        """Note equality can be either with old representation or new one."""
        if isinstance(value, tuple):
            return self.bounding_box() == value
        elif isinstance(value, ModelBoundingBox):
            return (self.intervals == value.intervals) and (
                self.ignored == value.ignored
            )
        else:
            return False

    def __setitem__(self, key, value):
        """Validate and store interval under key (input index or input name)."""
        index = self._get_index(key)
        if index in self._ignored:
            self._ignored.remove(index)

        self._intervals[index] = _Interval.validate(value)

    def __delitem__(self, key):
        """Delete stored interval."""
        index = self._get_index(key)
        if index in self._ignored:
            raise RuntimeError(f"Cannot delete ignored input: {key}!")
        del self._intervals[index]
        self._ignored.append(index)

    def _validate_dict(self, bounding_box: dict):
        """Validate passing dictionary of intervals and setting them."""
        for key, value in bounding_box.items():
            self[key] = value

    @property
    def _available_input_index(self):
        model_input_index = [self._get_index(_input) for _input in self._model.inputs]

        return [_input for _input in model_input_index if _input not in self._ignored]

    def _validate_sequence(self, bounding_box, order: str | None = None):
        """
        Validate passing tuple of tuples representation (or related) and setting them.
        """
        order = self._get_order(order)
        if order == "C":
            # If bounding_box is C/python ordered, it needs to be reversed
            # to be in Fortran/mathematical/input order.
            bounding_box = bounding_box[::-1]

        for index, value in enumerate(bounding_box):
            self[self._available_input_index[index]] = value

    @property
    def _n_inputs(self) -> int:
        n_inputs = self._model.n_inputs - len(self._ignored)
        if n_inputs > 0:
            return n_inputs
        else:
            return 0

    def _validate_iterable(self, bounding_box, order: str | None = None):
        """Validate and set any iterable representation."""
        if len(bounding_box) != self._n_inputs:
            raise ValueError(
                f"Found {len(bounding_box)} intervals, "
                f"but must have exactly {self._n_inputs}."
            )

        if isinstance(bounding_box, dict):
            self._validate_dict(bounding_box)
        else:
            self._validate_sequence(bounding_box, order)

    def _validate(self, bounding_box, order: str | None = None):
        """Validate and set any representation."""
        if self._n_inputs == 1 and not isinstance(bounding_box, dict):
            self[self._available_input_index[0]] = bounding_box
        else:
            self._validate_iterable(bounding_box, order)

    @classmethod
    def validate(
        cls,
        model,
        bounding_box,
        ignored: list | None = None,
        order: str = "C",
        _preserve_ignore: bool = False,
        **kwargs,
    ):
        """
        Construct a valid bounding box for a model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The model for which this will be a bounding_box
        bounding_box : dict, tuple
            A possible representation of the bounding box
        order : optional, str
            The order that a tuple representation will be assumed to be
                Default: 'C'
        """
        if isinstance(bounding_box, ModelBoundingBox):
            order = bounding_box.order
            if _preserve_ignore:
                ignored = bounding_box.ignored
            bounding_box = bounding_box.named_intervals

        new = cls({}, model, ignored=ignored, order=order)
        new._validate(bounding_box)

        return new

    def fix_inputs(self, model, fixed_inputs: dict, _keep_ignored=False):
        """
        Fix the bounding_box for a `fix_inputs` compound model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The new model for which this will be a bounding_box
        fixed_inputs : dict
            Dictionary of inputs which have been fixed by this bounding box.
        keep_ignored : bool
            Keep the ignored inputs of the bounding box (internal argument only)
        """
        new = self.copy()

        for _input in fixed_inputs.keys():
            del new[_input]

        if _keep_ignored:
            ignored = new.ignored
        else:
            ignored = None

        return ModelBoundingBox.validate(
            model, new.named_intervals, ignored=ignored, order=new._order
        )

    @property
    def dimension(self):
        return len(self)

    def domain(self, resolution, order: str | None = None):
        inputs = self._model.inputs
        order = self._get_order(order)
        if order == "C":
            inputs = inputs[::-1]

        return [self[input_name].domain(resolution) for input_name in inputs]

    def _outside(self, input_shape, inputs):
        """
        Get all the input positions which are outside the bounding_box,
        so that the corresponding outputs can be filled with the fill
        value (default NaN).

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs

        Returns
        -------
        outside_index : bool-numpy array
            True  -> position outside bounding_box
            False -> position inside  bounding_box
        all_out : bool
            if all of the inputs are outside the bounding_box
        """
        all_out = False

        outside_index = np.zeros(input_shape, dtype=bool)
        for index, _input in enumerate(inputs):
            _input = np.asanyarray(_input)

            outside = np.broadcast_to(self[index].outside(_input), input_shape)
            outside_index[outside] = True

            if outside_index.all():
                all_out = True
                break

        return outside_index, all_out

    def _valid_index(self, input_shape, inputs):
        """
        Get the indices of all the inputs inside the bounding_box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs

        Returns
        -------
        valid_index : numpy array
            array of all indices inside the bounding box
        all_out : bool
            if all of the inputs are outside the bounding_box
        """
        outside_index, all_out = self._outside(input_shape, inputs)

        valid_index = np.atleast_1d(np.logical_not(outside_index)).nonzero()
        if len(valid_index[0]) == 0:
            all_out = True

        return valid_index, all_out

    def prepare_inputs(self, input_shape, inputs) -> tuple[Any, Any, Any]:
        """
        Get prepare the inputs with respect to the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs

        Returns
        -------
        valid_inputs : list
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : array_like
            array of all indices inside the bounding box
        all_out: bool
            if all of the inputs are outside the bounding_box
        """
        valid_index, all_out = self._valid_index(input_shape, inputs)

        valid_inputs = []
        if not all_out:
            for _input in inputs:
                if input_shape:
                    valid_input = np.broadcast_to(np.atleast_1d(_input), input_shape)[
                        valid_index
                    ]
                    if np.isscalar(_input):
                        valid_input = valid_input.item(0)
                    valid_inputs.append(valid_input)
                else:
                    valid_inputs.append(_input)

        return tuple(valid_inputs), valid_index, all_out


class _BaseSelectorArgument(NamedTuple):
    index: int
    ignore: bool


class _SelectorArgument(_BaseSelectorArgument):
    """
    Contains a single CompoundBoundingBox slicing input.

    Parameters
    ----------
    index : int
        The index of the input in the input list

    ignore : bool
        Whether or not this input will be ignored by the bounding box.

    Methods
    -------
    validate :
        Returns a valid SelectorArgument for a given model.

    get_selector :
        Returns the value of the input for use in finding the correct
        bounding_box.

    get_fixed_value :
        Gets the slicing value from a fix_inputs set of values.
    """

    def __new__(cls, index, ignore):
        self = super().__new__(cls, index, ignore)

        return self

    @classmethod
    def validate(cls, model, argument, ignored: bool = True):
        """
        Construct a valid selector argument for a CompoundBoundingBox.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The model for which this will be an argument for.
        argument : int or str
            A representation of which evaluation input to use
        ignored : optional, bool
            Whether or not to ignore this argument in the ModelBoundingBox.

        Returns
        -------
        Validated selector_argument
        """
        return cls(get_index(model, argument), ignored)

    def get_selector(self, *inputs):
        """
        Get the selector value corresponding to this argument.

        Parameters
        ----------
        *inputs :
            All the processed model evaluation inputs.
        """
        _selector = inputs[self.index]
        if isiterable(_selector):
            if len(_selector) == 1:
                return _selector[0]
            else:
                return tuple(_selector)
        return _selector

    def name(self, model) -> str:
        """
        Get the name of the input described by this selector argument.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return get_name(model, self.index)

    def pretty_repr(self, model):
        """
        Get a pretty-print representation of this object.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return f"Argument(name='{self.name(model)}', ignore={self.ignore})"

    def get_fixed_value(self, model, values: dict):
        """
        Gets the value fixed input corresponding to this argument.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.

        values : dict
            Dictionary of fixed inputs.
        """
        if self.index in values:
            return values[self.index]
        else:
            if self.name(model) in values:
                return values[self.name(model)]
            else:
                raise RuntimeError(
                    f"{self.pretty_repr(model)} was not found in {values}"
                )

    def is_argument(self, model, argument) -> bool:
        """
        Determine if passed argument is described by this selector argument.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.

        argument : int or str
            A representation of which evaluation input is being used
        """
        return self.index == get_index(model, argument)

    def named_tuple(self, model):
        """
        Get a tuple representation of this argument using the input
        name from the model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return (self.name(model), self.ignore)


class _SelectorArguments(tuple):
    """
    Contains the CompoundBoundingBox slicing description.

    Parameters
    ----------
    input_ :
        The SelectorArgument values

    Methods
    -------
    validate :
        Returns a valid SelectorArguments for its model.

    get_selector :
        Returns the selector a set of inputs corresponds to.

    is_selector :
        Determines if a selector is correctly formatted for this CompoundBoundingBox.

    get_fixed_value :
        Gets the selector from a fix_inputs set of values.
    """

    _kept_ignore = None

    def __new__(cls, input_: tuple[_SelectorArgument], kept_ignore: list | None = None):
        self = super().__new__(cls, input_)

        if kept_ignore is None:
            self._kept_ignore = []
        else:
            self._kept_ignore = kept_ignore

        return self

    def pretty_repr(self, model):
        """
        Get a pretty-print representation of this object.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.
        """
        parts = ["SelectorArguments("]
        for argument in self:
            parts.append(f"    {argument.pretty_repr(model)}")
        parts.append(")")

        return "\n".join(parts)

    @property
    def ignore(self):
        """Get the list of ignored inputs."""
        ignore = [argument.index for argument in self if argument.ignore]
        ignore.extend(self._kept_ignore)

        return ignore

    @property
    def kept_ignore(self):
        """The arguments to persist in ignoring."""
        return self._kept_ignore

    @classmethod
    def validate(cls, model, arguments, kept_ignore: list | None = None):
        """
        Construct a valid Selector description for a CompoundBoundingBox.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        arguments :
            The individual argument information

        kept_ignore :
            Arguments to persist as ignored
        """
        inputs = []
        for argument in arguments:
            _input = _SelectorArgument.validate(model, *argument)
            if _input.index in [this.index for this in inputs]:
                raise ValueError(
                    f"Input: '{get_name(model, _input.index)}' has been repeated."
                )
            inputs.append(_input)

        if len(inputs) == 0:
            raise ValueError("There must be at least one selector argument.")

        if isinstance(arguments, _SelectorArguments):
            if kept_ignore is None:
                kept_ignore = []

            kept_ignore.extend(arguments.kept_ignore)

        return cls(tuple(inputs), kept_ignore)

    def get_selector(self, *inputs):
        """
        Get the selector corresponding to these inputs.

        Parameters
        ----------
        *inputs :
            All the processed model evaluation inputs.
        """
        return tuple(argument.get_selector(*inputs) for argument in self)

    def is_selector(self, _selector):
        """
        Determine if this is a reasonable selector.

        Parameters
        ----------
        _selector : tuple
            The selector to check
        """
        return isinstance(_selector, tuple) and len(_selector) == len(self)

    def get_fixed_values(self, model, values: dict):
        """
        Gets the value fixed input corresponding to this argument.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        values : dict
            Dictionary of fixed inputs.
        """
        return tuple(argument.get_fixed_value(model, values) for argument in self)

    def is_argument(self, model, argument) -> bool:
        """
        Determine if passed argument is one of the selector arguments.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which evaluation input is being used
        """
        return any(selector_arg.is_argument(model, argument) for selector_arg in self)

    def selector_index(self, model, argument):
        """
        Get the index of the argument passed in the selector tuples.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which argument is being used
        """
        for index, selector_arg in enumerate(self):
            if selector_arg.is_argument(model, argument):
                return index
        raise ValueError(f"{argument} does not correspond to any selector argument.")

    def reduce(self, model, argument):
        """
        Reduce the selector arguments by the argument given.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which argument is being used
        """
        arguments = list(self)
        kept_ignore = [arguments.pop(self.selector_index(model, argument)).index]
        kept_ignore.extend(self._kept_ignore)

        return _SelectorArguments.validate(model, tuple(arguments), kept_ignore)

    def add_ignore(self, model, argument):
        """
        Add argument to the kept_ignore list.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which argument is being used
        """
        if self.is_argument(model, argument):
            raise ValueError(
                f"{argument}: is a selector argument and cannot be ignored."
            )

        kept_ignore = [get_index(model, argument)]

        return _SelectorArguments.validate(model, self, kept_ignore)

    def named_tuple(self, model):
        """
        Get a tuple of selector argument tuples using input names.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.
        """
        return tuple(selector_arg.named_tuple(model) for selector_arg in self)


class CompoundBoundingBox(_BoundingDomain):
    """
    A model's compound bounding box.

    Parameters
    ----------
    bounding_boxes : dict
        A dictionary containing all the ModelBoundingBoxes that are possible
            keys   -> _selector (extracted from model inputs)
            values -> ModelBoundingBox

    model : `~astropy.modeling.Model`
        The Model this compound bounding_box is for.

    selector_args : _SelectorArguments
        A description of how to extract the selectors from model inputs.

    create_selector : optional
        A method which takes in the selector and the model to return a
        valid bounding corresponding to that selector. This can be used
        to construct new bounding_boxes for previously undefined selectors.
        These new boxes are then stored for future lookups.

    order : optional, str
        The ordering that is assumed for the tuple representation of the
        bounding_boxes.
    """

    def __init__(
        self,
        bounding_boxes: dict[Any, ModelBoundingBox],
        model,
        selector_args: _SelectorArguments,
        create_selector: Callable | None = None,
        ignored: list[int] | None = None,
        order: str = "C",
    ):
        super().__init__(model, ignored, order)

        self._create_selector = create_selector
        self._selector_args = _SelectorArguments.validate(model, selector_args)

        self._bounding_boxes = {}
        self._validate(bounding_boxes)

    def copy(self):
        bounding_boxes = {
            selector: bbox.copy(self.selector_args.ignore)
            for selector, bbox in self._bounding_boxes.items()
        }

        return CompoundBoundingBox(
            bounding_boxes,
            self._model,
            selector_args=self._selector_args,
            create_selector=copy.deepcopy(self._create_selector),
            order=self._order,
        )

    def __repr__(self):
        parts = ["CompoundBoundingBox(", "    bounding_boxes={"]
        # bounding_boxes
        for _selector, bbox in self._bounding_boxes.items():
            bbox_repr = bbox.__repr__().split("\n")
            parts.append(f"        {_selector} = {bbox_repr.pop(0)}")
            for part in bbox_repr:
                parts.append(f"            {part}")
        parts.append("    }")

        # selector_args
        selector_args_repr = self.selector_args.pretty_repr(self._model).split("\n")
        parts.append(f"    selector_args = {selector_args_repr.pop(0)}")
        for part in selector_args_repr:
            parts.append(f"        {part}")
        parts.append(")")

        return "\n".join(parts)

    @property
    def bounding_boxes(self) -> dict[Any, ModelBoundingBox]:
        return self._bounding_boxes

    @property
    def selector_args(self) -> _SelectorArguments:
        return self._selector_args

    @selector_args.setter
    def selector_args(self, value):
        self._selector_args = _SelectorArguments.validate(self._model, value)

        warnings.warn(
            "Overriding selector_args may cause problems you should re-validate "
            "the compound bounding box before use!",
            RuntimeWarning,
        )

    @property
    def named_selector_tuple(self) -> tuple:
        return self._selector_args.named_tuple(self._model)

    @property
    def create_selector(self):
        return self._create_selector

    @staticmethod
    def _get_selector_key(key):
        if isiterable(key):
            return tuple(key)
        else:
            return (key,)

    def __setitem__(self, key, value):
        _selector = self._get_selector_key(key)
        if not self.selector_args.is_selector(_selector):
            raise ValueError(f"{_selector} is not a selector!")

        ignored = self.selector_args.ignore + self.ignored
        self._bounding_boxes[_selector] = ModelBoundingBox.validate(
            self._model, value, ignored, order=self._order
        )

    def _validate(self, bounding_boxes: dict):
        for _selector, bounding_box in bounding_boxes.items():
            self[_selector] = bounding_box

    def __eq__(self, value):
        if isinstance(value, CompoundBoundingBox):
            return (
                self.bounding_boxes == value.bounding_boxes
                and self.selector_args == value.selector_args
                and self.create_selector == value.create_selector
            )
        else:
            return False

    @classmethod
    def validate(
        cls,
        model,
        bounding_box: dict,
        selector_args=None,
        create_selector=None,
        ignored: list | None = None,
        order: str = "C",
        _preserve_ignore: bool = False,
        **kwarg,
    ):
        """
        Construct a valid compound bounding box for a model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The model for which this will be a bounding_box
        bounding_box : dict
            Dictionary of possible bounding_box representations
        selector_args : optional
            Description of the selector arguments
        create_selector : optional, callable
            Method for generating new selectors
        order : optional, str
            The order that a tuple representation will be assumed to be
                Default: 'C'
        """
        if isinstance(bounding_box, CompoundBoundingBox):
            if selector_args is None:
                selector_args = bounding_box.selector_args
            if create_selector is None:
                create_selector = bounding_box.create_selector
            order = bounding_box.order
            if _preserve_ignore:
                ignored = bounding_box.ignored
            bounding_box = bounding_box.bounding_boxes

        if selector_args is None:
            raise ValueError(
                "Selector arguments must be provided "
                "(can be passed as part of bounding_box argument)"
            )

        return cls(
            bounding_box,
            model,
            selector_args,
            create_selector=create_selector,
            ignored=ignored,
            order=order,
        )

    def __contains__(self, key):
        return key in self._bounding_boxes

    def _create_bounding_box(self, _selector):
        self[_selector] = self._create_selector(_selector, model=self._model)

        return self[_selector]

    def __getitem__(self, key):
        _selector = self._get_selector_key(key)
        if _selector in self:
            return self._bounding_boxes[_selector]
        elif self._create_selector is not None:
            return self._create_bounding_box(_selector)
        else:
            raise RuntimeError(f"No bounding box is defined for selector: {_selector}.")

    def _select_bounding_box(self, inputs) -> ModelBoundingBox:
        _selector = self.selector_args.get_selector(*inputs)

        return self[_selector]

    def prepare_inputs(self, input_shape, inputs) -> tuple[Any, Any, Any]:
        """
        Get prepare the inputs with respect to the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs

        Returns
        -------
        valid_inputs : list
            The inputs reduced to just those inputs which are all inside
            their respective bounding box intervals
        valid_index : array_like
            array of all indices inside the bounding box
        all_out: bool
            if all of the inputs are outside the bounding_box
        """
        bounding_box = self._select_bounding_box(inputs)
        return bounding_box.prepare_inputs(input_shape, inputs)

    def _matching_bounding_boxes(self, argument, value) -> dict[Any, ModelBoundingBox]:
        selector_index = self.selector_args.selector_index(self._model, argument)
        matching = {}
        for selector_key, bbox in self._bounding_boxes.items():
            if selector_key[selector_index] == value:
                new_selector_key = list(selector_key)
                new_selector_key.pop(selector_index)

                if bbox.has_interval(argument):
                    new_bbox = bbox.fix_inputs(
                        self._model, {argument: value}, _keep_ignored=True
                    )
                else:
                    new_bbox = bbox.copy()

                matching[tuple(new_selector_key)] = new_bbox

        if len(matching) == 0:
            raise ValueError(
                f"Attempting to fix input {argument}, but there are no "
                f"bounding boxes for argument value {value}."
            )

        return matching

    def _fix_input_selector_arg(self, argument, value):
        matching_bounding_boxes = self._matching_bounding_boxes(argument, value)

        if len(self.selector_args) == 1:
            return matching_bounding_boxes[()]
        else:
            return CompoundBoundingBox(
                matching_bounding_boxes,
                self._model,
                self.selector_args.reduce(self._model, argument),
            )

    def _fix_input_bbox_arg(self, argument, value):
        bounding_boxes = {}
        for selector_key, bbox in self._bounding_boxes.items():
            bounding_boxes[selector_key] = bbox.fix_inputs(
                self._model, {argument: value}, _keep_ignored=True
            )

        return CompoundBoundingBox(
            bounding_boxes,
            self._model,
            self.selector_args.add_ignore(self._model, argument),
        )

    def fix_inputs(self, model, fixed_inputs: dict):
        """
        Fix the bounding_box for a `fix_inputs` compound model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The new model for which this will be a bounding_box
        fixed_inputs : dict
            Dictionary of inputs which have been fixed by this bounding box.
        """
        fixed_input_keys = list(fixed_inputs.keys())
        argument = fixed_input_keys.pop()
        value = fixed_inputs[argument]

        if self.selector_args.is_argument(self._model, argument):
            bbox = self._fix_input_selector_arg(argument, value)
        else:
            bbox = self._fix_input_bbox_arg(argument, value)

        if len(fixed_input_keys) > 0:
            new_fixed_inputs = fixed_inputs.copy()
            del new_fixed_inputs[argument]

            bbox = bbox.fix_inputs(model, new_fixed_inputs)

        if isinstance(bbox, CompoundBoundingBox):
            selector_args = bbox.named_selector_tuple
            bbox_dict = bbox
        elif isinstance(bbox, ModelBoundingBox):
            selector_args = None
            bbox_dict = bbox.named_intervals

        return bbox.__class__.validate(
            model, bbox_dict, order=bbox.order, selector_args=selector_args
        )
