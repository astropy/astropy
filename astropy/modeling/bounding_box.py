# Licensed under a 3-clause BSD style license - see LICENSE.rst


"""
This module is to contain an improved bounding box
"""

import abc
import copy
import warnings
from collections import namedtuple
from typing import Any, Callable, Dict, List, Tuple

import numpy as np

from astropy.units import Quantity
from astropy.utils import isiterable

__all__ = ['ModelBoundingBox', 'CompoundBoundingBox']


_BaseInterval = namedtuple('_BaseInterval', "lower upper")


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
        Contructs a valid interval

    outside :
        Determine which parts of an input array are outside the interval.

    domain :
        Contructs a discretization of the points inside the interval.
    """

    def __repr__(self):
        return f"Interval(lower={self.lower}, upper={self.upper})"

    def copy(self):
        return copy.deepcopy(self)

    @staticmethod
    def _validate_shape(interval):
        """Validate the shape of an interval representation"""
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
            valid_shape = (len(shape) > 0) and (shape[0] == 2) and \
                all(isinstance(b, np.ndarray) for b in interval)

        if not isiterable(interval) or not valid_shape:
            raise ValueError(MESSAGE)

    @classmethod
    def _validate_bounds(cls, lower, upper):
        """Validate the bounds are reasonable and construct an interval from them."""
        if (np.asanyarray(lower) > np.asanyarray(upper)).all():
            warnings.warn(f"Invalid interval: upper bound {upper} "
                          f"is strictly less than lower bound {lower}.", RuntimeWarning)

        return cls(lower, upper)

    @classmethod
    def validate(cls, interval):
        """
        Construct and validate an interval

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
            raise IndexError(f"Integer key: {key} must be non-negative and < {len(model.inputs)}.")
    else:
        raise ValueError(f"Key value: {key} must be string or integer.")

    return index


def get_name(model, key):
    """Get the input name corresponding to the input index"""

    return model.inputs[get_index(model, key)]


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

    def __init__(self, model = None, ignored: List[str] = None, order: str = 'C'):
        """Note this init should be called last in sub-classes"""

        self._order = self._get_order(order)

        if ignored is None:
            self._ignored = []
        else:
            self._ignored = ignored

        self.model = model

    @property
    def model(self):
        if self._model is None:
            raise RuntimeError("Method requires a model to function, please attach to a model")
        else:
            return self._model

    @model.setter
    def model(self, model):
        self.verify(model)

    @property
    def _has_model(self) -> bool:
        return self._model is not None

    @property
    def order(self) -> str:
        return self._order

    @property
    def ignored(self) -> List[str]:
        return self._ignored

    def _get_order(self, order: str = None) -> str:
        """
        Get if bounding_box is C/python ordered or Fortran/mathematically
        ordered
        """
        if order is None:
            order = self._order

        if order not in ('C', 'F'):
            raise ValueError("order must be either 'C' (C/python order) or "
                             f"'F' (Fortran/mathematical order), got: {order}.")

        return order

    def _get_index(self, key) -> int:
        """
        Get the input index corresponding to the given key.
            Can pass in either:
                the string name of the input or
                the input index itself.
        """

        return get_index(self.model, key)

    def _get_name(self, key):
        """Get the input name corresponding to the input index"""
        return get_name(self.model, key)

    @property
    def ignored_inputs(self) -> List[int]:
        return [self._get_index(name) for name in self._ignored]

    def verify(self, model, _external_ignored: List[str] = None):
        """
        Fully integrate the domain with the model and verify its functionality.

        Parameters
        ----------
        model :
            The astropy model to integrate functionality with
        _external_ignored : list
            Convenience for passing externally ignored model inputs.
            For compound bounding_box and fix_inputs support.
        """

        self._model = model
        if self._has_model:
            self._ignored = [self._get_name(key) for key in self._ignored]

        self._verify(_external_ignored)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError(
            "This bounding box is fixed by the model and does not have "
            "adjustable parameters.")

    @abc.abstractclassmethod
    def _verify(self, _external_ignored: List[str] = None):
        """
        Subclass verification method

        _external_ignored : list
            Convenience for passing externally ignored model inputs.
        """

        pass # pragma: no cover

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

        pass # pragma: no cover

    @abc.abstractmethod
    def prepare_inputs(self, input_shape, inputs, ignored: List[str] = []) -> Tuple[Any, Any, Any]:
        """
        Get prepare the inputs with respect to the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs
        ignored : List
            List of inputs to ignore by name

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

        pass # pragma: no cover

    @staticmethod
    def _base_output(input_shape, fill_value):
        """
        Create a baseline output, assuming that the entire input is outside
        the bounding box

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
        Create output if all inputs are outside the domain

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

        return [self._base_output(input_shape, fill_value)
                for _ in range(self._model.n_outputs)], None

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
            outputs.append(self._modify_output(valid_output, valid_index, input_shape, fill_value))

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

        return self._prepare_outputs(valid_outputs, valid_index, input_shape, fill_value)

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
            return getattr(valid_outputs, 'unit', None)

    def _evaluate_model(self, evaluate: Callable, valid_inputs, valid_index,
                        input_shape, fill_value, with_units: bool):
        """
        Evaluate the model using the given evaluate routine

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

        return self.prepare_outputs(valid_outputs, valid_index,
                                    input_shape, fill_value), valid_outputs_unit

    def _evaluate(self, evaluate: Callable, inputs, input_shape,
                  fill_value, with_units: bool):
        """
        Perform model evaluation steps:
            prepare_inputs -> evaluate -> prepare_outputs

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
            return self._evaluate_model(evaluate, valid_inputs, valid_index,
                                        input_shape, fill_value, with_units)

    @staticmethod
    def _set_outputs_unit(outputs, valid_outputs_unit):
        """
        Set the units on the outputs
            prepare_inputs -> evaluate -> prepare_outputs -> set output units

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
            return Quantity(outputs, valid_outputs_unit, copy=False)

        return outputs

    def evaluate(self, evaluate: Callable, inputs, fill_value):
        """
        Perform full model evaluation steps:
            prepare_inputs -> evaluate -> prepare_outputs -> set output units

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
        outputs, valid_outputs_unit = self._evaluate(evaluate, inputs, input_shape,
                                                     fill_value, self._model.bbox_with_units)
        return tuple(self._set_outputs_unit(outputs, valid_outputs_unit))


class ModelBoundingBox(_BoundingDomain):
    """
    A model's bounding box

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

    def __init__(self, intervals: Dict[str, _Interval], model = None,
                 ignored: List[str] = None, order: str = 'C'):

        # HACK, prevents a huge number of errors in legacy bounding box related tests.
        #       should be comparing to () but there is a bug in __eq__.
        if isinstance(intervals, tuple) and len(intervals) == 0:
            intervals = {}

        self._intervals = intervals
        super().__init__(model, ignored, order)

    def _pop_intervals(self):
        intervals = self._intervals
        self._intervals = {}

        return intervals

    def _verify_intervals_dict(self, _external_ignored: List[str] = None):
        if _external_ignored is None:
            _external_ignored = []

        intervals = self._pop_intervals()

        for key, value in intervals.items():
            name = self._get_name(key)
            if name in _external_ignored:
                raise ValueError(f"Interval: {name} is being externally ignored!")

            self[name] = value

        if len(self.intervals) != 0 and len(self.intervals) != len(names := self._interval_names(_external_ignored)):
            raise ValueError(f"Given {len(self.intervals)}, need {len(names)} intervals!")

    def _interval_names(self, _external_ignored: List[str] = None):
        if _external_ignored is None:
            _external_ignored = []

        return [name for name in self.model.inputs
                if name not in self._ignored
                and name not in _external_ignored]

    def _verify_intervals_sequence(self, _external_ignored: List[str] = None):
        names = self._interval_names(_external_ignored)
        if len(names) == 0 and len(self.intervals) > 0:
            raise ValueError("All intervals have been ignored!")

        intervals = self._pop_intervals()
        if not isiterable(intervals) or len(intervals) <= 1:
            raise ValueError(f"The intervals: {intervals}, do not contain enough information to construct a bounding_box!")

        if len(names) == 1: # Handle the 1D tuple case.
            self[names[0]] = intervals
        elif len(names) == len(intervals):
            if self._order == 'C':
                intervals = intervals[::-1]

            for value in intervals:
                self[names.pop(0)] = value
        else:
            raise ValueError(f"Given {len(intervals)}, need {len(names)} intervals!")

    def _verify_intervals(self, _external_ignored: List[str] = None):
        ignored = self._ignored.copy()

        if isinstance(self._intervals, dict):
            self._verify_intervals_dict(_external_ignored)
        else:
            self._verify_intervals_sequence(_external_ignored)

        if ignored != self._ignored or any(name in self._intervals for name in ignored):
            raise ValueError("At least one interval is being ignored")

    def _verify(self, _external_ignored: List[str] = None):
        if self._has_model:
            self._verify_intervals(_external_ignored)

    def copy(self, _external_ignored=None):
        intervals = {name: interval.copy()
                     for name, interval in self._intervals.items()}

        return ModelBoundingBox.validate(self._model, intervals, ignored=self.ignored.copy(),
                                         order=self.order, _external_ignored=_external_ignored)

    @property
    def intervals(self) -> Dict[str, _Interval]:
        """Return bounding_box labeled using input positions"""
        return self._intervals

    @property
    def indexed_intervals(self) -> Dict[int, _Interval]:
        """Return bounding_box labeled using input names"""
        return {self._get_index(name): bbox for name, bbox in self._intervals.items()}

    def __repr__(self):
        parts = [
            'ModelBoundingBox(',
            '    intervals={'
        ]

        for name, interval in self.intervals.items():
            parts.append(f"        {name}: {interval}")

        parts.append('    }')
        if len(self._ignored) > 0:
            parts.append(f"    ignored={self.ignored}")

        parts.append(f'    model={self._model.__class__.__name__}(inputs={self._model.inputs})')
        parts.append(f"    order='{self._order}'")
        parts.append(')')

        return '\n'.join(parts)

    def __len__(self):
        return len(self._intervals)

    def __contains__(self, key):
        try:
            return self._get_name(key) in self._intervals or self._ignored
        except (IndexError, ValueError):
            return False

    def has_interval(self, key):
        return self._get_name(key) in self._intervals

    def __getitem__(self, key):
        """Get bounding_box entries by either input name or input index"""
        name = self._get_name(key)
        if name in self._ignored:
            return _ignored_interval
        else:
            return self._intervals[self._get_name(key)]

    def bounding_box(self, order: str = None):
        """
        Return the old tuple of tuples representation of the bounding_box
            order='C' corresponds to the old bounding_box ordering
            order='F' corresponds to the gwcs bounding_box ordering.
        """
        if len(self._intervals) == 1:
            return tuple(list(self._intervals.values())[0])
        else:
            order = self._get_order(order)
            inputs = self._model.inputs
            if order == 'C':
                inputs = inputs[::-1]

            bbox = tuple([tuple(self[input_name]) for input_name in inputs])
            if len(bbox) == 1:
                bbox = bbox[0]

            return bbox

    def __eq__(self, value):
        """Note equality can be either with old representation or new one."""
        if isinstance(value, tuple):
            if len(self._intervals) == 0:
                return tuple(_ignored_interval) == value
            else:
                return self.bounding_box() == value
        elif isinstance(value, ModelBoundingBox):
            return (self.intervals == value.intervals) and (self.ignored == value.ignored)
        else:
            return False

    def __setitem__(self, key, value):
        """Validate and store interval under key (input index or input name)."""

        if self._has_model:
            name = self._get_name(key)
            if name in self._ignored:
                self._ignored.remove(name)
        else:
            name = key

        self._intervals[name] = _Interval.validate(value)

    def __delitem__(self, key):
        """Delete stored interval"""
        if self._has_model:
            name = self._get_name(key)
            if name in self._ignored:
                raise RuntimeError(f"Cannot delete ignored input: {key}!")
            self._ignored.append(name)
        else:
            name = key

        del self._intervals[name]

    @classmethod
    def validate(cls, model, bounding_box,
                 ignored: List[str] = None, order: str = 'C', _external_ignored: List[str] = None, **kwargs):
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
            ignored = bounding_box.ignored
            bounding_box = bounding_box.intervals

        new = cls(bounding_box, ignored=ignored, order=order)
        new.verify(model, _external_ignored)

        return new

    def fix_inputs(self, model, fixed_inputs: dict,
                   _external_ignored: List[str] = None,
                   _compound_ignore=False):
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
        if _external_ignored is None:
            _external_ignored = []

        new = self.copy(_external_ignored)

        for _input in fixed_inputs.keys():
            del new[_input]
            if _compound_ignore:
                new._ignored.remove(_input)
                _external_ignored.append(_input)

        return ModelBoundingBox.validate(model, new.intervals,
                                         order=new._order,
                                         _external_ignored=_external_ignored)

    @property
    def dimension(self):
        return len(self)

    def domain(self, resolution, order: str = None):
        inputs = self._model.inputs
        order = self._get_order(order)
        if order == 'C':
            inputs = inputs[::-1]

        return [self[input_name].domain(resolution) for input_name in inputs]

    def _get_interval(self, index: int, ignored: List[str]) -> _Interval:
        """
        Get the interval for the given input index accounting for externally
        ignored intervals

        Parameters
        ----------
        index : int
            The index of the input
        ignored : list
            List of inputs to ignore by name

        Returns
        -------
        an interval
        """
        if self._get_name(index) in ignored:
            return _ignored_interval
        else:
            return self[index]

    def _outside(self,  input_shape, inputs, ignored: List[str]):
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
        ignored : List
            List of inputs to ignore by name

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

            outside = np.broadcast_to(self._get_interval(index, ignored).outside(_input),
                                      input_shape)
            outside_index[outside] = True

            if outside_index.all():
                all_out = True
                break

        return outside_index, all_out

    def _valid_index(self, input_shape, inputs, ignored: List[str]):
        """
        Get the indices of all the inputs inside the bounding_box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs
        ignored : List
            List of inputs to ignore by name

        Returns
        -------
        valid_index : numpy array
            array of all indices inside the bounding box
        all_out : bool
            if all of the inputs are outside the bounding_box
        """
        outside_index, all_out = self._outside(input_shape, inputs, ignored)

        valid_index = np.atleast_1d(np.logical_not(outside_index)).nonzero()
        if len(valid_index[0]) == 0:
            all_out = True

        return valid_index, all_out

    def prepare_inputs(self, input_shape, inputs, ignored: List[str] = None) -> Tuple[Any, Any, Any]:
        """
        Get prepare the inputs with respect to the bounding box.

        Parameters
        ----------
        input_shape : tuple
            The shape that all inputs have be reshaped/broadcasted into
        inputs : list
            List of all the model inputs
        ignored : List
            List of inputs to ignore by name

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
        if ignored is None:
            ignored = []

        valid_index, all_out = self._valid_index(input_shape, inputs, ignored)

        valid_inputs = []
        if not all_out:
            for _input in inputs:
                if input_shape:
                    valid_input = np.broadcast_to(np.atleast_1d(_input), input_shape)[valid_index]
                    if np.isscalar(_input):
                        valid_input = valid_input.item(0)
                    valid_inputs.append(valid_input)
                else:
                    valid_inputs.append(_input)

        return tuple(valid_inputs), valid_index, all_out


_BaseSelectorArgument = namedtuple('_BaseSelectorArgument', "name ignore")


class _SelectorArgument(_BaseSelectorArgument):
    """
    Contains a single CompoundBoundingBox slicing input.

    Parameters
    ----------
    name : str
        The name of the input in the input list

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
        return cls(get_name(model, argument), ignored)

    def index(self, model) -> int:
        """
        Get the name of the input described by this selector argument

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return get_index(model, self.name)

    def get_selector(self, model, *inputs):
        """
        Get the selector value corresponding to this argument

        Parameters
        ----------
        *inputs :
            All the processed model evaluation inputs.
        """
        _selector = inputs[self.index(model)]
        if isiterable(_selector):
            if len(_selector) == 1:
                return _selector[0]
            else:
                return tuple(_selector)
        return _selector

    def __repr__(self):
        """
        Get a pretty-print representation of this object

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return f"Argument(name='{self.name}', ignore={self.ignore})"

    def get_fixed_value(self, model, values: dict):
        """
        Gets the value fixed input corresponding to this argument

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.

        values : dict
            Dictionary of fixed inputs.
        """
        if self.name in values:
            return values[self.name]
        else:
            if self.index(model) in values:
                return values[self.index(model)]
            else:
                raise RuntimeError(f"{self} was not found in {values}")

    def is_argument(self, model, argument) -> bool:
        """
        Determine if passed argument is described by this selector argument

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.

        argument : int or str
            A representation of which evaluation input is being used
        """

        return self.name == get_name(model, argument)

    def index_tuple(self, model):
        """
        Get a tuple representation of this argument using the input
        name from the model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model this selector argument is for.
        """
        return (self.index(model), self.ignore)


class _SelectorArguments(tuple):
    """
    Contains the CompoundBoundingBox slicing description

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

    def __new__(cls, input_: Tuple[_SelectorArgument]):
        self = super().__new__(cls, input_)

        return self

    def __repr__(self):
        """
        Get a pretty-print representation of this object

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.
        """
        parts = ['SelectorArguments(']
        for argument in self:
            parts.append(
                f"    {argument}"
            )
        parts.append(')')

        return '\n'.join(parts)

    @property
    def ignored(self) -> List[str]:
        """Get the list of ignored inputs"""
        return [argument.name for argument in self if argument.ignore]

    @classmethod
    def validate(cls, model, arguments):
        """
        Construct a valid Selector description for a CompoundBoundingBox.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        arguments :
            The individual argument informations
        """
        inputs = []
        for argument in arguments:
            _input = _SelectorArgument.validate(model, *argument)
            if _input.name in [this.name for this in inputs]:
                raise ValueError(f"Input: '{_input.name}' has been repeated.")
            inputs.append(_input)

        if len(inputs) == 0:
            raise ValueError("There must be at least one selector argument.")

        return cls(tuple(inputs))

    def get_selector(self, model, *inputs):
        """
        Get the selector corresponding to these inputs

        Parameters
        ----------
        *inputs :
            All the processed model evaluation inputs.
        """
        return tuple([argument.get_selector(model, *inputs) for argument in self])

    def is_selector(self, _selector):
        """
        Determine if this is a reasonable selector

        Parameters
        ----------
        _selector : tuple
            The selector to check
        """
        return isinstance(_selector, tuple) and len(_selector) == len(self)

    def get_fixed_values(self, model, values: dict):
        """
        Gets the value fixed input corresponding to this argument

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        values : dict
            Dictionary of fixed inputs.
        """
        return tuple([argument.get_fixed_value(model, values) for argument in self])

    def is_argument(self, model, argument) -> bool:
        """
        Determine if passed argument is one of the selector arguments

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which evaluation input is being used
        """

        for selector_arg in self:
            if selector_arg.is_argument(model, argument):
                return True
        else:
            return False

    def selector_index(self, model, argument):
        """
        Get the index of the argument passed in the selector tuples

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
        else:
            raise ValueError(f"{argument} does not correspond to any selector argument.")

    def reduce(self, model, argument):
        """
        Reduce the selector arguments by the argument given

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.

        argument : int or str
            A representation of which argument is being used
        """

        arguments = list(self)
        arguments.pop(self.selector_index(model, argument))

        return _SelectorArguments.validate(model, tuple(arguments))

    def index_tuple(self, model):
        """
        Get a tuple of selector argument tuples using input names

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The Model these selector arguments are for.
        """
        return tuple([selector_arg.index_tuple(model) for selector_arg in self])


class CompoundBoundingBox(_BoundingDomain):
    """
    A model's compound bounding box

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
    def __init__(self, bounding_boxes: Dict[Any, ModelBoundingBox], model = None,
                 selector_args: _SelectorArguments = None, create_selector: Callable = None,
                 ignored: List[str] = None, order: str = 'C'):

        self._create_selector = create_selector

        self._selector_args = selector_args
        self._bounding_boxes = bounding_boxes

        super().__init__(model, ignored, order)

    def _verify_selector_args(self):
        if self._selector_args is not None:
            self._selector_args = _SelectorArguments.validate(self.model, self._selector_args)

    def _pop_bounding_boxes(self):
        bounding_boxes = self._bounding_boxes
        self._bounding_boxes = {}

        return bounding_boxes

    def _verify_bounding_boxes(self, _external_ignored: List[str] = None):
        if self._selector_args is not None:
            bounding_boxes = self._pop_bounding_boxes()
            for selector, bounding_box in bounding_boxes.items():
                self.__setitem__(selector, bounding_box, _external_ignored)

    def _verify(self, _external_ignored: List[str] = None):
        if self._has_model:
            self._verify_selector_args()
            self._verify_bounding_boxes(_external_ignored)

    def copy(self):
        bounding_boxes = {selector: bbox.copy(self.ignored)
                          for selector, bbox in self._bounding_boxes.items()}

        return CompoundBoundingBox(bounding_boxes, self._model,
                                   selector_args=self._selector_args,
                                   create_selector=copy.deepcopy(self._create_selector),
                                   order=self._order)

    def __repr__(self):
        parts = ['CompoundBoundingBox(',
                 '    bounding_boxes={']
        # bounding_boxes
        for _selector, bbox in self._bounding_boxes.items():
            bbox_repr = bbox.__repr__().split('\n')
            parts.append(f"        {_selector} = {bbox_repr.pop(0)}")
            for part in bbox_repr:
                parts.append(f"            {part}")
        parts.append('    }')

        if len(self._ignored) > 0:
            parts.append(f"    ignored={self._ignored}")

        # selector_args
        selector_args_repr = self.selector_args.__repr__().split('\n')
        parts.append(f"    selector_args = {selector_args_repr.pop(0)}")
        for part in selector_args_repr:
            parts.append(f"        {part}")
        parts.append(')')

        return '\n'.join(parts)

    @property
    def ignored(self) -> List[str]:
        return super().ignored + self.selector_args.ignored

    @property
    def bounding_boxes(self) -> Dict[Any, ModelBoundingBox]:
        return self._bounding_boxes

    @property
    def selector_args(self) -> _SelectorArguments:
        if self._selector_args is None:
            raise RuntimeError("selector_args must be specified for a fully functional compound bounding box!")

        return self._selector_args

    @selector_args.setter
    def selector_args(self, value):
        if self._selector_args is not None:
            warnings.warn("Overriding selector_args may cause problems you should re-validate "
                          "the compound bounding box before use!", RuntimeWarning)
        self._selector_args = _SelectorArguments.validate(self.model, value)

    @property
    def create_selector(self):
        return self._create_selector

    @staticmethod
    def _get_selector_key(key):
        if isiterable(key):
            return tuple(key)
        else:
            return (key,)

    def __setitem__(self, key, value, _external_ignored: List[str] = None):
        if _external_ignored is None:
            _external_ignored = []
        _external_ignored.extend(self.ignored)

        selector = self._get_selector_key(key)
        if not self.selector_args.is_selector(selector):
            raise ValueError(f"{selector} is not a selector!")

        self._bounding_boxes[selector] = ModelBoundingBox.validate(self.model, value,
                                                                   _external_ignored=_external_ignored,
                                                                   order=self._order)

    def __eq__(self, value):
        if isinstance(value, CompoundBoundingBox):
            return (self.bounding_boxes == value.bounding_boxes) and \
                (self._selector_args == value._selector_args) and \
                (self.create_selector == value.create_selector)
        else:
            return False

    @classmethod
    def validate(cls, model, bounding_box: dict, selector_args=None, create_selector=None,
                 ignored: list = None, order: str = 'C',
                 _external_ignored: List[str] = None, **kwarg):
        """
        Construct a valid compound bounding box for a model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The model for which this will be a bounding_box
        bounding_box : dict
            Dictionary of possible bounding_box respresentations
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
            ignored = bounding_box._ignored
            bounding_box = bounding_box.bounding_boxes

        if selector_args is None:
            raise ValueError("Selector arguments must be provided (can be passed as part of bounding_box argument)!")

        new = cls(bounding_box, selector_args=selector_args,
                  create_selector=create_selector, ignored=ignored, order=order)
        new.verify(model, _external_ignored)

        return new

    def __contains__(self, key):
        return key in self._bounding_boxes

    def _create_bounding_box(self, _selector):
        self[_selector] = self._create_selector(_selector, model=self._model)

        return self[_selector]

    def __getitem__(self, key):
        selector = self._get_selector_key(key)
        if selector in self:
            return self._bounding_boxes[selector]
        elif self._create_selector is not None:
            return self._create_bounding_box(selector)
        else:
            raise RuntimeError(f"No bounding box is defined for selector: {selector}.")

    def _select_bounding_box(self, inputs) -> ModelBoundingBox:
        selector = self.selector_args.get_selector(self._model, *inputs)

        return self[selector]

    def prepare_inputs(self, input_shape, inputs, ignored: List[str] = None) -> Tuple[Any, Any, Any]:
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
        if ignored is None:
            ignored = []

        ignored.extend(self.ignored)
        bounding_box = self._select_bounding_box(inputs)

        return bounding_box.prepare_inputs(input_shape, inputs, ignored)

    def _fix_inputs_matching_bounding_boxes(self, argument, value,
                                            _external_ignored: List[str] = None) -> Dict[Any, ModelBoundingBox]:
        """
        Fix inputs for single input when fixed-input is a selector arg.
            - Returns bounding_boxes for a compound bounding box, if
              all selector arguments are removed, will be dict with a
              single entry, with key ().
        """
        if _external_ignored is None:
            _external_ignored = []
        _external_ignored.extend(self.ignored)

        selector_index = self.selector_args.selector_index(self._model, argument)
        matching = {}
        for selector_key, bbox in self._bounding_boxes.items():
            if selector_key[selector_index] == value:
                new_selector_key = list(selector_key)
                new_selector_key.pop(selector_index)

                if bbox.has_interval(argument):
                    new_bbox = bbox.fix_inputs(self._model, {argument: value},
                                               _compound_ignore=True,
                                               _external_ignored=_external_ignored)
                else:
                    new_bbox = bbox.copy(_external_ignored)

                matching[tuple(new_selector_key)] = new_bbox

        if len(matching) == 0:
            raise ValueError(f"Attempting to fix input {argument}, but there are no "
                             f"bounding boxes for argument value {value}.")

        return matching

    def _fix_input_selector_arg(self, argument, value,
                                _external_ignored: List[str] = None):
        """
        Fix inputs for single input when fixed-input is a selector arg.
            - Drops from compound bounding box to a normal bounding box
              if the selector arguments become redundant.
        """

        matching_bounding_boxes = self._fix_inputs_matching_bounding_boxes(argument, value, _external_ignored)
        if len(self.selector_args) == 1:
            return matching_bounding_boxes[()]
        else:
            return CompoundBoundingBox.validate(self.model, matching_bounding_boxes,
                                                selector_args=self.selector_args.reduce(self.model, argument),
                                                _external_ignored=[argument])

    def _fix_input_bbox_arg(self, argument, value,
                            _external_ignored: List[str] = None):
        """
        Fix inputs for single input when fixed-input is not part of the selector_args.
        """
        if _external_ignored is None:
            _external_ignored = []
        _external_ignored.extend(self.ignored)

        bounding_boxes = {}
        for selector_key, bbox in self._bounding_boxes.items():
            bounding_boxes[selector_key] = bbox.fix_inputs(self.model, {argument: value},
                                                           _external_ignored=_external_ignored.copy(),
                                                           _compound_ignore=True)

        return CompoundBoundingBox.validate(self.model, bounding_boxes,
                                            selector_args=self.selector_args,
                                            _external_ignored=[argument])

    def fix_inputs(self, model, fixed_inputs: dict,
                   _external_ignored: List[str] = None):
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
            bbox = self._fix_input_selector_arg(argument, value, _external_ignored)
        else:
            bbox = self._fix_input_bbox_arg(argument, value, _external_ignored)

        if len(fixed_input_keys) > 0:
            new_fixed_inputs = fixed_inputs.copy()
            del new_fixed_inputs[argument]

            bbox = bbox.fix_inputs(model, new_fixed_inputs,
                                   _external_ignored=[argument])

        if isinstance(bbox, CompoundBoundingBox):
            selector_args = bbox.selector_args
            bbox_dict = bbox
        elif isinstance(bbox, ModelBoundingBox):
            selector_args = None
            bbox_dict = bbox.intervals

        return bbox.__class__.validate(model, bbox_dict,
                                       order=bbox.order,
                                       selector_args=selector_args,
                                       _external_ignored=self.ignored)
