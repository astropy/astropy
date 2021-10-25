# Licensed under a 3-clause BSD style license - see LICENSE.rst


"""
This module is to contain an improved bounding box
"""

from collections import namedtuple
from typing import Dict

from astropy.utils import isiterable

import warnings
import numpy as np


__all__ = ['Interval', 'BoundingBox']


_BaseInterval = namedtuple('_BaseInterval', "lower upper")


class Interval(_BaseInterval):
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


class BoundingBox(object):
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

    order : optional, str
        The ordering that is assumed for the tuple representation of this
        bounding_box. Options: 'C': C/Python order, e.g. z, y, x.
        (default), 'F': Fortran/mathematical notation order, e.g. x, y, z.

    Methods
    -------
    validate :
        Constructs a valid bounding_box from any of the allowed
        respresentations of a bounding_box.

    bounding_box :
        Contructs a tuple respresentation

    domain :
        Contructs a discretization of the points inside the bounding_box

    prepare_inputs :
        Generates the necessary input information so that model can
        be evaluated only for input points entirely inside bounding_box.

    prepare_outputs :
        Fills the output values in for any input points outside the
        bounding_box.
    """
    def __init__(self, intervals: Dict[int, Interval], model, order: str = 'C'):
        self._model = model
        self._order = order

        self._intervals = {}
        if intervals != () and intervals != {}:
            self._validate(intervals, order=order)

    @property
    def intervals(self) -> Dict[int, Interval]:
        """Return bounding_box labeled using input positions"""
        return self._intervals

    def _get_name(self, index: int):
        """Get the input name corresponding to the input index"""
        return self._model.inputs[index]

    @property
    def named_intervals(self) -> Dict[str, Interval]:
        """Return bounding_box labeled using input names"""
        return {self._get_name(index): bbox for index, bbox in self._intervals.items()}

    def __repr__(self):
        parts = [
            'BoundingBox(',
            '    intervals={'
        ]

        for name, interval in self.named_intervals.items():
            parts.append(f"        {name}: {interval}")

        parts.append('    }')
        parts.append(f'    model={self._model.__class__.__name__}(inputs={self._model.inputs})')
        parts.append(f"    order='{self._order}'")
        parts.append(')')

        return '\n'.join(parts)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError(
            "This bounding box is fixed by the model and does not have "
            "adjustable parameters.")

    def _get_index(self, key):
        """
        Get the input index corresponding to the given key.
            Can pass in either:
                the string name of the input or
                the input index itself.
        """
        if isinstance(key, str):
            if key in self._model.inputs:
                index = self._model.inputs.index(key)
            else:
                raise ValueError(f"'{key}' is not one of the inputs: {self._model.inputs}.")
        elif np.issubdtype(type(key), np.integer):
            if key < len(self._model.inputs):
                index = key
            else:
                raise IndexError(f"Integer key: {key} must be < {len(self._model.inputs)}.")
        else:
            raise ValueError(f"Key value: {key} must be string or integer.")

        return index

    def __len__(self):
        return len(self._intervals)

    def __contains__(self, key):
        try:
            return self._get_index(key) in self._intervals
        except (IndexError, ValueError):
            return False

    def __getitem__(self, key):
        """Get bounding_box entries by either input name or input index"""
        return self._intervals[self._get_index(key)]

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

    def bounding_box(self, order: str = None):
        """
        Return the old tuple of tuples representation of the bounding_box
            order='C' corresponds to the old bounding_box ordering
            order='F' corresponds to the gwcs bounding_box ordering.
        """
        if len(self._intervals) == 1:
            return self._intervals[0]
        else:
            order = self._get_order(order)
            inputs = self._model.inputs
            if order == 'C':
                inputs = inputs[::-1]

            return tuple([self[input_name] for input_name in inputs])

    def __eq__(self, value):
        """Note equality can be either with old representation or new one."""
        if isinstance(value, tuple):
            return self.bounding_box() == value
        elif isinstance(value, BoundingBox):
            return self.intervals == value.intervals
        else:
            return False

    def __setitem__(self, key, value):
        """Validate and store interval under key (input index or input name)."""
        self._intervals[self._get_index(key)] = Interval.validate(value)

    def __delitem__(self, key):
        """Delete stored interval"""
        del self._intervals[self._get_index(key)]

    def _validate_dict(self, bounding_box: dict):
        """Validate passing dictionary of intervals and setting them."""
        for key, value in bounding_box.items():
            self[key] = value

    def _validate_sequence(self, bounding_box, order: str = None):
        """Validate passing tuple of tuples representation (or related) and setting them."""
        order = self._get_order(order)
        if order == 'C':
            # If bounding_box is C/python ordered, it needs to be reversed
            # to be in Fortran/mathematical/input order.
            bounding_box = bounding_box[::-1]

        for index, value in enumerate(bounding_box):
            self[index] = value

    def _validate_iterable(self, bounding_box, order: str = None):
        """Validate and set any iterable representation"""
        if len(bounding_box) != self._model.n_inputs:
            raise ValueError(f"Found {len(bounding_box)} intervals, "
                             f"but must have exactly {self._model.n_inputs}.")

        if isinstance(bounding_box, dict):
            self._validate_dict(bounding_box)
        else:
            self._validate_sequence(bounding_box, order)

    def _validate(self, bounding_box, order: str = None):
        """Validate and set any representation"""
        if self._model.n_inputs == 1 and not isinstance(bounding_box, dict):
            self[0] = bounding_box
        else:
            self._validate_iterable(bounding_box, order)

    @classmethod
    def validate(cls, model, bounding_box,
                 order: str = 'C'):
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
        if isinstance(bounding_box, BoundingBox):
            bounding_box = bounding_box.intervals

        new = cls({}, model, order)
        new._validate(bounding_box)

        return new

    def copy(self):
        return BoundingBox(self._intervals.copy(), self._model, self._order)

    def fix_inputs(self, model, fixed_inputs: list):
        """
        Fix the bounding_box for a `fix_inputs` compound model.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The new model for which this will be a bounding_box
        fixed_inputs : list
            List if inputs that have been fixed in this bounding_box
        """

        new = self.copy()

        for _input in fixed_inputs:
            del new[_input]

        return BoundingBox.validate(model, new.named_intervals, new._order)

    @property
    def dimension(self):
        return len(self)

    def domain(self, resolution, order: str = None):
        inputs = self._model.inputs
        order = self._get_order(order)
        if order == 'C':
            inputs = inputs[::-1]

        return [self[input_name].domain(resolution) for input_name in inputs]

    def _outside(self,  input_shape, inputs):
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
            outside = self[index].outside(_input)

            if _input.shape:
                outside_index[outside] = True
            else:
                outside_index |= outside
                if outside_index:
                    all_out = True

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

    def prepare_inputs(self, input_shape, inputs):
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
                    valid_inputs.append(np.atleast_1d(_input)[valid_index])
                else:
                    valid_inputs.append(_input)

        return valid_inputs, valid_index, all_out

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
