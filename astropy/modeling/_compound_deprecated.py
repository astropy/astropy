# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains implementations of the old compound model interfaces that
are now deprecated.

The classes exported by this module should be imported from `astropy.modeling`
rather than from `astropy.modeling.core` or from this module directly.
"""

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

import functools
import operator
import warnings

import numpy as np

from ..utils import deprecated, indent
from ..utils.compat.odict import OrderedDict
from ..utils.exceptions import AstropyDeprecationWarning
from .core import Model


__all__ = ['LabeledInput', 'SerialCompositeModel', 'SummedCompositeModel']


@deprecated('1.0', alternative=':ref:`compound-models` as described in the '
                               'Astropy documentation')
class LabeledInput(OrderedDict):
    """
    Used by `SerialCompositeModel` and `SummedCompositeModel` to choose input
    data using labels.

    This is a container assigning labels (names) to all input data arrays to a
    composite model.

    Parameters
    ----------
    data : list
        List of all input data
    labels : list of strings
        names matching each coordinate in data

    Examples
    --------
    >>> y, x = np.mgrid[:5, :5]
    >>> l = np.arange(10)
    >>> labeled_input = LabeledInput([x, y, l], ['x', 'y', 'pixel'])
    >>> labeled_input.x
    array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]])
    >>> labeled_input['x']
    array([[0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4],
           [0, 1, 2, 3, 4]])
    """

    def __init__(self, data, labels):
        if len(labels) != len(data):
            raise TypeError("Number of labels and data doesn't match")

        super(LabeledInput, self).__init__(zip(labels, data))

    def __getattr__(self, label):
        try:
            return self[label]
        except KeyError:
            raise AttributeError(label)

    def __setattr__(self, label, data):
        if label.startswith('_'):
            super(LabeledInput, self).__setattr__(label, data)
        else:
            self[label] = data

    def __delattr__(self, label):
        try:
            del self[label]
        except KeyError:
            raise AttributeError(label)

    @property
    def labels(self):
        return tuple(self.keys())

    def add(self, label=None, value=None, **kw):
        """
        Add input data to a LabeledInput object

        Parameters
        --------------
        label : str
            coordinate label
        value : numerical type
            coordinate value
        kw : dictionary
            if given this is a dictionary of ``{label: value}`` pairs
        """

        if ((label is None and value is not None) or
                (label is not None and value is None)):
            raise TypeError("Expected label and value to be defined")

        kw[label] = value

        self.update(kw)

    def copy(self):
        return LabeledInput(self.values(), self.labels)


class _LabeledInputMapping(Model):
    def __init__(self, labeled_input, inmap, outmap):
        self._labeled_input = labeled_input
        self._inmap = tuple(inmap)
        self._outmap = tuple(outmap)
        super(_LabeledInputMapping, self).__init__()

    def __repr__(self):
        return '<{0}>'.format(self.name)

    @property
    def inputs(self):
        return self._outmap

    @property
    def outputs(self):
        return self._inmap

    @property
    def name(self):
        return '{0}({1} -> {2})'.format(self.__class__.__name__,
                                        self._outmap, self._inmap)

    def evaluate(self, *inputs):
        for idx, label in enumerate(self._outmap):
            self._labeled_input[label] = inputs[idx]

        result = tuple(self._labeled_input[label] for label in self._inmap)

        if len(result) == 1:
            return result[0]
        else:
            return result


class _CompositeModel(Model):
    """Base class for all composite models."""

    _operator = None
    fittable = False

    def __init__(self, transforms, n_inputs, n_outputs, inmap=None,
                 outmap=None):
        self._transforms = transforms
        param_names = []
        for tr in self._transforms:
            param_names.extend(tr.param_names)
        super(_CompositeModel, self).__init__()
        self.param_names = param_names
        self._n_inputs = n_inputs
        self._n_outputs = n_outputs
        self._basic_transform = None

        self._inmap = inmap
        self._outmap = outmap

    def __repr__(self):
        return '<{0}([\n{1}\n])>'.format(
            self.__class__.__name__,
            indent(',\n'.join(repr(tr) for tr in self._transforms),
                   width=4))

    def __str__(self):
        parts = ['Model: {0}'.format(self.__class__.__name__)]
        for tr in self._transforms:
            parts.append(indent(str(tr), width=4))
        return '\n'.join(parts)

    @property
    def inputs(self):
        return self._transforms[0].inputs

    @property
    def outputs(self):
        return self._transforms[-1].outputs

    @property
    def n_inputs(self):
        return self._n_inputs

    @n_inputs.setter
    def n_inputs(self, val):
        warnings.warn(
            'Setting n_inputs on {0} objects is undefined and should not '
            'be used.'.format(self.__class__.__name__),
            AstropyDeprecationWarning)
        self._n_inputs = val

    @property
    def n_outputs(self):
        return self._n_outputs

    @n_outputs.setter
    def n_outputs(self, val):
        warnings.warn(
            'Setting n_outputs on {0} objects is undefined and should not '
            'be used.'.format(self.__class__.__name__),
            AstropyDeprecationWarning)
        self._n_outputs = val

    def invert(self):
        raise NotImplementedError("Subclasses should implement this")

    @property
    def parameters(self):
        raise NotImplementedError(
            "Composite models do not currently support the .parameters "
            "array.")

    def evaluate(self, *inputs):
        """
        Specialized `Model.evaluate` implementation that allows `LabeledInput`
        inputs to be handled when calling this model.

        This ignores any passed in parameter values, as _CompositeModels can't
        be fitted anyways.
        """

        # Drop parameter arguments
        inputs = inputs[:self.n_inputs]

        if len(inputs) == 1 and isinstance(inputs[0], LabeledInput):
            labeled_input = inputs[0].copy()
            transform = self._make_labeled_transform(labeled_input)
            inputs = [labeled_input[label] for label in self._inmap[0]]
            result = transform(*inputs)

            if self._transforms[-1].n_outputs == 1:
                labeled_input[self._outmap[-1][0]] = result
            else:
                for label, output in zip(self._outmap[-1], result):
                    labeled_input[label] = output

            return labeled_input
        else:
            if self._basic_transform is None:
                transform = self._transforms[0]
                for t in self._transforms[1:]:
                    transform = self._operator(transform, t)

                self._basic_transform = transform

            return self._basic_transform(*inputs)

    def __call__(self, *inputs):
        """
        Specialized `Model.__call__` implementation that allows
        `LabeledInput` inputs to be handled when calling this model.
        """

        return self.evaluate(*inputs)

    def _param_sets(self, raw=False):
        all_params = tuple(m._param_sets(raw=raw) for m in self._transforms)
        return np.vstack(all_params)

    def _make_labeled_transform(self, labeled_input):
        """
        Build up a transformation graph that incorporates the instructions
        encoded in the `LabeledInput` object.

        This requires use of the ``_inmap`` and ``_outmap`` attributes set
        when instantiating this `_CompositeModel`.
        """

        if self._inmap is None:
            raise TypeError("Parameter 'inmap' must be provided when "
                            "input is a labeled object.")
        if self._outmap is None:
            raise TypeError("Parameter 'outmap' must be provided when "
                            "input is a labeled object")

        transforms = [self._transforms[0]]
        previous_outmap = self._outmap[0]
        for model, inmap, outmap in zip(self._transforms[1:], self._inmap[1:],
                                        self._outmap[1:]):
            mapping = _LabeledInputMapping(labeled_input, inmap,
                                           previous_outmap)
            transforms.append(mapping | model)
            previous_outmap = outmap

        return functools.reduce(self._operator, transforms)


@deprecated('1.0', alternative=':ref:`compound-models` as described in the '
                               'Astropy documentation')
class SerialCompositeModel(_CompositeModel):
    """
    Composite model that evaluates models in series.

    Parameters
    ----------
    transforms : list
        a list of transforms in the order to be executed
    inmap : list of lists or None
        labels in an input instance of LabeledInput
        if None, the number of input coordinates is exactly what
        the transforms expect
    outmap : list or None
        labels in an input instance of LabeledInput
        if None, the number of output coordinates is exactly what
        the transforms expect
    n_inputs : int
        dimension of input space (e.g. 2 for a spatial model)
    n_outputs : int
        dimension of output

    Notes
    -----
    Output values of one model are used as input values of another.
    Obviously the order of the models matters.

    Examples
    --------
    Apply a 2D rotation followed by a shift in x and y::

        >>> import numpy as np
        >>> from astropy.modeling import models, LabeledInput, SerialCompositeModel
        >>> y, x = np.mgrid[:5, :5]
        >>> rotation = models.Rotation2D(angle=23.5)
        >>> offset_x = models.Shift(-4.23)
        >>> offset_y = models.Shift(2)
        >>> labeled_input = LabeledInput([x, y], ["x", "y"])
        >>> transform = SerialCompositeModel([rotation, offset_x, offset_y],
        ...                                  inmap=[['x', 'y'], ['x'], ['y']],
        ...                                  outmap=[['x', 'y'], ['x'], ['y']])
        >>> result = transform(labeled_input)
    """

    _operator = operator.or_

    def __init__(self, transforms, inmap=None, outmap=None, n_inputs=None,
                 n_outputs=None):
        if n_inputs is None:
            n_inputs = max([tr.n_inputs for tr in transforms])
            # the output dimension is equal to the output dim of the last
            # transform
            n_outputs = transforms[-1].n_outputs
        else:
            if n_outputs is None:
                raise TypeError("Expected n_inputs and n_outputs")

        if transforms and inmap and outmap:
            if not (len(transforms) == len(inmap) == len(outmap)):
                raise ValueError("Expected sequences of transform, "
                                 "inmap and outmap to have the same length")

        super(SerialCompositeModel, self).__init__(
                transforms, n_inputs, n_outputs, inmap=inmap, outmap=outmap)

    def inverse(self):
        try:
            transforms = []
            for transform in self._transforms[::-1]:
                transforms.append(transform.inverse)
        except NotImplementedError:
            raise NotImplementedError(
                "An analytical inverse has not been implemented for "
                "{0} models.".format(transform.__class__.__name__))
        if self._inmap is not None:
            inmap = self._inmap[::-1]
            outmap = self._outmap[::-1]
        else:
            inmap = None
            outmap = None
        return SerialCompositeModel(transforms, inmap, outmap)


@deprecated('1.0', alternative=':ref:`compound-models` as described in the '
                               'Astropy documentation')
class SummedCompositeModel(_CompositeModel):
    """
    Composite model that evaluates models in parallel.

    Parameters
    --------------
    transforms : list
        transforms to be executed in parallel
    inmap : list or None
        labels in an input instance of LabeledInput
        if None, the number of input coordinates is exactly what the
        transforms expect
    outmap : list or None

    Notes
    -----
    Evaluate each model separately and add the results to the input_data.
    """

    _operator = operator.add

    def __init__(self, transforms, inmap=None, outmap=None):
        n_inputs = transforms[0].n_inputs
        n_outputs = n_inputs
        for transform in transforms:
            if not (transform.n_inputs == transform.n_outputs == n_inputs):
                raise ValueError("A SummedCompositeModel expects n_inputs = "
                                 "n_outputs for all transforms")

        super(SummedCompositeModel, self).__init__(transforms, n_inputs,
                                                   n_outputs, inmap=inmap,
                                                   outmap=outmap)
