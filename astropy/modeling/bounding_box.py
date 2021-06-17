# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections import UserDict
from .utils import _BoundingBox

import numpy as np


class BoundingBox(_BoundingBox):
    @classmethod
    def validate(cls, model, bounding_box, slice_arg=None):
        """
        Validate a given bounding box sequence against the given model (which
        may be either a subclass of `~astropy.modeling.Model` or an instance
        thereof, so long as the ``.inputs`` attribute is defined.

        Currently this just checks that the bounding_box is either a 2-tuple
        of lower and upper bounds for 1-D models, or an N-tuple of 2-tuples
        for N-D models.

        This also returns a normalized version of the bounding_box input to
        ensure it is always an N-tuple (even for the 1-D case).
        """

        if isinstance(bounding_box, dict) or isinstance(bounding_box, CompoundBoundingBox):
            return CompoundBoundingBox.validate(model, bounding_box, slice_arg=slice_arg)

        nd = model.n_inputs
        if slice_arg is not None:
            if isinstance(slice_arg, tuple):
                nd -= len(slice_arg)
            else:
                nd -= 1

        return cls._validate(model, bounding_box, nd)

    def reverse(self, axes_ind):
        bbox = [self[ind] for ind in axes_ind][::-1]

        return _BoundingBox(bbox, self._model)


class CompoundBoundingBox(UserDict):
    def __init__(self, bounding_box,  model=None, slice_arg=None, remove_slice_arg=False):
        super().__init__(bounding_box)
        self._model = model

        if self._model is None:
            self._slice_arg = slice_arg
        else:
            self.set_slice_arg(slice_arg)

        self._remove_slice_arg = remove_slice_arg

    @property
    def slice_arg(self):
        return self._slice_arg

    @property
    def remove_slice_arg(self):
        return self._remove_slice_arg

    def _get_arg_index(self, slice_arg):
        if np.issubdtype(type(slice_arg), np.integer):
            arg_index = slice_arg
        else:
            if slice_arg in self._model.inputs:
                arg_index = self._model.inputs.index(slice_arg)
            else:
                raise ValueError(f'{slice_arg} is not an input of of your model inputs {self._model.inputs}')

        if arg_index < self._model.n_inputs:
            return arg_index
        else:
            raise ValueError(f'{arg_index} is out of model argument bounds')

    @classmethod
    def validate(cls, model, bounding_box, slice_arg=None, remove_slice_arg=None):
        if isinstance(bounding_box, CompoundBoundingBox) and slice_arg is None:
            slice_arg = bounding_box.slice_arg

        if remove_slice_arg is None:
            if isinstance(bounding_box, CompoundBoundingBox):
                remove_slice_arg = bounding_box.remove_slice_arg
            else:
                remove_slice_arg = False

        new_box = cls({}, model, slice_arg, remove_slice_arg)

        if not remove_slice_arg:
            slice_arg = None

        for slice_index, slice_box in bounding_box.items():
            new_box[slice_index] = BoundingBox.validate(model, slice_box, slice_arg)

        return new_box

    def set_slice_arg(self, slice_arg):
        if slice_arg is None:
            self._slice_arg = slice_arg
        elif isinstance(slice_arg, tuple):
            self._slice_arg = tuple([self._get_arg_index(arg) for arg in slice_arg])
        else:
            self._slice_arg = self._get_arg_index(slice_arg)

    def _get_slice_index(self, inputs, slice_arg):
        if inputs is None:
            raise RuntimeError('Inputs must not be None in order to lookup slice')

        slice_index = inputs[self._get_arg_index(slice_arg)]
        if isinstance(slice_index, np.ndarray):
            slice_index = slice_index.item()

        return slice_index

    def get_bounding_box(self, inputs=None, slice_index=True):
        if isinstance(slice_index, bool) and slice:
            if self._slice_arg is None:
                return None
            else:
                if isinstance(self._slice_arg, tuple):
                    slice_index = tuple([self._get_slice_index(inputs, slice_arg)
                                         for slice_arg in self._slice_arg])
                else:
                    slice_index = self._get_slice_index(inputs, self._slice_arg)

        if slice_index in self:
            return self[slice_index]
        else:
            raise RuntimeError(f"No bounding_box is defined for slice: {slice_index}!")

    def reverse(self, axes_ind):
        # bbox = {slice_index: slice_box.reverse(axes_ind)
        #         for slice_index, slice_box in self.items()}
        bbox = {}
        for slice_index, slice_box in self.items():
            bbox[slice_index] = [slice_box[ind] for ind in axes_ind][::-1]

        return CompoundBoundingBox(bbox, self._model, self._slice_arg, self._remove_slice_arg)

    def py_order(self, axes_order):
        bbox = {}
        for slice_index, slice_box in self.items():
            bbox[slice_index] = tuple(slice_box[::-1][i] for i in axes_order)

        return CompoundBoundingBox(bbox, self._model, self._slice_arg, self._remove_slice_arg)
