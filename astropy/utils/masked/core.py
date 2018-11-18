# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.
"""
from functools import reduce
import operator

import numpy as np


__all__ = ['Masked']


REDUCTION_FILL_VALUES = {
    np.add: 0,
    np.subtract: 0,
    np.multiply: 1}


class Masked:
    """A scalar value or array of values with associated mask.

    This object will take its exact type from whatever the contents are.
    In general this is expected to be an `~astropy.units.Quantity` or
    `numpy.ndarray`, although anything compatible with `numpy.asanyarray` is
    possible.

    Parameters
    ----------
    data : array-like
        The data for which a mask is to be added.
    """
    _generated_subclasses = {}

    def __new__(cls, data, mask=None):
        data_mask = getattr(data, 'mask', None)

        if data_mask is None:
            data = np.asanyarray(data)
            if mask is None:
                mask = False
        else:
            data = data.data
            if mask is None:
                mask = data_mask

        if data.dtype.names:
            raise NotImplementedError("cannot deal with structured dtype.")

        data_cls = type(data)
        if not issubclass(data_cls, Masked):
            # Make sure first letter is uppercase, but note that we can't use
            # str.capitalize since that converts the rest of the name to lowercase.

            new_cls = cls._generated_subclasses.get(data_cls, None)
            if new_cls is None:
                new_name = (cls.__name__ +
                            data_cls.__name__[0].upper() + data_cls.__name__[1:])
                new_cls = type(new_name, (cls, data_cls),
                               {'_data_cls': data_cls})
                cls._generated_subclasses[data_cls] = new_cls

        self = data.view(new_cls)
        self._data_cls = data_cls
        self.mask = mask
        return self

    @property
    def data(self):
        return self.view(self._data_cls)

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        self._mask = np.zeros(self.shape, dtype='?')
        self._mask[...] = mask

    def filled(self, fill_value):
        result = self.data.copy()
        result[self.mask] = fill_value
        return result

    def __getitem__(self, item):
        result = super().__getitem__(item)
        mask = self.mask[item]
        if isinstance(result, Masked):
            result._mask = mask
            return result
        else:
            return self.__class__(result, mask=mask)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.pop('out', None)
        if out is not None:
            kwargs['out'] = tuple((out_.data if hasattr(out_, 'mask')
                                   else out_) for out_ in out)
            if ufunc.nout == 1:
                out = out[0]

        if method == '__call__':
            converted = []
            masks = []
            for input_ in inputs:
                mask = getattr(input_, 'mask', None)
                if mask is not None:
                    masks.append(mask)
                    converted.append(input_.data)
                else:
                    converted.append(input_)

            result = getattr(ufunc, method)(*converted, **kwargs)
            if masks:
                mask = reduce(operator.or_, masks)
            else:
                mask = False

        elif method == 'reduce':

            fill_value = REDUCTION_FILL_VALUES.get(ufunc, None)
            if fill_value is None:
                return NotImplemented

            mask = getattr(inputs[0], 'mask', False)
            converted = inputs[0].filled(fill_value)
            result = getattr(ufunc, method)(converted, **kwargs)
            if mask is not False:
                mask = np.logical_and.reduce(mask, **kwargs)

        elif method in {'accumulate', 'reduceat'}:
            return NotImplemented

        if mask is False or result is None or result is NotImplemented:
            return result

        return self._masked_result(result, mask, out)

    def _masked_result(self, result, mask, out):
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            return tuple(self._masked_result(result_, mask, out_)
                         for (result_, out_) in zip(result, out))

        if out is None:
            return Masked(result, mask)

        assert isinstance(out, Masked)
        out._mask = mask
        return out

    # TODO: improve (greatly) repr and str!!
    def __repr__(self):
        reprarr = repr(self.data)
        if reprarr.endswith('>'):
            firstspace = reprarr.find(' ')
            reprarr = reprarr[firstspace+1:-1]  # :-1] removes the ending '>'
            return '<{} {} with mask={}>'.format(self.__class__.__name__,
                                                 reprarr, self.mask)
        else:  # numpy array-like
            firstparen = reprarr.find('(')
            reprarr = reprarr[firstparen:]
            return '{}{} with mask={}'.format(self.__class__.__name__,
                                              reprarr, self.mask)
            return reprarr

    def __str__(self):
        datastr = str(self.data)
        toadd = ' with mask={}'.format(self.mask)
        return datastr + toadd
