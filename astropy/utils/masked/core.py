# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.
"""
import numpy as np


__all__ = ['Masked']


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

        data_cls = type(data)
        if not issubclass(data_cls, Masked):
            # Make sure first letter is uppercase, but note that we can't use
            # str.capitalize since that converts the rest of the name to lowercase.
            new_name = cls.__name__ + data_cls.__name__[0].upper() + data_cls.__name__[1:]
            if new_name in cls._generated_subclasses:
                new_cls = cls._generated_subclasses[new_name]
            else:
                new_cls = type(new_name, (cls, data_cls),
                               {'_data_cls': data_cls})
                cls._generated_subclasses[new_name] = new_cls

        self = data.view(new_cls)
        self._data_cls = data_cls
        self._mask = np.zeros(data.shape, dtype='?')
        self._mask[...] = mask
        return self

    @property
    def data(self):
        return self.view(self._data_cls)

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, mask):
        self._mask[...] = mask

    def __getitem__(self, item):
        result = super().__getitem__(item)
        mask = self.mask[item]
        if isinstance(result, Masked):
            result._mask = mask
            return result
        else:
            return self.__class__(result, mask=mask)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        converted = []
        if method in {'reduce', 'accumulate', 'reduceat'}:
            raise NotImplementedError

        masks = []
        for input_ in inputs:
            mask = getattr(input_, 'mask', None)
            if mask is not None:
                masks.append(mask)
                converted.append(input_.data)
            else:
                converted.append(input_)

        results = getattr(ufunc, method)(*converted, **kwargs)
        if masks:
            mask = masks[0]
            for mask_ in masks[1:]:
                mask |= mask_
        else:
            mask = False

        if not isinstance(results, tuple):
            return Masked(results, mask)
        else:
            return [Masked(result, mask) for result in results]

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
