# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.
"""
from functools import reduce
import operator

import numpy as np


__all__ = ['Masked']


# Note: at present, cannot directly use ShapedLikeNDArray, because of
# metaclass problems.
# TODO: split this off in misc.py
class NDArrayShapeMethods:
    def __getitem__(self, item):
        try:
            return self._apply('__getitem__', item)
        except IndexError:
            if self.shape == ():
                raise TypeError('scalar {0!r} object is not subscriptable.'
                                .format(self.__class__.__name__))
            else:
                raise

    def copy(self, *args, **kwargs):
        """Return an instance containing copies of the internal data.

        Parameters are as for :meth:`~numpy.ndarray.copy`.
        """
        return self._apply('copy', *args, **kwargs)

    def reshape(self, *args, **kwargs):
        """Returns an instance containing the same data with a new shape.

        Parameters are as for :meth:`~numpy.ndarray.reshape`.  Note that it is
        not always possible to change the shape of an array without copying the
        data (see :func:`~numpy.reshape` documentation). If you want an error
        to be raise if the data is copied, you should assign the new shape to
        the shape attribute (note: this may not be implemented for all classes
        using ``ShapedLikeNDArray``).
        """
        return self._apply('reshape', *args, **kwargs)

    def ravel(self, *args, **kwargs):
        """Return an instance with the array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.ravel`. Note that it is
        not always possible to unravel an array without copying the data.
        If you want an error to be raise if the data is copied, you should
        should assign shape ``(-1,)`` to the shape attribute.
        """
        return self._apply('ravel', *args, **kwargs)

    def flatten(self, *args, **kwargs):
        """Return a copy with the array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.flatten`.
        """
        return self._apply('flatten', *args, **kwargs)

    def transpose(self, *args, **kwargs):
        """Return an instance with the data transposed.

        Parameters are as for :meth:`~numpy.ndarray.transpose`.  All internal
        data are views of the data of the original.
        """
        return self._apply('transpose', *args, **kwargs)

    def swapaxes(self, *args, **kwargs):
        """Return an instance with the given axes interchanged.

        Parameters are as for :meth:`~numpy.ndarray.swapaxes`:
        ``axis1, axis2``.  All internal data are views of the data of the
        original.
        """
        return self._apply('swapaxes', *args, **kwargs)

    def diagonal(self, *args, **kwargs):
        """Return an instance with the specified diagonals.

        Parameters are as for :meth:`~numpy.ndarray.diagonal`.  All internal
        data are views of the data of the original.
        """
        return self._apply('diagonal', *args, **kwargs)

    def squeeze(self, *args, **kwargs):
        """Return an instance with single-dimensional shape entries removed

        Parameters are as for :meth:`~numpy.ndarray.squeeze`.  All internal
        data are views of the data of the original.
        """
        return self._apply('squeeze', *args, **kwargs)

    def take(self, indices, axis=None, mode='raise'):
        """Return a new instance formed from the elements at the given indices.

        Parameters are as for :meth:`~numpy.ndarray.take`, except that,
        obviously, no output array can be given.
        """
        return self._apply('take', indices, axis=axis, mode=mode)


class Masked:
    """A scalar value or array of values with associated mask.

    This object will take its exact type from whatever the contents are.

    Parameters
    ----------
    data : array-like
        The data for which a mask is to be added.  The result will be a
        a subclass of the type of ``data``.
    mask : array-like of bool, optional
        The initial mask to assign.  If not given, taken from the data.

    """
    _generated_subclasses = {}

    _REDUCTION_HELPERS = {
        np.minimum: lambda ma: np.max(ma.data),
        np.maximum: lambda ma: np.min(ma.data)}
    """For ufuncs without an identity but where reduction makes sense."""

    def __new__(cls, data, mask=None):
        data, data_mask = cls._data_mask(data)
        data = np.asanyarray(data)

        if mask is None:
            mask = False if data_mask is None else data_mask

        if data.dtype.names:
            raise NotImplementedError("cannot deal with structured dtype.")

        data_cls = type(data)
        masked_cls = cls._generated_subclasses.get(data_cls, None)
        if masked_cls is None:
            # Make sure first letter is uppercase, but note that we can't use
            # str.capitalize since that converts the rest of the name to lowercase.
            new_name = (cls.__name__ +
                        data_cls.__name__[0].upper() + data_cls.__name__[1:])
            masked_cls = type(new_name,
                              (cls, NDArrayShapeMethods, data_cls), {})
            cls._generated_subclasses[data_cls] = masked_cls

        self = data.view(masked_cls)
        self._data = data
        self.mask = mask
        return self

    @staticmethod
    def _data_mask(data):
        mask = getattr(data, 'mask', None)
        if mask is not None:
            data = getattr(data, 'data', data)

        return data, mask

    @property
    def mask(self):
        """The mask.

        If set, replace the original mask, with whatever it is set with,
        using a view if no broadcasting or type conversion is required.
        """
        return self._mask

    @mask.setter
    def mask(self, mask):
        ma = np.asanyarray(mask, dtype='?')
        if ma.shape != self.shape:
            # This will fail (correctly) if not broadcastable.
            self._mask = np.empty(self.shape, dtype='?')
            self._mask[...] = ma
        elif ma is mask:
            # Use a view so that shape setting does not propagate.
            self._mask = mask.view()
        else:
            self._mask = ma

    # TODO: should this have a more unique name? E.g., unmask?
    # Then, fill_value=None could indicate to not fill.
    def filled(self, fill_value):
        result = self.data.copy()
        result[self.mask] = fill_value
        return result

    # TODO: should this have a more unique name? E.g., unmasked?
    # (could just be unmasked = property(unmask) with above)
    @property
    def data(self):
        return self._data

    def _apply(self, method, *args, **kwargs):
        # For use with NDArrayShapeMethods.
        if callable(method):
            data = method(self.data, *args, **kwargs)
            mask = method(self.mask, *args, **kwargs)
        else:
            data = getattr(self.data, method)(*args, **kwargs)
            mask = getattr(self.mask, method)(*args, **kwargs)

        return self.__class__(data, mask=mask)

    def __setitem__(self, item, value):
        value, mask = self._data_mask(value)
        self._data[item] = value
        self._mask[item] = mask

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.pop('out', None)
        if out is not None:
            kwargs['out'] = tuple((out_.data if isinstance(out_, Masked)
                                   else out_) for out_ in out)
            if ufunc.nout == 1:
                out = out[0]

        if method == '__call__':
            converted = []
            masks = []
            for input_ in inputs:
                if isinstance(input_, Masked):
                    masks.append(input_.mask)
                    converted.append(input_.data)
                else:
                    converted.append(input_)

            result = getattr(ufunc, method)(*converted, **kwargs)
            if masks:
                mask = reduce(operator.or_, masks)
            else:
                mask = False

        elif method == 'reduce':
            if isinstance(inputs[0], Masked):
                initial = kwargs.pop('initial', None)
                # Calculate mask *without* a possible initial kwarg
                # (but passing on axis and where arguments, etc.).
                mask = np.logical_and.reduce(inputs[0].mask, **kwargs)
                converted = inputs[0].data
                # But use it for real reduction (or provide it, as needed).
                if initial is None:
                    if ufunc.identity is None:
                        helper = self._REDUCTION_HELPERS.get(ufunc, None)
                        if helper:
                            kwargs['initial'] = helper(inputs[0])
                        else:
                            return NotImplemented
                else:
                    kwargs['initial'] = initial

                # Use where to take into account the mask.
                where = ~self._mask
                if 'where' in kwargs:
                    # In-place in our inverted mask to ensure shape is OK.
                    where &= kwargs['where']

                kwargs['where'] = where

            elif 'out' in kwargs and isinstance(kwargs['out'][0], Masked):
                mask = False
                converted = inputs[0]
            else:
                return NotImplemented

            result = getattr(ufunc, method)(converted, **kwargs)

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
