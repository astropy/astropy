# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.
"""
import functools

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


class Masked(NDArrayShapeMethods):
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
    _data_cls = np.ndarray
    _mask = None

    def __new__(cls, data, mask=None, copy=False):
        data = np.array(data, subok=True, copy=copy)
        data, data_mask = cls._data_mask(data)
        if mask is None:
            mask = False if data_mask is None else data_mask

        if data.dtype.names:
            raise NotImplementedError("cannot deal with structured dtype.")

        subclass = cls._get_subclass(data.__class__)
        self = data.view(subclass)
        self._set_mask(mask, copy=copy)
        return self

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        subclass = cls.__mro__[1]
        # If not explicitly defined for this class, use defaults.
        # (TODO: metaclass better in this case?)
        if '_data_cls' not in cls.__dict__:
            cls._data_cls = subclass
        cls._generated_subclasses[subclass] = cls

    @classmethod
    def _get_subclass(cls, subclass):
        if subclass is np.ndarray:
            return MaskedNDArray
        elif subclass is None:
            return cls
        elif issubclass(subclass, Masked):
            return subclass

        if not issubclass(subclass, np.ndarray):
            raise ValueError('can only pass in an ndarray subtype.')

        masked_subclass = cls._generated_subclasses.get(subclass)
        if masked_subclass is None:
            # Create (and therefore register) new Measurement subclass.
            new_name = 'Masked' + subclass.__name__
            # Walk through MRO and find closest generated class
            # (e.g., MaskedQuantity for Angle).
            for mro_item in subclass.__mro__:
                base_cls = cls._generated_subclasses.get(mro_item)
                if base_cls is not None:
                    break
            else:
                base_cls = MaskedNDArray

            masked_subclass = type(new_name, (subclass, base_cls), {})

        return masked_subclass

    def view(self, dtype=None, type=None):
        if type is None and issubclass(dtype, np.ndarray):
            type = dtype
            dtype = None
        elif dtype is not None:
            raise ValueError('{} cannot be viewed with new dtype.'
                             .format(self.__class__))
        return super().view(self._get_subclass(type))

    def __array_finalize__(self, obj):
        if obj is None or type(obj) is np.ndarray:
            return None

        super_array_finalize = super().__array_finalize__
        if super_array_finalize:
            super_array_finalize(obj)

        if self._mask is None:
            self._mask = getattr(obj, '_mask', None)

    @staticmethod
    def _data_mask(data):
        mask = getattr(data, 'mask', None)
        if mask is not None:
            if isinstance(data, MaskedNDArray):
                data = data.unmasked
            elif hasattr(data, 'filled'):
                data = data.filled()

        return data, mask

    def _get_mask(self):
        """The mask.

        If set, replace the original mask, with whatever it is set with,
        using a view if no broadcasting or type conversion is required.
        """
        return self._mask

    def _set_mask(self, mask, copy=False):
        ma = np.asanyarray(mask, dtype='?')
        if ma.shape != self.shape:
            # This will fail (correctly) if not broadcastable.
            self._mask = np.empty(self.shape, dtype='?')
            self._mask[...] = ma
        elif ma is mask:
            # Even if not copying use a view so that shape setting
            # does not propagate.
            self._mask = mask.copy() if copy else mask.view()
        else:
            self._mask = ma

    mask = property(_get_mask, _set_mask)

    def unmask(self, fill_value=None):
        result = super().view(self._data_cls)
        if fill_value is None:
            return result
        else:
            result = result.copy()
            result[self.mask] = fill_value
            return result

    unmasked = property(unmask)

    def _apply(self, method, *args, **kwargs):
        # For use with NDArrayShapeMethods.
        if callable(method):
            data = method(self.unmasked, *args, **kwargs)
            mask = method(self._mask, *args, **kwargs)
        else:
            data = getattr(self.unmasked, method)(*args, **kwargs)
            mask = getattr(self._mask, method)(*args, **kwargs)

        result = data[...].view(self.__class__)
        result._mask = mask
        return result

    def __setitem__(self, item, value):
        value, mask = self._data_mask(value)
        super().__setitem__(item, value)
        self._mask[item] = mask

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.pop('out', None)
        if out is not None:
            if any(out_ is not None and not isinstance(out_, Masked)
                   for out_ in out):
                return NotImplemented
            kwargs['out'] = tuple((None if out_ is None else out_._data)
                                  for out_ in out)
            if ufunc.nout == 1:
                out = out[0]

        if method == '__call__':
            converted = []
            masks = []
            for input_ in inputs:
                if isinstance(input_, Masked):
                    masks.append(input_.mask)
                    converted.append(input_.unmasked)
                else:
                    converted.append(input_)

            result = getattr(ufunc, method)(*converted, **kwargs)
            if masks:
                mask = functools.reduce(np.logical_or, masks)
            else:
                mask = False

        elif method == 'reduce':
            if isinstance(inputs[0], Masked):
                # By default, we simply propagate masks, since for
                # things like np.sum, it makes no sense to do otherwise.
                # Individual methods need to override as needed.
                mask = np.logical_or.reduce(inputs[0].mask, **kwargs)
                converted = inputs[0].unmasked

            elif 'out' in kwargs and isinstance(kwargs['out'][0], Masked):
                mask = False
                converted = inputs[0]
            else:
                return NotImplemented

            result = getattr(ufunc, method)(converted, **kwargs)

        elif method in {'accumulate', 'reduceat', 'at'}:
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

    def _masked_function(self, function, axis=None, out=None, keepdims=False,
                         where=True, **kwargs):
        if out is not None:
            if not isinstance(out, Masked):
                raise TypeError('output should be a Masked instance')
            out_mask = out.mask
            out_data = out.unmasked
        else:
            out_mask = out_data = None

        # Output mask is set if *all* elements are masked.  Calculate *without*
        # a possible initial kwarg (but passing on axis, where, etc.).
        mask = np.logical_and.reduce(self.mask, axis=axis, out=out_mask,
                                     keepdims=keepdims, where=where)
        where_data = ~self.mask
        if where is not True:
            # In-place in our inverted mask to ensure shape is OK.
            where_data &= where

        result = function(self.unmasked, axis=axis, out=out_data,
                          keepdims=keepdims, where=where_data, **kwargs)
        if out is None:
            out = Masked(result, mask)

        return out

    def min(self, axis=None, out=None, keepdims=False, initial=None,
            where=True):
        if initial is None:
            initial = self.unmasked.max()
        return self._masked_function(np.min, axis=axis, out=out,
                                     keepdims=False, initial=initial,
                                     where=where)

    def max(self, axis=None, out=None, keepdims=False, initial=None,
            where=True):
        if initial is None:
            initial = self.unmasked.min()
        return self._masked_function(np.max, axis=axis, out=out,
                                     keepdims=False, initial=initial,
                                     where=where)

    def mean(self, axis=None, dtype=None, out=None, keepdims=False):
        result = self._masked_function(np.sum, axis=axis, out=out,
                                       keepdims=keepdims, dtype=dtype)
        n = np.add.reduce(~self.mask, axis=axis, keepdims=keepdims)
        result /= n
        return result

    # TODO: improve (greatly) repr and str!!
    def __repr__(self):
        reprarr = repr(self._data)
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
        datastr = str(self._data)
        toadd = ' with mask={}'.format(self.mask)
        return datastr + toadd


class MaskedNDArray(Masked, np.ndarray):
    _data_cls = np.ndarray
