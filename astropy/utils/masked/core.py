# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.
"""
import builtins
import functools

import numpy as np

from astropy.utils.shapes import NDArrayShapeMethods

from .function_helpers import APPLY_TO_BOTH_FUNCTIONS


__all__ = ['Masked']


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
        if type is None and (isinstance(dtype, builtins.type) and
                             issubclass(dtype, np.ndarray)):
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
            kwargs['out'] = tuple((None if out_ is None else out_.unmasked)
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
                if ufunc.signature:
                    if ufunc.nin > 1:
                        return NotImplemented
                    # TODO: need to check for axes, keepdims too...
                    mask = np.logical_or.reduce(masks[0],
                                                axis=kwargs.get('axis', -1))
                else:
                    mask = functools.reduce(np.logical_or, masks)
            else:
                mask = False

        elif method == 'reduce':
            if isinstance(inputs[0], Masked):
                # By default, we simply propagate masks, since for
                # things like np.sum, it makes no sense to do otherwise.
                # Individual methods need to override as needed.
                # TODO take care of 'out' too?
                axis = kwargs.get('axis', None)
                keepdims = kwargs.get('keepdims', False)
                where = kwargs.get('where', True)
                mask = np.logical_or.reduce(inputs[0].mask, where=where,
                                            axis=axis, keepdims=keepdims)
                if where is not True:
                    # Mask also whole rows that were not selected by where.
                    mask |= np.logical_and.reduce(inputs[0].mask, where=where,
                                                  axis=axis, keepdims=keepdims)
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

    def __array_function__(self, function, types, args, kwargs):
        if function in APPLY_TO_BOTH_FUNCTIONS:
            helper = APPLY_TO_BOTH_FUNCTIONS[function]
            data, masks, args, kwargs, out = helper(*args, **kwargs)
            mask = function(masks, *args, **kwargs)
            if out is not None:
                if not isinstance(out, Masked):
                    return NotImplemented
                kwargs['out'] = out.unmasked
            result = function(data, *args, **kwargs)
            return self._masked_result(result, mask, out)

        if function is np.array2string:
            # Complete hack.
            if self.shape == ():
                return str(self)

            kwargs.setdefault('formatter', {'all': self._to_string})

        return super().__array_function__(function, types, args, kwargs)

    def _masked_result(self, result, mask, out):
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            return tuple(self._masked_result(result_, mask, out_)
                         for (result_, out_) in zip(result, out))

        if out is None:
            return Masked(result, mask)

        assert isinstance(out, Masked)
        out._mask[...] = mask
        return out

    def _reduce_defaults(self, kwargs, initial_func=None):
        if 'where' not in kwargs:
            kwargs['where'] = ~self.mask
        if initial_func is not None and 'initial' not in kwargs:
            kwargs['initial'] = initial_func(self.unmasked)
        return kwargs

    def min(self, axis=None, out=None, **kwargs):
        return super().min(axis=axis, out=out,
                           **self._reduce_defaults(kwargs, np.max))

    def max(self, axis=None, out=None, **kwargs):
        return super().max(axis=axis, out=out,
                           **self._reduce_defaults(kwargs, np.min))

    def mean(self, axis=None, dtype=None, out=None, keepdims=False):
        # Cast bool, unsigned int, and int to float64 by default,
        # and do float16 at higher precision.
        is_float16_result = False
        if dtype is None:
            if issubclass(self.dtype.type, (np.integer, np.bool_)):
                dtype = np.dtype('f8')
            elif issubclass(self.dtype.type, np.float16):
                dtype = np.dtype('f4')
                is_float16_result = out is None

        result = self.sum(axis=axis, dtype=dtype, out=out,
                          keepdims=keepdims, where=~self.mask)
        n = np.add.reduce(~self.mask, axis=axis, keepdims=keepdims)
        result /= n
        if is_float16_result:
            result = result.astype(self.dtype)
        return result

    def var(self, axis=None, dtype=None, out=None, ddof=0, keepdims=False):
        n = np.add.reduce(~self.mask, axis=axis, keepdims=keepdims)[...]

        # Cast bool, unsigned int, and int to float64 by default.
        if dtype is None and issubclass(self.dtype.type,
                                        (np.integer, np.bool_)):
            dtype = np.dtype('f8')
        mean = self.mean(axis=axis, dtype=dtype, keepdims=True)

        x = self - mean
        x *= x.conjugate()  # Conjugate just returns x if not complex.

        result = x.sum(axis=axis, dtype=dtype, out=out,
                       keepdims=keepdims, where=~x.mask)
        n -= ddof
        n = np.maximum(n, 0, out=n)
        result /= n
        result._mask |= (n == 0)
        return result

    def std(self, axis=None, dtype=None, out=None, ddof=0, keepdims=False):
        result = self.var(axis=axis, dtype=dtype, out=out, ddof=ddof,
                          keepdims=keepdims)
        return np.sqrt(result, out=result)

    def __bool__(self):
        # First get result from array itself; this will error if not a scalar.
        result = super().__bool__()
        return result and not self.mask

    def any(self, axis=None, out=None, keepdims=False):
        return np.logical_or.reduce(self, axis=axis, out=out,
                                    keepdims=keepdims, where=~self.mask)

    def all(self, axis=None, out=None, keepdims=False):
        return np.logical_and.reduce(self, axis=axis, out=out,
                                     keepdims=keepdims, where=~self.mask)

    def _to_string(self, a):
        # This exists only to work around a numpy annoyance that array scalars
        # always get turned into plain ndarray items and thus loose the mask.
        if a.shape == () and self.shape == () and a == self.unmasked:
            a = self

        if a.shape == ():
            string = str(a.unmasked)
            if a.mask:
                # Strikethrough would be neat, but it doesn't show in konsole.
                # return ''.join(s+'\u0336' for s in string)
                return ' ' * (len(string)-3) + '\u2014' * 3
            else:
                return string

    # TODO: improve (greatly) str and repr!!
    def __str__(self):
        with np.printoptions(formatter={'all': self._to_string}):
            return super().__str__()

    def __repr__(self):
        with np.printoptions(formatter={'all': self._to_string}):
            return super().__repr__()


class MaskedIterator:
    """
    Flat iterator object to iterate over Masked Arrays.

    A `~astropy.utils.masked.MaskedIterator` iterator is returned by ``m.flat``
    for any masked array ``m``.  It allows iterating over the array as if it
    were a 1-D array, either in a for-loop or by calling its `next` method.

    Iteration is done in C-contiguous style, with the last index varying the
    fastest. The iterator can also be indexed using basic slicing or
    advanced indexing.

    Notes
    -----
    The design of `~astropy.utils.masked.MaskedIterator` follows that of
    `~numpy.ma.core.MaskedIterator`.  It is not exported by the
    `~astropy.utils.masked` module.  Instead of instantiating directly,
    use the ``flat`` method in the masked array instance.
    """

    def __init__(self, m):
        self._masked = m
        self._dataiter = m.unmasked.flat
        self._maskiter = m.mask.flat

    def __iter__(self):
        return self

    def __getitem__(self, indx):
        out = self._dataiter.__getitem__(indx)
        mask = self._maskiter.__getitem__(indx)
        # For single elements, ndarray.flat.__getitem__ returns scalars; these
        # need a new view as a Masked array.
        if not isinstance(out, np.ndarray):
            out = out[...]
            mask = mask[...]

        out = out.view(self._masked.__class__)
        out._mask = mask
        return out

    def __setitem__(self, index, value):
        data, mask = self._masked._data_mask(value)
        self._dataiter[index] = data
        self._maskiter[index] = mask

    def __next__(self):
        """
        Return the next value, or raise StopIteration.
        """
        out = next(self._dataiter)[...]
        mask = next(self._maskiter)[...]
        out = out.view(self._masked.__class__)
        out._mask = mask
        return out

    next = __next__


class MaskedNDArray(Masked, np.ndarray):
    _data_cls = np.ndarray

    @property
    def flat(self):
        """A 1-D iterator over the Quantity array.

        This returns a ``QuantityIterator`` instance, which behaves the same
        as the `~numpy.flatiter` instance returned by `~numpy.ndarray.flat`,
        and is similar to, but not a subclass of, Python's built-in iterator
        object.
        """
        return MaskedIterator(self)
