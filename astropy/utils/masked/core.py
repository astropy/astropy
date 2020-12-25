# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.

The design uses `Masked` as a factory class which automatically
generates new subclasses for any data class that is itself a
subclass of a predefined masked class, with `MaskedNDArray`
providing such a predefined class for `~numpy.ndarray`.

Generally, any new predefined class has to provide a
``from_unmasked(data, mask, copy=False)`` method as well
as an ``unmasked`` property or attribute that returns just
the data.  The `Masked` class itself provides a base ``mask``
property, which can be overridden if needed.
"""
import builtins
import functools

import numpy as np

from astropy.utils.shapes import NDArrayShapeMethods

from .function_helpers import (MASKED_SAFE_FUNCTIONS,
                               APPLY_TO_BOTH_FUNCTIONS,
                               DISPATCHED_FUNCTIONS)


__all__ = ['Masked', 'MaskedNDArray']


get__doc__ = functools.partial(format, """
Masked version of {0.__name__}.

Except for the ability to pass in a ``mask``, parameters are
as for `{0.__module__}.{0.__name__}`.
""")


class Masked(NDArrayShapeMethods):
    """A scalar value or array of values with associated mask.

    The resulting instance will take its exact type from whatever the
    contents are, with the type generated on the fly as needed.

    Parameters
    ----------
    data : array-like
        The data for which a mask is to be added.  The result will be a
        a subclass of the type of ``data``.
    mask : array-like of bool, optional
        The initial mask to assign.  If not given, taken from the data.
    copy : bool
        Whether the data and mask should be copied. Default: `False`.

    """

    _base_classes = {}
    """Explicitly defined masked classes keyed by their unmasked counterparts.

    For subclasses of these unmasked classes, masked counterparts can be generated.
    """

    _masked_classes = {}
    """Masked classes keyed by their unmasked data counterparts."""

    def __new__(cls, *args, **kwargs):
        if cls is Masked:
            # Initializing with Masked itself means we're in "factory mode".
            if not kwargs and len(args) == 1 and isinstance(args[0], type):
                # Create a new masked class.
                return cls._get_masked_cls(args[0])
            else:
                return cls._get_masked_instance(*args, **kwargs)
        else:
            # Otherwise we're a subclass and should just pass information on.
            return super().__new__(cls, *args, **kwargs)

    def __init_subclass__(cls, base_cls=None, data_cls=None, **kwargs):
        """Register a Masked subclass.

        Parameters
        ----------
        base_cls : type, optional
            If given, it is taken to mean that ``cls`` can be used as
            a base for masked versions of all subclasses of ``base_cls``,
            so it is registered as such in ``_base_classes``.
        data_cls : type, optional
            If given, ``cls`` should will be registered as the masked version of
            ``data_cls``.  Will set the private ``cls._data_cls`` attribute,
            and auto-generate a docstring if not present already.
        **kwargs
            Passed on for possible further initialization by superclasses.

        """
        if base_cls is not None:
            Masked._base_classes[base_cls] = cls

        if data_cls is not None:
            cls._data_cls = data_cls
            cls._masked_classes[data_cls] = cls
            if '__doc__' not in cls.__dict__:
                cls.__doc__ = get__doc__(data_cls)

        super().__init_subclass__(**kwargs)

    @classmethod
    def _get_masked_instance(cls, data, mask=None, copy=False):
        data, data_mask = cls._data_mask(data)
        if data is None:
            raise NotImplementedError("cannot initialize with np.ma.masked.")
        if mask is None:
            mask = False if data_mask is None else data_mask

        masked_cls = cls._get_masked_cls(data.__class__)
        return masked_cls.from_unmasked(data, mask, copy)

    @classmethod
    def _get_masked_cls(cls, data_cls):
        """Get the masked wrapper for a given data class.

        If the data class does not exist yet but is a subclass of any of the
        registered base data classes, it is automatically generated.
        """
        if issubclass(data_cls, Masked):
            return data_cls

        masked_cls = cls._masked_classes.get(data_cls)
        if masked_cls is None:
            # Walk through MRO and find closest base data class.
            # Note: right now, will basically always be ndarray, but
            # one could imagine needing some special care for one subclass,
            # which would then get its own entry.  E.g., if MaskedAngle
            # defined something special, then MaskedLongitude should depend
            # on it.
            for mro_item in data_cls.__mro__:
                base_cls = cls._base_classes.get(mro_item)
                if base_cls is not None:
                    break
            else:
                # Just hope that MaskedNDArray can handle it.
                # TODO: this covers the case where a user puts in a list or so,
                # but for those one could just explicitly do something like
                # _masked_classes[list] = MaskedNDArray.
                return MaskedNDArray

            # Create (and therefore register) new Masked subclass for the
            # given data_cls.
            masked_cls = type('Masked' + data_cls.__name__,
                              (data_cls, base_cls), {}, data_cls=data_cls)

        return masked_cls

    @classmethod
    def _data_mask(cls, data):
        """Split data into unmasked and mask, if present.

        Parameters
        ----------
        data : array-like
            Possibly masked item, judged by whether it has a ``mask`` attribute.
            If so, checks for being an instance of `~astropy.utils.masked.Masked`
            or `~numpy.ma.MaskedArray`, and gets unmasked data appropriately.

        Returns
        -------
        unmasked, mask : array-like
            Unmasked will be `None` for `~numpy.ma.masked`.

        """
        mask = getattr(data, 'mask', None)
        if mask is not None:
            if isinstance(data, Masked):
                data = data.unmasked
            elif data is np.ma.masked:
                data = None
            elif isinstance(data, np.ma.MaskedArray):
                data = data.data

        return data, mask

    def _get_mask(self):
        """The mask.

        If set, replace the original mask, with whatever it is set with,
        using a view if no broadcasting or type conversion is required.
        """
        return self._mask

    def _set_mask(self, mask, copy=False):
        mask_dtype = (np.ma.make_mask_descr(self.dtype)
                      if self.dtype.names else np.dtype('?'))
        ma = np.asanyarray(mask, dtype=mask_dtype)
        if ma.shape != self.shape:
            # This will fail (correctly) if not broadcastable.
            self._mask = np.empty(self.shape, dtype=mask_dtype)
            self._mask[...] = ma
        elif ma is mask:
            # Even if not copying use a view so that shape setting
            # does not propagate.
            self._mask = mask.copy() if copy else mask.view()
        else:
            self._mask = ma

    mask = property(_get_mask, _set_mask)

    # Note: any subclass needs to define an unmasked attribute or property.
    def unmask(self, fill_value=None):
        """Get the underlying data, possibly filling masked values."""
        unmasked = self.unmasked
        if fill_value is None:
            return unmasked
        else:
            unmasked = unmasked.copy()
            if self.dtype.names:
                np.ma.core._recursive_filled(unmasked, self.mask, fill_value)
            else:
                unmasked[self.mask] = fill_value
            return unmasked

    def _apply(self, method, *args, **kwargs):
        # Required method for NDArrayShapeMethods, to help provide __getitem__
        # and shape-changing methods.
        if callable(method):
            data = method(self.unmasked, *args, **kwargs)
            mask = method(self.mask, *args, **kwargs)
        else:
            data = getattr(self.unmasked, method)(*args, **kwargs)
            mask = getattr(self.mask, method)(*args, **kwargs)

        return self.from_unmasked(data, mask, copy=False)

    def __setitem__(self, item, value):
        value, mask = self._data_mask(value)
        if value is not None:
            self.unmasked[item] = value
        self.mask[item] = mask


def _comparison_method(op):
    """
    Create a comparison operator for MaskedNDArray.

    Needed since for string dtypes the base operators bypass __array_ufunc__
    and hence return unmasked results.
    """
    def _compare(self, other):
        other_data, other_mask = self._data_mask(other)
        result = getattr(self.unmasked, op)(other_data)
        if result is NotImplemented:
            return NotImplemented
        mask = self.mask | (other_mask if other_mask is not None else False)
        return self._masked_result(result, mask, None)

    return _compare


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

        return self._masked.from_unmasked(out, mask, copy=False)

    def __setitem__(self, index, value):
        data, mask = self._masked._data_mask(value)
        if data is not None:
            self._dataiter[index] = data
        self._maskiter[index] = mask

    def __next__(self):
        """
        Return the next value, or raise StopIteration.
        """
        out = next(self._dataiter)[...]
        mask = next(self._maskiter)[...]
        return self._masked.from_unmasked(out, mask, copy=False)

    next = __next__


class MaskedNDArray(Masked, np.ndarray, base_cls=np.ndarray, data_cls=np.ndarray):
    _mask = None

    def __new__(cls, *args, mask=False, **kwargs):
        """Get data class instance from arguments and then set mask."""
        self = super().__new__(cls, *args, **kwargs)
        self.mask = mask
        return self

    def __init_subclass__(cls, **kwargs):
        # For all subclasses we should set a default __new__ that passes on
        # arguments other than mask to the data class, and then sets the mask.
        if '__new__' not in cls.__dict__:
            def __new__(newcls, *args, mask=False, **kwargs):
                """Get data class instance from arguments and then set mask."""
                # Need to explicitly mention classes outside of class definition.
                self = super(cls, newcls).__new__(newcls, *args, **kwargs)
                self.mask = mask
                return self
            cls.__new__ = __new__

        super().__init_subclass__(cls, **kwargs)

    # The two required pieces.
    @classmethod
    def from_unmasked(cls, data, mask=None, copy=False):
        """Create an instance from unmasked data and a mask."""
        data = np.array(data, subok=True, copy=copy)
        self = data.view(cls)
        self._set_mask(mask, copy=copy)
        return self

    @property
    def unmasked(self):
        """The unmasked values."""
        return super().view(self._data_cls)

    @classmethod
    def _get_masked_cls(cls, data_cls):
        # Short-cuts
        if data_cls is np.ndarray:
            return MaskedNDArray
        elif data_cls is None:  # for .view()
            return cls

        if not issubclass(data_cls, np.ndarray):
            raise ValueError('can only pass in an ndarray subtype.')

        return super()._get_masked_cls(data_cls)

    @property
    def flat(self):
        """A 1-D iterator over the Masked array.

        This returns a ``MaskedIterator`` instance, which behaves the same
        as the `~numpy.flatiter` instance returned by `~numpy.ndarray.flat`,
        and is similar to Python's built-in iterator, except that it also
        allows assignment.
        """
        return MaskedIterator(self)

    def view(self, dtype=None, type=None):
        """New view of the masked array.

        Like `numpy.ndarray.view`, but always returning a masked array subclass.
        """
        if type is None and (isinstance(dtype, builtins.type) and
                             issubclass(dtype, np.ndarray)):
            type = dtype
            dtype = None
        elif dtype is not None:
            raise NotImplementedError('{} cannot be viewed with new dtype.'
                                      .format(self.__class__))
        return super().view(self._get_masked_cls(type))

    def __array_finalize__(self, obj):
        if obj is None:
            return None

        # Logically, this should come from ndarray and hence be None, but
        # just in case someone creates a new mixin, we check.
        super_array_finalize = super().__array_finalize__
        if super_array_finalize:
            super_array_finalize(obj)

        if self._mask is None:
            # Got here after, e.g., a view of another masked class.
            # Get its mask, or initialize ours.
            self._set_mask(getattr(obj, '_mask', False))

    _eq_simple = _comparison_method('__eq__')
    _ne_simple = _comparison_method('__ne__')
    __lt__ = _comparison_method('__lt__')
    __le__ = _comparison_method('__le__')
    __gt__ = _comparison_method('__gt__')
    __ge__ = _comparison_method('__ge__')

    def __eq__(self, other):
        if not self.dtype.names:
            return self._eq_simple(other)

        # For structured arrays, we treat this as a reduction over the fields,
        # where masked fields are skipped and thus do not influence the result.
        other = np.asanyarray(other, dtype=self.dtype)
        result = np.stack([self[field] == other[field]
                           for field in self.dtype.names], axis=-1)
        return result.all(axis=-1)

    def __ne__(self, other):
        if not self.dtype.names:
            return self._ne_simple(other)

        # For structured arrays, we treat this as a reduction over the fields,
        # where masked fields are skipped and thus do not influence the result.
        other = np.asanyarray(other, dtype=self.dtype)
        result = np.stack([self[field] != other[field]
                           for field in self.dtype.names], axis=-1)
        return result.any(axis=-1)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.pop('out', None)
        if out is not None:
            converted_out = []
            for out_ in out:
                if out_ is None:
                    converted_out.append(None)
                elif isinstance(out_, Masked):
                    converted_out.append(out_.unmasked)
                else:
                    return NotImplemented
            kwargs['out'] = tuple(converted_out)
            if ufunc.nout == 1:
                out = out[0]

        if ufunc.signature:
            # We're dealing with a gufunc. For now, only deal with
            # np.matmul and gufunc with a single axis.  Also ignore
            # axes keyword for now...  TODO: generally, the mask can
            # be generated purely based on the signature.
            assert 'axes' not in kwargs
            converted = []
            masks = []
            for input_ in inputs:
                if isinstance(input_, Masked):
                    masks.append(input_.mask)
                    converted.append(input_.unmasked)
                else:
                    masks.append(np.zeros(getattr(input_, 'shape', ()), bool))
                    converted.append(input_)

            result = ufunc(*converted, **kwargs)
            if ufunc is np.matmul:
                # np.matmul is tricky and its signature cannot be parsed by
                # _parse_gufunc_signature.
                mask0, mask1 = masks
                if mask1.ndim > 1:
                    mask0 = np.logical_or.reduce(mask0, axis=-1, keepdims=True)
                    mask1 = np.logical_or.reduce(mask1, axis=-2, keepdims=True)
                else:
                    mask0 = np.logical_or.reduce(mask0, axis=-1)
                    mask1 = np.logical_or.reduce(mask1)

                mask = np.logical_or(mask0, mask1)

            else:
                # Parse signature with private numpy function. Note it
                # cannot handle spaces in tuples, so remove those.
                in_sig, out_sig = np.lib.function_base._parse_gufunc_signature(
                    ufunc.signature.replace(' ', ''))
                axis = kwargs.get('axis', None)
                keepdims = kwargs.get('keepdims', False)
                in_masks = []
                for sig, mask in zip(in_sig, masks):
                    if sig:
                        # Input has core dimensions.
                        if axis is None:
                            mask = np.logical_or.reduce(
                                mask, axis=tuple(range(-len(sig), 0)))
                        else:
                            mask = np.logical_or.reduce(
                                mask, axis=axis, keepdims=keepdims)
                    in_masks.append(mask)

                mask = functools.reduce(np.logical_or, in_masks)
                out_masks = []
                for os in out_sig:
                    if os:
                        # Output has core dimensions.
                        if axis is None:
                            out_mask = mask[(Ellipsis,)+(np.newaxis,)*len(os)]
                        else:
                            out_mask = np.expand_dims(mask, axis)
                    else:
                        out_mask = mask
                    out_masks.append(out_mask)

                mask = out_masks if len(out_masks) > 1 else out_masks[0]

        elif method == '__call__':
            # Regular ufunc call.
            converted = []
            masks = []
            for input_ in inputs:
                if isinstance(input_, Masked):
                    masks.append(input_.mask)
                    converted.append(input_.unmasked)
                else:
                    converted.append(input_)

            result = ufunc(*converted, **kwargs)
            if masks:
                mask = functools.reduce(np.logical_or, masks)
            else:
                mask = False

        elif method == 'reduce':
            # Reductions like np.add.reduce (sum).
            if isinstance(inputs[0], Masked):
                # By default, we simply propagate masks, since for
                # things like np.sum, it makes no sense to do otherwise.
                # Individual methods need to override as needed.
                # TODO: take care of 'out' too?
                axis = kwargs.get('axis', None)
                keepdims = kwargs.get('keepdims', False)
                where = kwargs.get('where', True)
                mask = np.logical_or.reduce(inputs[0].mask, where=where,
                                            axis=axis, keepdims=keepdims)
                if where is not True:
                    # Mask also whole rows that were not selected by where,
                    # so would have been left as unmasked above.
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
            # TODO: implement things like np.add.accumulate (used for cumsum).
            return NotImplemented

        if mask is False or result is None or result is NotImplemented:
            return result

        return self._masked_result(result, mask, out)

    def __array_function__(self, function, types, args, kwargs):
        # TODO: go through functions systematically to see which ones
        # work and/or can be supported.
        if function in MASKED_SAFE_FUNCTIONS:
            return super().__array_function__(function, types, args, kwargs)

        elif function in APPLY_TO_BOTH_FUNCTIONS:
            helper = APPLY_TO_BOTH_FUNCTIONS[function]
            data_args, mask_args, kwargs, out = helper(*args, **kwargs)
            if mask_args is not None:
                mask = function(*mask_args, **kwargs)
            if out is not None:
                if isinstance(out, Masked):
                    if mask is None:
                        return NotImplemented
                    kwargs['out'] = out.unmasked
                elif mask is not None:
                    return NotImplemented
                kwargs['out'] = out
            result = function(*data_args, **kwargs)

        elif function in DISPATCHED_FUNCTIONS:
            dispatched_function = DISPATCHED_FUNCTIONS[function]
            result, mask, out = dispatched_function(*args, **kwargs)

        else:
            # By default, just pass it through for now.
            if function is np.array2string:
                # Complete hack.
                if self.shape == ():
                    return str(self)

                kwargs.setdefault('formatter', {'all': self._to_string})

            return super().__array_function__(function, types, args, kwargs)

        if mask is None:
            return result
        else:
            return self._masked_result(result, mask, out)

    def _masked_result(self, result, mask, out):
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            if not isinstance(mask, (list, tuple)):
                mask = [mask] * len(result)
            return tuple(self._masked_result(result_, mask_, out_)
                         for (result_, mask_, out_) in zip(result, mask, out))

        if out is None:
            # Note that we cannot count on result being the same class as
            # 'self' (e.g., comparison of quantity results in an ndarray, most
            # operations on Longitude and Latitude result in Angle or
            # Quantity), so use Masked to determine the appropriate class.
            return Masked(result, mask)

        # TODO: remove this sanity check once test cases are more complete.
        assert isinstance(out, Masked)
        # If we have an output, the result was writtin in-place, so we should
        # also write the mask in-place.
        out._mask[...] = mask
        return out

    # Below are ndarray methods that need to be overridden as masked elements
    # need to be skipped and/or an initial value needs to be set.
    def _reduce_defaults(self, kwargs, initial_func=None):
        """Get default where and initial for masked reductions.

        Generally, the default should be to skip all masked elements.  For
        reductions such as np.minimum.reduce, we also need an initial value,
        which can be determined using ``initial_func``.

        """
        if 'where' not in kwargs:
            kwargs['where'] = ~self.mask
        if initial_func is not None and 'initial' not in kwargs:
            kwargs['initial'] = initial_func(self.unmasked)
        return kwargs

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        # Unfortunately, cannot override the call to diagonal inside trace, so
        # duplicate implementation in numpy/core/src/multiarray/calculation.c.
        diagonal = self.diagonal(offset=offset, axis1=axis1, axis2=axis2)
        return diagonal.sum(-1, dtype=dtype, out=out)

    def min(self, axis=None, out=None, **kwargs):
        return super().min(axis=axis, out=out,
                           **self._reduce_defaults(kwargs, np.max))

    def max(self, axis=None, out=None, **kwargs):
        return super().max(axis=axis, out=out,
                           **self._reduce_defaults(kwargs, np.min))

    def argmin(self, axis=None, fill_value=None, out=None):
        if fill_value is None:
            fill_value = np.max(self.unmasked)
        return self.unmask(fill_value=fill_value).argmin(axis=axis, out=out)

    def argmax(self, axis=None, fill_value=None, out=None):
        if fill_value is None:
            fill_value = np.min(self.unmasked)
        return self.unmask(fill_value=fill_value).argmax(axis=axis, out=out)

    def argsort(self, axis=-1, kind=None, order=None):
        """Returns the indices that would sort an array.

        Perform an indirect sort along the given axis on both the array
        and the mask, with masked items being sorted to the end.

        Parameters
        ----------
        axis : int or None, optional
            Axis along which to sort.  The default is -1 (the last axis).
            If None, the flattened array is used.
        kind : str or None, ignored.
            The kind of sort.  Present only to allow subclasses to work.
        order : str or list of str, not yet implemented
            Way to sort structured arrays.

        Returns
        -------
        index_array : ndarray, int
            Array of indices that sorts along the specified ``axis``.  Use
            ``np.take_along_axis(self, index_array, axis=axis)`` to obtain
            the sorted array.

        """
        if order is not None:
            raise NotImplementedError("structured arrays cannot yet "
                                      "be sorted.")
        data, mask = self.unmasked, self.mask
        if axis is None:
            data = data.ravel()
            mask = mask.ravel()
            axis = -1
        return np.lexsort((data, mask), axis=axis)

    def sort(self, axis=-1, kind=None, order=None):
        # TODO: probably possible to do this faster than going through argsort!
        indices = self.argsort(axis, kind=kind, order=order)
        self[:] = np.take_along_axis(self, indices, axis=axis)

    def mean(self, axis=None, dtype=None, out=None, keepdims=False):
        # Implementation based on that in numpy/core/_methods.py
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
        # Simplified implementation based on that in numpy/core/_methods.py
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

    # Following overrides needed since somehow the ndarray implementation
    # does not actually call these.
    def __str__(self):
        return np.array_str(self)

    def __repr__(self):
        return np.array_repr(self)
