# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Built-in mask mixin class.

The design uses `Masked` as a factory class which automatically
generates new subclasses for any data class that is itself a
subclass of a predefined masked class, with `MaskedNDArray`
providing such a predefined class for `~numpy.ndarray`.

Generally, any new predefined class should override the
``from_unmasked(data, mask, copy=False)`` class method that
creates an instance from unmasked data and a mask, as well as
the ``unmasked`` property that returns just the data.
The `Masked` class itself provides a base ``mask`` property,
which can also be overridden if needed.

"""
import builtins

import numpy as np

from astropy.utils.compat import NUMPY_LT_2_0
from astropy.utils.data_info import ParentDtypeInfo
from astropy.utils.shapes import NDArrayShapeMethods

from .function_helpers import (
    APPLY_TO_BOTH_FUNCTIONS,
    DISPATCHED_FUNCTIONS,
    MASKED_SAFE_FUNCTIONS,
    UNSUPPORTED_FUNCTIONS,
)

__all__ = ["Masked", "MaskedNDArray"]


get__doc__ = """Masked version of {0.__name__}.

Except for the ability to pass in a ``mask``, parameters are
as for `{0.__module__}.{0.__name__}`.
""".format


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
            if cls.__doc__ is None:
                cls.__doc__ = get__doc__(data_cls)

        super().__init_subclass__(**kwargs)

    # This base implementation just uses the class initializer.
    # Subclasses can override this in case the class does not work
    # with this signature, or to provide a faster implementation.
    @classmethod
    def from_unmasked(cls, data, mask=None, copy=False):
        """Create an instance from unmasked data and a mask."""
        return cls(data, mask=mask, copy=copy)

    @classmethod
    def _get_masked_instance(cls, data, mask=None, copy=False):
        data, data_mask = cls._get_data_and_mask(data)
        if mask is None:
            mask = False if data_mask is None else data_mask

        masked_cls = cls._get_masked_cls(data.__class__)
        return masked_cls.from_unmasked(data, mask, copy)

    @classmethod
    def _get_masked_cls(cls, data_cls):
        """Get the masked wrapper for a given data class.

        If the data class does not exist yet but is a subclass of any of the
        registered base data classes, it is automatically generated
        (except we skip `~numpy.ma.MaskedArray` subclasses, since then the
        masking mechanisms would interfere).
        """
        if issubclass(data_cls, (Masked, np.ma.MaskedArray)):
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
            masked_cls = type(
                "Masked" + data_cls.__name__,
                (data_cls, base_cls),
                {},
                data_cls=data_cls,
            )

        return masked_cls

    @classmethod
    def _get_data_and_mask(cls, data, allow_ma_masked=False):
        """Split data into unmasked and mask, if present.

        Parameters
        ----------
        data : array-like
            Possibly masked item, judged by whether it has a ``mask`` attribute.
            If so, checks for being an instance of `~astropy.utils.masked.Masked`
            or `~numpy.ma.MaskedArray`, and gets unmasked data appropriately.
        allow_ma_masked : bool, optional
            Whether or not to process `~numpy.ma.masked`, i.e., an item that
            implies no data but the presence of a mask.

        Returns
        -------
        unmasked, mask : array-like
            Unmasked will be `None` for `~numpy.ma.masked`.

        Raises
        ------
        ValueError
            If `~numpy.ma.masked` is passed in and ``allow_ma_masked`` is not set.

        """
        mask = getattr(data, "mask", None)
        if mask is not None:
            try:
                data = data.unmasked
            except AttributeError:
                if not isinstance(data, np.ma.MaskedArray):
                    raise
                if data is np.ma.masked:
                    if allow_ma_masked:
                        data = None
                    else:
                        raise ValueError("cannot handle np.ma.masked here.") from None
                else:
                    data = data.data

        return data, mask

    @classmethod
    def _get_data_and_masks(cls, *args):
        data_masks = [cls._get_data_and_mask(arg) for arg in args]
        return (
            tuple(data for data, _ in data_masks),
            tuple(mask for _, mask in data_masks),
        )

    def _get_mask(self):
        """The mask.

        If set, replace the original mask, with whatever it is set with,
        using a view if no broadcasting or type conversion is required.
        """
        return self._mask

    def _set_mask(self, mask, copy=False):
        self_dtype = getattr(self, "dtype", None)
        mask_dtype = (
            np.ma.make_mask_descr(self_dtype)
            if self_dtype and self_dtype.names
            else np.dtype("?")
        )
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

    # Note: subclass should generally override the unmasked property.
    # This one assumes the unmasked data is stored in a private attribute.
    @property
    def unmasked(self):
        """The unmasked values.

        See Also
        --------
        astropy.utils.masked.Masked.filled
        """
        return self._unmasked

    def filled(self, fill_value):
        """Get a copy of the underlying data, with masked values filled in.

        Parameters
        ----------
        fill_value : object
            Value to replace masked values with.

        See Also
        --------
        astropy.utils.masked.Masked.unmasked
        """
        unmasked = self.unmasked.copy()
        if self.mask.dtype.names:
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

        result = self.from_unmasked(data, mask, copy=False)
        if "info" in self.__dict__:
            result.info = self.info

        return result

    def __setitem__(self, item, value):
        value, mask = self._get_data_and_mask(value, allow_ma_masked=True)
        if value is not None:
            self.unmasked[item] = value
        self.mask[item] = mask


class MaskedInfoBase:
    mask_val = np.ma.masked

    def __init__(self, bound=False):
        super().__init__(bound)

        # If bound to a data object instance then create the dict of attributes
        # which stores the info attribute values.
        if bound:
            # Specify how to serialize this object depending on context.
            self.serialize_method = {
                "fits": "null_value",
                "ecsv": "null_value",
                "hdf5": "data_mask",
                "parquet": "data_mask",
                None: "null_value",
            }


class MaskedNDArrayInfo(MaskedInfoBase, ParentDtypeInfo):
    """
    Container for meta information like name, description, format.
    """

    # Add `serialize_method` attribute to the attrs that MaskedNDArrayInfo knows
    # about.  This allows customization of the way that MaskedColumn objects
    # get written to file depending on format.  The default is to use whatever
    # the writer would normally do, which in the case of FITS or ECSV is to use
    # a NULL value within the data itself.  If serialize_method is 'data_mask'
    # then the mask is explicitly written out as a separate column if there
    # are any masked values.  This is the same as for MaskedColumn.
    attr_names = ParentDtypeInfo.attr_names | {"serialize_method"}

    # When `serialize_method` is 'data_mask', and data and mask are being written
    # as separate columns, use column names <name> and <name>.mask (instead
    # of default encoding as <name>.data and <name>.mask).
    _represent_as_dict_primary_data = "data"

    def _represent_as_dict(self):
        out = super()._represent_as_dict()

        masked_array = self._parent

        # If the serialize method for this context (e.g. 'fits' or 'ecsv') is
        # 'data_mask', that means to serialize using an explicit mask column.
        method = self.serialize_method[self._serialize_context]

        if method == "data_mask":
            out["data"] = masked_array.unmasked

            if np.any(masked_array.mask):
                # Only if there are actually masked elements do we add the ``mask`` column
                out["mask"] = masked_array.mask

        elif method == "null_value":
            out["data"] = np.ma.MaskedArray(
                masked_array.unmasked, mask=masked_array.mask
            )

        else:
            raise ValueError(
                'serialize method must be either "data_mask" or "null_value"'
            )

        return out

    def _construct_from_dict(self, map):
        # Override usual handling, since MaskedNDArray takes shape and buffer
        # as input, which is less useful here.
        # The map can contain either a MaskedColumn or a Column and a mask.
        # Extract the mask for the former case.
        map.setdefault("mask", getattr(map["data"], "mask", False))
        return self._parent_cls.from_unmasked(**map)


class MaskedArraySubclassInfo(MaskedInfoBase):
    """Mixin class to create a subclasses such as MaskedQuantityInfo."""

    # This is used below in __init_subclass__, which also inserts a
    # 'serialize_method' attribute in attr_names.

    def _represent_as_dict(self):
        # Use the data_cls as the class name for serialization,
        # so that we do not have to store all possible masked classes
        # in astropy.table.serialize.__construct_mixin_classes.
        out = super()._represent_as_dict()
        data_cls = self._parent._data_cls
        out.setdefault("__class__", data_cls.__module__ + "." + data_cls.__name__)
        return out


def _comparison_method(op):
    """
    Create a comparison operator for MaskedNDArray.

    Needed since for string dtypes the base operators bypass __array_ufunc__
    and hence return unmasked results.
    """

    def _compare(self, other):
        other_data, other_mask = self._get_data_and_mask(other)
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
        data, mask = self._masked._get_data_and_mask(value, allow_ma_masked=True)
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

    info = MaskedNDArrayInfo()

    def __new__(cls, *args, mask=None, **kwargs):
        """Get data class instance from arguments and then set mask."""
        self = super().__new__(cls, *args, **kwargs)
        if mask is not None:
            self.mask = mask
        elif self._mask is None:
            self.mask = False
        return self

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(cls, **kwargs)
        # For all subclasses we should set a default __new__ that passes on
        # arguments other than mask to the data class, and then sets the mask.
        if "__new__" not in cls.__dict__:

            def __new__(newcls, *args, mask=None, **kwargs):
                """Get data class instance from arguments and then set mask."""
                # Need to explicitly mention classes outside of class definition.
                self = super(cls, newcls).__new__(newcls, *args, **kwargs)
                if mask is not None:
                    self.mask = mask
                elif self._mask is None:
                    self.mask = False
                return self

            cls.__new__ = __new__

        if "info" not in cls.__dict__ and hasattr(cls._data_cls, "info"):
            data_info = cls._data_cls.info
            attr_names = data_info.attr_names | {"serialize_method"}
            new_info = type(
                cls.__name__ + "Info",
                (MaskedArraySubclassInfo, data_info.__class__),
                dict(attr_names=attr_names),
            )
            cls.info = new_info()

    # The two pieces typically overridden.
    @classmethod
    def from_unmasked(cls, data, mask=None, copy=False):
        # Note: have to override since __new__ would use ndarray.__new__
        # which expects the shape as its first argument, not an array.

        if not NUMPY_LT_2_0 and copy is False:
            # strict backward compatibility
            # see https://github.com/astropy/astropy/pull/16142
            copy = None

        data = np.array(data, subok=True, copy=copy)
        self = data.view(cls)
        self._set_mask(mask, copy=copy)
        return self

    @property
    def unmasked(self):
        return super().view(self._data_cls)

    @classmethod
    def _get_masked_cls(cls, data_cls):
        # Short-cuts
        if data_cls is np.ndarray:
            return MaskedNDArray
        elif data_cls is None:  # for .view()
            return cls

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

    @property
    def _baseclass(self):
        """Work-around for MaskedArray initialization.

        Allows the base class to be inferred correctly when a masked instance
        is used to initialize (or viewed as) a `~numpy.ma.MaskedArray`.

        """
        return self._data_cls

    def view(self, dtype=None, type=None):
        """New view of the masked array.

        Like `numpy.ndarray.view`, but always returning a masked array subclass.
        """
        if type is None and (
            isinstance(dtype, builtins.type) and issubclass(dtype, np.ndarray)
        ):
            return super().view(self._get_masked_cls(dtype))

        if dtype is None:
            return super().view(self._get_masked_cls(type))

        dtype = np.dtype(dtype)
        if not (
            dtype.itemsize == self.dtype.itemsize
            and (dtype.names is None or len(dtype.names) == len(self.dtype.names))
        ):
            raise NotImplementedError(
                f"{self.__class__} cannot be viewed with a dtype with a "
                "with a different number of fields or size."
            )

        return super().view(dtype, self._get_masked_cls(type))

    def __array_finalize__(self, obj):
        # If we're a new object or viewing an ndarray, nothing has to be done.
        if obj is None or obj.__class__ is np.ndarray:
            return

        # Logically, this should come from ndarray and hence be None, but
        # just in case someone creates a new mixin, we check.
        super_array_finalize = super().__array_finalize__
        if super_array_finalize:  # pragma: no cover
            super_array_finalize(obj)

        if self._mask is None:
            # Got here after, e.g., a view of another masked class.
            # Get its mask, or initialize ours.
            self._set_mask(getattr(obj, "_mask", False))

        if "info" in obj.__dict__:
            self.info = obj.info

    @property
    def shape(self):
        """The shape of the data and the mask.

        Usually used to get the current shape of an array, but may also be
        used to reshape the array in-place by assigning a tuple of array
        dimensions to it.  As with `numpy.reshape`, one of the new shape
        dimensions can be -1, in which case its value is inferred from the
        size of the array and the remaining dimensions.

        Raises
        ------
        AttributeError
            If a copy is required, of either the data or the mask.

        """
        # Redefinition to allow defining a setter and add a docstring.
        return super().shape

    @shape.setter
    def shape(self, shape):
        old_shape = self.shape
        self._mask.shape = shape
        # Reshape array proper in try/except just in case some broadcasting
        # or so causes it to fail.
        try:
            super(MaskedNDArray, type(self)).shape.__set__(self, shape)
        except Exception as exc:
            self._mask.shape = old_shape
            # Given that the mask reshaping succeeded, the only logical
            # reason for an exception is something like a broadcast error in
            # in __array_finalize__, or a different memory ordering between
            # mask and data.  For those, give a more useful error message;
            # otherwise just raise the error.
            if "could not broadcast" in exc.args[0]:
                raise AttributeError(
                    "Incompatible shape for in-place modification. "
                    "Use `.reshape()` to make a copy with the desired "
                    "shape."
                ) from None
            else:  # pragma: no cover
                raise

    _eq_simple = _comparison_method("__eq__")
    _ne_simple = _comparison_method("__ne__")
    __lt__ = _comparison_method("__lt__")
    __le__ = _comparison_method("__le__")
    __gt__ = _comparison_method("__gt__")
    __ge__ = _comparison_method("__ge__")

    def __eq__(self, other):
        if not self.dtype.names:
            return self._eq_simple(other)

        # For structured arrays, we treat this as a reduction over the fields,
        # where masked fields are skipped and thus do not influence the result.
        other = np.asanyarray(other, dtype=self.dtype)
        result = np.stack(
            [self[field] == other[field] for field in self.dtype.names], axis=-1
        )
        return result.all(axis=-1)

    def __ne__(self, other):
        if not self.dtype.names:
            return self._ne_simple(other)

        # For structured arrays, we treat this as a reduction over the fields,
        # where masked fields are skipped and thus do not influence the result.
        other = np.asanyarray(other, dtype=self.dtype)
        result = np.stack(
            [self[field] != other[field] for field in self.dtype.names], axis=-1
        )
        return result.any(axis=-1)

    def _combine_masks(self, masks, out=None, where=True, copy=True):
        """Combine masks, possibly storing it in some output.

        Parameters
        ----------
        masks : tuple of array of bool or None
            Input masks.  Any that are `None` or `False` are ignored.
            Should broadcast to each other.
        out : output mask array, optional
            Possible output array to hold the result.
        where : array of bool, optional
            Which elements of the output array to fill.
        copy : bool optional
            Whether to ensure a copy is made. Only relevant if a single
            input mask is not `None`, and ``out`` is not given.
        """
        masks = [m for m in masks if m is not None and m is not False]
        if not masks:
            return False
        if len(masks) == 1:
            if out is None:
                return masks[0].copy() if copy else masks[0]
            else:
                np.copyto(out, masks[0], where=where)
                return out

        # [...] at the end to ensure we have an array, not a scalar, and
        # thus can be used for in-place changes in the loop.
        out = np.logical_or(masks[0], masks[1], out=out, where=where)[...]
        for mask in masks[2:]:
            np.logical_or(out, mask, out=out, where=where)
        return out

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        out = kwargs.pop("out", None)
        out_unmasked = None
        out_mask = None
        if out is not None:
            out_unmasked, out_masks = self._get_data_and_masks(*out)
            for d, m in zip(out_unmasked, out_masks):
                if m is None:
                    # TODO: allow writing to unmasked output if nothing is masked?
                    if d is not None:
                        raise TypeError("cannot write to unmasked output")
                elif out_mask is None:
                    out_mask = m

        # TODO: where is only needed for __call__ and reduce;
        # this is very fast, but still worth separating out?
        where = kwargs.pop("where", True)
        if where is True:
            where_unmasked = True
            where_mask = None
        else:
            where_unmasked, where_mask = self._get_data_and_mask(where)

        unmasked, masks = self._get_data_and_masks(*inputs)

        if ufunc.signature:
            # We're dealing with a gufunc. For now, only deal with
            # np.matmul and gufuncs for which the mask of any output always
            # depends on all core dimension values of all inputs.
            # Also ignore axes keyword for now...
            # TODO: in principle, it should be possible to generate the mask
            # purely based on the signature.
            if "axes" in kwargs:
                raise NotImplementedError(
                    "Masked does not yet support gufunc calls with 'axes'."
                )
            if ufunc is np.matmul:
                # np.matmul is tricky and its signature cannot be parsed by
                # _parse_gufunc_signature.
                unmasked = np.atleast_1d(*unmasked)
                mask0, mask1 = masks
                masks = []
                is_mat1 = unmasked[1].ndim >= 2
                if mask0 is not None:
                    masks.append(np.logical_or.reduce(mask0, axis=-1, keepdims=is_mat1))

                if mask1 is not None:
                    masks.append(
                        np.logical_or.reduce(mask1, axis=-2, keepdims=True)
                        if is_mat1
                        else np.logical_or.reduce(mask1)
                    )

                mask = self._combine_masks(masks, out=out_mask, copy=False)

            else:
                # Parse signature with private numpy function. Note it
                # cannot handle spaces in tuples, so remove those.
                if NUMPY_LT_2_0:
                    in_sig, out_sig = np.lib.function_base._parse_gufunc_signature(
                        ufunc.signature.replace(" ", "")
                    )
                else:
                    (
                        in_sig,
                        out_sig,
                    ) = np.lib._function_base_impl._parse_gufunc_signature(
                        ufunc.signature.replace(" ", "")
                    )
                axis = kwargs.get("axis")
                keepdims = kwargs.get("keepdims", False)
                in_masks = []
                for sig, mask in zip(in_sig, masks):
                    if mask is not None:
                        if sig:
                            # Input has core dimensions.  Assume that if any
                            # value in those is masked, the output will be
                            # masked too (TODO: for multiple core dimensions
                            # this may be too strong).
                            in_axis = (
                                tuple(range(-1, -1 - len(sig), -1))
                                if axis is None
                                else axis
                            )
                            mask = np.logical_or.reduce(
                                mask, axis=in_axis, keepdims=keepdims
                            )
                        in_masks.append(mask)

                mask = self._combine_masks(in_masks)
                result_masks = []
                for os in out_sig:
                    if os:
                        # Output has core dimensions.  Assume all those
                        # get the same mask.
                        out_axis = (
                            tuple(range(-1, -1 - len(os), -1)) if axis is None else axis
                        )
                        result_mask = np.expand_dims(mask, out_axis)
                    else:
                        result_mask = mask
                    result_masks.append(result_mask)

                mask = result_masks if len(result_masks) > 1 else result_masks[0]

        elif method == "__call__":
            # Regular ufunc call.
            # Combine the masks from the input, possibly selecting elements.
            mask = self._combine_masks(masks, out=out_mask, where=where_unmasked)
            # If relevant, also mask output elements for which where was masked.
            if where_mask is not None:
                mask |= where_mask

        elif method == "outer":
            # Must have two arguments; adjust masks as will be done for data.
            m0, m1 = masks
            if m0 is not None and m0.ndim > 0:
                m0 = m0[(...,) + (np.newaxis,) * np.ndim(unmasked[1])]
            mask = self._combine_masks((m0, m1), out=out_mask)

        elif method in {"reduce", "accumulate"}:
            # Reductions like np.add.reduce (sum).
            # Treat any masked where as if the input element was masked.
            mask = self._combine_masks((masks[0], where_mask), copy=False)
            if mask is not False:
                # By default, we simply propagate masks, since for
                # things like np.sum, it makes no sense to do otherwise.
                # Individual methods need to override as needed.
                if method == "reduce":
                    axis = kwargs.get("axis", None)
                    keepdims = kwargs.get("keepdims", False)
                    mask = np.logical_or.reduce(
                        mask,
                        where=where_unmasked,
                        axis=axis,
                        keepdims=keepdims,
                        out=out_mask,
                    )
                    if where_unmasked is not True:
                        # Mask also whole rows in which no elements were selected;
                        # those will have been left as unmasked above.
                        mask |= ~np.logical_or.reduce(
                            where_unmasked, axis=axis, keepdims=keepdims
                        )

                else:
                    # Accumulate
                    axis = kwargs.get("axis", 0)
                    mask = np.logical_or.accumulate(mask, axis=axis, out=out_mask)

            elif out is None:
                # Can only get here if neither input nor output was masked, but
                # perhaps where was masked (possible in "not NUMPY_LT_1_25" and
                # in NUMPY_LT_1_21 (latter also allowed axis).
                # We don't support this.
                return NotImplemented

        elif method in {"reduceat", "at"}:  # pragma: no cover
            raise NotImplementedError(
                "masked instances cannot yet deal with 'reduceat' or 'at'."
            )

        if out_unmasked is not None:
            kwargs["out"] = out_unmasked
        if where_unmasked is not True:
            kwargs["where"] = where_unmasked
        result = getattr(ufunc, method)(*unmasked, **kwargs)

        if result is None:  # pragma: no cover
            # This happens for the "at" method.
            return result

        if out is not None and len(out) == 1:
            out = out[0]
        return self._masked_result(result, mask, out)

    def __array_function__(self, function, types, args, kwargs):
        # TODO: go through functions systematically to see which ones
        # work and/or can be supported.
        if function in MASKED_SAFE_FUNCTIONS:
            return super().__array_function__(function, types, args, kwargs)

        elif function in APPLY_TO_BOTH_FUNCTIONS:
            helper = APPLY_TO_BOTH_FUNCTIONS[function]
            try:
                helper_result = helper(*args, **kwargs)
            except NotImplementedError:
                return self._not_implemented_or_raise(function, types)

            data_args, mask_args, kwargs, out = helper_result
            if out is not None:
                if not isinstance(out, Masked):
                    return self._not_implemented_or_raise(function, types)
                function(*mask_args, out=out.mask, **kwargs)
                function(*data_args, out=out.unmasked, **kwargs)
                return out

            mask = function(*mask_args, **kwargs)
            result = function(*data_args, **kwargs)

        elif function in DISPATCHED_FUNCTIONS:
            dispatched_function = DISPATCHED_FUNCTIONS[function]
            try:
                dispatched_result = dispatched_function(*args, **kwargs)
            except NotImplementedError:
                return self._not_implemented_or_raise(function, types)

            if dispatched_result is None:
                return None

            result, mask, out = dispatched_result

        elif function in UNSUPPORTED_FUNCTIONS:
            return NotImplemented

        else:  # pragma: no cover
            # By default, just pass it through for now.
            return super().__array_function__(function, types, args, kwargs)

        if mask is None:
            return result
        else:
            return self._masked_result(result, mask, out)

    def _not_implemented_or_raise(self, function, types):
        # Our function helper or dispatcher found that the function does not
        # work with Masked.  In principle, there may be another class that
        # knows what to do with us, for which we should return NotImplemented.
        # But if there is ndarray (or a non-Masked subclass of it) around,
        # it quite likely coerces, so we should just break.
        if any(issubclass(t, np.ndarray) and not issubclass(t, Masked) for t in types):
            raise TypeError(
                f"the MaskedNDArray implementation cannot handle {function} "
                "with the given arguments."
            ) from None
        else:
            return NotImplemented

    def _masked_result(self, result, mask, out):
        if isinstance(result, tuple):
            if out is None:
                out = (None,) * len(result)
            if not isinstance(mask, (list, tuple)):
                mask = (mask,) * len(result)
            return tuple(
                self._masked_result(result_, mask_, out_)
                for (result_, mask_, out_) in zip(result, mask, out)
            )

        if out is None:
            # Note that we cannot count on result being the same class as
            # 'self' (e.g., comparison of quantity results in an ndarray, most
            # operations on Longitude and Latitude result in Angle or
            # Quantity), so use Masked to determine the appropriate class.
            return Masked(result, mask)

        # TODO: remove this sanity check once test cases are more complete.
        assert isinstance(out, Masked)
        # If we have an output, the result was written in-place, so we should
        # also write the mask in-place (if not done already in the code).
        if out._mask is not mask:
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
        if "where" not in kwargs:
            kwargs["where"] = ~self.mask
        if initial_func is not None and "initial" not in kwargs:
            kwargs["initial"] = initial_func(self.unmasked)
        return kwargs

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        # Unfortunately, cannot override the call to diagonal inside trace, so
        # duplicate implementation in numpy/core/src/multiarray/calculation.c.
        diagonal = self.diagonal(offset=offset, axis1=axis1, axis2=axis2)
        return diagonal.sum(-1, dtype=dtype, out=out)

    def min(self, axis=None, out=None, **kwargs):
        return super().min(
            axis=axis, out=out, **self._reduce_defaults(kwargs, np.nanmax)
        )

    def max(self, axis=None, out=None, **kwargs):
        return super().max(
            axis=axis, out=out, **self._reduce_defaults(kwargs, np.nanmin)
        )

    def ptp(self, axis=None, out=None, **kwargs):
        result = self.max(axis=axis, out=out, **kwargs)
        result -= self.min(axis=axis, **kwargs)
        return result

    def nonzero(self):
        unmasked_nonzero = self.unmasked.nonzero()
        if self.ndim >= 1:
            not_masked = ~self.mask[unmasked_nonzero]
            return tuple(u[not_masked] for u in unmasked_nonzero)
        else:
            return unmasked_nonzero if not self.mask else np.nonzero(0)

    def compress(self, condition, axis=None, out=None):
        if out is not None:
            raise NotImplementedError("cannot yet give output")
        return self._apply("compress", condition, axis=axis)

    def repeat(self, repeats, axis=None):
        return self._apply("repeat", repeats, axis=axis)

    def choose(self, choices, out=None, mode="raise"):
        # Let __array_function__ take care since choices can be masked too.
        return np.choose(self, choices, out=out, mode=mode)

    def argmin(self, axis=None, out=None, *, keepdims=False):
        # TODO: should this return a masked integer array, with masks
        # if all elements were masked?
        at_min = self == self.min(axis=axis, keepdims=True)
        return at_min.filled(False).argmax(axis=axis, out=out, keepdims=keepdims)

    def argmax(self, axis=None, out=None, *, keepdims=False):
        at_max = self == self.max(axis=axis, keepdims=True)
        return at_max.filled(False).argmax(axis=axis, out=out, keepdims=keepdims)

    def argsort(self, axis=-1, kind=None, order=None, *, stable=None):
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
        order : str or list of str.
            For an array with fields defined, the fields to compare first,
            second, etc.  A single field can be specified as a string, and not
            all fields need be specified, but unspecified fields will still be
            used, in dtype order, to break ties.
        stable: bool, keyword-only, ignored
            Sort stability. Present only to allow subclasses to work.

        Returns
        -------
        index_array : ndarray, int
            Array of indices that sorts along the specified ``axis``.  Use
            ``np.take_along_axis(self, index_array, axis=axis)`` to obtain
            the sorted array.

        """
        if axis is None:
            data = self.ravel()
            axis = -1
        else:
            data = self

        if self.dtype.names:
            # As done inside the argsort implementation in multiarray/methods.c.
            if order is None:
                order = self.dtype.names
            elif NUMPY_LT_2_0:
                order = np.core._internal._newnames(self.dtype, order)
            else:
                order = np._core._internal._newnames(self.dtype, order)

            keys = tuple(data[name] for name in order[::-1])

        elif order is not None:
            raise ValueError("Cannot specify order when the array has no fields.")

        else:
            keys = (data,)

        return np.lexsort(keys, axis=axis)

    def sort(self, axis=-1, kind=None, order=None, *, stable=False):
        """Sort an array in-place. Refer to `numpy.sort` for full documentation.

        Notes
        -----
        Masked items will be sorted to the end. The implementation
        is via `numpy.lexsort` and thus ignores the ``kind`` and ``stable`` arguments;
        they are present only so that subclasses can pass them on.
        """
        # TODO: probably possible to do this faster than going through argsort!
        argsort_kwargs = dict(kind=kind, order=order)
        if not NUMPY_LT_2_0:
            argsort_kwargs["stable"] = stable
        indices = self.argsort(axis, **argsort_kwargs)
        self[:] = np.take_along_axis(self, indices, axis=axis)

    def argpartition(self, kth, axis=-1, kind="introselect", order=None):
        # TODO: should be possible to do this faster than with a full argsort!
        return self.argsort(axis=axis, order=order)

    def partition(self, kth, axis=-1, kind="introselect", order=None):
        # TODO: should be possible to do this faster than with a full argsort!
        return self.sort(axis=axis, order=None)

    def cumsum(self, axis=None, dtype=None, out=None):
        if axis is None:
            self = self.ravel()
            axis = 0
        return np.add.accumulate(self, axis=axis, dtype=dtype, out=out)

    def cumprod(self, axis=None, dtype=None, out=None):
        if axis is None:
            self = self.ravel()
            axis = 0
        return np.multiply.accumulate(self, axis=axis, dtype=dtype, out=out)

    def clip(self, min=None, max=None, out=None, **kwargs):
        """Return an array whose values are limited to ``[min, max]``.

        Like `~numpy.clip`, but any masked values in ``min`` and ``max``
        are ignored for clipping.  The mask of the input array is propagated.
        """
        # TODO: implement this at the ufunc level.
        dmin, mmin = self._get_data_and_mask(min)
        dmax, mmax = self._get_data_and_mask(max)
        if mmin is None and mmax is None:
            # Fast path for unmasked max, min.
            return super().clip(min, max, out=out, **kwargs)

        masked_out = np.positive(self, out=out)
        out = masked_out.unmasked
        if dmin is not None:
            np.maximum(out, dmin, out=out, where=True if mmin is None else ~mmin)
        if dmax is not None:
            np.minimum(out, dmax, out=out, where=True if mmax is None else ~mmax)
        return masked_out

    def mean(self, axis=None, dtype=None, out=None, keepdims=False, *, where=True):
        # Implementation based on that in numpy/core/_methods.py
        # Cast bool, unsigned int, and int to float64 by default,
        # and do float16 at higher precision.
        is_float16_result = False
        if dtype is None:
            if issubclass(self.dtype.type, (np.integer, np.bool_)):
                dtype = np.dtype("f8")
            elif issubclass(self.dtype.type, np.float16):
                dtype = np.dtype("f4")
                is_float16_result = out is None

        where = ~self.mask & where

        result = self.sum(
            axis=axis, dtype=dtype, out=out, keepdims=keepdims, where=where
        )
        n = np.add.reduce(where, axis=axis, keepdims=keepdims)

        # catch the case when an axis is fully masked to prevent div by zero:
        n = np.add.reduce(where, axis=axis, keepdims=keepdims)
        neq0 = n == 0
        n += neq0
        result /= n

        # correct fully-masked slice results to what is expected for 0/0 division
        result.unmasked[neq0] = np.nan

        if is_float16_result:
            result = result.astype(self.dtype)
        return result

    def var(
        self, axis=None, dtype=None, out=None, ddof=0, keepdims=False, *, where=True
    ):
        where_final = ~self.mask & where

        # Simplified implementation based on that in numpy/core/_methods.py
        n = np.add.reduce(where_final, axis=axis, keepdims=keepdims)[...]

        # Cast bool, unsigned int, and int to float64 by default.
        if dtype is None and issubclass(self.dtype.type, (np.integer, np.bool_)):
            dtype = np.dtype("f8")
        mean = self.mean(axis=axis, dtype=dtype, keepdims=True, where=where)

        x = self - mean
        x *= x.conjugate()  # Conjugate just returns x if not complex.

        result = x.sum(
            axis=axis, dtype=dtype, out=out, keepdims=keepdims, where=where_final
        )
        n -= ddof
        n = np.maximum(n, 0, out=n)
        result /= n
        result._mask |= n == 0
        return result

    def std(
        self, axis=None, dtype=None, out=None, ddof=0, keepdims=False, *, where=True
    ):
        result = self.var(
            axis=axis, dtype=dtype, out=out, ddof=ddof, keepdims=keepdims, where=where
        )
        return np.sqrt(result, out=result)

    def __bool__(self):
        # First get result from array itself; this will error if not a scalar.
        result = super().__bool__()
        return result and not self.mask

    def any(self, axis=None, out=None, keepdims=False, *, where=True):
        return np.logical_or.reduce(
            self, axis=axis, out=out, keepdims=keepdims, where=~self.mask & where
        )

    def all(self, axis=None, out=None, keepdims=False, *, where=True):
        return np.logical_and.reduce(
            self, axis=axis, out=out, keepdims=keepdims, where=~self.mask & where
        )

    # Following overrides needed since somehow the ndarray implementation
    # does not actually call these.
    def __str__(self):
        return np.array_str(self)

    def __repr__(self):
        return np.array_repr(self)

    def __format__(self, format_spec):
        string = super().__format__(format_spec)
        if self.shape == () and self.mask:
            n = min(3, max(1, len(string)))
            return " " * (len(string) - n) + "\u2014" * n
        else:
            return string


class MaskedRecarray(np.recarray, MaskedNDArray, data_cls=np.recarray):
    # Explicit definition since we need to override some methods.

    def __array_finalize__(self, obj):
        # recarray.__array_finalize__ does not do super, so we do it
        # explicitly.
        super().__array_finalize__(obj)
        super(np.recarray, self).__array_finalize__(obj)

    # __getattribute__, __setattr__, and field use these somewhat
    # obscrure ndarray methods.  TODO: override in MaskedNDArray?
    def getfield(self, dtype, offset=0):
        for field, info in self.dtype.fields.items():
            if offset == info[1] and dtype == info[0]:
                return self[field]

        raise NotImplementedError("can only get existing field from structured dtype.")

    def setfield(self, val, dtype, offset=0):
        for field, info in self.dtype.fields.items():
            if offset == info[1] and dtype == info[0]:
                self[field] = val
                return

        raise NotImplementedError("can only set existing field from structured dtype.")
