# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""The ShapedLikeNDArray mixin class and shape-related functions."""

from __future__ import annotations

import abc
import numbers
from itertools import zip_longest
from math import prod
from typing import TYPE_CHECKING

import numpy as np

from astropy.utils.compat import NUMPY_LT_2_0
from astropy.utils.decorators import deprecated

if NUMPY_LT_2_0:
    import numpy.core as np_core
    from numpy.core.multiarray import normalize_axis_index
else:
    import numpy._core as np_core
    from numpy.lib.array_utils import normalize_axis_index

__all__ = [
    "IncompatibleShapeError",
    "NDArrayShapeMethods",
    "ShapedLikeNDArray",
    "check_broadcast",
    "simplify_basic_index",
    "unbroadcast",
]

if TYPE_CHECKING:
    from collections.abc import Sequence
    from types import EllipsisType
    from typing import Self, TypeVar

    from numpy.typing import NDArray

    DT = TypeVar("DT", bound=np.generic)


class NDArrayShapeMethods:
    """Mixin class to provide shape-changing methods.

    The class proper is assumed to have some underlying data, which are arrays
    or array-like structures. It must define a ``shape`` property, which gives
    the shape of those data, as well as an ``_apply`` method that creates a new
    instance in which a `~numpy.ndarray` method has been applied to those.

    Furthermore, for consistency with `~numpy.ndarray`, it is recommended to
    define a setter for the ``shape`` property, which, like the
    `~numpy.ndarray.shape` property allows in-place reshaping the internal data
    (and, unlike the ``reshape`` method raises an exception if this is not
    possible).

    This class only provides the shape-changing methods and is meant in
    particular for `~numpy.ndarray` subclasses that need to keep track of
    other arrays.  For other classes, `~astropy.utils.shapes.ShapedLikeNDArray`
    is recommended.

    """

    # Note to developers: if new methods are added here, be sure to check that
    # they work properly with the classes that use this, such as Time and
    # BaseRepresentation, i.e., look at their ``_apply`` methods and add
    # relevant tests.  This is particularly important for methods that imply
    # copies rather than views of data (see the special-case treatment of
    # 'flatten' in Time).

    def __getitem__(self, item):
        return self._apply("__getitem__", item)

    def copy(self, *args, **kwargs):
        """Return an instance containing copies of the internal data.

        Parameters are as for :meth:`~numpy.ndarray.copy`.
        """
        return self._apply("copy", *args, **kwargs)

    def reshape(self, *args, **kwargs):
        """Returns an instance containing the same data with a new shape.

        Parameters are as for :meth:`~numpy.ndarray.reshape`.  Note that it is
        not always possible to change the shape of an array without copying the
        data (see :func:`~numpy.reshape` documentation). If you want an error
        to be raise if the data is copied, you should assign the new shape to
        the shape attribute (note: this may not be implemented for all classes
        using ``NDArrayShapeMethods``).
        """
        return self._apply("reshape", *args, **kwargs)

    def ravel(self, *args, **kwargs):
        """Return an instance with the array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.ravel`. Note that it is
        not always possible to unravel an array without copying the data.
        If you want an error to be raise if the data is copied, you should
        should assign shape ``(-1,)`` to the shape attribute.
        """
        return self._apply("ravel", *args, **kwargs)

    def flatten(self, *args, **kwargs):
        """Return a copy with the array collapsed into one dimension.

        Parameters are as for :meth:`~numpy.ndarray.flatten`.
        """
        return self._apply("flatten", *args, **kwargs)

    def transpose(self, *args, **kwargs):
        """Return an instance with the data transposed.

        Parameters are as for :meth:`~numpy.ndarray.transpose`.  All internal
        data are views of the data of the original.
        """
        return self._apply("transpose", *args, **kwargs)

    @property
    def T(self) -> Self:
        """Return an instance with the data transposed.

        Parameters are as for :attr:`~numpy.ndarray.T`.  All internal
        data are views of the data of the original.
        """
        if self.ndim < 2:
            return self
        else:
            return self.transpose()

    def swapaxes(self, *args, **kwargs):
        """Return an instance with the given axes interchanged.

        Parameters are as for :meth:`~numpy.ndarray.swapaxes`:
        ``axis1, axis2``.  All internal data are views of the data of the
        original.
        """
        return self._apply("swapaxes", *args, **kwargs)

    def diagonal(self, *args, **kwargs):
        """Return an instance with the specified diagonals.

        Parameters are as for :meth:`~numpy.ndarray.diagonal`.  All internal
        data are views of the data of the original.
        """
        return self._apply("diagonal", *args, **kwargs)

    def squeeze(self, *args, **kwargs):
        """Return an instance with single-dimensional shape entries removed.

        Parameters are as for :meth:`~numpy.ndarray.squeeze`.  All internal
        data are views of the data of the original.
        """
        return self._apply("squeeze", *args, **kwargs)

    def take(self, indices, axis=None, out=None, mode="raise"):
        """Return a new instance formed from the elements at the given indices.

        Parameters are as for :meth:`~numpy.ndarray.take`, except that,
        obviously, no output array can be given.
        """
        if out is not None:
            return NotImplementedError("cannot pass 'out' argument to 'take.")

        return self._apply("take", indices, axis=axis, mode=mode)


class ShapedLikeNDArray(NDArrayShapeMethods, metaclass=abc.ABCMeta):
    """Mixin class to provide shape-changing methods.

    The class proper is assumed to have some underlying data, which are arrays
    or array-like structures. It must define a ``shape`` property, which gives
    the shape of those data, as well as an ``_apply`` method that creates a new
    instance in which a `~numpy.ndarray` method has been applied to those.

    Furthermore, for consistency with `~numpy.ndarray`, it is recommended to
    define a setter for the ``shape`` property, which, like the
    `~numpy.ndarray.shape` property allows in-place reshaping the internal data
    (and, unlike the ``reshape`` method raises an exception if this is not
    possible).

    This class also defines default implementations for ``ndim`` and ``size``
    properties, calculating those from the ``shape``.  These can be overridden
    by subclasses if there are faster ways to obtain those numbers.

    """

    # Note to developers: if new methods are added here, be sure to check that
    # they work properly with the classes that use this, such as Time and
    # BaseRepresentation, i.e., look at their ``_apply`` methods and add
    # relevant tests.  This is particularly important for methods that imply
    # copies rather than views of data (see the special-case treatment of
    # 'flatten' in Time).

    @property
    @abc.abstractmethod
    def shape(self) -> tuple[int, ...]:
        """The shape of the underlying data."""

    @abc.abstractmethod
    def _apply(method, *args, **kwargs):
        """Create a new instance, with ``method`` applied to underlying data.

        The method is any of the shape-changing methods for `~numpy.ndarray`
        (``reshape``, ``swapaxes``, etc.), as well as those picking particular
        elements (``__getitem__``, ``take``, etc.). It will be applied to the
        underlying arrays (e.g., ``jd1`` and ``jd2`` in `~astropy.time.Time`),
        with the results used to create a new instance.

        Parameters
        ----------
        method : str
            Method to be applied to the instance's internal data arrays.
        args : tuple
            Any positional arguments for ``method``.
        kwargs : dict
            Any keyword arguments for ``method``.

        """

    @property
    def ndim(self) -> int:
        """The number of dimensions of the instance and underlying arrays."""
        return len(self.shape)

    @property
    def size(self) -> int:
        """The size of the object, as calculated from its shape."""
        return prod(self.shape)

    @property
    def isscalar(self) -> bool:
        return self.shape == ()

    def __len__(self) -> int:
        if self.isscalar:
            raise TypeError(f"Scalar {self.__class__.__name__!r} object has no len()")
        return self.shape[0]

    def __bool__(self) -> bool:
        """Any instance should evaluate to True, except when it is empty."""
        return self.size > 0

    def __getitem__(self, item):
        try:
            return self._apply("__getitem__", item)
        except IndexError:
            if self.isscalar:
                raise TypeError(
                    f"scalar {self.__class.__name__!r} object is not subscriptable."
                )
            else:
                raise

    def __iter__(self):
        if self.isscalar:
            raise TypeError(
                f"scalar {self.__class__.__name__!r} object is not iterable."
            )

        # We cannot just write a generator here, since then the above error
        # would only be raised once we try to use the iterator, rather than
        # upon its definition using iter(self).
        def self_iter():
            for idx in range(len(self)):
                yield self[idx]

        return self_iter()

    # Functions that change shape or essentially do indexing.
    _APPLICABLE_FUNCTIONS = {
        np.moveaxis,
        np.rollaxis,
        np.atleast_1d,
        np.atleast_2d,
        np.atleast_3d,
        np.expand_dims,
        np.broadcast_to,
        np.flip,
        np.fliplr,
        np.flipud,
        np.rot90,
        np.roll,
        np.delete,
    }

    # Functions that themselves defer to a method. Those are all
    # defined in np.core.fromnumeric, but exclude alen as well as
    # sort and partition, which make copies before calling the method.
    _METHOD_FUNCTIONS = {
        getattr(np, name): {
            "amax": "max",
            "amin": "min",
            "around": "round",
            "round_": "round",
            "alltrue": "all",
            "sometrue": "any",
        }.get(name, name)
        for name in np_core.fromnumeric.__all__
        if name not in ["alen", "sort", "partition"]
    }
    # Add np.copy, which we may as well let defer to our method.
    _METHOD_FUNCTIONS[np.copy] = "copy"

    # Could be made to work with a bit of effort:
    # np.where, np.compress, np.extract,
    # np.diag_indices_from, np.triu_indices_from, np.tril_indices_from
    # np.tile, np.repeat (need .repeat method)
    # TODO: create a proper implementation.
    # Furthermore, some arithmetic functions such as np.mean, np.median,
    # could work for Time, and many more for TimeDelta, so those should
    # override __array_function__.
    def __array_function__(self, function, types, args, kwargs):
        """Wrap numpy functions that make sense."""
        if function in self._APPLICABLE_FUNCTIONS:
            if function is np.broadcast_to:
                # Ensure that any ndarray subclasses used are
                # properly propagated.
                kwargs.setdefault("subok", True)
            elif (
                function in {np.atleast_1d, np.atleast_2d, np.atleast_3d}
                and len(args) > 1
            ):
                return tuple(function(arg, **kwargs) for arg in args)

            if self is not args[0]:
                return NotImplemented

            return self._apply(function, *args[1:], **kwargs)

        # For functions that defer to methods, use the corresponding
        # method/attribute if we have it.  Otherwise, fall through.
        if self is args[0] and function in self._METHOD_FUNCTIONS:
            method = getattr(self, self._METHOD_FUNCTIONS[function], None)
            if method is not None:
                if callable(method):
                    return method(*args[1:], **kwargs)
                else:
                    # For np.shape, etc., just return the attribute.
                    return method

        # Fall-back, just pass the arguments on since perhaps the function
        # works already (see above).
        return function.__wrapped__(*args, **kwargs)


class IncompatibleShapeError(ValueError):
    def __init__(
        self,
        shape_a: tuple[int, ...],
        shape_a_idx: int,
        shape_b: tuple[int, ...],
        shape_b_idx: int,
    ) -> None:
        super().__init__(shape_a, shape_a_idx, shape_b, shape_b_idx)


@deprecated("7.0", alternative="np.broadcast_shapes")
def check_broadcast(*shapes: tuple[int, ...]) -> tuple[int, ...]:
    """
    Determines whether two or more Numpy arrays can be broadcast with each
    other based on their shape tuple alone.

    Parameters
    ----------
    *shapes : tuple
        All shapes to include in the comparison.  If only one shape is given it
        is passed through unmodified.  If no shapes are given returns an empty
        `tuple`.

    Returns
    -------
    broadcast : `tuple`
        If all shapes are mutually broadcastable, returns a tuple of the full
        broadcast shape.
    """
    if len(shapes) == 0:
        return ()
    elif len(shapes) == 1:
        return shapes[0]

    reversed_shapes = (reversed(shape) for shape in shapes)

    full_shape = []

    for dims in zip_longest(*reversed_shapes, fillvalue=1):
        max_dim = 1
        max_dim_idx = None
        for idx, dim in enumerate(dims):
            if dim == 1:
                continue

            if max_dim == 1:
                # The first dimension of size greater than 1
                max_dim = dim
                max_dim_idx = idx
            elif dim != max_dim:
                raise IncompatibleShapeError(
                    shapes[max_dim_idx], max_dim_idx, shapes[idx], idx
                )

        full_shape.append(max_dim)

    return tuple(full_shape[::-1])


def unbroadcast(array: NDArray[DT]) -> NDArray[DT]:
    """
    Given an array, return a new array that is the smallest subset of the
    original array that can be re-broadcasted back to the original array.

    See https://stackoverflow.com/questions/40845769/un-broadcasting-numpy-arrays
    for more details.
    """
    if array.ndim == 0:
        return array

    array = array[
        tuple((slice(0, 1) if stride == 0 else slice(None)) for stride in array.strides)
    ]

    # Remove leading ones, which are not needed in numpy broadcasting.
    first_not_unity = next(
        (i for (i, s) in enumerate(array.shape) if s > 1), array.ndim
    )

    return array.reshape(array.shape[first_not_unity:])


def simplify_basic_index(
    basic_index: int | slice | Sequence[int | slice | EllipsisType | None],
    *,
    shape: Sequence[int],
) -> tuple[int | slice, ...]:
    """
    Given a Numpy basic index, return a tuple of integers and slice objects
    with no default values (`None`) if possible.

    If one of the dimensions has a slice and the step is negative and the stop
    value of the slice was originally `None`, the new stop value of the slice
    may still be set to `None`.

    For more information on valid basic indices, see
    https://numpy.org/doc/stable/user/basics.indexing.html#basic-indexing

    Parameters
    ----------
    basic_index
        A valid Numpy basic index
    shape
        The shape of the array being indexed
    """
    ndim = len(shape)

    if not isinstance(basic_index, (tuple, list)):  # We just have a single int
        basic_index = (basic_index,)

    new_index = list(basic_index)

    if Ellipsis in new_index:
        if new_index.count(Ellipsis) > 1:
            raise IndexError("an index can only have a single ellipsis ('...')")

        # Replace the Ellipsis with the correct number of slice(None)s
        e_ind = new_index.index(Ellipsis)
        new_index.remove(Ellipsis)
        n_e = ndim - len(new_index)
        for i in range(n_e):
            ind = e_ind + i
            new_index.insert(ind, slice(0, shape[ind], 1))

    if len(new_index) > ndim:
        raise ValueError(
            f"The dimensionality of the basic index {basic_index} can not be greater "
            f"than the dimensionality ({ndim}) of the data."
        )

    for i in range(ndim):
        if i < len(new_index):
            slc = new_index[i]
            if isinstance(slc, slice):
                indices = list(slc.indices(shape[i]))
                # The following case is the only one where slice(*indices) does
                # not give the 'correct' answer because it will set stop to -1
                # which means the last element in the array.
                if indices[1] == -1:
                    indices[1] = None
                new_index[i] = slice(*indices)
            elif isinstance(slc, numbers.Integral):
                new_index[i] = normalize_axis_index(int(slc), shape[i])
            else:
                raise ValueError(f"Unexpected index element in basic index: {slc}")
        else:
            new_index.append(slice(0, shape[i], 1))

    return tuple(new_index)
