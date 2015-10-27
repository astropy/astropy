# coding: utf-8
# Licensed like the corresponding numpy file; see licenses/NUMPY_LICENSE.rst
"""
Utilities that manipulate strides to achieve desirable effects.

An explanation of strides can be found in the "ndarray.rst" file in the
NumPy reference guide.

Notes
-----
The version provided here ensures broadcast_arrays passes on subclasses
if one sets ``subok=True``; see https://github.com/numpy/numpy/pull/4622

"""
from __future__ import division, absolute_import, print_function

import numpy as np

__all__ = ['broadcast_arrays', 'broadcast_to', 'GE1P10']


def GE1P10(module=np):
    return hasattr(module, 'broadcast_to')


if GE1P10():
    from numpy.lib.stride_tricks import broadcast_arrays, broadcast_to

else:
    from numpy.lib.stride_tricks import DummyArray

    def _maybe_view_as_subclass(original_array, new_array):
        if type(original_array) is not type(new_array):
            # if input was an ndarray subclass and subclasses were OK,
            # then view the result as that subclass.
            new_array = new_array.view(type=type(original_array))
            # Since we have done something akin to a view from original_array, we
            # should let the subclass finalize (if it has it implemented, i.e., is
            # not None).
            if new_array.__array_finalize__:
                new_array.__array_finalize__(original_array)
        return new_array


    def as_strided(x, shape=None, strides=None, subok=False):
        """ Make an ndarray from the given array with the given shape and strides.
        """
        # first convert input to array, possibly keeping subclass
        x = np.array(x, copy=False, subok=subok)
        interface = dict(x.__array_interface__)
        if shape is not None:
            interface['shape'] = tuple(shape)
        if strides is not None:
            interface['strides'] = tuple(strides)
        array = np.asarray(DummyArray(interface, base=x))

        if array.dtype.fields is None and x.dtype.fields is not None:
            # This should only happen if x.dtype is [('', 'Vx')]
            array.dtype = x.dtype

        return _maybe_view_as_subclass(x, array)


    def _broadcast_to(array, shape, subok, readonly):
        shape = tuple(shape) if np.iterable(shape) else (shape,)
        array = np.array(array, copy=False, subok=subok)
        if not shape and array.shape:
            raise ValueError('cannot broadcast a non-scalar to a scalar array')
        if any(size < 0 for size in shape):
            raise ValueError('all elements of broadcast shape must be non-'
                             'negative')
        broadcast = np.nditer(
            (array,), flags=['multi_index', 'refs_ok', 'zerosize_ok'],
            op_flags=['readonly'], itershape=shape, order='C').itviews[0]
        result = _maybe_view_as_subclass(array, broadcast)
        if not readonly and array.flags.writeable:
            result.flags.writeable = True
        return result


    def broadcast_to(array, shape, subok=False):
        """Broadcast an array to a new shape.

        Parameters
        ----------
        array : array_like
            The array to broadcast.
        shape : tuple
            The shape of the desired array.
        subok : bool, optional
            If True, then sub-classes will be passed-through, otherwise
            the returned array will be forced to be a base-class array (default).

        Returns
        -------
        broadcast : array
            A readonly view on the original array with the given shape. It is
            typically not contiguous. Furthermore, more than one element of a
            broadcasted array may refer to a single memory location.

        Raises
        ------
        ValueError
            If the array is not compatible with the new shape according to NumPy's
            broadcasting rules.

        Examples
        --------
        >>> x = np.array([1, 2, 3])
        >>> np.broadcast_to(x, (3, 3))
        array([[1, 2, 3],
               [1, 2, 3],
               [1, 2, 3]])
        """
        return _broadcast_to(array, shape, subok=subok, readonly=True)


    def _broadcast_shape(*args):
        """Returns the shape of the arrays that would result from broadcasting the
        supplied arrays against each other.
        """
        if not args:
            raise ValueError('must provide at least one argument')
        if len(args) == 1:
            # a single argument does not work with np.broadcast
            return np.asarray(args[0]).shape
        # use the old-iterator because np.nditer does not handle size 0 arrays
        # consistently
        b = np.broadcast(*args[:32])
        # unfortunately, it cannot handle 32 or more arguments directly
        for pos in range(32, len(args), 31):
            b = np.broadcast(b, *args[pos:(pos + 31)])
        return b.shape


    def broadcast_arrays(*args, **kwargs):
        """
        Broadcast any number of arrays against each other.

        Parameters
        ----------
        `*args` : array_likes
            The arrays to broadcast.

        subok : bool, optional
            If True, then sub-classes will be passed-through, otherwise
            the returned arrays will be forced to be a base-class array (default).

        Returns
        -------
        broadcasted : list of arrays
            These arrays are views on the original arrays.  They are typically
            not contiguous.  Furthermore, more than one element of a
            broadcasted array may refer to a single memory location.  If you
            need to write to the arrays, make copies first.

        Examples
        --------
        >>> x = np.array([[1,2,3]])
        >>> y = np.array([[1],[2],[3]])
        >>> np.broadcast_arrays(x, y)
        [array([[1, 2, 3],
               [1, 2, 3],
               [1, 2, 3]]), array([[1, 1, 1],
               [2, 2, 2],
               [3, 3, 3]])]

        Here is a useful idiom for getting contiguous copies instead of
        non-contiguous views.

        >>> [np.array(a) for a in np.broadcast_arrays(x, y)]
        [array([[1, 2, 3],
               [1, 2, 3],
               [1, 2, 3]]), array([[1, 1, 1],
               [2, 2, 2],
               [3, 3, 3]])]

        """
        # nditer is not used here to avoid the limit of 32 arrays.
        # Otherwise, something like the following one-liner would suffice:
        # return np.nditer(args, flags=['multi_index', 'zerosize_ok'],
        #                  order='C').itviews

        subok = kwargs.pop('subok', False)
        if kwargs:
            raise TypeError('broadcast_arrays() got an unexpected keyword '
                            'argument {}'.format(kwargs.pop()))
        args = [np.array(_m, copy=False, subok=subok) for _m in args]

        shape = _broadcast_shape(*args)

        if all(array.shape == shape for array in args):
            # Common case where nothing needs to be broadcasted.
            return args

        # TODO: consider making the results of broadcast_arrays readonly to match
        # broadcast_to. This will require a deprecation cycle.
        return [_broadcast_to(array, shape, subok=subok, readonly=False)
                for array in args]
