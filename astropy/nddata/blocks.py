# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""

import numpy as np

from .decorators import support_nddata

__all__ = ['reshape_as_blocks', 'block_reduce', 'block_replicate']


def _process_block_inputs(data, block_size):
    data = np.asanyarray(data)
    block_size = np.atleast_1d(block_size)

    if np.any(block_size <= 0):
        raise ValueError('block_size elements must be strictly positive')

    if data.ndim > 1 and len(block_size) == 1:
        block_size = np.repeat(block_size, data.ndim)

    if len(block_size) != data.ndim:
        raise ValueError('block_size must be a scalar or have the same '
                         'length as the number of data dimensions')

    block_size_int = block_size.astype(int)
    if np.any(block_size_int != block_size):  # e.g., 2.0 is OK, 2.1 is not
        raise ValueError('block_size elements must be integers')

    return data, block_size_int


def reshape_as_blocks(data, block_size):
    """
    Reshape a data array into blocks.

    This is useful to efficiently apply functions on block subsets of
    the data instead of using loops.  The reshaped array is a view of
    the input data array.

    .. versionadded:: 4.1

    Parameters
    ----------
    data : `~numpy.ndarray`
        The input data array.

    block_size : int or array_like (int)
        The integer block size along each axis.  If ``block_size`` is a
        scalar and ``data`` has more than one dimension, then
        ``block_size`` will be used for for every axis.  Each dimension
        of ``block_size`` must divide evenly into the corresponding
        dimension of ``data``.

    Returns
    -------
    output : `~numpy.ndarray`
        The reshaped array as a view of the input ``data`` array.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata import reshape_as_blocks
    >>> data = np.arange(16).reshape(4, 4)
    >>> data
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])
    >>> reshape_as_blocks(data, (2, 2))
    array([[[[ 0,  1],
             [ 4,  5]],
            [[ 2,  3],
             [ 6,  7]]],
           [[[ 8,  9],
             [12, 13]],
            [[10, 11],
             [14, 15]]]])
    """

    data, block_size = _process_block_inputs(data, block_size)

    if np.any(np.mod(data.shape, block_size) != 0):
        raise ValueError('Each dimension of block_size must divide evenly '
                         'into the corresponding dimension of data')

    nblocks = np.array(data.shape) // block_size
    new_shape = tuple(k for ij in zip(nblocks, block_size) for k in ij)
    nblocks_idx = tuple(range(0, len(new_shape), 2))  # even indices
    block_idx = tuple(range(1, len(new_shape), 2))  # odd indices

    return data.reshape(new_shape).transpose(nblocks_idx + block_idx)


@support_nddata
def block_reduce(data, block_size, func=np.sum):
    """
    Downsample a data array by applying a function to local blocks.

    If ``data`` is not perfectly divisible by ``block_size`` along a
    given axis then the data will be trimmed (from the end) along that
    axis.

    Parameters
    ----------
    data : array_like
        The data to be resampled.

    block_size : int or array_like (int)
        The integer block size along each axis.  If ``block_size`` is a
        scalar and ``data`` has more than one dimension, then
        ``block_size`` will be used for for every axis.

    func : callable, optional
        The method to use to downsample the data.  Must be a callable
        that takes in a `~numpy.ndarray` along with an ``axis`` keyword,
        which defines the axis or axes along which the function is
        applied.  The ``axis`` keyword must accept multiple axes as a
        tuple.  The default is `~numpy.sum`, which provides block
        summation (and conserves the data sum).

    Returns
    -------
    output : array_like
        The resampled data.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata import block_reduce
    >>> data = np.arange(16).reshape(4, 4)
    >>> block_reduce(data, 2)  # doctest: +FLOAT_CMP
    array([[10, 18],
           [42, 50]])

    >>> block_reduce(data, 2, func=np.mean)  # doctest: +FLOAT_CMP
    array([[  2.5,   4.5],
           [ 10.5,  12.5]])
    """

    data, block_size = _process_block_inputs(data, block_size)
    nblocks = np.array(data.shape) // block_size
    size_init = nblocks * block_size  # evenly-divisible size

    # trim data if necessary
    for axis in range(data.ndim):
        if data.shape[axis] != size_init[axis]:
            data = data.swapaxes(0, axis)
            data = data[:size_init[axis]]
            data = data.swapaxes(0, axis)

    reshaped = reshape_as_blocks(data, block_size)
    axis = tuple(range(data.ndim, reshaped.ndim))

    return func(reshaped, axis=axis)


@support_nddata
def block_replicate(data, block_size, conserve_sum=True):
    """
    Upsample a data array by block replication.

    Parameters
    ----------
    data : array_like
        The data to be block replicated.

    block_size : int or array_like (int)
        The integer block size along each axis.  If ``block_size`` is a
        scalar and ``data`` has more than one dimension, then
        ``block_size`` will be used for for every axis.

    conserve_sum : bool, optional
        If `True` (the default) then the sum of the output
        block-replicated data will equal the sum of the input ``data``.

    Returns
    -------
    output : array_like
        The block-replicated data.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata import block_replicate
    >>> data = np.array([[0., 1.], [2., 3.]])
    >>> block_replicate(data, 2)  # doctest: +FLOAT_CMP
    array([[0.  , 0.  , 0.25, 0.25],
           [0.  , 0.  , 0.25, 0.25],
           [0.5 , 0.5 , 0.75, 0.75],
           [0.5 , 0.5 , 0.75, 0.75]])

    >>> block_replicate(data, 2, conserve_sum=False)  # doctest: +FLOAT_CMP
    array([[0., 0., 1., 1.],
           [0., 0., 1., 1.],
           [2., 2., 3., 3.],
           [2., 2., 3., 3.]])
    """

    data, block_size = _process_block_inputs(data, block_size)
    for i in range(data.ndim):
        data = np.repeat(data, block_size[i], axis=i)

    if conserve_sum:
        data = data / float(np.prod(block_size))

    return data
