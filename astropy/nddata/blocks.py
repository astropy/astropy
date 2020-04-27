# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module includes helper functions for array operations.
"""

import numpy as np

from .decorators import support_nddata

__all__ = ['block_reduce', 'block_replicate']


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
        which defines the axis along which the function is applied.  The
        default is `~numpy.sum`, which provides block summation (and
        conserves the data sum).

    Returns
    -------
    output : array_like
        The resampled data.

    Examples
    --------
    >>> import numpy as np
    >>> from astropy.nddata import block_reduce
    >>> data = np.arange(16).reshape(4, 4)
    >>> block_reduce(data, 2)    # doctest: +SKIP
    array([[10, 18],
           [42, 50]])

    >>> block_reduce(data, 2, func=np.mean)    # doctest: +SKIP
    array([[  2.5,   4.5],
           [ 10.5,  12.5]])
    """

    from skimage.measure import block_reduce

    data = np.asanyarray(data)

    block_size = np.atleast_1d(block_size)
    if data.ndim > 1 and len(block_size) == 1:
        block_size = np.repeat(block_size, data.ndim)

    if len(block_size) != data.ndim:
        raise ValueError('`block_size` must be a scalar or have the same '
                         'length as `data.shape`')

    block_size = np.array([int(i) for i in block_size])
    size_resampled = np.array(data.shape) // block_size
    size_init = size_resampled * block_size

    # trim data if necessary
    for i in range(data.ndim):
        if data.shape[i] != size_init[i]:
            data = data.swapaxes(0, i)
            data = data[:size_init[i]]
            data = data.swapaxes(0, i)

    return block_reduce(data, tuple(block_size), func=func)


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

    data = np.asanyarray(data)

    block_size = np.atleast_1d(block_size)
    if data.ndim > 1 and len(block_size) == 1:
        block_size = np.repeat(block_size, data.ndim)

    if len(block_size) != data.ndim:
        raise ValueError('`block_size` must be a scalar or have the same '
                         'length as `data.shape`')

    for i in range(data.ndim):
        data = np.repeat(data, block_size[i], axis=i)

    if conserve_sum:
        data = data / float(np.prod(block_size))

    return data
