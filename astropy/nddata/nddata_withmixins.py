# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module implements a class based on NDData with all Mixins.
"""


from .mixins.ndarithmetic import NDArithmeticMixin
from .mixins.ndio import NDIOMixin
from .mixins.ndslicing import NDSlicingMixin
from .nddata import NDData

__all__ = ["NDDataRef"]


class NDDataRef(NDArithmeticMixin, NDIOMixin, NDSlicingMixin, NDData):
    """Implements `NDData` with all Mixins.

    This class implements a `NDData`-like container that supports reading and
    writing as implemented in the ``astropy.io.registry`` and also slicing
    (indexing) and simple arithmetic (add, subtract, divide and multiply).

    Notes
    -----
    A key distinction from `NDDataArray` is that this class does not attempt
    to provide anything that was not defined in any of the parent classes.

    See Also
    --------
    NDData
    NDArithmeticMixin
    NDSlicingMixin
    NDIOMixin

    Examples
    --------
    The mixins allow operation that are not possible with `NDData` or
    `NDDataBase`, i.e. simple arithmetic::

        >>> from astropy.nddata import NDDataRef, StdDevUncertainty
        >>> import numpy as np

        >>> data = np.ones((3,3), dtype=float)
        >>> ndd1 = NDDataRef(data, uncertainty=StdDevUncertainty(data))
        >>> ndd2 = NDDataRef(data, uncertainty=StdDevUncertainty(data))

        >>> ndd3 = ndd1.add(ndd2)
        >>> ndd3.data  # doctest: +FLOAT_CMP
        array([[2., 2., 2.],
               [2., 2., 2.],
               [2., 2., 2.]])
        >>> ndd3.uncertainty.array  # doctest: +FLOAT_CMP
        array([[1.41421356, 1.41421356, 1.41421356],
               [1.41421356, 1.41421356, 1.41421356],
               [1.41421356, 1.41421356, 1.41421356]])

    see `NDArithmeticMixin` for a complete list of all supported arithmetic
    operations.

    But also slicing (indexing) is possible::

        >>> ndd4 = ndd3[1,:]
        >>> ndd4.data  # doctest: +FLOAT_CMP
        array([2., 2., 2.])
        >>> ndd4.uncertainty.array  # doctest: +FLOAT_CMP
        array([1.41421356, 1.41421356, 1.41421356])

    See `NDSlicingMixin` for a description how slicing works (which attributes)
    are sliced.
    """

    pass
