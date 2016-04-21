# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module implements a class based on NDData with all Mixins.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .nddata import NDData

from .mixins.ndslicing import NDSlicingMixin
from .mixins.ndarithmetic import NDArithmeticMixin
from .mixins.ndio import NDIOMixin

__all__ = ['NDDataAllMixins']


class NDDataAllMixins(NDArithmeticMixin, NDIOMixin, NDSlicingMixin, NDData):
    """Implements `NDData` with all Mixins.

    This class implements a `NDData`-like container that supports reading and
    writing as implemented in the ``astropy.io.registry`` and also slicing
    (indexing) and simple arithmetics (add, subtract, divide and multiply).

    Notes
    -----
    A key distinction from `NDDataArray` is that this class does not attempt
    to provide anything that was not defined in any of the parent classes.

    See also
    --------
    NDData
    NDArithmeticMixin
    NDSlicingMixin
    NDIOMixin

    Examples
    --------
    Simple arithmetics::

        >>> from astropy.nddata import NDDataAllMixins, StdDevUncertainty
        >>> import numpy as np

        >>> data = np.ones((3,3))
        >>> ndd1 = NDDataAllMixins(data, uncertainty=StdDevUncertainty(data))
        >>> ndd2 = NDDataAllMixins(data, uncertainty=StdDevUncertainty(data))

        >>> ndd3 = ndd1.add(ndd2)
        >>> ndd3.data
        array([[ 2.,  2.,  2.],
               [ 2.,  2.,  2.],
               [ 2.,  2.,  2.]])
        >>> ndd3.uncertainty.array
        array([[ 1.41421356,  1.41421356,  1.41421356],
               [ 1.41421356,  1.41421356,  1.41421356],
               [ 1.41421356,  1.41421356,  1.41421356]])

    Slicing (Indexing)::

        >>> ndd4 = ndd3[1,:]
        >>> ndd4.data
        array([ 2.,  2.,  2.])
        >>> ndd4.uncertainty.array
        array([ 1.41421356,  1.41421356,  1.41421356])
    """
    pass
