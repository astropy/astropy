# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module contains a class equivalent to pre-1.0 NDData.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from .. units import UnitsError

from .nddata import NDData
from .slicing import NDSlicingMixin
from .arithmetic import NDArithmeticMixin
from .nduncertainty import NDUncertainty

__all__ = ['NDDataArray']


class NDDataArray(NDArithmeticMixin, NDSlicingMixin, NDData):
    """
    An ``NDData`` object with arithmetic. This class is functionally equivalent
    to ``NDData`` in astropy  versions prior to 1.0.
    """

    def __init__(self, *arg, **kwd):
        # Initialize with the parent...
        super(NDDataArray, self).__init__(*arg, **kwd)

        # ...then reset uncertainty to force it to go through the
        # setter logic below. In base NDData all that is done is to
        # set self._uncertainty to whatever uncertainty is passed in.
        self.uncertainty = self._uncertainty

        # Same thing for mask...
        self.mask = self._mask

    # Implement uncertainty as NDUncertainty to support propagation of
    # uncertainties in arithmetic operations
    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            if isinstance(value, NDUncertainty):
                class_name = self.__class__.__name__
                if self.unit and value._unit:
                    try:
                        scaling = (1 * value._unit).to(self.unit)
                    except UnitsError:
                        raise UnitsError('Cannot convert unit of uncertainty '
                                         'to unit of '
                                         '{0} object.'.format(class_name))
                    value.array *= scaling
                elif not self.unit and value._unit:
                    # Raise an error if uncertainty has unit and data does not
                    raise ValueError("Cannot assign an uncertainty with unit "
                                     "to {0} without "
                                     "a unit".format(class_name))
                self._uncertainty = value
                self._uncertainty.parent_nddata = self
            else:
                raise TypeError("Uncertainty must be an instance of "
                                "a NDUncertainty object")
        else:
            self._uncertainty = value

    # Implement mask in a way that converts nicely to a numpy masked array
    @property
    def mask(self):
        if self._mask is np.ma.nomask:
            return None
        else:
            return self._mask

    @mask.setter
    def mask(self, value):
        # Check that value is not either type of null mask.
        if (value is not None) and (value is not np.ma.nomask):
            mask = np.array(value, dtype=np.bool_, copy=False)
            if mask.shape != self.data.shape:
                raise ValueError("dimensions of mask do not match data")
            else:
                self._mask = mask
        else:
            # internal representation should be one numpy understands
            self._mask = np.ma.nomask

    def __array__(self):
        """
        This allows code that requests a Numpy array to use an NDData
        object as a Numpy array.
        """
        if self.mask is not None:
            return np.ma.masked_array(self.data, self.mask)
        else:
            return np.array(self.data)

    def __array_prepare__(self, array, context=None):
        """
        This ensures that a masked array is returned if self is masked.
        """
        if self.mask is not None:
            return np.ma.masked_array(array, self.mask)
        else:
            return array
