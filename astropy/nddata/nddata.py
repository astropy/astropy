# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import numpy as np

from .nddata_base import NDDataBase
from .nduncertainty import NDUncertainty
from ..units import Unit, Quantity
from .. import log
from ..utils.compat.odict import OrderedDict
from ..extern import six

__all__ = ['NDData']

__doctest_skip__ = ['NDData']

class NDData(NDDataBase):
    """
    A basic class for array-based data.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units,
    and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray`, `~numpy.ndarray`-like, or `NDData`
        The actual data contained in this `NDData` object. If *possible*, data
        will not be copied `data`, so you should copy the ``data`` before
        passing it in if that's the desired behavior.

    uncertainty : any type, optional
        Uncertainty on the data. The uncertainty *should* have a string
        attribute named ``uncertainty_type``, but there is otherwise no
        restriction.

    mask : any type, optional
        Mask for the data.

    wcs : undefined, optional
        WCS-object containing the world coordinate system for the data.

    meta : `dict`-like object, optional
        Metadata for this object. Must be dict-like but no further restriction
        is placed on meta.

    unit : `~astropy.units.UnitBase` instance or str, optional
        The units of the data.

    Notes
    -----
    The data in a `NDData` object should be accessed through the data
    attribute.

    For example::

        >>> from astropy.nddata import NDData
        >>> x = NDData([1,2,3])
        >>> x.data
        array([1, 2, 3])
    """

    def __init__(self, data, uncertainty=None, mask=None, wcs=None,
                 meta=None, unit=None):

        super(NDData, self).__init__()

        # Check if data is any type from which to collect some not explicitly
        # passed parameters.
        if isinstance(data, NDData):  # don't use self.__class__ (issue #4137)
            # Of course we need to check the data because subclasses with other
            # init-logic might try to init this class.
            if unit is not None and data.unit is not None:
                log.info("Overwriting NDData's current "
                         "unit with specified unit")
            elif data.unit is not None:
                unit = data.unit

            if uncertainty is not None and data.uncertainty is not None:
                log.info("Overwriting NDData's current "
                         "uncertainty with specified uncertainty")
            elif data.uncertainty is not None:
                uncertainty = data.uncertainty

            if mask is not None and data.mask is not None:
                log.info("Overwriting NDData's current "
                         "mask with specified mask")
            elif data.mask is not None:
                mask = data.mask

            if wcs is not None and data.wcs is not None:
                log.info("Overwriting NDData's current "
                         "wcs with specified wcs")
            elif data.wcs is not None:
                wcs = data.wcs

            if meta is not None and data.meta is not None:
                log.info("Overwriting NDData's current "
                         "meta with specified meta")
            elif data.meta is not None:
                meta = data.meta

            data = data.data

        else:
            # We actually need the data to have a mask _and_ data attribute
            if hasattr(data, 'mask') and hasattr(data, 'data'):
                if mask is not None:
                    log.info("Overwriting Masked Objects's current "
                             "mask with specified mask")
                else:
                    mask = data.mask
                # Just save the data for further processing, we could handle
                # a masked Quantity here or something else entirely.
                data = data.data

            if isinstance(data, Quantity):
                if unit is not None:
                    log.info("Overwriting Quantity's current "
                             "unit with specified unit")
                else:
                    unit = data.unit
                data = data.value

        # Quick check on the parameters if they match the requirements
        if (not hasattr(data, 'shape') or
            not hasattr(data, '__getitem__') or
            not hasattr(data, '__array__')):
            # Data doesn't look like a numpy array, try converting it to
            # one.
            self._data = np.array(data, subok=True, copy=False)
            # Quick check to see if what we got out looks like an array
            # rather than an object (since numpy will convert a
            # non-numerical input to an array of objects).
            if self._data.dtype == 'O':
                raise TypeError("Could not convert data to numpy array.")
        else:
            self._data = data  # np.array(data, subok=True, copy=False)

        self._mask = mask

        self._wcs = wcs

        if meta is None:
            self._meta = OrderedDict()
        elif not isinstance(meta, collections.Mapping):
            raise TypeError("meta attribute must be dict-like")
        else:
            self._meta = meta

        if unit is not None:
            self._unit = Unit(unit)
        else:
            self._unit = None
        # This must come after self's unit has been set so that the unit
        # of the uncertainty, if any, can be converted to the unit of the
        # unit of self.
        self.uncertainty = uncertainty

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        prefix = self.__class__.__name__ + '('
        body = np.array2string(self.data, separator=', ', prefix=prefix)
        return ''.join([prefix, body, ')'])

    @property
    def data(self):
        """
        `~numpy.ndarray`: the data
        """
        return self._data

    @property
    def mask(self):
        """
        any type: Mask for the data, if any.

        Using `~numpy.ndarray`-like masks containing ``bools`` is *recommended*
        if the `NDData` is used for any `~numpy.ma.MaskedArray` operations.
        Valid values should have a corresponding mask value of ``False`` while
        invalid ones should be ``True`` (following the numpy convention).
        """
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def unit(self):
        """
        `~astropy.units.Unit`: Unit for the data, if any.
        """
        return self._unit

    @property
    def wcs(self):
        """
        any type: WCS for the data, if any.

        Even though nothing is enforced using `~astropy.wcs.WCS` as `WCS`
        attribute is encouraged.
        """
        return self._wcs

    @property
    def meta(self):
        """
        `dict`-like: Meta information, if any.
        """
        return self._meta

    @property
    def uncertainty(self):
        """
        any type: Uncertainty in the data, if any.

        Uncertainty should have an attribute ``uncertainty_type`` that is
        a string. `~astropy.nddata.NDUncertainty`-subclasses are recommended,
        if `~astropy.nddata.NDArithmeticMixin` is used.
        """
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            if (not hasattr(value, 'uncertainty_type') or
                    not isinstance(value.uncertainty_type, six.string_types)):
                log.info('Uncertainty should have attribute uncertainty_type '
                         ' whose type is string.')
            elif isinstance(value, NDUncertainty):
                # If it is a subclass of NDUncertainty we must set the
                # parent_nddata attribute.
                value.parent_nddata = self
        self._uncertainty = value
