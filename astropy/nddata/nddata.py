# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import numpy as np
from copy import deepcopy

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

    The key distinction from raw `numpy.ndarray` is the presence of
    additional metadata such as uncertainties, a mask, units,
    and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray`, `~numpy.ndarray`-like, or `NDData`
        The actual `data` contained in this `NDData` object.

    uncertainty : any type, optional
        Uncertainty on the data. The `uncertainty` *should* have a string
        attribute named ``uncertainty_type``, but there is otherwise no
        restriction.

    mask : any type, optional
        Mask for the data.

    wcs : undefined, optional
        WCS-object containing the world coordinate system for the data.

    meta : `dict`-like object, optional
        Metadata for this object. Must be `dict`-like but no further
        restriction is placed on meta.

    unit : `~astropy.units.UnitBase` instance or `str`, optional
        The units of the data.

    copy : `bool`
        ``True`` if the passed parameters should be copied for the new instance
        or ``False`` if not. This affects every parameter even the data!
        Default is False.

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
                 meta=None, unit=None, copy=False):

        # Rather pointless since the NDDataBase does not implement any setting
        # but before this PR (#4234) the NDDataBase did call the uncertainty
        # setter and if anyone wants to alter this behaviour again this call
        # to the superclass NDDataBase should be in here.
        super(NDData, self).__init__()

        # Check if data is any type from which to collect some implicitly
        # passed parameters.
        if isinstance(data, NDData):  # don't use self.__class__ (issue #4137)
            # Of course we need to check the data because subclasses with other
            # init-logic might try to init this class. We could skip these
            # tests if we compared for self.__class__ but that has other
            # drawbacks.

            # Comparing if there is an explicit and an implicit unit parameter.
            # If that is the case use the explicit one and issue a warning
            # that there might be a conflict. In case there is no explicit
            # unit just overwrite the unit parameter with the NDData.unit
            # and proceed as if that one was given as parameter. Same for the
            # other parameters.
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
            if hasattr(data, 'mask') and hasattr(data, 'data'):
                # We actually need the data to have a mask _and_ data attribute
                if mask is not None:
                    log.info("Overwriting Masked Objects's current "
                             "mask with specified mask")
                else:
                    mask = data.mask

                # Just save the data for further processing, we could have
                # a masked Quantity here or something else entirely. Better to
                # check it first.
                data = data.data

            if isinstance(data, Quantity):
                if unit is not None:
                    log.info("Overwriting Quantity's current "
                             "unit with specified unit")
                else:
                    unit = data.unit
                data = data.value

        # Quick check on the parameters if they match the requirements
        # Could be completely moved inside setters...
        if (not hasattr(data, 'shape') or
            not hasattr(data, '__getitem__') or
            not hasattr(data, '__array__')):
            # Data doesn't look like a numpy array, try converting it to
            # one.
            data = np.array(data, subok=True, copy=False)
            # Quick check to see if what we got out looks like an array
            # rather than an object (since numpy will convert a
            # non-numerical/string inputs to an array of objects).
            if data.dtype == 'O':
                raise TypeError("Could not convert data to numpy array.")

        if meta is None:
            meta = OrderedDict()
        elif not isinstance(meta, collections.Mapping):
            raise TypeError("meta attribute must be dict-like")

        if unit is not None:
            unit = Unit(unit)

        if copy:
            data = deepcopy(data)
            mask = deepcopy(mask)
            wcs = deepcopy(wcs)
            meta = deepcopy(meta)
            uncertainty = deepcopy(uncertainty)
            # Actually - copying the unit is unnecessary but better safe
            # than sorry :-)
            unit = deepcopy(unit)

        # Store the attributes
        self._data = data
        self._mask = mask
        self._wcs = wcs
        self._meta = meta
        self._unit = unit
        # Call the setter for uncertainty to further check the uncertainty
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
        `~numpy.ndarray`: the data.
        """
        return self._data

    @property
    def mask(self):
        """
        any type: Mask for the data, if any.

        Using `~numpy.ndarray`-like masks containing ``bools`` is *recommended*
        if the `NDData` is used for any `~numpy.ma.MaskedArray` operations.
        Valid values should have a corresponding mask value of ``False`` while
        invalid ones should be ``True`` (following the `numpy` convention).
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
        any type: wcs for the data, if any.

        Even though nothing is enforced, using `~astropy.wcs.WCS` as ``wcs``
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

        Uncertainty *should* have an attribute ``uncertainty_type`` that is
        a string. `~astropy.nddata.NDUncertainty`-subclasses are recommended,
        if `~astropy.nddata.NDArithmeticMixin` with uncertainty propagation
        is used.
        """
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            # Check for the uncertainty_type attribute and that it is a string
            # and issue a warning if it has not.
            if (not hasattr(value, 'uncertainty_type') or
                    not isinstance(value.uncertainty_type, six.string_types)):
                log.info('Uncertainty should have attribute uncertainty_type '
                         'whose type is string.')
            if isinstance(value, NDUncertainty):
                # If it is a subclass of NDUncertainty we must set the
                # parent_nddata attribute. (#4152)
                value.parent_nddata = self
        self._uncertainty = value
