# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import numpy as np
from copy import deepcopy

from .nddata_base import NDDataBase
from .nduncertainty import NDUncertainty, UnknownUncertainty
from .. import log
from ..units import Unit, Quantity
from ..utils.compat.odict import OrderedDict
from ..extern import six

__all__ = ['NDData']

__doctest_skip__ = ['NDData']


class NDData(NDDataBase):
    """
    A container for `~numpy.ndarray`-based datasets, using the
    `~astropy.nddata.NDDataBase` interface.

    The key distinction from raw `numpy.ndarray` is the presence of
    additional metadata such as uncertainty, mask, unit, a coordinate system
    and/or a `dict` containg further meta information. This class *only*
    provides a container for *storing* such datasets. For further functionality
    take a look at the ``See also`` section.

    Parameters
    -----------
    data : `~numpy.ndarray`-like, `NDData`-like or list
        The dataset.

    uncertainty : any type, optional
        Uncertainty in the dataset.
        Should have an attribute
        ``uncertainty_type`` that defines what kind of uncertainty is stored,
        such as ``std`` for standard deviation or
        ``var`` for variance. A metaclass defining such an interface is
        `~astropy.nddata.NDUncertainty` but isn't mandatory.
        Defaults to ``None``.

        TODO: *Must* or *Should* have this ``uncertainty_type`` ? (mwcraig)
        I implemented UnknownUncertainty to work around that so what is
        saved in _uncertainty has an uncertainty_type but it's unknown :-)

    mask : any type, optional
        Mask for the dataset.
        Masks should follow the ``numpy`` convention that valid data points are
        marked by ``False`` and invalid ones with ``True``.
        Defaults to ``None``.

    wcs : any type, optional
        A world coordinate system (WCS) for the dataset.
        Default is ``None``.

    meta : `dict`-like object, optional
        Meta information about the dataset. If no meta is provided an empty
        `collections.OrderedDict` is created.

    unit : `~astropy.units.Unit`-like or str, optional
        Unit for the dataset.
        Default is ``None``.

    copy : `bool`, optional
        Save the attributes as copy or as reference. ``True`` copies every
        attribute before saving it while ``False`` tries to save every
        parameter as reference. Note however that it is not always possible to
        save each input as reference.
        Default is ``False``.

    Raises
    ------
    TypeError:
        In case any parameter does not fulfill the classes restrictions.

    Notes
    -----
    Each attribute can be accessed through the homonymous instance attribute,
    such as data in a `NDData` object can be accessed through the ``data``
    attribute.

    For example::

        >>> from astropy.nddata import NDData
        >>> nd = NDData([1,2,3])
        >>> nd.data
        array([1, 2, 3])

    Given an implicit parameter and an explicit one during initialization, for
    example the ``data`` is a `~astropy.units.Quantity` and the unit parameter
    is not None. Then the implicit parameter is replaced (without conversion)
    by the explicit one and a warning is issued::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> q = np.array([1,2,3,4]) * u.m
        >>> nd2 = NDData(q, unit=u.cm)  # doctest: +SKIP
        INFO: Overwriting Quantity's current unit with specified unit [astropy.nddata.nddata]
        >>> nd2.data  # doctest: +SKIP
        array([ 1.,  2.,  3.,  4.])
        >>> nd2.unit  # doctest: +SKIP
        Unit("cm")

    See also
    --------

    NDIOMixin
    NDSlicingMixin
    NDArithmeticMixin
    NDDataArray
    """

    def __init__(self, data, uncertainty=None, mask=None, wcs=None,
                 meta=None, unit=None, copy=False):

        # Rather pointless since the NDDataBase does not implement any setting
        # but before (#4234) the NDDataBase did call the uncertainty
        # setter. But if anyone wants to alter this behaviour again the call
        # to the superclass NDDataBase should be in here.
        super(NDData, self).__init__()

        # Check if data is any type from which to collect some implicitly
        # passed parameters.
        if isinstance(data, NDData):  # don't use self.__class__ (issue #4137)
            # Of course we need to check the data because subclasses with other
            # init-logic might be passed in here. We could skip these
            # tests if we compared for self.__class__ but that has other
            # drawbacks.

            # Comparing if there is an explicit and an implicit unit parameter.
            # If that is the case use the explicit one and issue a warning
            # that there might be a conflict. In case there is no explicit
            # unit just overwrite the unit parameter with the NDData.unit
            # and proceed as if that one was given as parameter. Same for the
            # other parameters.
            if (unit is not None and data.unit is not None and
                    unit != data.unit):
                # TODO: Clarify warning that unit is replaced instead of
                # converted?
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
                # Seperating data and mask
                if mask is not None:
                    log.info("Overwriting Masked Objects's current "
                             "mask with specified mask")
                else:
                    mask = data.mask

                # Just save the data for further processing, we could be given
                # a masked Quantity or something else entirely. Better to check
                # it first.
                data = data.data

            if isinstance(data, Quantity):
                if unit is not None and unit != data.unit:
                    log.info("Overwriting Quantity's current "
                             "unit with specified unit")
                else:
                    unit = data.unit
                data = data.value

        # Quick check on the parameters if they match the requirements.
        if (not hasattr(data, 'shape') or not hasattr(data, '__getitem__') or
                not hasattr(data, '__array__')):
            # Data doesn't look like a numpy array, try converting it to
            # one.
            data = np.array(data, subok=True, copy=False)

        # Another quick check to see if what we got looks like an array
        # rather than an object (since numpy will convert a
        # non-numerical/non-string inputs to an array of objects).
        if data.dtype == 'O':
            raise TypeError("Could not convert data to numpy array.")

        # Check if meta is a dict and create an empty one if no meta was given
        if meta is None:
            meta = OrderedDict()
        elif not isinstance(meta, collections.Mapping):
            raise TypeError("Meta attribute must be dict-like")

        if unit is not None:
            unit = Unit(unit)

        if copy:
            # Data might have been copied before but no way of validating
            # without another variable.
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
        `~numpy.ndarray`-like : The stored dataset.
        """
        return self._data

    @property
    def mask(self):
        """
        any type : Mask for the dataset, if any.

        Masks should follow the ``numpy`` convention that valid data points are
        marked by ``False`` and invalid ones with ``True``.
        """
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def unit(self):
        """
        `~astropy.units.Unit` : Unit for the dataset, if any.
        """
        return self._unit

    @property
    def wcs(self):
        """
        any type : A world coordinate system (WCS) for the dataset, if any.
        """
        return self._wcs

    @property
    def meta(self):
        """
        `dict`-like : Meta information about the dataset, if any.
        """
        return self._meta

    @property
    def uncertainty(self):
        """
        any type : Uncertainty in the dataset, if any.

        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, such as ``std`` for standard deviation or
        ``var`` for variance. A metaclass defining such an interface is
        `~astropy.nddata.NDUncertainty` but isn't mandatory.

        TODO: *Must* or *Should* have this ``uncertainty_type`` ? (mwcraig)
        I implemented UnknownUncertainty to work around that so what is
        saved in _uncertainty has an uncertainty_type but it's unknown :-)
        """
        if isinstance(self._uncertainty, UnknownUncertainty):
            return self._uncertainty.array
        else:
            return self._uncertainty

    @property
    def uncertainty_type(self):
        """
        str: The kind of uncertainty that is saved.

        The uncertainty_type is using the convention that standard deviation
        is abbreviated by ``std``, variance by ``var``, ...

        TODO: This is new and only a shortcut to the uncertainty type since
        the property hides the really saved uncertainty if it is an
        ``UnknownUncertainty``.
        """
        if self._uncertainty is None:
            return None
        else:
            return self._uncertainty.uncertainty_type

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            # There is one requirements on the uncertainty: That
            # it has an attribute 'uncertainty_type' which is a string.
            # If it does not match this requirement convert it to an unknown
            # uncertainty.
            if (not hasattr(value, 'uncertainty_type') or
                    not isinstance(value.uncertainty_type, six.string_types)):
                log.info('Uncertainty should have attribute uncertainty_type '
                         'whose type is string.')
                value = UnknownUncertainty(value, copy=False)

            # If it is a subclass of NDUncertainty we must set the
            # parent_nddata attribute. (#4152)
            if isinstance(value, NDUncertainty):
                value.parent_nddata = self
        self._uncertainty = value
