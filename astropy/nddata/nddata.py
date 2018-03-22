# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.


import numpy as np
from copy import deepcopy

from .nddata_base import NDDataBase
from .nduncertainty import NDUncertainty, UnknownUncertainty
from .. import log
from ..units import Unit, Quantity
from ..utils.metadata import MetaData

__all__ = ['NDData']

_meta_doc = """`dict`-like : Additional meta information about the dataset."""


class NDData(NDDataBase):
    """
    A container for `numpy.ndarray`-based datasets, using the
    `~astropy.nddata.NDDataBase` interface.

    The key distinction from raw `numpy.ndarray` is the presence of
    additional metadata such as uncertainty, mask, unit, a coordinate system
    and/or a dictionary containing further meta information. This class *only*
    provides a container for *storing* such datasets. For further functionality
    take a look at the ``See also`` section.

    Parameters
    -----------
    data : `numpy.ndarray`-like or `NDData`-like
        The dataset.

    uncertainty : any type, optional
        Uncertainty in the dataset.
        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, for example ``"std"`` for standard deviation or
        ``"var"`` for variance. A metaclass defining such an interface is
        `NDUncertainty` - but isn't mandatory. If the uncertainty has no such
        attribute the uncertainty is stored as `UnknownUncertainty`.
        Defaults to ``None``.

    mask : any type, optional
        Mask for the dataset. Masks should follow the ``numpy`` convention that
        **valid** data points are marked by ``False`` and **invalid** ones with
        ``True``.
        Defaults to ``None``.

    wcs : any type, optional
        World coordinate system (WCS) for the dataset.
        Default is ``None``.

    meta : `dict`-like object, optional
        Additional meta information about the dataset. If no meta is provided
        an empty `collections.OrderedDict` is created.
        Default is ``None``.

    unit : `~astropy.units.Unit`-like or str, optional
        Unit for the dataset. Strings that can be converted to a
        `~astropy.units.Unit` are allowed.
        Default is ``None``.

    copy : `bool`, optional
        Indicates whether to save the arguments as copy. ``True`` copies
        every attribute before saving it while ``False`` tries to save every
        parameter as reference.
        Note however that it is not always possible to save the input as
        reference.
        Default is ``False``.

        .. versionadded:: 1.2

    Raises
    ------
    TypeError
        In case ``data`` or ``meta`` don't meet the restrictions.

    Notes
    -----
    Each attribute can be accessed through the homonymous instance attribute:
    ``data`` in a `NDData` object can be accessed through the `data`
    attribute::

        >>> from astropy.nddata import NDData
        >>> nd = NDData([1,2,3])
        >>> nd.data
        array([1, 2, 3])

    Given a conflicting implicit and an explicit parameter during
    initialization, for example the ``data`` is a `~astropy.units.Quantity` and
    the unit parameter is not ``None``, then the implicit parameter is replaced
    (without conversion) by the explicit one and a warning is issued::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> q = np.array([1,2,3,4]) * u.m
        >>> nd2 = NDData(q, unit=u.cm)
        INFO: overwriting Quantity's current unit with specified unit. [astropy.nddata.nddata]
        >>> nd2.data  # doctest: +FLOAT_CMP
        array([1., 2., 3., 4.])
        >>> nd2.unit
        Unit("cm")

    See also
    --------
    NDDataRef
    NDDataArray
    """

    # Instead of a custom property use the MetaData descriptor also used for
    # Tables. It will check if the meta is dict-like or raise an exception.
    meta = MetaData(doc=_meta_doc, copy=False)

    def __init__(self, data, uncertainty=None, mask=None, wcs=None,
                 meta=None, unit=None, copy=False):

        # Rather pointless since the NDDataBase does not implement any setting
        # but before the NDDataBase did call the uncertainty
        # setter. But if anyone wants to alter this behaviour again the call
        # to the superclass NDDataBase should be in here.
        super().__init__()

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
                log.info("overwriting NDData's current "
                         "unit with specified unit.")
            elif data.unit is not None:
                unit = data.unit

            if uncertainty is not None and data.uncertainty is not None:
                log.info("overwriting NDData's current "
                         "uncertainty with specified uncertainty.")
            elif data.uncertainty is not None:
                uncertainty = data.uncertainty

            if mask is not None and data.mask is not None:
                log.info("overwriting NDData's current "
                         "mask with specified mask.")
            elif data.mask is not None:
                mask = data.mask

            if wcs is not None and data.wcs is not None:
                log.info("overwriting NDData's current "
                         "wcs with specified wcs.")
            elif data.wcs is not None:
                wcs = data.wcs

            if meta is not None and data.meta is not None:
                log.info("overwriting NDData's current "
                         "meta with specified meta.")
            elif data.meta is not None:
                meta = data.meta

            data = data.data

        else:
            if hasattr(data, 'mask') and hasattr(data, 'data'):
                # Separating data and mask
                if mask is not None:
                    log.info("overwriting Masked Objects's current "
                             "mask with specified mask.")
                else:
                    mask = data.mask

                # Just save the data for further processing, we could be given
                # a masked Quantity or something else entirely. Better to check
                # it first.
                data = data.data

            if isinstance(data, Quantity):
                if unit is not None and unit != data.unit:
                    log.info("overwriting Quantity's current "
                             "unit with specified unit.")
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
            raise TypeError("could not convert data to numpy array.")

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
        self.mask = mask
        self._wcs = wcs
        self.meta = meta  # TODO: Make this call the setter sometime
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
    def uncertainty(self):
        """
        any type : Uncertainty in the dataset, if any.

        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, such as ``'std'`` for standard deviation or
        ``'var'`` for variance. A metaclass defining such an interface is
        `~astropy.nddata.NDUncertainty` but isn't mandatory.
        """
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            # There is one requirements on the uncertainty: That
            # it has an attribute 'uncertainty_type'.
            # If it does not match this requirement convert it to an unknown
            # uncertainty.
            if not hasattr(value, 'uncertainty_type'):
                log.info('uncertainty should have attribute uncertainty_type.')
                value = UnknownUncertainty(value, copy=False)

            # If it is a subclass of NDUncertainty we must set the
            # parent_nddata attribute. (#4152)
            if isinstance(value, NDUncertainty):
                # In case the uncertainty already has a parent create a new
                # instance because we need to assume that we don't want to
                # steal the uncertainty from another NDData object
                if value._parent_nddata is not None:
                    value = value.__class__(value, copy=False)
                # Then link it to this NDData instance (internally this needs
                # to be saved as weakref but that's done by NDUncertainty
                # setter).
                value.parent_nddata = self
        self._uncertainty = value
