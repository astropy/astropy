# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module contains a class equivalent to pre-1.0 NDData.


import numpy as np

from astropy import log
from astropy.units import Unit, UnitConversionError, UnitsError  # noqa: F401

from .flag_collection import FlagCollection
from .mixins.ndarithmetic import NDArithmeticMixin
from .mixins.ndio import NDIOMixin
from .mixins.ndslicing import NDSlicingMixin
from .nddata import NDData
from .nduncertainty import NDUncertainty

__all__ = ["NDDataArray"]


class NDDataArray(NDArithmeticMixin, NDSlicingMixin, NDIOMixin, NDData):
    """
    An ``NDData`` object with arithmetic. This class is functionally equivalent
    to ``NDData`` in astropy  versions prior to 1.0.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units, flags,
    and/or a coordinate system.

    See also: https://docs.astropy.org/en/stable/nddata/

    Parameters
    ----------
    data : ndarray or `NDData`
        The actual data contained in this `NDData` object. Not that this
        will always be copies by *reference* , so you should make copy
        the ``data`` before passing it in if that's the  desired behavior.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    mask : array-like, optional
        Mask for the data, given as a boolean Numpy array or any object that
        can be converted to a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where
        the data is *valid* and ``True`` when it is not (like Numpy
        masked arrays). If ``data`` is a numpy masked array, providing
        ``mask`` here will causes the mask from the masked array to be
        ignored.

    flags : array-like or `~astropy.nddata.FlagCollection`, optional
        Flags giving information about each pixel. These can be specified
        either as a Numpy array of any type (or an object which can be converted
        to a Numpy array) with a shape matching that of the
        data, or as a `~astropy.nddata.FlagCollection` instance which has a
        shape matching that of the data.

    wcs : None, optional
        WCS-object containing the world coordinate system for the data.

        .. warning::
            This is not yet defined because the discussion of how best to
            represent this class's WCS system generically is still under
            consideration. For now just leave it as None

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object.  e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    unit : `~astropy.units.UnitBase` instance or str, optional
        The units of the data.


    Raises
    ------
    ValueError :
        If the `uncertainty` or `mask` inputs cannot be broadcast (e.g., match
        shape) onto ``data``.
    """

    def __init__(self, data, *args, flags=None, **kwargs):
        # Initialize with the parent...
        super().__init__(data, *args, **kwargs)

        # ...then reset uncertainty to force it to go through the
        # setter logic below. In base NDData all that is done is to
        # set self._uncertainty to whatever uncertainty is passed in.
        self.uncertainty = self._uncertainty

        # Same thing for mask.
        self.mask = self._mask

        # Initial flags because it is no longer handled in NDData
        # or NDDataBase.
        if isinstance(data, NDDataArray):
            if flags is None:
                flags = data.flags
            else:
                log.info(
                    "Overwriting NDDataArrays's current flags with specified flags"
                )
        self.flags = flags

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
                if not self.unit and value._unit:
                    # Raise an error if uncertainty has unit and data does not
                    raise ValueError(
                        "Cannot assign an uncertainty with unit "
                        "to {} without "
                        "a unit".format(class_name)
                    )
                self._uncertainty = value
                self._uncertainty.parent_nddata = self
            else:
                raise TypeError(
                    "Uncertainty must be an instance of a NDUncertainty object"
                )
        else:
            self._uncertainty = value

    # Override unit so that we can add a setter.
    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        from . import conf

        try:
            if self._unit is not None and conf.warn_setting_unit_directly:
                log.info(
                    "Setting the unit directly changes the unit without "
                    "updating the data or uncertainty. Use the "
                    ".convert_unit_to() method to change the unit and "
                    "scale values appropriately."
                )
        except AttributeError:
            # raised if self._unit has not been set yet, in which case the
            # warning is irrelevant
            pass

        if value is None:
            self._unit = None
        else:
            self._unit = Unit(value)

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
                raise ValueError(
                    f"dimensions of mask {mask.shape} and data {self.data.shape} do not match"
                )
            else:
                self._mask = mask
        else:
            # internal representation should be one numpy understands
            self._mask = np.ma.nomask

    @property
    def shape(self):
        """
        shape tuple of this object's data.
        """
        return self.data.shape

    @property
    def size(self):
        """
        integer size of this object's data.
        """
        return self.data.size

    @property
    def dtype(self):
        """
        `numpy.dtype` of this object's data.
        """
        return self.data.dtype

    @property
    def ndim(self):
        """
        integer dimensions of this object's data.
        """
        return self.data.ndim

    @property
    def flags(self):
        return self._flags

    @flags.setter
    def flags(self, value):
        if value is not None:
            if isinstance(value, FlagCollection):
                if value.shape != self.shape:
                    raise ValueError("dimensions of FlagCollection does not match data")
                else:
                    self._flags = value
            else:
                flags = np.array(value, copy=False)
                if flags.shape != self.shape:
                    raise ValueError("dimensions of flags do not match data")
                else:
                    self._flags = flags
        else:
            self._flags = value

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

    def convert_unit_to(self, unit, equivalencies=[]):
        """
        Returns a new `NDData` object whose values have been converted
        to a new unit.

        Parameters
        ----------
        unit : `astropy.units.UnitBase` instance or str
            The unit to convert to.

        equivalencies : list of tuple
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`astropy:unit_equivalencies`.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Raises
        ------
        `~astropy.units.UnitsError`
            If units are inconsistent.

        """
        if self.unit is None:
            raise ValueError("No unit specified on source data")
        data = self.unit.to(unit, self.data, equivalencies=equivalencies)
        if self.uncertainty is not None:
            uncertainty_values = self.unit.to(
                unit, self.uncertainty.array, equivalencies=equivalencies
            )
            # should work for any uncertainty class
            uncertainty = self.uncertainty.__class__(uncertainty_values)
        else:
            uncertainty = None
        if self.mask is not None:
            new_mask = self.mask.copy()
        else:
            new_mask = None
        # Call __class__ in case we are dealing with an inherited type
        result = self.__class__(
            data,
            uncertainty=uncertainty,
            mask=new_mask,
            wcs=self.wcs,
            meta=self.meta,
            unit=unit,
        )

        return result
