# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

__all__ = ['NDData']

import numpy as np

from ..units import Unit
from .. import log

from .flag_collection import FlagCollection
from .nduncertainty import IncompatibleUncertaintiesException, NDUncertainty
from ..utils.compat.odict import OrderedDict
from ..io import registry as io_registry
from ..config import ConfigurationItem

WARN_UNSUPPORTED_CORRELATED = ConfigurationItem(
    'warn_unsupported_correlated', True,
    'Whether to issue a warning if NDData arithmetic is performed with '
    'uncertainties and the uncertainties do not support the propagation '
    'of correlated uncertainties.'
)


class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units, flags,
    and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray` or `~astropy.nddata.NDData`
        The actual data contained in this `NDData` object. Not that this
        will always be copies by *reference* , so you should make copy
        the `data` before passing it in if that's the  desired behavior.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    mask : `~numpy.ndarray`, optional
        Mask for the data, given as a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where
        the data is *valid* and ``True`` when it is not (like Numpy
        masked arrays). If `data` is a numpy masked array, providing
        `mask` here will causes the mask from the masked array to be
        ignored.

    flags : `~numpy.ndarray` or `~astropy.nddata.FlagCollection`, optional
        Flags giving information about each pixel. These can be specified
        either as a Numpy array of any type with a shape matching that of the
        data, or as a `~astropy.nddata.FlagCollection` instance which has a
        shape matching that of the data.

    wcs : undefined, optional
        WCS-object containing the world coordinate system for the data.

        .. warning::
            This is not yet defind because the discussion of how best to
            represent this class's WCS system generically is still under
            consideration. For now just leave it as None

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute
        of this particular object.  e.g., creation date, unique identifier,
        simulation parameters, exposure time, telescope name, etc.

    unit : `astropy.units.UnitBase` instance or str, optional
        The units of the data.


    Raises
    ------
    ValueError
        If the `uncertainty` or `mask` inputs cannot be broadcast (e.g., match
        shape) onto `data`.

    Notes
    -----
    `NDData` objects can be easily converted to a regular Numpy array
    using `numpy.asarray`

    For example::

        >>> from astropy.nddata import NDData
        >>> import numpy as np
        >>> x = NDData([1,2,3])
        >>> np.asarray(x)
        array([1, 2, 3])

    If the `NDData` object has a `mask`, `numpy.asarray` will return a
    Numpy masked array.

    This is useful, for example, when plotting a 2D image using
    matplotlib::

        >>> from astropy.nddata import NDData
        >>> from matplotlib import pyplot as plt
        >>> x = NDData([[1,2,3], [4,5,6]])
        >>> plt.imshow(x)

    """

    def __init__(self, data, uncertainty=None, mask=None, flags=None, wcs=None,
                 meta=None, unit=None):

        if isinstance(data, self.__class__):
            self.data = np.array(data.data, subok=True, copy=False)

            if uncertainty is not None:
                self.uncertainty = uncertainty
                log.info("Overwriting NDData's current uncertainty being overwritten with specified uncertainty")

            if mask is not None:
                self.mask = mask
                log.info("Overwriting NDData's current mask being overwritten with specified mask")

            if flags is not None:
                self.flags = flags
                log.info("Overwriting NDData's current flags being overwritten with specified flag")

            if wcs is not None:
                self.wcs = wcs
                log.info("Overwriting NDData's current wcs being overwritten with specified wcs")

            if meta is not None:
                self.meta = meta
                log.info("Overwriting NDData's current meta being overwritten with specified meta")

            if unit is not None:
                raise ValueError('To convert to different unit please use .to')
        else:
            if hasattr(data, 'mask'):
                self.data = np.array(data.data, subok=True, copy=False)

                if mask is not None:
                    self.mask = mask
                    log.info("NDData was created with a masked array, and a "
                             "mask was explictly provided to NDData. The explicitly "
                             "passed-in mask will be used and the masked array's "
                             "mask will be ignored.")
                else:
                    self.mask = data.mask
            else:
                self.data = np.array(data, subok=True, copy=False)
                self.mask = mask

            self.uncertainty = uncertainty
            self.flags = flags
            self.wcs = wcs
            self.meta = meta
            self.unit = unit

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        if value is not None:
            if isinstance(value, np.ndarray):
                if value.dtype != np.bool_:
                    raise TypeError("mask must be a boolean Numpy array")
                else:
                    if value.shape != self.shape:
                        raise ValueError("dimensions of mask do not match data")
                    else:
                        self._mask = value
            else:
                raise TypeError("mask must be a Numpy array")
        else:
            self._mask = value

    @property
    def flags(self):
        return self._flags

    @flags.setter
    def flags(self, value):
        if value is not None:
            if isinstance(value, np.ndarray):
                if value.shape != self.shape:
                    raise ValueError("dimensions of flags do not match data")
                else:
                    self._flags = value
            elif isinstance(value, FlagCollection):
                if value.shape != self.shape:
                    raise ValueError("dimensions of FlagCollection does not match data")
                else:
                    self._flags = value
            else:
                raise TypeError("flags should be a Numpy array or a FlagCollection instance")
        else:
            self._flags = value

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            if isinstance(value, NDUncertainty):
                self._uncertainty = value
                self._uncertainty.parent_nddata = self
            else:
                raise TypeError("Uncertainty must be an instance of a NDUncertainty object")
        else:
            self._uncertainty = value

    @property
    def meta(self):
        return self._meta

    @meta.setter
    def meta(self, value):
        if value is None:
            self._meta = OrderedDict()
        else:
            try:
                self._meta = OrderedDict(value)
            except ValueError:
                raise TypeError('NDData meta attribute must be dict-like')

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        if value is None:
            self._unit = None
        else:
            self._unit = Unit(value)

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
        integer dimensions of this object's data
        """
        return self.data.ndim

    def __array__(self):
        """
        This allows code that requests a Numpy array to use an NDData
        object as a Numpy array.
        """
        if self.mask is not None:
            return np.ma.masked_array(self.data, self.mask)
        else:
            return self.data

    def __getitem__(self, item):

        new_data = self.data[item]

        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[item]
        else:
            new_uncertainty = None

        if self.mask is not None:
            new_mask = self.mask[item]
        else:
            new_mask = None

        if self.flags is not None:
            if isinstance(self.flags, np.ndarray):
                new_flags = self.flags[item]
            elif isinstance(self.flags, FlagCollection):
                raise NotImplementedError('Slicing complex Flags is currently not implemented')
        else:
            new_flags = None

        if self.wcs is not None:
            raise NotImplementedError('Slicing for WCS is not currently implemented')
        else:
            new_wcs = None

        return self.__class__(new_data, uncertainty=new_uncertainty, mask=new_mask, flags=new_flags, wcs=new_wcs,
                              meta=self.meta, unit=self.unit)

    def _arithmetic(self, operand, propagate_uncertainties, name, operation):
        """
        {name} another dataset (`operand`) to this dataset.

        Parameters
        ----------
        operand : `~astropy.nddata.NDData`
            The second operand in the operation a {operator} b
        propagate_uncertainties : bool
            Whether to propagate uncertainties following the propagation rules
            defined by the class used for the `uncertainty` attribute.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Notes
        -----
        This method requires the datasets to have identical WCS properties,
        equivalent units, and identical shapes. Flags and meta-data get set to
        None in the resulting dataset. The unit in the result is the same as
        the unit in `self`. Uncertainties are propagated, although correlated
        errors are not supported by any of the built-in uncertainty classes.
        If uncertainties are assumed to be correlated, a warning is issued by
        default (though this can be disabled via the
        `WARN_UNSUPPORTED_CORRELATED` configuration item). Values masked in
        either dataset before the operation are masked in the resulting
        dataset.
        """
        if self.wcs != operand.wcs:
            raise ValueError("WCS properties do not match")

        if not (self.unit is None and operand.unit is None):
            if (self.unit is None or operand.unit is None
                    or not self.unit.is_equivalent(operand.unit)):
                raise ValueError("operand units do not match")

        if self.shape != operand.shape:
            raise ValueError("operand shapes do not match")

        if self.unit is not None:
            operand_data = operand.unit.to(self.unit, operand.data)
        else:
            operand_data = operand.data
        data = operation(self.data, operand_data)
        result = self.__class__(data)  # in case we are dealing with an inherited type

        if propagate_uncertainties is None:
            result.uncertainty = None
        elif self.uncertainty is None and operand.uncertainty is None:
            result.uncertainty = None
        elif self.uncertainty is None:
            result.uncertainty = operand.uncertainty
        elif operand.uncertainty is None:
            result.uncertainty = self.uncertainty
        else:  # both self and operand have uncertainties
            if WARN_UNSUPPORTED_CORRELATED() and \
               (not self.uncertainty.support_correlated or
               not operand.uncertainty.support_correlated):
                log.info("The uncertainty classes used do not support the "
                         "propagation of correlated errors, so uncertainties"
                         " will be propagated assuming they are uncorrelated")
            try:
                method = getattr(self.uncertainty, propagate_uncertainties)
                result.uncertainty = method(operand, result.data)
            except IncompatibleUncertaintiesException:
                raise IncompatibleUncertaintiesException(
                    "Cannot propagate uncertainties of type {0:s} with uncertainties of "
                    "type {1:s} for {2:s}".format(
                        self.uncertainty.__class__.__name__,
                        operand.uncertainty.__class__.__name__,
                        name))

        if self.mask is None and operand.mask is None:
            result.mask = None
        elif self.mask is None:
            result.mask = operand.mask
        elif operand.mask is None:
            result.mask = self.mask
        else:  # combine masks as for Numpy masked arrays
            result.mask = self.mask & operand.mask

        result.flags = None
        result.wcs = self.wcs
        result.meta = None
        result.unit = self.unit

        return result

    def add(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_add"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "addition", np.add)
    add.__doc__ = _arithmetic.__doc__.format(name="Add", operator="+")

    def subtract(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_subtract"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "subtraction", np.subtract)
    subtract.__doc__ = _arithmetic.__doc__.format(name="Subtract", operator="-")

    def multiply(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_multiply"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "multiplication", np.multiply)
    multiply.__doc__ = _arithmetic.__doc__.format(name="Multiply", operator="*")

    def divide(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_divide"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "division", np.divide)
    divide.__doc__ = _arithmetic.__doc__.format(name="Divide", operator="/")

    def convert_unit_to(self, unit, equivalencies=[]):
        """
        Returns a new `NDData` object whose values have been converted
        to a new unit.

        Parameters
        ----------
        unit : `astropy.units.UnitBase` instance or str
            The unit to convert to.

        equivalencies : list of equivalence pairs, optional
           A list of equivalence pairs to try if the units are not
           directly convertible.  See :ref:`unit_equivalencies`.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Raises
        ------
        UnitsException
            If units are inconsistent.
        """
        if self.unit is None:
            raise ValueError("No unit specified on source data")
        data = self.unit.to(unit, self.data, equivalencies=equivalencies)
        result = self.__class__(data)  # in case we are dealing with an inherited type

        result.uncertainty = self.uncertainty
        result.mask = self.mask
        result.flags = None
        result.wcs = self.wcs
        result.meta = self.meta
        result.unit = unit

        return result

    read = classmethod(io_registry.read)
    write = io_registry.write
