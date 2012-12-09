# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the base NDData class.

__all__ = ['NDData']

import numpy as np
import copy


from ..units import Unit, Quantity
from ..logger import log

from .flag_collection import FlagCollection
from .nduncertainty import IncompatibleUncertaintiesException, NDUncertainty
from ..utils.compat.odict import OrderedDict

from ..config import ConfigurationItem



def merge_meta(meta1, meta2, operation_name):
    """
    Merging meta-data for `NDData`-objects. Currently it will build a new dictionary using the operation name
    as a dictionary key.

    Parameters
    ----------

    meta1 : `~OrderedDict`
        first meta dictionary

    meta2 : `~OrderedDict`
        second meta dictionary


    Returns
    -------
        `~OrderedDict` containing the merged dictionaries

    """

    if (meta1 is None) and (meta2 is None):
        new_meta = None
    elif (meta1 is None) and (meta2 is not None):
        new_meta = meta2
    elif (meta1 is not None) and (meta2 is None):
        new_meta = meta1

    elif (meta1 is not None) and (meta2 is not None):
        new_meta = OrderedDict()
        new_meta[operation_name] = (meta1, meta2)




    return new_meta



class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of additional
    metadata such as uncertainties, a mask, units, flags, and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray` or `~astropy.nddata.NDData`
        The actual data contained in this `NDData` object.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    mask : `~numpy.ndarray`, optional
        Mask for the data, given as a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where the
        data is *valid* and ``True`` when it is not (as for Numpy masked
        arrays).

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

    units : `astropy.units.UnitBase` instance or str, optional
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
            self.data = np.array(data.data, subok=True)

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
                self.data = np.array(data.data, subok=True)
                self.mask = data.mask
            else:
                self.data = np.array(data, subok=True)
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
        return self._units

    @unit.setter
    def unit(self, value):
        if value is None:
            self._units = None
        else:
            self._units = Unit(value)

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
            meta=self.meta, units=self.unit, copy=False)

    def _arithmetic(self, operand, operation_name):
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
        operation_name2operation = {'add':'__add__', 'subtract':'__sub__', 'multiply':'__mul__', 'divide':'__div__'}

        operation = operation_name2operation[operation_name]

        if isinstance(operand, Quantity):
            new_data = getattr(self.__array__(), operation)(other.to(self.unit).value)
            return self.__class__(new_data, unit=self.unit)

        elif isinstance(operand, self.__class__):
            if self.shape != operand.shape:
                raise ValueError("operand shapes do not match")

            if self.wcs != operand.wcs:
                raise ValueError("WCS properties do not match")

            if (self.flags is not None) or (operand.flags is not None):
                log.warn('Currently flag arithmetic is not implemented. The result of this operation will have None'
                         'as a flag')

            if (self.uncertainty is not None) or (operand.uncertainty is not None):
                log.warn('Currently uncertainty propagation is not implemented. The result of this operation '
                         'will have None as an uncertainty')

            if (self.unit is not None) and (operand.unit is not None):

                if operation_name in ('add', 'subtract'):
                    operand_data = operand.unit.to(self.unit, operand.data)
                    new_unit = self.unit
                    new_data = getattr(self.__array__(), operation)(operand_data)

                elif operation_name in ('multiply', 'divide'):
                    new_unit = self.unit * operand.unit
                    new_data = getattr(self.__array__(), operation)(operand.data)

            elif (self.unit is None) and (operand.unit is None):
                new_unit = None
                new_data = getattr(self.__array__(), operation)(operand.data)

            else:
                raise TypeError('Both units must be None or an actual Unit. For a dimensionless unit please use Unit(1)')


            return self.__class__(new_data, unit=self.unit, wcs=self.wcs,
                                    meta=merge_meta(self.meta, operand.meta, operation_name))

        else:
            raise TypeError('Cannot {operation_name} of type {self_type} with {other_type}'.format(
                                    operation_name=operation_name, self_type=type(self), operand_type=type(operand)))



    def __add__(self, other):
        return self._arithmetic(other, 'add')

    def __sub__(self, other):
        return self._arithmetic(other, 'subtract')

    def __mul__(self, other):
        return self._arithmetic(other, 'multiply')

    def __div__(self, other):
        return self._arithmetic(other, 'divide')




    def convert_units_to(self, unit):
        """
        Returns a new `NDData` object whose values have been converted
        to a new unit.

        Parameters
        ----------
        unit : `astropy.units.UnitBase` instance or str
            The unit to convert to.

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
            raise ValueError("No units specified on source data")
        data = self.unit.to(unit, self.data)
        result = self.__class__(data)  # in case we are dealing with an inherited type

        result.uncertainty = self.uncertainty
        result.mask = self.mask
        result.flags = None
        result.wcs = self.wcs
        result.meta = self.meta
        result.unit = unit

        return result
