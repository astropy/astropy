# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the base NDData class.

__all__ = ['NDData']

import numpy as np

from .flag_collection import FlagCollection
from .nderror import IncompatibleErrors, NDError


class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of additional
    metadata such as errors, a mask, units, flags, and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray`
        The actual data contained in this `NDData` object.

    error : `~astropy.nddata.NDError`, optional
        Errors on the data.

    mask : `~numpy.ndarray`, optional
        Mask for the data, given as a boolean Numpy array with a shape
        matching that of the data. The values should be ``False`` where the
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

    units : undefined, optional
        The units of the data.

        .. warning::
            The units scheme is under development. For now, just supply a
            string when relevant - the units system will likely be compatible
            with providing strings to initialize itself.

    copy : bool, optional
        If True, the array will be *copied* from the provided `data`, otherwise
        it will be referenced if possible (see `numpy.array` :attr:`copy`
        argument for details).

    Raises
    ------
    ValueError
        If the `error` or `mask` inputs cannot be broadcast (e.g., match
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

    def __init__(self, data, error=None, mask=None, flags=None, wcs=None,
                 meta=None, units=None, copy=True):

        self.data = np.array(data, subok=True, copy=copy)

        self.error = error
        self.mask = mask
        self.flags = flags
        self.wcs = wcs

        if meta is None:
            self.meta = {}
        else:
            self.meta = dict(meta)  # makes a *copy* of the passed-in meta

        self.units = units

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        if value is not None:
            if isinstance(value, np.ndarray):
                if value.dtype != np.bool_:
                    raise TypeError("mask should be a boolean Numpy array")
                else:
                    if value.shape != self.shape:
                        raise ValueError("dimensions of mask do not match data")
                    else:
                        self._mask = value
            else:
                raise TypeError("mask should be a Numpy array")
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
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        if value is not None:
            if isinstance(value, NDError):
                self._error = value
                self._error.parent = self
            else:
                raise TypeError("error should be an instance of a NDError object")
        else:
            self._error = value

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

    # At the moment, the units are required to be the same for addition and
    # subtraction. Once the units framework is implemented, we no longer have
    # to have this requirement, and can also add methods for multiplication
    # and division.

    def add(self, operand, propagate_errors=True):
        """
        Add another dataset (`operand`) to this dataset.

        Parameters
        ----------
        operand : `~astropy.nddata.NDData`
            The second operand in the operation a + b
        propagate_errors : bool
            Whether to propagate errors

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Notes
        -----
        This method requires the datasets to have identical WCS properties,
        identical units, and identical shapes. Flags and meta-data get set to
        None in the resulting dataset. Errors are propagated, although this
        feature is experimental, and errors are assumed to be uncorrelated.
        Values masked in either dataset before the operation are masked in
        the resulting dataset.
        """

        if self.wcs != operand.wcs:
            raise ValueError("WCS properties do not match")

        if self.units != operand.units:
            raise ValueError("operand units do not match")

        if self.shape != operand.shape:
            raise ValueError("operand shapes do not match")

        result = self.__class__(self.data + operand.data)  # in case we are dealing with an inherited type

        if not propagate_errors:
            result.error = None
        elif self.error is None and operand.error is None:
            result.error = None
        elif self.error is None:
            result.error = operand.error
        elif operand.error is None:
            result.error = self.error
        else:  # both self and operand have errors
            try:
                result.error = self.error.propagate_add(operand, result.data)
            except IncompatibleErrors:
                raise IncompatibleErrors("Cannot propagate errors of type {0:s} with errors of type {1:s} for addition".format(self.error.__class__.__name__, operand.error.__class__.__name__))

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
        result.units = self.units

        return result

    def subtract(self, operand, propagate_errors=True):
        """
        Subtract another dataset (`operand`) from this dataset.

        Parameters
        ----------
        operand : `~astropy.nddata.NDData`
            The second operand in the operation a - b
        propagate_errors : bool
            Whether to propagate errors

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Notes
        -----
        This method requires the datasets to have identical WCS properties,
        identical units, and identical shapes. Flags and meta-data get set to
        None in the resulting dataset. Errors are propagated, although this
        feature is experimental, and errors are assumed to be uncorrelated.
        Values masked in either dataset before the operation are masked in
        the resulting dataset.
        """

        if self.wcs != operand.wcs:
            raise ValueError("WCS properties do not match")

        if self.units != operand.units:
            raise ValueError("operand units do not match")

        if self.shape != operand.shape:
            raise ValueError("operand shapes do not match")

        result = self.__class__(self.data - operand.data)  # in case we are dealing with an inherited type

        if not propagate_errors:
            result.error = None
        elif self.error is None and operand.error is None:
            result.error = None
        elif self.error is None:
            result.error = operand.error
        elif operand.error is None:
            result.error = self.error
        else:  # both self and operand have errors
            try:
                result.error = self.error.propagate_subtract(operand, result.data)
            except IncompatibleErrors:
                raise IncompatibleErrors("Cannot propagate errors of type {0:s} with errors of type {1:s} for subtractition".format(self.error.__class__.__name__, operand.error.__class__.__name__))

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
        result.units = self.units

        return result
