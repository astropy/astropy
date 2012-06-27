# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the base NDData class.

__all__ = ['NDData']

import numpy as np


class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as error arrays, bad pixel masks, units
    and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray`
        The actual data contained in this `NDData` object.

    error : `~numpy.ndarray`, optional
        Error of the data. This should be interpreted as a 1-sigma error (e.g,
        square root of the variance), under the assumption of Gaussian errors.
        Must be a shape that can be broadcast onto `data`.

         .. warning::
             The physical interpretation of the `error` array may change in the
             future, as it has not been intensively discussed. For now assume
             the above description holds, using an `error` property if
             necessary, but feel free to use the most convinient internal
             representation in subclasses

    mask : `~numpy.ndarray`, optional
        Masking of the data; Should be ``False`` where the data is *valid* and
        ``True`` when it is not (as for Numpy masked arrays).

    flags : `~numpy.ndarray`, optional
        Flags giving information about each pixel. While the flags he values
        inside the flags array. Flags can be any valid Numpy array type.

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

        if error is None:
            self.error = None
        else:
            self.error = np.array(error, subok=True, copy=copy)

        if mask is None:
            self.mask = None
        else:
            self.mask = np.array(mask, subok=True, copy=copy)

        if flags is None:
            self.flags = None
        else:
            self.flags = np.array(flags, subok=True, copy=copy)

        self.wcs = wcs
        self.units = units

        if meta is None:
            self.meta = {}
        else:
            self.meta = dict(meta)  # makes a *copy* of the passed-in meta

    def _get_mask(self):
        return self._mask

    def _set_mask(self, value):
        if value is not None:
            if isinstance(value, np.ndarray):
                if value.dtype != np.bool_:
                    raise TypeError("`mask` should be a boolean Numpy array")
                else:
                    if value.shape != self.shape:
                        raise ValueError("dimensions of `mask` do not match data")
                    else:
                        self._mask = value
            else:
                raise TypeError("`mask` should be a Numpy array")
        else:
            self._mask = value

    mask = property(_get_mask, _set_mask)

    def _get_flags(self):
        return self._flags

    def _set_flags(self, value):
        if value is not None:
            if isinstance(value, np.ndarray):
                try:
                    np.broadcast(self.data, value)
                except ValueError:
                    raise ValueError("dimensions of `flags` do not match data")
                else:
                    self._flags = value
            else:
                raise TypeError("`flags` should be a Numpy array")
        else:
            self._flags = value

    flags = property(_get_flags, _set_flags)

    def _get_error(self):
        return self._error

    def _set_error(self, value):
        if value is not None:
            try:
                np.broadcast(self.data, value)
            except ValueError:
                raise ValueError("dimensions of `error` do not match data")
            else:
                self._error = value
        else:
            self._error = value

    error = property(_get_error, _set_error)

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

    @property
    def boolmask(self):
        """
        The mask as a boolean array (or None if the mask is None).

        This mask is True where the data is *valid*, and False where the data
        should be *masked*.  This is the opposite of the convention used for
        `mask`, but allows simple retrieval of the unmasked data points as
        ``ndd.data[ndd.boolmask]``.
        """
        if self.mask is None:
            return None
        else:
            return ~self.mask

    def __array__(self):
        """
        This allows code that requests a Numpy array to use an NDData
        object as a Numpy array.
        """
        if self.mask is not None:
            return np.ma.masked_array(self.data, self.mask)
        else:
            return self.data
