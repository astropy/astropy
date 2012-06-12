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
        Masking of the data; Should be False/0 (or the empty string) where the
        data is *valid*.  All other values indicate that the value should be
        masked. Must be a shape that can be broadcast onto `data`.

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

    validate : bool, optional
        If False, no type or shape-checking or array conversion will occur.
        Note that if `validate` is False, :attr:`copy` will be ignored.

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
    def __init__(self, data, error=None, mask=None, wcs=None, meta=None,
                 units=None, copy=True, validate=True):
        if validate:
            self.data = np.array(data, subok=True, copy=copy)

            if error is None:
                self.error = None
            else:
                self.error = np.array(error, subok=True, copy=copy)

            if mask is None:
                self.mask = None
            else:
                self.mask = np.array(mask, subok=True, copy=copy)

            self._validate_mask_and_error()

            self.wcs = wcs
            self.units = units
            if meta is None:
                self.meta = {}
            else:
                self.meta = dict(meta)  # makes a *copy* of the passed-in meta
        else:
            self.data = data
            self.error = error
            self.mask = mask
            self.wcs = wcs
            self.meta = meta
            self.units = units

    def _validate_mask_and_error(self):
        """
        Raises ValueError if they don't match (using ~numpy.broadcast)
        """

        try:
            if self.mask is not None:
                np.broadcast(self.data, self.mask)
            maskmatch = True
        except ValueError:
            maskmatch = False

        try:
            if self.error is not None:
                np.broadcast(self.data, self.error)
            errmatch = True
        except ValueError:
            errmatch = False

        if not errmatch and not maskmatch:
            raise ValueError('NDData error and mask do not match data')
        elif not errmatch:
            raise ValueError('NDData error does not match data')
        elif not maskmatch:
            raise ValueError('NDData mask does not match data')

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
            dtchar = self.mask.dtype.char
            if dtchar == 'U':
                return self.mask == u''
            elif dtchar == 'S':
                return self.mask == b''
            else:
                return ~self.mask.astype(bool)

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
