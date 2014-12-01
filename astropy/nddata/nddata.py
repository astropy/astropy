# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from ..units import Unit, Quantity
from .. import log

from ..io import registry as io_registry
from ..config import ConfigAlias
from ..utils.metadata import MetaData
from ..extern import six

__all__ = ['NDData']


__doctest_skip__ = ['NDData']


WARN_UNSUPPORTED_CORRELATED = ConfigAlias(
    '0.4', 'WARN_UNSUPPORTED_CORRELATED', 'warn_unsupported_correlated',
    'astropy.nddata.nddata', 'astropy.nddata')


class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units,
    and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray` or `NDData`
        The actual data contained in this `NDData` object. Not that this
        will always be copies by *reference* , so you should make copy
        the ``data`` before passing it in if that's the  desired behavior.

    uncertainty : `~astropy.nddata.NDUncertainty`, optional
        Uncertainties on the data.

    mask : `~numpy.ndarray`-like, optional
        Mask for the data, given as a boolean Numpy array or any object that
        can be converted to a boolean Numpy array with a shape
        matching that of the data. The values must be ``False`` where
        the data is *valid* and ``True`` when it is not (like Numpy
        masked arrays). If ``data`` is a numpy masked array, providing
        ``mask`` here will causes the mask from the masked array to be
        ignored.

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

    unit : `~astropy.units.UnitBase` instance or str, optional
        The units of the data.


    Raises
    ------
    ValueError :
        If the `uncertainty` or `mask` inputs cannot be broadcast (e.g., match
        shape) onto ``data``.

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

    meta = MetaData()

    def __init__(self, data, uncertainty=None, mask=None, wcs=None,
                 meta=None, unit=None):

        if isinstance(data, self.__class__):
            self._data = np.array(data.data, subok=True, copy=False)
            self.uncertainty = data.uncertainty
            self._mask = data.mask
            self._wcs = data.wcs
            self.meta = data.meta
            self._unit = data.unit

            if uncertainty is not None:
                self._uncertainty = uncertainty
                log.info("Overwriting NDData's current uncertainty being"
                         " overwritten with specified uncertainty")

            if mask is not None:
                self._mask = mask
                log.info("Overwriting NDData's current mask with specified mask")

            if wcs is not None:
                self._wcs = wcs
                log.info("Overwriting NDData's current wcs with specified wcs")

            if meta is not None:
                self.meta = meta
                log.info("Overwriting NDData's current meta with specified meta")

            if unit is not None:
                raise ValueError('To convert to different unit please use .to')
        else:
            if not hasattr(data, 'shape') or not hasattr(data, '__getitem__'):
                raise TypeError('Data must have shape attribute '
                                'and be slicable.')
            if hasattr(data, 'mask'):
                self._data = np.array(data.data, subok=True, copy=False)

                if mask is not None:
                    self._mask = mask
                    log.info("NDData was created with a masked array, and a "
                             "mask was explictly provided to NDData. The explicitly "
                             "passed-in mask will be used and the masked array's "
                             "mask will be ignored.")
                else:
                    self._mask = data.mask
            elif isinstance(data, Quantity):
                self._data = np.array(data.value, subok=True, copy=False)
                self._mask = mask
            else:
                self._data = data #np.array(data, subok=True, copy=False)
                self._mask = mask

            self._wcs = wcs
            self.meta = meta
            if isinstance(data, Quantity):
                if unit is not None:
                    raise ValueError("Cannot use the unit argument when data "
                                     "is a Quantity")
                else:
                    self._unit = data.unit
            else:
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
        return self._data

    @property
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            try:
                not_good_uncertainty = not isinstance(value.uncertainty_type,
                                                      six.string_types)
            except AttributeError:
                not_good_uncertainty = True
            finally:
                if not_good_uncertainty:
                    raise TypeError('Uncertainty must have attribute '
                                    'uncertainty_type whose type is string.')
        self._uncertainty = value

    @property
    def unit(self):
        return self._unit

    @property
    def wcs(self):
        return self._wcs

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

    read = classmethod(io_registry.read)
    write = io_registry.write
