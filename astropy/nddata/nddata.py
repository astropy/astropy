# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from ..units import Unit, Quantity
from .. import log

from ..io import registry as io_registry
from ..config import ConfigAlias
from ..utils.metadata import MetaData

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
            self.data = np.array(data.data, subok=True, copy=False)
            self.uncertainty = data.uncertainty
            self.mask = data.mask
            self.wcs = data.wcs
            self.meta = data.meta
            self.unit = data.unit

            if uncertainty is not None:
                self.uncertainty = uncertainty
                log.info("Overwriting NDData's current uncertainty being"
                         " overwritten with specified uncertainty")

            if mask is not None:
                self.mask = mask
                log.info("Overwriting NDData's current mask with specified mask")

            if wcs is not None:
                self.wcs = wcs
                log.info("Overwriting NDData's current wcs with specified wcs")

            if meta is not None:
                self.meta = meta
                log.info("Overwriting NDData's current meta with specified meta")

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
            elif isinstance(data, Quantity):
                self.data = np.array(data.value, subok=True, copy=False)
                self.mask = mask
            else:
                self.data = np.array(data, subok=True, copy=False)
                self.mask = mask

            self.wcs = wcs
            self.meta = meta
            if isinstance(data, Quantity):
                if unit is not None:
                    raise ValueError("Cannot use the unit argument when data "
                                     "is a Quantity")
                else:
                    self.unit = data.unit
            else:
                self.unit = unit
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
            if mask.shape != self.shape:
                raise ValueError("dimensions of mask do not match data")
            else:
                self._mask = mask
        else:
            # internal representation should be one numpy understands
            self._mask = np.ma.nomask

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        self._uncertainty = value

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        from . import conf

        try:
            if self._unit is not None and conf.warn_setting_unit_directly:
                log.info('Setting the unit directly changes the unit without '
                         'updating the data or uncertainty. Use the '
                         '.convert_unit_to() method to change the unit and '
                         'scale values appropriately.')
        except AttributeError:
            # raised if self._unit has not been set yet, in which case the
            # warning is irrelevant
            pass

        if value is None:
            self._unit = None
        else:
            self._unit = Unit(value)

    @property
    def wcs(self):
        return self._wcs

    @wcs.setter
    def wcs(self, value):
        self._wcs = value

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
            return np.array(self.data)

    def __array_prepare__(self, array, context=None):
        """
        This ensures that a masked array is returned if self is masked.
        """
        if self.mask is not None:
            return np.ma.masked_array(array, self.mask)
        else:
            return array

    def __getitem__(self, item):

        new_data = self.data[item]

        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[item]
        else:
            new_uncertainty = None

        if self.mask is not None:
            new_mask = self.mask[item]
            # mask setter expects an array, always
            if new_mask.shape == ():
                new_mask = np.array(new_mask)
        else:
            new_mask = None

        if self.wcs is not None:
            raise NotImplementedError('Slicing for WCS is not currently implemented')
        else:
            new_wcs = None

        return self.__class__(new_data, uncertainty=new_uncertainty,
                              mask=new_mask, wcs=new_wcs,
                              meta=self.meta, unit=self.unit)

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
        UnitsError
            If units are inconsistent.

        """
        if self.unit is None:
            raise ValueError("No unit specified on source data")
        data = self.unit.to(unit, self.data, equivalencies=equivalencies)
        if self.uncertainty is not None:
            uncertainty_values = self.unit.to(unit, self.uncertainty.array,
                                              equivalencies=equivalencies)
            # should work for any uncertainty class
            uncertainty = self.uncertainty.__class__(uncertainty_values)
        else:
            uncertainty = None
        if self.mask is not None:
            new_mask = self.mask.copy()
        else:
            new_mask = None
        # Call __class__ in case we are dealing with an inherited type
        result = self.__class__(data, uncertainty=uncertainty,
                                mask=new_mask,
                                wcs=self.wcs,
                                meta=self.meta, unit=unit)

        return result

    read = classmethod(io_registry.read)
    write = io_registry.write
