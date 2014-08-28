# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.

from __future__ import absolute_import, division, print_function, unicode_literals

from copy import deepcopy

import numpy as np

from ..units import Unit, Quantity, UnitsError, dimensionless_unscaled
from .. import log

from .flag_collection import FlagCollection
from .nduncertainty import IncompatibleUncertaintiesException, NDUncertainty
from ..io import registry as io_registry
from ..config import ConfigAlias
from ..utils.metadata import MetaData

__all__ = ['NDData', 'NDArithmetic', 'NDDataArithmetic']


__doctest_skip__ = ['NDData']


WARN_UNSUPPORTED_CORRELATED = ConfigAlias(
    '0.4', 'WARN_UNSUPPORTED_CORRELATED', 'warn_unsupported_correlated',
    'astropy.nddata.nddata', 'astropy.nddata')


class NDData(object):
    """A Superclass for array-based data in Astropy.

    The key distinction from raw numpy arrays is the presence of
    additional metadata such as uncertainties, a mask, units, flags,
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

    flags : `~numpy.ndarray`-like or `~astropy.nddata.FlagCollection`, optional
        Flags giving information about each pixel. These can be specified
        either as a Numpy array of any type (or an object which can be converted
        to a Numpy array) with a shape matching that of the
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

    def __init__(self, data, uncertainty=None, mask=None, flags=None, wcs=None,
                 meta=None, unit=None):

        if isinstance(data, self.__class__):
            self.data = np.array(data.data, subok=True, copy=False)
            self.uncertainty = data.uncertainty
            self.mask = data.mask
            self.flags = data.flags
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

            if flags is not None:
                self.flags = flags
                log.info("Overwriting NDData's current flags with specified flag")

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

            self.flags = flags
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

    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            if isinstance(value, NDUncertainty):
                class_name = self.__class__.__name__
                if self.unit and value._unit:
                    try:
                        scaling = value._unit.to(self.unit)
                    except UnitsError:
                        raise UnitsError('Cannot convert unit of uncertainty '
                                         'to unit of '
                                         '{0} object.'.format(class_name))
                    value.array *= scaling
                elif not self.unit and value._unit:
                    # Raise an error if uncertainty has unit and data does not
                    raise ValueError("Cannot assign an uncertainty with unit "
                                     "to {0} without "
                                     "a unit".format(class_name))
                self._uncertainty = value
                self._uncertainty.parent_nddata = self
            else:
                raise TypeError("Uncertainty must be an instance of a NDUncertainty object")
        else:
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

        if self.flags is not None:
            if isinstance(self.flags, np.ndarray):
                new_flags = self.flags[item]
                # flags setter expects an array, always
                if new_flags.shape == ():
                    new_flags = np.array(new_flags)
            elif isinstance(self.flags, FlagCollection):
                raise NotImplementedError('Slicing complex Flags is currently not implemented')
        else:
            new_flags = None

        if self.wcs is not None:
            raise NotImplementedError('Slicing for WCS is not currently implemented')
        else:
            new_wcs = None

        return self.__class__(new_data, uncertainty=new_uncertainty,
                              mask=new_mask, flags=new_flags, wcs=new_wcs,
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

        Notes
        -----
        Flags are set to None in the result.
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
                                mask=new_mask, flags=self.flags,
                                wcs=self.wcs,
                                meta=self.meta, unit=unit)

        return result

    read = classmethod(io_registry.read)
    write = io_registry.write


class NDArithmetic(object):
    """
    Mixin class to add arithmetic to an NDData object.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArithmetic` for an example.
    """
    def _arithmetic(self, operand, propagate_uncertainties, name, operation):
        """
        {name} another dataset (``operand``) to this dataset.

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
        This method requires the datasets to have identical WCS
        properties, equivalent units, and identical shapes. Flags and
        meta-data get set to None in the resulting dataset. The unit
        in the result is the same as the unit in ``self``. Uncertainties
        are propagated, although correlated errors are not supported
        by any of the built-in uncertainty classes.  If uncertainties
        are assumed to be correlated, a warning is issued by default
        (though this can be disabled via the
        ``astropy.nddata.conf.warn_unsupported_correlated``
        configuration item). Values masked in either dataset before
        the operation are masked in the resulting dataset.
        """
        from . import conf

        if self.wcs != operand.wcs:
            raise ValueError("WCS properties do not match")

        # get a sensible placeholder if .unit is None
        self_unit = self.unit or dimensionless_unscaled
        operand_unit = operand.unit or dimensionless_unscaled

        # This check could be rolled into the calculation of the result
        # but checking now avoids a potentially expensive calculation that
        # would fail anyway.
        try:
            # Quantity is designed to work with numpy ufuncs, but plain
            # Unit is not, so convert units to quantities
            result_unit = operation(1 * self_unit, 1 * operand_unit).unit
        except UnitsError:
            # current API raises ValueError in this case, not UnitError
            raise ValueError("operand units do not match")

        if self.shape != operand.shape:
            raise ValueError("operand shapes do not match")

        # Instead of manually scaling the operand data just let Quantity
        # handle it.
        # Order of the arguments is important here if the operation is
        # addition or subtraction and the units of the operands are different
        # but compatible. NDData follows the convention that Quantity follows
        # in that case, with the units of the first operand (i.e. self)
        # determining the units of the result.
        data = operation(self.data * self_unit, operand.data * operand_unit)

        result_unit = data.unit
        # If neither self nor operand had units then should return a result
        # that has no unit. A check that the result_unit is dimensionless
        # should not be necessary, but also does no harm.
        if self.unit is None and operand.unit is None:
            if result_unit is dimensionless_unscaled:
                result_unit = None
            else:
                raise ValueError("arithmetic result was not unitless even "
                                 "though operands were unitless")
        data = data.value
        new_wcs = deepcopy(self.wcs)

        # Call __class__ in case we are dealing with an inherited type
        result = self.__class__(data, uncertainty=None,
                                mask=None, flags=None, wcs=new_wcs,
                                meta=None, unit=result_unit)

        # Prepare to scale uncertainty if it is needed
        if operand.uncertainty:
            operand_uncert_value = operand.uncertainty.array

        # By this point the arithmetic has succeeded, so the input units were
        # consistent with each other given the operation.
        #
        # If the operation is addition or subtraction then need to ensure that
        # the uncertainty of the operand is the same units as the result
        # (which will be the same as self.unit).

        # The data ought to also be scaled in this case -- for addition of
        # a StdDevUncertainty this isn't really necessary but other
        # uncertainties when added/subtracted may depend on both the operand
        # uncertainty and the operand data.

        # Since the .unit.to methods create a copy, avoid the conversion
        # unless it is necessary.
        if (operation in [np.add, np.subtract] and
                self.unit != operand.unit):
            operand_data = operand.unit.to(self.unit, operand.data)
            if operand.uncertainty:
                operand_uncert_value = operand.unit.to(self.unit,
                                                       operand_uncert_value)
        else:
            operand_data = operand.data

        if operand.uncertainty:
            # Create a copy here in case this is returned as the uncertainty
            # of the result.
            operand_uncertainty = \
                operand.uncertainty.__class__(operand_uncert_value, copy=True)
        else:
            operand_uncertainty = None

        if propagate_uncertainties is None:
            result.uncertainty = None
        elif self.uncertainty is None and operand.uncertainty is None:
            result.uncertainty = None
        elif self.uncertainty is None:
            result.uncertainty = operand_uncertainty
        elif operand.uncertainty is None:
            result.uncertainty = self.uncertainty.__class__(self.uncertainty,
                                                            copy=True)
        else:  # both self and operand have uncertainties
            if (conf.warn_unsupported_correlated and
                (not self.uncertainty.support_correlated or
                 not operand.uncertainty.support_correlated)):
                log.info("The uncertainty classes used do not support the "
                         "propagation of correlated errors, so uncertainties"
                         " will be propagated assuming they are uncorrelated")
            operand_scaled = operand.__class__(operand_data,
                                               uncertainty=operand_uncertainty,
                                               unit=operand.unit,
                                               wcs=operand.wcs,
                                               mask=operand.mask,
                                               flags=operand.flags,
                                               meta=operand.meta)
            try:
                method = getattr(self.uncertainty, propagate_uncertainties)
                result.uncertainty = method(operand_scaled, result.data)
            except IncompatibleUncertaintiesException:
                raise IncompatibleUncertaintiesException(
                    "Cannot propagate uncertainties of type {0:s} with "
                    "uncertainties of type {1:s} for {2:s}".format(
                        self.uncertainty.__class__.__name__,
                        operand.uncertainty.__class__.__name__,
                        name))

        if self.mask is None and operand.mask is None:
            result.mask = None
        elif self.mask is None:
            result.mask = operand.mask.copy()
        elif operand.mask is None:
            result.mask = self.mask.copy()
        else:  # combine masks as for Numpy masked arrays
            result.mask = self.mask | operand.mask  # copy implied by operator

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


class NDDataArithmetic(NDArithmetic, NDData):
    """
    An ``NDData`` object with arithmetic. This class is functionally equivalent
    to ``NDData`` in astropy  versions prior to 1.0.
    """
    pass
