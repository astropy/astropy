# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import absolute_import, division, print_function, unicode_literals

from copy import deepcopy

import numpy as np

from .. units import dimensionless_unscaled, UnitsError
from .. import log
from .nddata import NDData
from .nduncertainty import IncompatibleUncertaintiesException, NDUncertainty

__all__ = ['NDArithmetic', 'NDDataArithmetic']


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
            defined by the class used for the ``uncertainty`` attribute.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        Notes
        -----
        This method requires the datasets to have identical WCS
        properties, equivalent units, and identical shapes.
        Meta-data get set to None in the resulting dataset. The unit
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
                                mask=None, wcs=new_wcs,
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


class NDSlicing(object):
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

        if self.wcs is not None:
            raise NotImplementedError('Slicing for WCS is not currently implemented')
        else:
            new_wcs = None

        return self.__class__(new_data, uncertainty=new_uncertainty,
                              mask=new_mask, wcs=new_wcs,
                              meta=self.meta, unit=self.unit)


class NDDataArithmetic(NDArithmetic, NDSlicing, NDData):
    """
    An ``NDData`` object with arithmetic. This class is functionally equivalent
    to ``NDData`` in astropy  versions prior to 1.0.
    """

    def __init__(self, *arg, **kwd):
        # Initialize with the parent...
        super(NDDataArithmetic, self).__init__(*arg, **kwd)

        # ...then reset uncertainty to force it to go through the
        # setter logic below. In base NDData all that is done is to
        # set self._uncertainty to whatever uncertainty is passed in.
        self.uncertainty = self._uncertainty

        # Same thing for mask...
        self.mask = self._mask

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
                if self.unit and value._unit:
                    try:
                        scaling = (1 * value._unit).to(self.unit)
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
            if mask.shape != self.shape:
                raise ValueError("dimensions of mask do not match data")
            else:
                self._mask = mask
        else:
            # internal representation should be one numpy understands
            self._mask = np.ma.nomask

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
