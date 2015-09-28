# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np

from ... units import dimensionless_unscaled, UnitsError
from ... import log
from ...extern.six import string_types
from ..nduncertainty import IncompatibleUncertaintiesException

__all__ = ['NDArithmeticMixin']


class NDArithmeticMixin(object):
    """
    Mixin class to add arithmetic to an NDData object.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArray` for an example.
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

        from .. import conf

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

        if self.data.shape != operand.data.shape:
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
                if hasattr(operand.uncertainty, '_unit') and \
                                operand.uncertainty._unit is not None:
                    operand_uncert_value *= operand.uncertainty._unit
                else:
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
        else:
            if self.uncertainty is None:
                self.uncertainty = operand.uncertainty.__class__(None)
            elif operand_uncertainty is None:
                operand_uncertainty = self.uncertainty.__class__(None)
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

    if isinstance(_arithmetic.__doc__, string_types):
        add.__doc__ = _arithmetic.__doc__.format(name="Add", operator="+")

    def subtract(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_subtract"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "subtraction", np.subtract)

    if isinstance(_arithmetic.__doc__, string_types):
        subtract.__doc__ = _arithmetic.__doc__.format(name="Subtract",
                                                      operator="-")

    def multiply(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_multiply"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "multiplication", np.multiply)

    if isinstance(_arithmetic.__doc__, string_types):
        multiply.__doc__ = _arithmetic.__doc__.format(name="Multiply",
                                                      operator="*")

    def divide(self, operand, propagate_uncertainties=True):
        if propagate_uncertainties:
            propagate_uncertainties = "propagate_divide"
        else:
            propagate_uncertainties = None
        return self._arithmetic(
            operand, propagate_uncertainties, "division", np.divide)

    if isinstance(_arithmetic.__doc__, string_types):
        divide.__doc__ = _arithmetic.__doc__.format(name="Divide",
                                                    operator="/")
