# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np

from ...units import dimensionless_unscaled, Quantity
from ... import log
from ...extern.six import string_types
from ..nduncertainty import NDUncertainty
from ...wcs import WCS
from ...config import ConfigAlias

WARN_UNSUPPORTED_CORRELATED = ConfigAlias(
    '0.4', 'WARN_UNSUPPORTED_CORRELATED', 'warn_unsupported_correlated',
    'astropy.nddata.nddata', 'astropy.nddata')

__all__ = ['NDArithmeticMixin']


class NDArithmeticMixin(object):
    """
    Mixin class to add arithmetic to an NDData object.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArray` for an example.
    """

    def arithmetic(self, operation, operand,
                   uncertainty_correlation=0, meta_kwds_operate=None,
                   meta_kwds_set=None, handle_uncertainties=True,
                   handle_mask=True, compare_wcs=True):
        """
        Base method which calculates the result of the arithmetic operation.

        This class determines the result of the arithmetic operation on the
        ``data`` with respect to their units and then forwards to other methods
        to calculate the other properties for the result (like uncertainty).

        Parameters
        ----------
        operation: `str`
            The operation that is performed on the `NDData`. Supported are
            ``addition``, ``subtraction``, ``multiplication`` and ``division``.

        operand: `NDData` instance or something that can be converted to one.
            The second operand in the operation.

        uncertainty_correlation: ``Number`` or `~numpy.ndarray`
            The correlation (rho) is defined between the uncertainties in
            ``sigma_AB = sigma_A * sigma_B * rho`` (Latex?). If
            ``handle_uncertainties`` is ``False`` this will be ignored. A value
            of ``0`` means uncorrelated (which is also the default)

        meta_kwds_operate: ``None`` or `list`
            If ``None`` no arithmetic meta operations are done. If this is a
            list then each element of the list will be interpreted as a
            keyword which should be arithmetically changed. In case both
            operands have not-empty meta properties the resulting meta keyword
            will be operand1 keyword operated upon operand2 keyword. If only
            one operand has a not-empty meta, the resulting keyword will be
            the keyword of the operand operated with the ``data`` (this is
            only allowed if the data is actually a scalar) of the other
            operand. (There will be more explanation and exampled sometime
            soon)

        meta_kwds_set: ``None`` or `dict`
            If ``None`` there are no keyword settings after the arithmetic
            meta keyword operations. If this is a `dict` each key/value pair
            in the ``meta_kwds_set`` will be added (or replaced if the keyword
            already exists) in the results meta.

        handle_uncertainties: `bool`
            If ``True`` it will be tried to propagate the uncertainties. If
            ``False`` the results uncertainty will be the uncertainty of the
            instance which called this method.

        handle_mask: `bool`
            If ``True`` the mask of the result will be processed by
            `_arithmetic_mask`. If ``False`` the result will not have a mask.

        check_wcs: `bool`
            If ``True`` the wcs of both operands will be tested for equality.
            If ``False`` this will not be tested. The result will contain the
            wcs information of the first operand except if the first operand
            has no wcs information and the second has, then the wcs information
            from the second operand will be kept.

        Returns
        -------
        result: `NDData` subclass
            The result as the same class as the first operand in the operation.

        Notes
        -----
        1. Due to the internal mechanics of `~astropy.units.Quantity` the
           resulting data might have units like ``km/m`` if you divided for
           example 100km by 5m. So this class has adopted this behaviour.
        2. Uncertainty propagation allows that the unit of the uncertainty
           differs from the unit of the data. In these cases uncertainty
           propagation *trys* to keep these units but sometimes this will
           not be possible. Maybe it is possible but that's not so easy right
           now and therefore if you use different units for uncertainty and
           data there might be cases where the uncertainty's unit may seem
           rather odd after an arithmetic operation.
        3. Internally all properties of the base NDData are determined and then
           the new subclass is created. So if you need to extend this class to
           handle other properties as well but keep the existing operations
           just call ``result=super('yourclass', self)._arithmetic`` and then
           append the additional properties via setters.
        4. tbc...
        """
        # Convert the operand to the same class this allows for arithmetic
        # operations with numbers, lists, numpy arrays, numpy masked arrays
        # astropy quantities, masked quantities and of other subclasses of
        # NDData
        operand = self.__class__(operand)

        # TODO: Check that both data have numeric dtypes otherwise we could
        # end up with some arithmetic operation on string types. (init only
        # enforces that the dtype is NOT objects)

        kwargs = {}
        # First check that the WCS allows the arithmetic operation
        kwargs['wcs'] = self._arithmetic_wcs(operand, compare_wcs)
        # Then calculate the resulting data (which can but not needs to be a
        # quantity)
        result = self._arithmetic_data(operation, operand)
        # Determine the other properties (if requested)
        if handle_uncertainties:
            kwargs['uncertainty'] = self._arithmetic_uncertainty(operation,
                                    operand, result, uncertainty_correlation)
        if handle_mask:
            kwargs['mask'] = self._arithmetic_mask(operand)
        kwargs['meta'] = self._arithmetic_meta(operation, operand,
                                    meta_kwds_operate, meta_kwds_set)
        # Wrap the individual results into a new instance of the same class.
        return self.__class__(result, **kwargs)

    def _arithmetic_data(self, operation, operand):
        """
        Calculate the resulting data

        Parameters
        ----------
        operation: `str`
            see `NDArithmeticMixin.arithmetic` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        decompose_unit: `bool`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        Returns
        -------
        result_data: `~numpy.ndarray` or `~astropy.units.Quantity`
            if both operands had no unit the resulting data is a simple numpy
            array, but if any of the operands had a unit the return is a
            Quantity.
        """

        # Find the right function for that operation... (numpy ufuncs)
        if operation == 'addition':
            operator = np.add
        elif operation == 'subtraction':
            operator = np.subtract
        elif operation == 'multiplication':
            operator = np.multiply
        elif operation == 'division':
            # Unfortunatly the future import does not implicitly call the
            # true division of numpy so explicitly call np.true_divide instead
            # of np.divide (until python2 support is dropped ... just joking)
            operator = np.true_divide
        else:
            raise ValueError('Unsupported operation')

        # Do the calculation with or without units
        if self.unit is None and operand.unit is None:
            result = operator(self.data, operand.data)
        elif self.unit is None:
            result = operator(self.data * dimensionless_unscaled,
                              operand.data * operand.unit)
        elif operand.unit is None:
            result = operator(self.data * self.unit,
                              operand.data * dimensionless_unscaled)
        else:
            result = operator(self.data * self.unit,
                              operand.data * operand.unit)

        return result

    def _arithmetic_uncertainty(self, operation, operand, result, correlation):
        """
        Calculate the resulting uncertainty.

        Parameters
        ----------
        operation: `str`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        result: `~astropy.units.Quantity` or `~numpy.ndarray`
            The result of :meth:`NDArithmeticMixin._arithmetic_data`.

        correlation: `Number` or `~numpy.ndarray`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        Returns
        -------
        result_uncertainty: `NDUncertainty` subclass instance or None
            The resulting uncertainty already saved in the same `NDUncertainty`
            subclass that ``self`` had (or ``operand`` if self had no
            uncertainty). ``None`` only if both had no uncertainty.
        """

        # Make sure these uncertainties are NDUncertainties so this kind of
        # propagation is possible.
        if (self.uncertainty is not None and
                        not isinstance(self.uncertainty, NDUncertainty)):
            raise TypeError("Uncertainty propagation is only defined for "
                            "subclasses of NDUncertainty.")
        if (operand.uncertainty is not None and
                        not isinstance(operand.uncertainty, NDUncertainty)):
            raise TypeError("Uncertainty propagation is only defined for "
                            "subclasses of NDUncertainty.")

        # Now do the uncertainty propagation
        if self.uncertainty is None and operand.uncertainty is None:
            # Neither has uncertainties so the result should have none.
            return None
        elif self.uncertainty is None:
            # Create a temporary uncertainty to allow uncertainty propagation
            # to yield the correct results. (issue #4152)
            self.uncertainty = operand.uncertainty.__class__(None)
            result_uncert = self.uncertainty.propagate(operation, operand,
                                                       result, correlation)
            # Delete the temporary uncertainty again.
            self.uncertainty = None
            return result_uncert

        elif operand.uncertainty is None:
            # As with self.uncertainty is None but the other way around.
            operand.uncertainty = self.uncertainty.__class__(None)
            result_uncert = self.uncertainty.propagate(operation, operand,
                                                       result, correlation)
            operand.uncertainty = None
            return result_uncert

        else:
            # Both have uncertainties so just propagate.
            return self.uncertainty.propagate(operation, operand, result,
                                              correlation)

    def _arithmetic_mask(self, operand):
        """
        Calculate the resulting mask

        This is implemented as the piecewise ``or`` operation if both have a
        mask.

        Parameters
        ----------
        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        Returns
        -------
        result_mask: `bool` or `~numpy.ndarray` or None
            The mask created by a piecewise ``or`` operation. If only one mask
            was present this mask is returned instead. If neither had a mask
            ``None`` is returned. This kind of mask combination is only
            possible if both operands have a mask and this mask is either
            a boolean or an `~numpy.ndarray` containing booleans.
        """

        # If only one mask is present we need not bother about any type checks
        if self.mask is None and operand.mask is None:
            return None
        elif self.mask is None:
            # Make a copy so there is no reference in the result.
            return deepcopy(operand.mask)
        elif operand.mask is None:
            return deepcopy(self.mask)
        else:
            # Make sure these masks are boolean or ndarrays of bools so this
            # kind of mask operation is possible.
            if self.mask is not None:
                if isinstance(self.mask, bool):
                    pass
                elif (isinstance(self.mask, np.ndarray) and
                                                self.mask.dtype == 'bool'):
                    pass
                else:
                    raise TypeError("Mask arithmetics is only defined for "
                                    "boolean or arrays of booleans.")
            if operand.mask is not None:
                if isinstance(operand.mask, bool):
                    pass
                elif (isinstance(operand.mask, np.ndarray) and
                                                operand.mask.dtype == 'bool'):
                    pass
                else:
                    raise TypeError("Mask arithmetics is only defined for "
                                    "boolean or arrays of booleans.")

            # Now lets calculate the resulting mask (operation enforces copy)
            return self.mask | operand.mask

    def _arithmetic_wcs(self, operand, compare_wcs):
        """
        Calculate the resulting wcs.

        There is actually no calculation involved but it is a good place to
        compare wcs information of both operands. This is currently not working
        properly with `~astropy.wcs.WCS` (which is the suggested class for
        storing as wcs property) but it will not break it neither.

        Parameters
        ----------
        operand: `NDData` instance or subclass
            The second operand wrapped in an instance of the same class as
            self.

        compare_wcs: `bool`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        Returns
        -------
        result_wcs: any type
            The WCS information of the first operand is copyied if it is set.
            If it was not set the second operand will be tryed. If that one
            had no WCS too ``None`` is returned.
        """

        # ok, not really arithmetics but we need to check which wcs makes sense
        # for the result and this is an ideal place to compare the two WCS,
        # too.

        if self.wcs is None and operand.wcs is None:
            # None has any wcs information so the result should have None, no
            # need to compare the WCS here.
            return None
        elif self.wcs is None:
            # Only the second operand had WCS information so take this and
            # issue a warning if comparison of WCS is requested.
            if compare_wcs:
                log.info("Only the second operand had some WCS information.")
            return deepcopy(operand.wcs)
        elif operand.wcs is None:
            # Same as last case but reversed
            if compare_wcs:
                log.info("Only the first operand had some WCS information.")
            return deepcopy(self.wcs)
        else:
            # Both have WCS information. If requested compare the WCS and
            # issue a warning if they are not equal
            # TODO: Currently (astropy 1.04) this is a bit tricky because no
            # two astropy.wcs.WCS are really identical even if they are copied
            # so there is no reason to raise an Exception here, but issue the
            # warning nevertheless. Really unfortunate... I currently worked
            # around it so the warning differs if both are astropy.wcs.WCS
            # instances.
            if compare_wcs and self.wcs != operand.wcs:
                if isinstance(self.wcs, WCS) and isinstance(operand.wcs, WCS):
                    log.info("Astropy.WCS does not allow equality checks.")
                else:
                    log.info("WCS is not equal.")
            return deepcopy(self.wcs)

    def _arithmetic_meta(self, operation, operand, meta_kwds_operate,
                         meta_kwds_set):
        """
        Calculate the resulting meta.

        There are two ways of altering the resulting meta, one is by
        applying the operation that was used on the data also on a keyword or
        by manually setting the keyword. The operation on keywords is done
        *before* setting the manually defined keywords.

        Parameters
        ----------
        operation: `str`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        meta_kwds_operate: ``None`` or `list`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        meta_kwds_set: ``None`` or `dict`
            see :meth:`NDArithmeticMixin.arithmetic` parameter description.

        Returns
        -------
        result_meta: `dict`-like
            The class of the result depends on which meta is set. If only the
            second operand has a not empty meta property the return has the
            same class as the second operands meta. In all other cases the
            result has the same class as the first operands meta.
        """

        # Again this is only partly arithmetics. But it should be best to
        # keep everything relating meta in one operation.

        # Starting by creating the basic meta for the result
        if len(self.meta) == 0 and len(operand.meta) > 0:
            # The only case in which we would like the meta from the second
            # operand is if the first operand has no meta and the second one
            # contains at least one meta element.
            # There is no use-case where merging the meta of both operands
            # makes sense. In case one might really set the results meta
            # information to match the one of the second operand is to put a
            # reference or copy inside the parameter that sets the keywords
            # in here.
            result_meta = deepcopy(operand.meta)
        else:
            result_meta = deepcopy(self.meta)

        # Now get to work, first we check if there are any keywords which
        # should be arithmetically operated upon.
        if meta_kwds_operate is not None and len(meta_kwds_operate) > 0:

            # Get the operator function for these operation (use numpy ufuncs)
            if operation == 'addition':
                operator = np.add
            elif operation == 'subtraction':
                operator = np.subtract
            elif operation == 'multiplication':
                operator = np.multiply
            elif operation == 'division':
                operator = np.true_divide
            else:
                raise ValueError('Unsupported operation')

            # There are two possibilities, either operate keyword of operand1
            # with the same keyword of operand2. Like adding exposuretimes in
            # addition. So if both operands have meta information assume
            # this kind is wanted. The other possibility is to operate
            # the keyword of the only meta (since if both have meta we try
            # the other way) with the number of the other operand. Like
            # exposure time is multiplied with 5 if the data is multiplied by
            # 5.
            if len(result_meta) == 0:
                # No meta present. Nothing to do but issueing a warning that
                # it was not possible.
                log.info("No meta is present so no operation on them is"
                         "possible.")
            elif len(self.meta) > 0 and len(operand.meta) > 0:
                # Both have meta information so we assume that we operate
                # keyword on keyword.
                for i in meta_kwds_operate:
                    result_meta[i] = operator(self.meta[i], operand.meta[i])
            elif len(self.meta) > 0:
                # We have one meta, now we need to operate the keywords with
                # the second operands data. But that's only well defined if
                # that is a scalar... (ok, at least if we don't want to allow
                # the meta to contain numpy arrays :-) )
                if operand.data.size != 1:
                    log.info("The second operand is not a scalar. Cannot"
                             "operate keywords with arrays")
                else:
                    for i in meta_kwds_operate:
                        result_meta[i] = float(
                                        operator(self.meta[i], operand.data))
            elif len(operand.meta) > 0:
                # Same thing but reversed
                if self.data.size != 1:
                    log.info("The first operand is not a scalar. Cannot"
                             "operate keywords with arrays")
                else:
                    for i in meta_kwds_operate:
                        result_meta[i] = float(
                                        operator(self.data, operand.meta[i]))

        # Lastly, we need to set the keywords that should be set in the results
        # meta
        if meta_kwds_set is not None and len(meta_kwds_set) > 0:
            for i in meta_kwds_set:
                result_meta[i] = meta_kwds_set[i]

        return result_meta


    def add(self, operand, **kwargs):
        return self.arithmetic("addition", operand, **kwargs)

    def subtract(self, operand, **kwargs):
        return self.arithmetic("subtraction", operand, **kwargs)

    def multiply(self, operand, **kwargs):
        return self.arithmetic("multiplication", operand, **kwargs)

    def divide(self, operand, **kwargs):
        return self.arithmetic("division", operand, **kwargs)


    @classmethod
    def addition(cls, operand1, operand2, **kwargs):
        # Convert the first operand to the implicit passed class (cls)
        # this allows for reverse operations.
        operand1 = cls(operand1)
        # Call the instance function for addition to proceed
        return operand1.add(operand2, **kwargs)

    @classmethod
    def subtraction(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.subtract(operand2, **kwargs)

    @classmethod
    def multiplication(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.multiply(operand2, **kwargs)

    @classmethod
    def division(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.divide(operand2, **kwargs)

    if isinstance(arithmetic.__doc__, string_types):
        doc = """
        {name} another dataset (``operand``) to this dataset.

        Parameters
        ----------
        operand : `~astropy.nddata.NDData`-like
            The second operand in the operation a {operator} b

        kwargs :
            anything that can be passed as optional parameter to
            :meth:`NDArithmeticMixin.arithmetic`.

        Returns
        -------
        result : `~astropy.nddata.NDData`
            The resulting dataset

        See also
        --------
        :meth:`NDArithmeticMixin.arithmetic`
        """
        add.__doc__ = doc.format(name="Add", operator="+")
        subtract.__doc__ = doc.format(name="Subtract", operator="-")
        multiply.__doc__ = doc.format(name="Multiply", operator="*")
        divide.__doc__ = doc.format(name="Divide", operator="/")
        del doc
