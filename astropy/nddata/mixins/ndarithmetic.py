# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np

from ...units import dimensionless_unscaled
from ... import log
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

    def _arithmetic(self, operation, operand,
                   uncertainty_correlation=0, meta_kwds_operate=None,
                   meta_kwds_set=None, propagate_uncertainties=True,
                   handle_mask=True, handle_meta=True, compare_wcs=True):
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
            see :meth:`NDArithmeticMixin.add`

        uncertainty_correlation: ``Number`` or `~numpy.ndarray`, optional
            see :meth:`NDArithmeticMixin.add`

        meta_kwds_operate: ``None`` or `list`, optional
            see :meth:`NDArithmeticMixin.add`

        meta_kwds_set: ``None`` or `dict`, optional
            see :meth:`NDArithmeticMixin.add`

        propagate_uncertainties: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_mask: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_mata: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        check_wcs: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        Returns
        -------
        result: `~numpy.ndarray` or `~astropy.units.Quantity`
            The resulting data as array (in case both operands were without
            unit) or as quantity if at least one had a unit.

        kwargs: `dict`
            kwargs to create a new instance of the same class as the first
            operand was. The calling method or function is responsible for
            creating the instance

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
        if compare_wcs is None:
            kwargs['wcs'] = None
        elif not compare_wcs:
            if self.wcs is None:
                kwargs['wcs'] = operand.wcs
            else:
                kwargs['wcs'] = self.wcs
        else:
            kwargs['wcs'] = self._arithmetic_wcs(operand, compare_wcs)
        # Then calculate the resulting data (which can but not needs to be a
        # quantity)
        result = self._arithmetic_data(operation, operand)
        # Determine the other properties
        if propagate_uncertainties is None:
            kwargs['uncertainty'] = None
        elif not propagate_uncertainties:
            if self.uncertainty is None:
                kwargs['uncertainty'] = operand.uncertainty
            else:
                kwargs['uncertainty'] = self.uncertainty
        else:
            kwargs['uncertainty'] = self._arithmetic_uncertainty(operation,
                                    operand, result, uncertainty_correlation)
        if handle_mask is None:
            kwargs['mask'] = None
        elif not handle_mask:
            if self.mask is None:
                kwargs['mask'] = operand.mask
            else:
                kwargs['mask'] = self.mask
        else:
            kwargs['mask'] = self._arithmetic_mask(operand)

        if handle_meta is None:
            kwargs['meta'] = None
        elif not handle_meta:
            if len(self.meta) == 0:
                kwargs['meta'] = operand.meta
            else:
                kwargs['meta'] = self.meta
        else:
            kwargs['meta'] = self._arithmetic_meta(operation, operand,
                                        meta_kwds_operate, meta_kwds_set)

        # Wrap the individual results into a new instance of the same class.
        return result, kwargs

    def _arithmetic_data(self, operation, operand):
        """
        Calculate the resulting data

        Parameters
        ----------
        operation: `str`
            see `NDArithmeticMixin._arithmetic` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

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
            see :meth:`NDArithmeticMixin.add` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        result: `~astropy.units.Quantity` or `~numpy.ndarray`
            The result of :meth:`NDArithmeticMixin._arithmetic_data`.

        correlation: `Number` or `~numpy.ndarray`
            see :meth:`NDArithmeticMixin.add` parameter description.

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
            see :meth:`NDArithmeticMixin.add` parameter description.

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
            see :meth:`NDArithmeticMixin.add` parameter description.

        operand: `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        meta_kwds_operate: ``None`` or `list`
            see :meth:`NDArithmeticMixin.add` parameter description.

        meta_kwds_set: ``None`` or `dict`
            see :meth:`NDArithmeticMixin.add` parameter description.

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
        result, kwargs = self._arithmetic("addition", operand, **kwargs)
        return self.__class__(result, **kwargs)

    def subtract(self, operand, **kwargs):
        result, kwargs = self._arithmetic("subtraction", operand, **kwargs)
        return self.__class__(result, **kwargs)

    def multiply(self, operand, **kwargs):
        result, kwargs = self._arithmetic("multiplication", operand, **kwargs)
        return self.__class__(result, **kwargs)

    def divide(self, operand, **kwargs):
        result, kwargs = self._arithmetic("division", operand, **kwargs)
        return self.__class__(result, **kwargs)


    @classmethod
    def addition(cls, operand1, operand2, **kwargs):
        """
        Like :meth:`NDArithmeticMixin.add` you can add two operands.

        This classmethod allows for inverse operations and handling different
        `NDData` subclasses.

        Parameters
        ----------
        operand1: `NDData`-like or convertable to `NDData`
            the first operand in the operation.

        operand2: `NDData`-like or convertable to `NDData`
            the second operand in the operation.

        kwargs:
            see :meth:`NDArithmeticMixin.add`

        Returns
        -------
        result: `NDData`-like
            The result of the operation.

        Notes
        -----
        This method allows performing arithmetic operations where the first
        operand *should* not define the results class or where the inverse
        operation wouldn't be possible otherwise.
        """
        # Convert the first operand to the implicit passed class (cls)
        # this allows for reverse operations.
        operand1 = cls(operand1)
        # Call the instance function for addition to proceed
        return operand1.add(operand2, **kwargs)

    @classmethod
    def subtraction(cls, operand1, operand2, **kwargs):
        """
        Like :meth:`NDArithmeticMixin.subtract` you can subtract two operands.

        This classmethod allows for inverse operations and handling different
        `NDData` subclasses.

        Parameters
        ----------
        operand1: `NDData`-like or convertable to `NDData`
            the first operand in the operation.

        operand2: `NDData`-like or convertable to `NDData`
            the second operand in the operation.

        kwargs:
            see :meth:`NDArithmeticMixin.subtract`

        Returns
        -------
        result: `NDData`-like
            The result of the operation.

        Notes
        -----
        This method allows performing arithmetic operations where the first
        operand *should* not define the results class or where the inverse
        operation wouldn't be possible otherwise.
        """
        operand1 = cls(operand1)
        return operand1.subtract(operand2, **kwargs)

    @classmethod
    def multiplication(cls, operand1, operand2, **kwargs):
        """
        Like :meth:`NDArithmeticMixin.multiply` you can multiply two operands.

        This classmethod allows for inverse operations and handling different
        `NDData` subclasses.

        Parameters
        ----------
        operand1: `NDData`-like or convertable to `NDData`
            the first operand in the operation.

        operand2: `NDData`-like or convertable to `NDData`
            the second operand in the operation.

        kwargs:
            see :meth:`NDArithmeticMixin.multiply`

        Returns
        -------
        result: `NDData`-like
            The result of the operation.

        Notes
        -----
        This method allows performing arithmetic operations where the first
        operand *should* not define the results class or where the inverse
        operation wouldn't be possible otherwise.
        """
        operand1 = cls(operand1)
        return operand1.multiply(operand2, **kwargs)

    @classmethod
    def division(cls, operand1, operand2, **kwargs):
        """
        Like :meth:`NDArithmeticMixin.divide` you can divide two operands.

        This classmethod allows for inverse operations and handling different
        `NDData` subclasses.

        Parameters
        ----------
        operand1: `NDData`-like or convertable to `NDData`
            the first operand in the operation.

        operand2: `NDData`-like or convertable to `NDData`
            the second operand in the operation.

        kwargs:
            see :meth:`NDArithmeticMixin.divide`

        Returns
        -------
        result: `NDData`-like
            The result of the operation.

        Notes
        -----
        This method allows performing arithmetic operations where the first
        operand *should* not define the results class or where the inverse
        operation wouldn't be possible otherwise.
        """
        operand1 = cls(operand1)
        return operand1.divide(operand2, **kwargs)

    if True:
        doc = """
        {name} another dataset (``operand``) to this dataset.

        Parameters
        ----------
        operand : `~astropy.nddata.NDData`-like
            The second operand in the operation a {operator} b

        operand: `NDData` instance or something that can be converted to one.
            The second operand in the operation.

        uncertainty_correlation: ``Number`` or `~numpy.ndarray`, optional
            The correlation (rho) is defined between the uncertainties in
            ``sigma_AB = sigma_A * sigma_B * rho`` (Latex?). If
            ``propagate_uncertainties`` is *not* ``True`` this will be ignored.
            A value of ``0`` means uncorrelated (which is also the default)

        meta_kwds_operate: ``None`` or `list`, optional
            If ``None`` no arithmetic meta operations are done. If this is a
            list then each element of the list will be interpreted as a
            keyword which should be arithmetically changed. In case both
            operands have not-empty meta properties the resulting meta keyword
            will be operand1 keyword operated upon operand2 keyword. If only
            one operand has a not-empty meta, the resulting keyword will be
            the keyword of the operand operated with the ``data`` (this is
            only allowed if the data is actually a scalar) of the other
            operand. (There will be more explanation and exampled sometime
            soon). This parameter will be ignored if ``handle_meta`` is *not*
            ``True``. Default is ``None``.

        meta_kwds_set: ``None`` or `dict`, optional
            If ``None`` there are no keyword settings after the arithmetic
            meta keyword operations. If this is a `dict` each key/value pair
            in the ``meta_kwds_set`` will be added (or replaced if the keyword
            already exists) in the results meta. This parameter will be ignored
            if ``handle_meta`` is *not* ``True``. Default is ``None``.

        propagate_uncertainties: `bool` or ``None``, optional
            If ``None`` the result will have no uncertainty. If ``False`` the
            result will have a copied version of the first operand (or the
            second if the first had no uncertainty). In case this is ``True``
            the result will have a correctly propagated uncertainty from the
            uncertainties of the operands. Default is ``True``.

        handle_mask: `bool` or ``None``, optional
            If ``None`` the result will have no mask. If ``False`` the result
            will have a copied version of the mask of the first operand (or
            if the first operand has no mask then the one from the second is
            taken). If ``True`` the masks of both operands will be taken into
            account and combined (by a bitwise ``or`` operation). Default is
            ``True``.

        handle_mata: `bool` or ``None``, optional
            If ``None`` the result will have no meta. If ``False`` the result
            will have a copied version of the meta of the first operand (or
            if the first operand has no meta then the one from the second is
            taken). If ``True`` the meta of both operands will be changed
            depending on ``meta_kwds_operate`` and ``meta_kwds_set``. Default
            is ``True``.

        check_wcs: `bool` or ``None``, optional
            If ``None`` the result will have no wcs and no comparison between
            the wcs of the operands is made. If ``False`` the result will have
            the wcs of the first operand (or second if the first had ``None``).
            If ``True`` the resulting mask will be like in the case of
            ``False`` but the wcs information is checked if it is the same
            for both operands (in some cases it would be bad to apply any
            arithmetic operation on datasets which have different wcs
            informations). If they do not match a warning is issues. Default is
            ``True``.

        Returns
        -------
        result : `~astropy.nddata.NDData`-like
            The resulting dataset

        Notes
        -----
        1. It is not tried to decompose the units, mainly due to the internal
           mechanics of `~astropy.units.Quantity` the resulting data might have
           units like ``km/m`` if you divided for example 100km by 5m. So this
           Mixin has adopted this behaviour.

        2. Uncertainty propagation allows that the unit of the uncertainty
           differs from the unit of the data. In these cases uncertainty
           propagation *trys* to keep these units but sometimes this will
           not be possible. Maybe it is possible but that's not so easy right
           now and therefore if you use different units for uncertainty and
           data there might be cases where the uncertainty's unit may seem
           rather odd after an arithmetic operation.

        See also
        --------
        {other}
        """
        add.__doc__ = doc.format(name="Add", operator="+",
                            other=":meth:`NDArithmeticMixin.addition`")
        subtract.__doc__ = doc.format(name="Subtract", operator="-",
                            other=":meth:`NDArithmeticMixin.subtraction`")
        multiply.__doc__ = doc.format(name="Multiply", operator="*",
                            other=":meth:`NDArithmeticMixin.multiplication`")
        divide.__doc__ = doc.format(name="Divide", operator="/",
                            other=":meth:`NDArithmeticMixin.division`")

        examples_for_main_documentation = """
        For example::

            >>> from astropy.nddata import NDData, NDArithmeticMixin, NDSlicingMixin
            >>> class NDDataWithMath(NDArithmeticMixin, NDData): pass
            >>> nd = NDDataWithMath([1,2,3], unit='meter')
            >>> nd_inv = nd.division(1, nd)
            >>> nd.__class__.__name__
            'NDDataWithMath'
            >>> nd.data
            array([ 1.        ,  0.5       ,  0.33333333])
            >>> nd.unit
            Unit("1 / m")

        This method also allows that the result of unrelated objects is tried
        with what this class implements as arithmetics and give a result that
        is the same class as the class that called the method.

            >>> nd2 = nd.multiplication([1,2],3)
            >>> nd2.__class__.__name__
            'NDDataWithMath'
            >>> nd2.data
            array([3, 6])

        Since this is a classmethod we don't need an instance to use them:

            >>> nd3 = NDDataWithMath.subtraction(5, NDData([1,2,3]))
            >>> nd3.__class__.__name__
            'NDDataWithMath'
            >>> nd3.data
            array([4, 3, 2])

        And it allows to handle arithmetics with different subclasses.

            >>> class NDDataWithMathAndSlicing(NDSlicingMixin, NDArithmeticMixin, NDData): pass
            >>> nd = NDDataWithMath([5,5,5])
            >>> nd2 = NDDataWithMathAndSlicing([3,2,5])
            >>> nd3 = nd2.addition(nd, nd2)
            >>> nd3.__class__.__name__
            'NDDataWithMathAndSlicing'
            >>> nd3.data
            array([ 8,  7, 10])
        """
        del doc
