# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np

from ..nduncertainty import NDUncertainty
from ... import log
from ...units import dimensionless_unscaled
#from ...utils.decorators import sharedmethod
from ...wcs import WCS
from ...config import ConfigAlias
from ...extern import six

WARN_UNSUPPORTED_CORRELATED = ConfigAlias(
    '0.4', 'WARN_UNSUPPORTED_CORRELATED', 'warn_unsupported_correlated',
    'astropy.nddata.nddata', 'astropy.nddata')

__all__ = ['NDArithmeticMixin']


# TODO: General question of implementing an Operator class #3861? (embray)
# Maybe a bit too big for these changes.


# TODO: Delete this and rebase if #4242 is merged.
def tmp_deco(docstring, *args, **kwargs):
    def set_docstring(func):
        if not isinstance(docstring, six.string_types):
            doc = docstring.__doc__
        elif docstring != 'self':
            doc = docstring
        else:
            doc = func.__doc__
            func.__doc__ = None
        if not doc:
            raise ValueError
        kwargs['original_doc'] = func.__doc__ or ''
        func.__doc__ = doc.format(*args, **kwargs)
        return func
    return set_docstring

# TODO: Maybe not nice to pollute globals but I felt the same way about
# polluting the class namespace and this way potential subclasses may or refuse
# to pull this docstring into their class as well.

# Docstring templates for add, subtract, multiply, divide methods.
_arit_doc = """
    Performs {name} by evaluating ``self`` {op} ``operand``.

    Parameters
    ----------
    operand: `NDData`-like instance or convertable to one.
        The second operand in the operation.

    uncertainty_correlation: ``Number`` or `~numpy.ndarray`, optional
        The correlation (rho) is defined between the uncertainties in
        ``sigma_AB = sigma_A * sigma_B * rho`` (Latex?). If
        ``propagate_uncertainties`` is *not* ``True`` this will be ignored.
        A value of ``0`` means uncorrelated (which is also the default).
        TODO: TeX of formula?

    meta_kwds_operate: ``None`` or `list`, optional
        If ``None`` no arithmetic meta operations are done. If this is a
        list then each element of the list will be interpreted as a
        keyword which should be arithmetically changed. In case both
        operands have not-empty meta properties the resulting meta keyword
        will be operand1-keyword operated upon operand2-keyword. If only
        one operand has a not-empty meta, the resulting keyword will be
        the keyword of the operand operated with the ``data`` (this is
        only allowed if the data is actually a scalar) of the other
        operand. This parameter will be ignored if ``handle_meta`` is *not*
        ``True``. Default is ``None``.
        TODO: Create explanations and examples for this

    meta_kwds_set: ``None`` or `dict`, optional
        If ``None`` there are no special keyword settings after the arithmetic
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

    handle_meta: `bool` or ``None``, optional
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
        If ``True`` the resulting wcs will be like in the case of
        ``False`` but the wcs information is checked if it is the same
        for both operands (in some cases it would be bad to apply any
        arithmetic operation on datasets which have different wcs
        informations). If they do not match a warning is raised. Default is
        ``True``.
        TODO: Raise Exception instead of warning? This can lead to problems
        concerning astropy.wcs.WCS which has not an __eq__ method and
        therefore falls back to ``id`` equality which will be never fulfilled
        except for reference-copies.

    Returns
    -------
    result : `~astropy.nddata.NDData`-like
        The resulting dataset

    See also
    --------
    :meth:`NDArithmeticMixin.ic_{name}`
    """

_arit_cls_doc = """
    Like :meth:`NDArithmeticMixin.{0}` you can {0} two operands.

    This method is avaiable as classmethod in order to allow arithmetic
    operations between arbitary objects as long as they are convertable to
    the class that called the method. Therefor the name prefix ``ic_`` which
    stands for ``interclass`` operations.

    Parameters
    ----------
    operand1: `NDData`-like or convertable to `NDData`
        the first operand in the operation.

    operand2: `NDData`-like or convertable to `NDData`
        the second operand in the operation.

    kwargs:
        see :meth:`NDArithmeticMixin.{0}`

    Returns
    -------
    result: `NDData`-like
        The result of the operation. The class of the result is the class from
        where the method was invoked.
    """


class NDArithmeticMixin(object):
    """
    Mixin class to add arithmetic to an NDData object.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArray` for an example.

    Notes
    -----
    This class only aims at covering the most common cases so there are certain
    restrictions on the saved attributes::

        - ``uncertainty`` : has to be something that has a `NDUncertainty`-like
          interface for uncertainty propagation
        - ``mask`` : has to be something that can be used by a bitwise ``or``
          operation.
        - ``wcs`` : has to implement a way of comparing with ``=`` to allow
          the operation.

    But there is a workaround that allows to disable handling a specific
    attribute and to simply set the results attribute to ``None`` or to
    copy the existing attribute (and neglecting the other).
    For example for uncertainties not representing an `NDUncertainty`-like
    interface you can alter the ``propagate_uncertainties`` parameter in
    :meth:`NDArithmeticMixin.add`. ``None`` means that the result will have no
    uncertainty, ``False`` means it takes the uncertainty of the first operand
    (if this does not exist from the second operand) as the results
    uncertainty. This behaviour is also explained in the help-page for the
    different arithmetic operations.

    Notes
    -----
    1. It is not tried to decompose the units, mainly due to the internal
       mechanics of `~astropy.units.Quantity`, so the resulting data might have
       units like ``km/m`` if you divided for example 100km by 5m. So this
       Mixin has adopted this behaviour.

    Examples
    --------

    TODO: Just here because the examples fit more into the main documentation
    than in here.

    For example::

        >>> from astropy.nddata import *
        >>> class NDDataWithMath(NDArithmeticMixin, NDData): pass
        >>> nd = NDDataWithMath([1,2,3], unit='meter')
        >>> nd_inv = nd.ic_division(1, nd)
        >>> nd.__class__.__name__
        'NDDataWithMath'
        >>> nd.data
        array([ 1.        ,  0.5       ,  0.33333333])
        >>> nd.unit
        Unit("1 / m")

    This method also allows that the result of unrelated objects is tried
    with what this class implements as arithmetics and give a result that
    is the same class as the class that called the method.

        >>> nd2 = nd.ic_multiplication([1,2],3)
        >>> nd2.__class__.__name__
        'NDDataWithMath'
        >>> nd2.data
        array([3, 6])

    Since this is a classmethod we don't need an instance to use them:

        >>> nd3 = NDDataWithMath.ic_subtraction(5, NDData([1,2,3]))
        >>> nd3.__class__.__name__
        'NDDataWithMath'
        >>> nd3.data
        array([4, 3, 2])

    And it allows to handle arithmetics with different subclasses.

        >>> class NDDataWithMathAndSlicing(NDSlicingMixin,NDArithmeticMixin, NDData): pass
        >>> nd = NDDataWithMath([5,5,5])
        >>> nd2 = NDDataWithMathAndSlicing([3,2,5])
        >>> nd3 = nd2.ic_addition(nd, nd2)
        >>> nd3.__class__.__name__
        'NDDataWithMathAndSlicing'
        >>> nd3.data
        array([ 8,  7, 10])
    """

    def _arithmetic(self, operation, operand,
                    uncertainty_correlation=0, meta_kwds_operate=None,
                    meta_kwds_set=None, propagate_uncertainties=True,
                    handle_mask=True, handle_meta=True, compare_wcs=True,
                    **kwds):
        """
        Base method which calculates the result of the arithmetic operation.

        This class determines the result of the arithmetic operation on the
        ``data`` with respect to their units and then forwards to other methods
        to calculate the other properties for the result (like uncertainty).

        TODO: Simplify some parameter names, descriptions or presence? (embray)
        TODO: Maybe define the parameter descriptions as strings before so that
        this parameters description does not need to reference other methods...
        Or simply call them kwargs and say a full description of the parameters
        is given in method add?

        Parameters
        ----------
        operation: str
            The operation that is performed on the `NDData`. Supported are
            ``addition``, ``subtraction``, ``multiplication`` and ``division``.

            TODO: Use UFUNC instead of string? (embray, mwcraig)
            I didn't want to do it since each attribute can interpret the
            function that is used by the name of the operation better than by
            the function itself. Suppose we have ``data`` that is a
            `~numpy.matrix` or ``Pandas.DataFrame`` we need to use different
            UFUNCS (maybe) and that explodes the number of if/elif checks in
            uncertainty (and maybe other attributes like ``flags``).

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

        handle_meta: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        check_wcs: `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        kwargs:
            Any other parameter that should be passed to the
            different :meth:`NDArithmeticMixin._arithmetic_data` (or wcs, ...)
            methods.
            # TODO: Add Example for a use-case with subclassing.

        Returns
        -------
        result: `~numpy.ndarray` or `~astropy.units.Quantity`
            The resulting data as array (in case both operands were without
            unit) or as quantity if at least one had a unit.

        kwargs: `dict`
            the kwargs should contain all the other attributes (besides data
            and unit) to create a new instance for the result. Creating the
            new instance is up to the calling method, for example
            :meth:`NDArithmeticMixin.add`.

            TODO: Reword this? (embray)
            I did formulate a new text, does this look better to you?

        """
        # Convert the operand to the same class this allows for arithmetic
        # operations with numbers, lists, numpy arrays, numpy masked arrays
        # astropy quantities, masked quantities and of other subclasses of
        # NDData
        operand = self.__class__(operand)

        # TODO: Check that both data have numeric dtypes otherwise we could
        # end up with some arithmetic operation on string types. (init only
        # enforces that the dtype is NOT objects)?

        kwargs = {}
        # First check that the WCS allows the arithmetic operation
        if compare_wcs is None:
            kwargs['wcs'] = None
        elif not compare_wcs:
            if self.wcs is None:
                kwargs['wcs'] = deepcopy(operand.wcs)
            else:
                kwargs['wcs'] = deepcopy(self.wcs)
        else:
            kwargs['wcs'] = self._arithmetic_wcs(operand, compare_wcs,
                                                 **kwds)
        # Then calculate the resulting data (which can but not needs to be a
        # quantity)
        result = self._arithmetic_data(operation, operand, **kwds)
        # Determine the other properties
        if propagate_uncertainties is None:
            kwargs['uncertainty'] = None
        elif not propagate_uncertainties:
            if self.uncertainty is None:
                kwargs['uncertainty'] = deepcopy(operand.uncertainty)
            else:
                kwargs['uncertainty'] = deepcopy(self.uncertainty)
        else:
            kwargs['uncertainty'] = self._arithmetic_uncertainty(
                operation, operand, result, uncertainty_correlation, **kwds)
        if handle_mask is None:
            kwargs['mask'] = None
        elif not handle_mask:
            if self.mask is None:
                kwargs['mask'] = deepcopy(operand.mask)
            else:
                kwargs['mask'] = deepcopy(self.mask)
        else:
            kwargs['mask'] = self._arithmetic_mask(operand, **kwds)

        if handle_meta is None:
            kwargs['meta'] = None
        elif not handle_meta:
            if len(self.meta) == 0:
                kwargs['meta'] = deepcopy(operand.meta)
            else:
                kwargs['meta'] = deepcopy(self.meta)
        else:
            kwargs['meta'] = self._arithmetic_meta(
                operation, operand, meta_kwds_operate, meta_kwds_set, **kwds)

        # Wrap the individual results into a new instance of the same class.
        return result, kwargs

    def _arithmetic_data(self, operation, operand, **kwds):
        """
        Calculate the resulting data

        Parameters
        ----------
        operation: str
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

    def _arithmetic_uncertainty(self, operation, operand, result, correlation,
                                **kwds):
        """
        Calculate the resulting uncertainty.

        Parameters
        ----------
        operation: str
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
        # TODO: There is no enforced requirement that actually forbids the
        # uncertainty to have negative entries but with correlation the
        # sign of the uncertainty DOES matter.
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

    def _arithmetic_mask(self, operand, **kwds):
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
                                    "booleans or arrays of booleans.")
            if operand.mask is not None:
                if isinstance(operand.mask, bool):
                    pass
                elif (isinstance(operand.mask, np.ndarray) and
                        operand.mask.dtype == 'bool'):
                    pass
                else:
                    raise TypeError("Mask arithmetics is only defined for "
                                    "booleans or arrays of booleans.")

            # Now lets calculate the resulting mask (operation enforces copy)
            return self.mask | operand.mask

    def _arithmetic_wcs(self, operand, compare_wcs, **kwds):
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
                         meta_kwds_set, **kwds):
        """
        Calculate the resulting meta.

        There are two ways of altering the resulting meta, one is by
        applying the operation that was used on the data also on a keyword or
        by manually setting the keyword. The operation on keywords is done
        *before* setting the manually defined keywords.

        Parameters
        ----------
        operation: str
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
            # the other way) with the number of the other operand. For example
            # exposure time is multiplied with 5 if the data is multiplied by
            # 5.
            if len(result_meta) == 0:
                # No meta present. Nothing to do but issueing a warning that
                # it was not possible.
                log.info("No meta is present so no operation on them is "
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
                    log.info("The second operand is not a scalar. Cannot "
                             "operate meta keywords with array.s")
                else:
                    for i in meta_kwds_operate:
                        result_meta[i] = float(
                                        operator(self.meta[i], operand.data))
            elif len(operand.meta) > 0:
                # Same thing but reversed
                if self.data.size != 1:
                    log.info("The first operand is not a scalar. Cannot"
                             "operate meta keywords with arrays.")
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

    @tmp_deco(_arit_doc, name='addition', op='+')
    def add(self, operand, **kwargs):
        result, kwargs = self._arithmetic("addition", operand, **kwargs)
        return self.__class__(result, **kwargs)

    @tmp_deco(_arit_doc, name="subtraction", op="-")
    def subtract(self, operand, **kwargs):
        result, kwargs = self._arithmetic("subtraction", operand, **kwargs)
        return self.__class__(result, **kwargs)

    @tmp_deco(_arit_doc, name="multiplication", op="*")
    def multiply(self, operand, **kwargs):
        result, kwargs = self._arithmetic("multiplication", operand, **kwargs)
        return self.__class__(result, **kwargs)

    @tmp_deco(_arit_doc, name="division", op="/")
    def divide(self, operand, **kwargs):
        result, kwargs = self._arithmetic("division", operand, **kwargs)
        return self.__class__(result, **kwargs)

    # TODO: Sharedmethod instead of classmethod? (embray)
    # http://docs.astropy.org/en/v1.0.5/api/astropy.utils.decorators.sharedmethod.html#astropy.utils.decorators.sharedmethod
    # But that only allows for one documentation and we need the first operand.
    # Also it doesn't allow sphinx to build two documentations one for the
    # classmethod and one of the instance method and since creating one
    # doc that applies to both is not-trivial I set this option back for now.
    # So I decided to temporarly alter the name that there are not too obvious
    # name clashes by adding a 'ic_' prefix that stands for 'interclass'
    # arithmetic operations

    @classmethod
    @tmp_deco(_arit_cls_doc, 'add')
    def ic_addition(cls, operand1, operand2, **kwargs):
        # Convert the first operand to the implicit passed class (cls)
        # this allows for reverse operations.
        operand1 = cls(operand1)
        # Call the instance function for addition to proceed
        return operand1.add(operand2, **kwargs)

    @classmethod
    @tmp_deco(_arit_cls_doc, 'subtract')
    def ic_subtraction(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.subtract(operand2, **kwargs)

    @classmethod
    @tmp_deco(_arit_cls_doc, 'multiply')
    def ic_multiplication(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.multiply(operand2, **kwargs)

    @classmethod
    @tmp_deco(_arit_cls_doc, 'divide')
    def ic_division(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.divide(operand2, **kwargs)
