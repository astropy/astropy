# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np

from ..nduncertainty import NDUncertainty
from ...units import dimensionless_unscaled
from ...extern import six

__all__ = ['NDArithmeticMixin']

try:
    from ...utils import format_doc
except ImportError:
    # TODO: Delete this and rebase if #4242 is merged.
    def format_doc(docstring, *args, **kwargs):
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

# Global so it doesn't pollute the class dict unnecessarily:

# Docstring templates for add, subtract, multiply, divide methods.
_arit_doc = """
    Performs {name} by evaluating ``self`` {op} ``operand``.

    For the inverse operation or operations between different subclasses take a
    look at
    :meth:`~astropy.nddata.NDArithmeticMixin.ic_{name}` .

    Parameters
    ----------
    operand : `NDData`-like instance or convertible to one.
        The second operand in the operation.

    propagate_uncertainties : `bool` or ``None``, optional
        If ``None`` the result will have no uncertainty. If ``False`` the
        result will have a copied version of the first operand that has an
        uncertainty. If ``True`` the result will have a correctly propagated
        uncertainty from the uncertainties of the operands but this assumes
        that the uncertainties are `StdDevUncertainty`. Default is ``True``.

    handle_mask : callable, ``False`` or ``None``, optional
        If ``None`` the result will have no mask. If ``False`` the
        result will have a copied version of the first operand that has a
        mask). If it is a callable then the specified callable must
        create the results ``mask`` and if necessary provide a copy.
        Default is `numpy.logical_or`.

    handle_meta : callable, ``False`` or ``None``, optional
        If ``None`` the result will have no meta. If ``False`` the
        result will have a copied version of the first operand that has a
        (not empty) meta. If it is a callable then the specified callable must
        create the results ``meta`` and if necessary provide a copy.
        Default is ``False``.

    compare_wcs : callable, ``False`` or ``None``, optional
        If ``None`` the result will have no wcs and no comparison between
        the wcs of the operands is made. If ``False`` the
        result will have a copied version of the first operand that has a
        wcs. If it is a callable then the specified callable must
        compare the ``wcs``. The resulting ``wcs`` will be like if ``False``
        was given otherwise it raises a ``ValueError`` if the comparison was
        not successful. Default is ``False``.

    uncertainty_correlation : number or `~numpy.ndarray`, optional
        The correlation (rho) is defined between the uncertainties in
        ``sigma_AB = sigma_A * sigma_B * rho`` (Latex?). If
        ``propagate_uncertainties`` is *not* ``True`` this will be ignored.
        A value of ``0`` means uncorrelated (which is also the default).
        TODO: TeX of formula?

    kwargs :
        Any other parameter that should be passed to the callables used.

    Returns
    -------
    result : `~astropy.nddata.NDData`-like
        The resulting dataset

    Notes
    -----
    If a ``callable`` is used for ``mask``, ``wcs`` or ``meta`` the
    callable must accept the corresponding attributes as first two
    parameters. If the callable also needs additional parameters these can be
    defined as ``kwargs`` and must start with ``"wcs_"`` (for wcs callable) or
    ``"meta_"`` (for meta callable). This startstring is removed before the
    callable is called.
    """

_arit_cls_doc = """
    Like :meth:`NDArithmeticMixin.{0}` you can {0} two operands.

    This method is avaiable as classmethod in order to allow arithmetic
    operations between arbitary objects as long as they are convertable to
    the class that called the method. Therefore the name prefix ``ic_`` which
    stands for ``interclass`` operations.

    Parameters
    ----------
    operand1 : `NDData`-like or convertable to `NDData`
        the first operand in the operation.

    operand2 : `NDData`-like or convertable to `NDData`
        the second operand in the operation.

    kwargs :
        see :meth:`NDArithmeticMixin.{0}`

    Returns
    -------
    result : `NDData`-like
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
    (if this does not exist from the second operand) as the result's
    uncertainty. This behaviour is also explained in the docstring for the
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
        >>> nd_inv.__class__.__name__
        'NDDataWithMath'
        >>> nd_inv.data
        array([ 1.        ,  0.5       ,  0.33333333])
        >>> nd_inv.unit
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
                    propagate_uncertainties=True, handle_mask=np.logical_or,
                    handle_meta=False, uncertainty_correlation=0,
                    compare_wcs=False, **kwds):
        """
        Base method which calculates the result of the arithmetic operation.

        This method determines the result of the arithmetic operation on the
        ``data`` including their units and then forwards to other methods
        to calculate the other properties for the result (like uncertainty).

        Parameters
        ----------
        operation : callable
            The operation that is performed on the `NDData`. Supported are
            `numpy.add`, `numpy.subtract`, `numpy.multiply` and
            `numpy.true_divide`.

        operand : `NDData` instance or something that can be converted to one.
            see :meth:`NDArithmeticMixin.add`

        propagate_uncertainties : `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_mask : callable, ``False`` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_meta : callable, ``False`` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        compare_wcs : callable, ``False`` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        uncertainty_correlation : ``Number`` or `~numpy.ndarray`, optional
            see :meth:`NDArithmeticMixin.add`

        kwargs :
            Any other parameter that should be passed to the
            different :meth:`NDArithmeticMixin._arithmetic_mask` (or wcs, ...)
            methods.

        Returns
        -------
        result : `~numpy.ndarray` or `~astropy.units.Quantity`
            The resulting data as array (in case both operands were without
            unit) or as quantity if at least one had a unit.

        kwargs : `dict`
            The kwargs should contain all the other attributes (besides data
            and unit) needed to create a new instance for the result. Creating
            the new instance is up to the calling method, for example
            :meth:`NDArithmeticMixin.add`.

        """
        # Convert the operand to the same class this allows for arithmetic
        # operations with numbers, lists, numpy arrays, numpy masked arrays
        # astropy quantities, masked quantities and of other subclasses of
        # NDData
        operand = self.__class__(operand)

        # Find the appropriate keywords for the appropriate method (not sure
        # if data and uncertainty are ever used ...)
        kwds2 = {'mask': {}, 'meta': {}, 'wcs': {},
                 'data': {}, 'uncertainty': {}}
        for i in kwds:
            splitted = i.split('_', 1)
            try:
                kwds2[splitted[0]][splitted[1]] = kwds[i]
            except KeyError:
                raise KeyError('Unknown prefix {0} for parameter {1}'
                               ''.format(splitted[0], i))

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
            kwargs['wcs'] = self._arithmetic_wcs(operation, operand,
                                                 compare_wcs, **kwds2['wcs'])

        # Then calculate the resulting data (which can but not needs to be a
        # quantity)
        result = self._arithmetic_data(operation, operand, **kwds2['data'])

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
                operation, operand, result, uncertainty_correlation,
                **kwds2['uncertainty'])

        if handle_mask is None:
            kwargs['mask'] = None
        elif not handle_mask:
            if self.mask is None:
                kwargs['mask'] = deepcopy(operand.mask)
            else:
                kwargs['mask'] = deepcopy(self.mask)
        else:
            kwargs['mask'] = self._arithmetic_mask(operation, operand,
                                                   handle_mask,
                                                   **kwds2['mask'])

        if handle_meta is None:
            kwargs['meta'] = None
        elif not handle_meta:
            if not self.meta:
                kwargs['meta'] = deepcopy(operand.meta)
            else:
                kwargs['meta'] = deepcopy(self.meta)
        else:
            kwargs['meta'] = self._arithmetic_meta(
                operation, operand, handle_meta, **kwds2['meta'])

        # Wrap the individual results into a new instance of the same class.
        return result, kwargs

    def _arithmetic_data(self, operation, operand, **kwds):
        """
        Calculate the resulting data

        Parameters
        ----------
        operation : callable
            see `NDArithmeticMixin._arithmetic` parameter description.

        operand : `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        kwds :
            Additional parameters.

        Returns
        -------
        result_data : `~numpy.ndarray` or `~astropy.units.Quantity`
            If both operands had no unit the resulting data is a simple numpy
            array, but if any of the operands had a unit the return is a
            Quantity.
        """

        # Do the calculation with or without units
        if self.unit is None and operand.unit is None:
            result = operation(self.data, operand.data)
        elif self.unit is None:
            result = operation(self.data * dimensionless_unscaled,
                               operand.data * operand.unit)
        elif operand.unit is None:
            result = operation(self.data * self.unit,
                               operand.data * dimensionless_unscaled)
        else:
            result = operation(self.data * self.unit,
                               operand.data * operand.unit)

        return result

    def _arithmetic_uncertainty(self, operation, operand, result, correlation,
                                **kwds):
        """
        Calculate the resulting uncertainty.

        Parameters
        ----------
        operation : callable
            see :meth:`NDArithmeticMixin._arithmetic` parameter description.

        operand : `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        result : `~astropy.units.Quantity` or `~numpy.ndarray`
            The result of :meth:`NDArithmeticMixin._arithmetic_data`.

        correlation : number or `~numpy.ndarray`
            see :meth:`NDArithmeticMixin.add` parameter description.

        kwds :
            Additional parameters.

        Returns
        -------
        result_uncertainty : `NDUncertainty` subclass instance or None
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

    def _arithmetic_mask(self, operation, operand, handle_mask, **kwds):
        """
        Calculate the resulting mask

        This is implemented as the piecewise ``or`` operation if both have a
        mask.

        Parameters
        ----------
        operation : callable
            see :meth:`NDArithmeticMixin._arithmetic` parameter description.
            By default, the ``operation`` will be ignored.

        operand : `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        handle_mask : callable
            see :meth:`NDArithmeticMixin.add`

        kwds :
            Additional parameters given to ``handle_mask``.

        Returns
        -------
        result_mask : any type
            If only one mask was present this mask is returned.
            If neither had a mask ``None`` is returned. Otherwise
            ``handle_mask`` must create (and copy) the returned mask.
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
            # Now lets calculate the resulting mask (operation enforces copy)
            return handle_mask(self.mask, operand.mask, **kwds)

    def _arithmetic_wcs(self, operation, operand, compare_wcs, **kwds):
        """
        Calculate the resulting wcs.

        There is actually no calculation involved but it is a good place to
        compare wcs information of both operands. This is currently not working
        properly with `~astropy.wcs.WCS` (which is the suggested class for
        storing as wcs property) but it will not break it neither.

        Parameters
        ----------
        operation : callable
            see :meth:`NDArithmeticMixin._arithmetic` parameter description.
            By default, the ``operation`` will be ignored.

        operand : `NDData` instance or subclass
            The second operand wrapped in an instance of the same class as
            self.

        compare_wcs : callable
            see :meth:`NDArithmeticMixin.add` parameter description.

        kwds :
            Additional parameters given to ``compare_wcs``.

        Raises
        ------
        ValueError
            If ``compare_wcs`` returns ``False``.

        Returns
        -------
        result_wcs : any type
            The ``wcs`` of the first operand is returned.
        """

        # ok, not really arithmetics but we need to check which wcs makes sense
        # for the result and this is an ideal place to compare the two WCS,
        # too.

        # I'll assume that the comparison returned None or False in case they
        # are not equal.
        if not compare_wcs(self.wcs, operand.wcs, **kwds):
            raise ValueError("WCS are not equal.")

        return self.wcs

    def _arithmetic_meta(self, operation, operand, handle_meta, **kwds):
        """
        Calculate the resulting meta.

        Parameters
        ----------
        operation : callable
            see :meth:`NDArithmeticMixin._arithmetic` parameter description.
            By default, the ``operation`` will be ignored.

        operand : `NDData`-like instance
            The second operand wrapped in an instance of the same class as
            self.

        handle_meta : callable
            see :meth:`NDArithmeticMixin.add`

        kwds :
            Additional parameters given to ``handle_meta``.

        Returns
        -------
        result_meta : any type
            The result of ``handle_meta``.
        """
        # Just return what handle_meta does with both of the metas.
        return handle_meta(self.meta, operand.meta, **kwds)

    @format_doc(_arit_doc, name='addition', op='+')
    def add(self, operand, **kwargs):
        result, kwargs = self._arithmetic(np.add, operand, **kwargs)
        return self.__class__(result, **kwargs)

    @format_doc(_arit_doc, name="subtraction", op="-")
    def subtract(self, operand, **kwargs):
        result, kwargs = self._arithmetic(np.subtract, operand, **kwargs)
        return self.__class__(result, **kwargs)

    @format_doc(_arit_doc, name="multiplication", op="*")
    def multiply(self, operand, **kwargs):
        result, kwargs = self._arithmetic(np.multiply, operand, **kwargs)
        return self.__class__(result, **kwargs)

    @format_doc(_arit_doc, name="division", op="/")
    def divide(self, operand, **kwargs):
        result, kwargs = self._arithmetic(np.true_divide, operand, **kwargs)
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
    @format_doc(_arit_cls_doc, 'add')
    def ic_addition(cls, operand1, operand2, **kwargs):
        # Convert the first operand to the implicit passed class (cls)
        # this allows for reverse operations.
        operand1 = cls(operand1)
        # Call the instance function for addition to proceed
        return operand1.add(operand2, **kwargs)

    @classmethod
    @format_doc(_arit_cls_doc, 'subtract')
    def ic_subtraction(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.subtract(operand2, **kwargs)

    @classmethod
    @format_doc(_arit_cls_doc, 'multiply')
    def ic_multiplication(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.multiply(operand2, **kwargs)

    @classmethod
    @format_doc(_arit_cls_doc, 'divide')
    def ic_division(cls, operand1, operand2, **kwargs):
        operand1 = cls(operand1)
        return operand1.divide(operand2, **kwargs)

    # TODO: Only temporary since uncertainty setter was changed during refactor
    # Remove after #4270 is merged.
    @property
    def uncertainty(self):
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if isinstance(value, NDUncertainty):
            value.parent_nddata = self
        self._uncertainty = value
