# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Arithmetic mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from copy import deepcopy

import numpy as np
import warnings

from ..nduncertainty import NDUncertainty
from ...units import dimensionless_unscaled
from ...utils import format_doc, sharedmethod
from ...utils.exceptions import AstropyDeprecationWarning

__all__ = ['NDArithmeticMixin']

# Global so it doesn't pollute the class dict unnecessarily:

# Docstring templates for add, subtract, multiply, divide methods.
_arit_doc = """
    Performs {name} by evaluating ``self`` {op} ``operand``.

    Parameters
    ----------
    operand, operand2 : `NDData`-like instance or convertible to one.
        If ``operand2`` is ``None`` or not given it will perform the operation
        ``self`` {op} ``operand``.
        If ``operand2`` is given it will perform ``operand`` {op} ``operand2``.
        If the method was called on a class rather than on the instance
        ``operand2`` must be given.

    propagate_uncertainties : `bool` or ``None``, optional
        If ``None`` the result will have no uncertainty. If ``False`` the
        result will have a copied version of the first operand that has an
        uncertainty. If ``True`` the result will have a correctly propagated
        uncertainty from the uncertainties of the operands but this assumes
        that the uncertainties are `NDUncertainty`-like. Default is ``True``.

        .. versionchanged:: 1.2
            This parameter must be given as keyword-parameter. Using it as
            positional parameter is deprecated.
            ``None`` was added as valid parameter value.

    handle_mask : callable, ``'first_found'`` or ``None``, optional
        If ``None`` the result will have no mask. If ``'first_found'`` the
        result will have a copied version of the first operand that has a
        mask). If it is a callable then the specified callable must
        create the results ``mask`` and if necessary provide a copy.
        Default is `numpy.logical_or`.

        .. versionadded:: 1.2

    handle_meta : callable, ``'first_found'`` or ``None``, optional
        If ``None`` the result will have no meta. If ``'first_found'`` the
        result will have a copied version of the first operand that has a
        (not empty) meta. If it is a callable then the specified callable must
        create the results ``meta`` and if necessary provide a copy.
        Default is ``None``.

        .. versionadded:: 1.2

    compare_wcs : callable, ``'first_found'`` or ``None``, optional
        If ``None`` the result will have no wcs and no comparison between
        the wcs of the operands is made. If ``'first_found'`` the
        result will have a copied version of the first operand that has a
        wcs. If it is a callable then the specified callable must
        compare the ``wcs``. The resulting ``wcs`` will be like if ``False``
        was given otherwise it raises a ``ValueError`` if the comparison was
        not successful. Default is ``'first_found'``.

        .. versionadded:: 1.2

    uncertainty_correlation : number or `~numpy.ndarray`, optional
        The correlation between the two operands is used for correct error
        propagation for correlated data as given in:
        https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
        Default is 0.

        .. versionadded:: 1.2


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

    ``"first_found"`` can also be abbreviated with ``"ff"``.
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
    Using this Mixin with `~astropy.nddata.NDData`:

        >>> from astropy.nddata import NDData, NDArithmeticMixin
        >>> class NDDataWithMath(NDArithmeticMixin, NDData):
        ...     pass

    Using it with one operand on an instance::

        >>> ndd = NDDataWithMath(100)
        >>> ndd.add(20)
        NDDataWithMath(120)

    Using it with two operand on an instance::

        >>> ndd = NDDataWithMath(5)
        >>> ndd.divide(1, ndd)
        NDDataWithMath(0.2)

    Using it as classmethod requires two operands::

        >>> NDDataWithMath.subtract(5, 4)
        NDDataWithMath(1)

    """

    def _arithmetic(self, operation, operand,
                    propagate_uncertainties=True, handle_mask=np.logical_or,
                    handle_meta=None, uncertainty_correlation=0,
                    compare_wcs='first_found', **kwds):
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

        operand : same type (class) as self
            see :meth:`NDArithmeticMixin.add`

        propagate_uncertainties : `bool` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_mask : callable, ``'first_found'`` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        handle_meta : callable, ``'first_found'`` or ``None``, optional
            see :meth:`NDArithmeticMixin.add`

        compare_wcs : callable, ``'first_found'`` or ``None``, optional
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
        elif compare_wcs in ['ff', 'first_found']:
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
        elif handle_mask in ['ff', 'first_found']:
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
        elif handle_meta in ['ff', 'first_found']:
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

    @sharedmethod
    @format_doc(_arit_doc, name='addition', op='+')
    def add(self, operand, operand2=None, **kwargs):
        return self._prepare_then_do_arithmetic(np.add, operand, operand2,
                                                **kwargs)

    @sharedmethod
    @format_doc(_arit_doc, name='subtraction', op='-')
    def subtract(self, operand, operand2=None, **kwargs):
        return self._prepare_then_do_arithmetic(np.subtract, operand, operand2,
                                                **kwargs)

    @sharedmethod
    @format_doc(_arit_doc, name="multiplication", op="*")
    def multiply(self, operand, operand2=None, **kwargs):
        return self._prepare_then_do_arithmetic(np.multiply, operand, operand2,
                                                **kwargs)

    @sharedmethod
    @format_doc(_arit_doc, name="division", op="/")
    def divide(self, operand, operand2=None, **kwargs):
        return self._prepare_then_do_arithmetic(np.true_divide, operand,
                                                operand2, **kwargs)

    @sharedmethod
    def _prepare_then_do_arithmetic(self_or_cls, operation, operand, operand2,
                                    **kwargs):
        """Intermediate method called by public arithmetics (i.e. ``add``)
        before the processing method (``_arithmetic``) is invoked.

        .. warning::
            Do not override this method in subclasses.

        This method checks if it was called as instance or as class method and
        then wraps the operands and the result from ``_arithmetics`` in the
        appropriate subclass.

        Parameters
        ----------
        self_or_cls : instance or class
            ``sharedmethod`` behaves like a normal method if called on the
            instance (then this parameter is ``self``) but like a classmethod
            when called on the class (then this parameter is ``cls``).

        operations : callable
            The operation (normally a numpy-ufunc) that represents the
            appropriate action.

        operand, operand2, kwargs :
            See for example ``add``.

        Result
        ------
        result : `~astropy.nddata.NDData`-like
            Depending how this method was called either ``self_or_cls``
            (called on class) or ``self_or_cls.__class__`` (called on instance)
            is the NDData-subclass that is used as wrapper for the result.
        """
        # DO NOT OVERRIDE THIS METHOD IN SUBCLASSES.

        # TODO: Remove this in astropy 1.3 or 1.4:

        # Before 1.2 propagate_uncertainties could be given as positional
        # keyword, this is now deprecated:
        if (isinstance(operand2, bool) and
                'propagate_uncertainties' not in kwargs):
            # No explicit propagate_uncertainties was given but the second
            # operand was given as boolean. I'll assume that most don't want
            # to do arithmetics with a boolean operand, print a deprecation
            # warning. If someone really wanted to do arithmetics with a
            # boolean he should have set propagate_uncertainties. :-/
            warnings.warn('propagate_uncertainties should be given as keyword '
                          'parameter, i.e. "propagate_uncertainties={0}".'
                          ''.format(operand2), AstropyDeprecationWarning)
            # Set the kwarg and reset operand2.
            kwargs['propagate_uncertainties'] = operand2
            operand2 = None

        # TODO: The following parts must remain here if the above part is
        # removed.

        if isinstance(self_or_cls, NDArithmeticMixin):
            # True means it was called on the instance, so self_or_cls is
            # a reference to self
            cls = self_or_cls.__class__

            if operand2 is None:
                # Only one operand was given. Set operand2 to operand and
                # operand to self so that we call the appropriate method of the
                # operand.
                operand2 = operand
                operand = self_or_cls
            else:
                # Convert the first operand to the class of this method.
                # This is important so that always the correct _arithmetics is
                # called later that method.
                operand = cls(operand)

        else:
            # It was used as classmethod so self_or_cls represents the cls
            cls = self_or_cls

            # It was called on the class so we expect two operands!
            if operand2 is None:
                raise TypeError("operand2 must be given when the method isn't "
                                "called on an instance.")

            # Convert to this class. See above comment why.
            operand = cls(operand)

        # At this point operand, operand2, kwargs and cls are determined.

        # Let's try to convert operand2 to the class of operand to allows for
        # arithmetic operations with numbers, lists, numpy arrays, numpy masked
        # arrays, astropy quantities, masked quantities and of other subclasses
        # of NDData.
        operand2 = cls(operand2)

        # Now call the _arithmetics method to do the arithmetics.
        result, init_kwds = operand._arithmetic(operation, operand2, **kwargs)

        # Return a new class based on the result
        return cls(result, **init_kwds)
