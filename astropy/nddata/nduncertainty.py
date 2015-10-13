# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from abc import ABCMeta, abstractproperty, abstractmethod
from copy import deepcopy

#from ..utils.compat import ignored
from ..units import Unit, Quantity
from ..extern import six
from .. import log

__all__ = ['MissingDataAssociationException',
           'IncompatibleUncertaintiesException', 'NDUncertainty',
           'StdDevUncertainty']


class IncompatibleUncertaintiesException(Exception):
    """
    This exception should be used to indicate cases in which uncertainties
    with two different classes can not be propagated.
    """


class MissingDataAssociationException(Exception):
    """
    This exception should be used to indicate that an uncertainty instance has
    not been associated with a parent `~astropy.nddata.NDData` object.
    """



# Make a placeholder for the different uncertainty propagation methods.
_propagate_doc = """
Propagate uncertainties for {operation}.

Parameters
----------
other_uncert : {instance} instance
    The data for the uncertainty of b in a {operator} b
result_data : `~numpy.ndarray` instance or `~astropy.units.Quantity`
    The data array that is the result of the {operation}.
correlation: `Number` or `~numpy.ndarray`
    Array or scalar representing the correlation. If the subclass does not
    support correlated uncertainties this will be replaced by 0 (uncorrelated).

Returns
-------
result_uncertainty : {instance} instance
    The resulting uncertainty

Raises
------
ValueError
    Raised if the uncertainty arrays or resulting data cannot be used together
    for the arithmetic operation (not broadcastable) or conflicting units.

Notes
-----
Handling units (especially if the differ from the parents unit) is done in here
but be aware that since there are no checks for the unit that incorrect
assigned units may break the uncertainty propagation even if the resulting data
(with units) can be computed.
"""


@six.add_metaclass(ABCMeta)
class NDUncertainty(object):
    """
    This is the base class for uncertainty classes used with NDData.

    Parameters
    ----------
    array: any type, optional
        The array or value (the parameter name is due to historical reasons) of
        the uncertainty. There is nothing enforced but `~numpy.ndarray`,
        `~astropy.units.Quantity` or `NDUncertainty` subclasses are
        recommended. If the `array` is `list`-like it will be cast to a
        `~numpy.ndarray`.

    unit: `~astropy.units.Unit` or `str`, optional
        The unit of the uncertainty ``array``. If input is a `str` this will be
        converted to an `~astropy.units.Unit`.

    copy: `bool`, optional
        Should the uncertainty be saved as reference (``False``) or as a copy
        (``True``). Defaults to True.

    Notes
    -----
    1. NDUncertainty is an abstract class and should *never* be instantiated
       directly.

    2. NDUncertainty takes a ``unit`` parameter and if the `array` is something
       with a unit (like a `~astropy.units.Quantity` of another
       `NDUncertainty`) the `unit` parameter is considered the True unit of
       the `NDUncertainty` instance and a warning is issued that the unit is
       overwritten. *No* conversion is done while overwriting the unit so the
       uncertainty ``array`` is not altered.

    3. NDUncertainty provides a :meth:``__getitem__`` method so slicing is in
       theory supported but that relies on the `array` being something that
       actually can be sliced.

    4. Subclasses need to define:

       - property ``uncertainty_type`` which should be a small string defining
         the kind of uncertainty. Recommened is ``std`` for standard deviation
         ``var`` for variance (following the ``numpy`` conventions).
       - :meth:`_propagate_add` and similar ones which takes the
         uncertainty of the other element and the resulting data (or quantity)
         and calculates the array (or quantity) which represents the resulting
         uncertainty.
       - Most of the time a subclasses will need to extend or override the
         :meth:`NDUncertainty._convert_uncertainty`. This method is responsible
         for checking that the other uncertainty has a class that can be used
         for uncertainty propagation. (`NDUncertainty` only checks that it is
         another instance or subclass of `NDUncertainty`). For subclasses that
         cannot propagate with arbitary other uncertainties they should check
         that the other uncertainty has the same class as they are *or* convert
         it to such a class. (Handling units should be part of the
         ``_propagate_*`` methods and *NOT* done in there).
       - The :meth:`NDUncertainty.propagate` method should
         only be overriden if one wants to allow other operations than
         ``addition``, ``subtraction``, ``multiplication`` or ``division``.
         This function is the common entry point of any arithmetic computation
         that tries to propagate uncertainties and is responsible for calling
         the appropriate ``_propagate_*`` method and then returns another
         instance of the same class with the results uncertainty.

    5. `NDUncertainty` and it's subclasses can only be used for uncertainty
       propagation if their ``parent_nddata`` attribute is a reference to the
       `NDData` object. This is needed since many kinds of propagation need the
       actual data besides the uncertainty.

    6. If a ``unit`` is implicit (``data`` had a unit) or explicitly passed
       to the ``__init__`` one should be sure that this unit is identical or
       convertable to the unit of the parent. Otherwise uncertainty propagation
       should fail.
    """

    supports_correlated = False

    def __init__(self, array=None, unit=None, copy=True):
        # It is possible to create an uncertainty from another Uncertainty or
        # a quantity or a list
        if isinstance(array, NDUncertainty):
            # Check that the uncertainty type is the same otherwise every
            # uncertainty could be converted to another uncertainty_type.
            if array.uncertainty_type != self.uncertainty_type:
                raise IncompatibleUncertaintiesException
            # Check if two units are given and take the explicit one then.
            if unit is not None and array.unit is not None:
                log.info("Overwriting Uncertainty's current "
                         "unit with specified unit")
            elif array.unit is not None:
                unit = array.unit
            array = array.array

        elif isinstance(array, Quantity):
            # Check if two units are given and take the explicit one then.
            if unit is not None and array.unit is not None:
                log.info("Overwriting Quantity's current "
                         "unit with specified unit")
            elif array.unit is not None:
                unit = array.unit
            array = array.value

        elif isinstance(array, list):
            copy = False # since lists get copied while casting them to ndarray
            array = np.array(array)

        if copy:
            array = deepcopy(array)
            unit = deepcopy(unit)

        self._array = array
        if unit is None:
            self._unit = None
        else:
            self._unit = Unit(unit)
        self.parent_nddata = None # no associated NDData - until it is set!

    def __getitem__(self, item):
        """
        Simple slicing is allowed but returns a reference *not* a copy. This
        assimilates numpy slicing.
        """
        return self.__class__(self.array[item], unit=self.unit, copy=False)

    @abstractproperty
    def uncertainty_type(self):
        """
        `str`: Short description which kind of uncertainty is saved.

        Defined as abstract property so subclasses *have* to override this
        and return a string.
        """
        return None

    @property
    def parent_nddata(self):
        """
        `NDData` reference: The `NDData` whose uncertainty this is.

        In case the reference is not set uncertainty propagation will not be
        possible since almost all kinds of propagation need the uncertain
        data besides the uncertainty.
        """
        message = "Uncertainty is not associated with an NDData object"
        try:
            if self._parent_nddata is None:
                raise MissingDataAssociationException(message)
            else:
                return self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException(message)

    @parent_nddata.setter
    def parent_nddata(self, value):
        self._parent_nddata = value

    @property
    def array(self):
        """
        any type: the uncertainty's value.
        """
        return self._array

    @property
    def unit(self):
        """
        `~astropy.units.Unit`: The unit of the uncertainty.

        Even though it is not enforced the unit should be convertable to the
        ``parent_nddata`` unit. Otherwise uncertainty propagation might give
        wrong results.

        If the unit is not set the unit of the parent will be returned.
        """
        if self._unit is None:
            if (self._parent_nddata is None or
                                    self.parent_nddata.unit is None):
                return None
            else:
                return self.parent_nddata.unit
        return self._unit

    def propagate(self, operation, other_nddata, result_data, correlation):
        """
        Calculate the resulting uncertainty of the arithmetic operation.

        Parameters
        ----------
        operation: `str`
            The operation that is performed on the `NDData`. Supported are
            ``addition``, ``subtraction``, ``multiplication`` and ``division``.

        other_nddata: `NDData` instance
            The second NDData in the arithmetic operation.

        result_data: `~numpy.ndarray` or `~astropy.units.Quantity`
            The result of the arithmetic operations. This saves some duplicate
            calculations.

        correlation: ``Number`` or `~numpy.ndarray`
            The correlation (rho) is defined between the uncertainties in
            sigma_AB = sigma_A * sigma_B * rho. A value of ``0`` means
            uncorrelated operands.

        Returns
        -------
        resulting_uncertainty: `NDUncertainty` instance
            Another instance of the same `NDUncertainty` subclass containing
            the uncertainty of the result.

        Raises
        ------
        ValueError:
            If the ``operation`` is not supported.

        Notes
        -----
        This method at first tries to convert the ``uncertainty`` of the
        other operand to a `NDUncertainty` subclass that is useable for
        uncertainty propagation (this check and optional conversion is done
        in :meth:`NDUncertainty._convert_uncertainty`) and then the
        appropriate ``_propagate_*`` method for calculating the resulting
        uncertainty is invoked. Afterwards this function wraps it into
        another instance of `NDUncertainty` and returns it.
        """
        # Check if the subclass supports correlation
        if not self.supports_correlated:
            if correlation != 0:
                log.info("This subclass does not support correlation")
                correlation = 0

        # Get the other uncertainty (and convert it to a matching one)
        other_uncert = self._convert_uncertainty(other_nddata.uncertainty)

        if operation == 'addition':
            result = self._propagate_add(other_uncert, result_data,
                                         correlation)
        elif operation == 'subtraction':
            result = self._propagate_subtract(other_uncert, result_data,
                                              correlation)
        elif operation == 'multiplication':
            result = self._propagate_multiply(other_uncert, result_data,
                                              correlation)
        elif operation == 'division':
            result = self._propagate_divide(other_uncert, result_data,
                                            correlation)
        else:
            raise ValueError('Unsupported operation')

        return self.__class__(result, copy=False)

    def _convert_uncertainty(other_uncert):
        """
        Checks that the uncertainties are compatible for propagation.

        Parameters
        ----------
        other_uncert: `NDUncertainty` subclass
            The other uncertainty

        Returns
        -------
        other_uncert: `NDUncertainty` subclass
            but converted to a compatible `NDUncertainty` subclass if
            necessary.

        Raises
        ------
        IncompatibleUncertaintiesException:
            If the other uncertainty cannot be converted to a compatible
            `NDUncertainty` subclass.

        Notes
        -----
        For most subclasses a simple check that the classes are the same would
        be sufficient but in case one wishes to create a conversion scheme this
        would be the right place to implement it.
        """
        if isinstance(other_uncert, NDUncertainty):
            return other_uncert
        else:
            raise IncompatibleUncertaintiesException

    @abstractmethod
    def _propagate_add(self, other_uncert, result_data, correlation):
        return None

    @abstractmethod
    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return None

    @abstractmethod
    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return None

    @abstractmethod
    def _propagate_divide(self, other_uncert, result_data, correlation):
        return None

    # TODO: Make a decorator that adds the docstrings instead.

    # Apply docstrings
    _propagate_add.__doc__ = _propagate_doc.format(operation='addition',
            operator='+', instance='NDUncertainty')
    _propagate_subtract.__doc__ = _propagate_doc.format(
            operation='subtraction', operator='-', instance='NDUncertainty')
    _propagate_multiply.__doc__ = _propagate_doc.format(
            operation='multiplication', operator='*', instance='NDUncertainty')
    _propagate_divide.__doc__ = _propagate_doc.format(operation='divison',
            operator='/', instance='NDUncertainty')


class StdDevUncertainty(NDUncertainty):
    """
    A class for standard deviation uncertainty.

    This class implements uncertainty propagation for ``addition``,
    ``subtraction``, ``multiplication`` and ``division`` but only with other
    instances of `StdDevUncertainty`. The class can fully handle if the
    uncertainty has a unit that differs from (but is convertable to) the
    parents `NDData` unit but converts the unit of the propagated uncertainty
    to the unit of the resulting data. Also support for correlation is possible
    but that requires that the correlation is an input it cannot handle
    correlation determination itself.
    """

    supports_correlated = True

    @property
    def uncertainty_type(self):
        """
        `str`: ``std``
            `StdDevUncertainty` implements standard deviation.
        """
        return 'std'

    def _convert_uncertainty(self, other_uncert):
        if isinstance(other_uncert, StdDevUncertainty):
            return other_uncert
        else:
            raise IncompatibleUncertaintiesException


    def _propagate_add(self, other_uncert, result_data, correlation):

        if self.array is None:
            # Formula sigma = dB

            if other_uncert.unit is not None and (
                        result_data.unit != other_uncert.unit):
                # If the other uncertainty has a unit and this unit differs
                # from the unit of the result convert it to the results unit
                return (other_uncert.array * other_uncert.unit).to(
                            result_data.unit).value
            else:
                # Copy the result because _propagate will not copy it but for
                # arithmetic operations users will expect copys.
                return deepcopy(other_uncert.array)

        elif other_uncert.array is None:
            # Formula sigma = dA

            if self.unit is not None and self.unit != self.parent_nddata.unit:
                # If the uncertainty has a different unit than the result we
                # need to convert it to the results unit.
                return (self.array * self.unit).to(result_data.unit).value
            else:
                # Copy the result because _propagate will not copy it but for
                # arithmetic operations users will expect copys.
                return deepcopy(self.array)

        else:
            # Formula sigma = sqrt(dA**2 + dB**2 + 2*rho*dA*dB)

            if self.unit != other_uncert.unit:
                # In case the two uncertainties (or data) have different units
                # we need to use quantity operations. The case where only on
                # has a unit and the other doesn't is not possible with
                # addition and would have raised an exception in the data
                # computation
                this = self.array * self.unit
                other = other_uncert.array * other_uncert.unit
            else:
                # Since both units are the same or None we can just use
                # numpy operations
                this = self.array
                other = other_uncert.array

            # Determine the result depending on the correlation
            if correlation != 0:
                corr = 2 * correlation * this * other
                result = np.sqrt(this**2 + other**2 + corr)
            else:
                result = np.sqrt(this**2 + other**2)

            if isinstance(result, Quantity):
                # In case we worked with quantities we need to return the
                # uncertainty that has the same unit as the resulting data
                if result.unit == result_data.unit:
                    return result.value
                else:
                    # Convert it to the data's unit and then drop the unit.
                    return result.to(result_data.unit).value
            else:
                return result


    def _propagate_subtract(self, other_uncert, result_data, correlation):
        # Since the formulas are equivalent to addition you should look at the
        # explanations provided in _propagate_add

        if self.array is None:
            if other_uncert.unit is not None and (
                        result_data.unit != other_uncert.unit):
                return (other_uncert.array * other_uncert.unit).to(
                            result_data.unit).value
            else:
                return deepcopy(other_uncert.array)
        elif other_uncert.array is None:
            if self.unit is not None and self.unit != self.parent_nddata.unit:
                return (self.array * self.unit).to(result_data.unit).value
            else:
                return deepcopy(self.array)
        else:
            # Formula sigma = sqrt(dA**2 + dB**2 - 2*rho*dA*dB)
            if self.unit != other_uncert.unit:
                this = self.array * self.unit
                other = other_uncert.array * other_uncert.unit
            else:
                this = self.array
                other = other_uncert.array
            if correlation != 0:
                corr = 2 * correlation * this * other
                # The only difference to addition is that the correlation is
                # subtracted.
                result = np.sqrt(this**2 + other**2 - corr)
            else:
                result = np.sqrt(this**2 + other**2)
            if isinstance(result, Quantity):
                if result.unit == result_data.unit:
                    return result.value
                else:
                    return result.to(result_data.unit).value
            else:
                return result


    def _propagate_multiply(self, other_uncert, result_data, correlation):

        # For multiplication we don't need the result as quantity
        if isinstance(result_data, Quantity):
            result_data = result_data.value

        if self.array is None:
            # Formula sigma = A * dB

            # We want the result to have the same unit as the result so we
            # only need to convert the unit of the other uncertainty if it is
            # different from it's datas unit.
            if other_uncert.unit != other_uncert.parent_nddata.unit:
                other = (other_uncert.array * other_uncert.unit).to(
                            other_uncert.parent_nddata.unit).value
            else:
                other = other_uncert.array
            return np.abs(self.parent_nddata.data * other)

        elif other_uncert.array is None:
            # Formula sigma = dA * B

            # Just the reversed case
            if self.unit != self.parent_nddata.unit:
                this = (self.array * self.unit).to(
                                            self.parent_nddata.unit).value
            else:
                this = self.array
            return np.abs(other_uncert.parent_nddata.data * this)

        else:
            # Formula sigma = AB * sqrt((dA/A)**2 + (dB/B)**2 + 2*dA/A*dB/B)

            # In this case we just need to catch the case if one uncertainty
            # has a unit that is not the same as the data's because we try to
            # work with dimensionless fractions. Because it allows for numpy
            # operations instead of the slower quantity operations
            if self.unit != self.parent_nddata.unit:
                left = (self.array * self.unit).to(
                    self.parent_nddata.unit).value / self.parent_nddata.data
            else:
                left = self.array / self.parent_nddata.data

            if other_uncert.unit != other_uncert.parent_nddata.unit:
                right = ((other_uncert.array * other_uncert.unit).to(
                    other_uncert.parent_nddata.unit).value /
                    other_uncert.parent_nddata.data)
            else:
                right = (other_uncert.array / other_uncert.parent_nddata.data)

            if correlation != 0:
                corr = 2 * correlation * left * right
                return np.abs(result_data) * np.sqrt(left**2 + right**2 + corr)
            else:
                return np.abs(result_data) * np.sqrt(left**2 + right**2)


    def _propagate_divide(self, other_uncert, result_data, correlation):

        # For division we don't need the result as quantity
        if isinstance(result_data, Quantity):
            result_data = result_data.value

        if self.array is None:
            # Formula sigma = (A / B) * (dB / B)

            # We need (db / B) to be dimensionless so we convert (if necessary)
            # dB to the same unit as B
            if other_uncert.unit != other_uncert.parent_nddata.unit:
                right = ((other_uncert.array * other_uncert.unit).to(
                    other_uncert.parent_nddata.unit).value /
                    other_uncert.parent_nddata.data)
            else:
                right = (other_uncert.array / other_uncert.parent_nddata.data)
            return np.abs(result_data * right)

        elif other_uncert.array is None:
            # Formula sigma = dA / B.

            # We need to convert dA to the unit of A to have a result that
            # matches the resulting data's unit.
            if self.unit != self.parent_nddata.unit:
                this = (self.array * self.unit).to(self.parent_nddata.unit)
            else:
                this = self.array
            return np.abs(this / other_uncert.parent_nddata.data)

        else:
            # Formula sigma = AB * sqrt((dA/A)**2 + (dB/B)**2 - 2*dA/A*dB/B)
            if self.unit != self.parent_nddata.unit:
                left = (self.array * self.unit).to(
                    self.parent_nddata.unit).value / self.parent_nddata.data
            else:
                left = self.array / self.parent_nddata.data
            if other_uncert.unit != other_uncert.parent_nddata.unit:
                right = ((other_uncert.array * other_uncert.unit).to(
                    other_uncert.parent_nddata.unit).value /
                    other_uncert.parent_nddata.data)
            else:
                right = (other_uncert.array / other_uncert.parent_nddata.data)
            if correlation != 0:
                corr = 2 * correlation * left * right
                # This differs from multiplication because the correlation
                # term needs to be subtracted
                return np.abs(result_data) * np.sqrt(left**2 + right**2 - corr)
            else:
                return np.abs(result_data) * np.sqrt(left**2 + right**2)

    # TODO: Make a decorator that adds the docstrings instead.

    # Apply docstrings
    _propagate_add.__doc__ = _propagate_doc.format(operation='addition',
        operator='+', instance='StdDevUncertainty')
    _propagate_subtract.__doc__ = _propagate_doc.format(
        operation='subtraction', operator='-', instance='StdDevUncertainty')
    _propagate_multiply.__doc__ = _propagate_doc.format(
        operation='multiplication', operator='*', instance='StdDevUncertainty')
    _propagate_divide.__doc__ = _propagate_doc.format(operation='divison',
        operator='/', instance='StdDevUncertainty')
    _convert_uncertainty.__doc__ = NDUncertainty._convert_uncertainty.__doc__
