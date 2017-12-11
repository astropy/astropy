# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np
from abc import ABCMeta, abstractmethod
from copy import deepcopy
import weakref

# from ..utils.compat import ignored
from .. import log
from ..units import Unit, Quantity

__all__ = ['MissingDataAssociationException',
           'IncompatibleUncertaintiesException', 'NDUncertainty',
           'StdDevUncertainty', 'UnknownUncertainty']


class IncompatibleUncertaintiesException(Exception):
    """This exception should be used to indicate cases in which uncertainties
    with two different classes can not be propagated.
    """


class MissingDataAssociationException(Exception):
    """This exception should be used to indicate that an uncertainty instance
    has not been associated with a parent `~astropy.nddata.NDData` object.
    """


class NDUncertainty(metaclass=ABCMeta):
    """This is the metaclass for uncertainty classes used with `NDData`.

    Parameters
    ----------
    array : any type, optional
        The array or value (the parameter name is due to historical reasons) of
        the uncertainty. `numpy.ndarray`, `~astropy.units.Quantity` or
        `NDUncertainty` subclasses are recommended.
        If the `array` is `list`-like or `numpy.ndarray`-like it will be cast
        to a plain `numpy.ndarray`.
        Default is ``None``.

    unit : `~astropy.units.Unit` or str, optional
        Unit for the uncertainty ``array``. Strings that can be converted to a
        `~astropy.units.Unit` are allowed.
        Default is ``None``.

    copy : `bool`, optional
        Indicates whether to save the `array` as a copy. ``True`` copies it
        before saving, while ``False`` tries to save every parameter as
        reference. Note however that it is not always possible to save the
        input as reference.
        Default is ``True``.

    Raises
    ------
    IncompatibleUncertaintiesException
        If given another `NDUncertainty`-like class as ``array`` if their
        ``uncertainty_type`` is different.
    """

    def __init__(self, array=None, copy=True, unit=None):
        if isinstance(array, NDUncertainty):
            # Given an NDUncertainty class or subclass check that the type
            # is the same.
            if array.uncertainty_type != self.uncertainty_type:
                raise IncompatibleUncertaintiesException
            # Check if two units are given and take the explicit one then.
            if (unit is not None and unit != array._unit):
                # TODO : Clarify it (see NDData.init for same problem)?
                log.info("overwriting Uncertainty's current "
                         "unit with specified unit.")
            elif array._unit is not None:
                unit = array.unit
            array = array.array

        elif isinstance(array, Quantity):
            # Check if two units are given and take the explicit one then.
            if (unit is not None and array.unit is not None and
                    unit != array.unit):
                log.info("overwriting Quantity's current "
                         "unit with specified unit.")
            elif array.unit is not None:
                unit = array.unit
            array = array.value

        if unit is None:
            self._unit = None
        else:
            self._unit = Unit(unit)

        if copy:
            array = deepcopy(array)
            unit = deepcopy(unit)

        self.array = array
        self.parent_nddata = None  # no associated NDData - until it is set!

    @property
    @abstractmethod
    def uncertainty_type(self):
        """`str` : Short description of the type of uncertainty.

        Defined as abstract property so subclasses *have* to override this.
        """
        return None

    @property
    def supports_correlated(self):
        """`bool` : Supports uncertainty propagation with correlated \
                 uncertainties?

        .. versionadded:: 1.2
        """
        return False

    @property
    def array(self):
        """`numpy.ndarray` : the uncertainty's value.
        """
        return self._array

    @array.setter
    def array(self, value):
        if isinstance(value, (list, np.ndarray)):
            value = np.array(value, subok=False, copy=False)
        self._array = value

    @property
    def unit(self):
        """`~astropy.units.Unit` : The unit of the uncertainty, if any.

        Even though it is not enforced the unit should be convertible to the
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

    @property
    def quantity(self):
        """
        This uncertainty as an `~astropy.units.Quantity` object.
        """
        return Quantity(self.array, self.unit, copy=False, dtype=self.array.dtype)

    @property
    def parent_nddata(self):
        """`NDData` : reference to `NDData` instance with this uncertainty.

        In case the reference is not set uncertainty propagation will not be
        possible since propagation might need the uncertain data besides the
        uncertainty.
        """
        message = "uncertainty is not associated with an NDData object"
        try:
            if self._parent_nddata is None:
                raise MissingDataAssociationException(message)
            else:
                # The NDData is saved as weak reference so we must call it
                # to get the object the reference points to.
                if isinstance(self._parent_nddata, weakref.ref):
                    return self._parent_nddata()
                else:
                    log.info("parent_nddata should be a weakref to an NDData "
                             "object.")
                    return self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException(message)

    @parent_nddata.setter
    def parent_nddata(self, value):
        if value is not None and not isinstance(value, weakref.ref):
            # Save a weak reference on the uncertainty that points to this
            # instance of NDData. Direct references should NOT be used:
            # https://github.com/astropy/astropy/pull/4799#discussion_r61236832
            value = weakref.ref(value)
        self._parent_nddata = value

    def __repr__(self):
        prefix = self.__class__.__name__ + '('
        try:
            body = np.array2string(self.array, separator=', ', prefix=prefix)
        except AttributeError:
            # In case it wasn't possible to use array2string
            body = str(self.array)
        return ''.join([prefix, body, ')'])

    def __getitem__(self, item):
        """Normal slicing on the array, keep the unit and return a reference.
        """
        return self.__class__(self.array[item], unit=self.unit, copy=False)

    def propagate(self, operation, other_nddata, result_data, correlation):
        """Calculate the resulting uncertainty given an operation on the data.

        .. versionadded:: 1.2

        Parameters
        ----------
        operation : callable
            The operation that is performed on the `NDData`. Supported are
            `numpy.add`, `numpy.subtract`, `numpy.multiply` and
            `numpy.true_divide` (or `numpy.divide`).

        other_nddata : `NDData` instance
            The second operand in the arithmetic operation.

        result_data : `~astropy.units.Quantity` or `numpy.ndarray`
            The result of the arithmetic operations on the data.

        correlation : `numpy.ndarray` or number
            The correlation (rho) is defined between the uncertainties in
            sigma_AB = sigma_A * sigma_B * rho. A value of ``0`` means
            uncorrelated operands.

        Returns
        -------
        resulting_uncertainty : `NDUncertainty` instance
            Another instance of the same `NDUncertainty` subclass containing
            the uncertainty of the result.

        Raises
        ------
        ValueError
            If the ``operation`` is not supported or if correlation is not zero
            but the subclass does not support correlated uncertainties.

        Notes
        -----
        First this method checks if a correlation is given and the subclass
        implements propagation with correlated uncertainties.
        Then the second uncertainty is converted (or an Exception is raised)
        to the same class in order to do the propagation.
        Then the appropriate propagation method is invoked and the result is
        returned.
        """
        # Check if the subclass supports correlation
        if not self.supports_correlated:
            if isinstance(correlation, np.ndarray) or correlation != 0:
                raise ValueError("{0} does not support uncertainty propagation"
                                 " with correlation."
                                 "".format(self.__class__.__name__))

        # Get the other uncertainty (and convert it to a matching one)
        other_uncert = self._convert_uncertainty(other_nddata.uncertainty)

        if operation.__name__ == 'add':
            result = self._propagate_add(other_uncert, result_data,
                                         correlation)
        elif operation.__name__ == 'subtract':
            result = self._propagate_subtract(other_uncert, result_data,
                                              correlation)
        elif operation.__name__ == 'multiply':
            result = self._propagate_multiply(other_uncert, result_data,
                                              correlation)
        elif operation.__name__ in ['true_divide', 'divide']:
            result = self._propagate_divide(other_uncert, result_data,
                                            correlation)
        else:
            raise ValueError('unsupported operation')

        return self.__class__(result, copy=False)

    def _convert_uncertainty(self, other_uncert):
        """Checks if the uncertainties are compatible for propagation.

        Checks if the other uncertainty is `NDUncertainty`-like and if so
        verify that the uncertainty_type is equal. If the latter is not the
        case try returning ``self.__class__(other_uncert)``.

        Parameters
        ----------
        other_uncert : `NDUncertainty` subclass
            The other uncertainty.

        Returns
        -------
        other_uncert : `NDUncertainty` subclass
            but converted to a compatible `NDUncertainty` subclass if
            possible and necessary.

        Raises
        ------
        IncompatibleUncertaintiesException:
            If the other uncertainty cannot be converted to a compatible
            `NDUncertainty` subclass.
        """
        if isinstance(other_uncert, NDUncertainty):
            if self.uncertainty_type == other_uncert.uncertainty_type:
                return other_uncert
            else:
                return self.__class__(other_uncert)
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


class UnknownUncertainty(NDUncertainty):
    """This class implements any unknown uncertainty type.

    The main purpose of having an unknown uncertainty class is to prevent
    uncertainty propagation.

    Parameters
    ----------
    args, kwargs :
        see `NDUncertainty`
    """

    @property
    def supports_correlated(self):
        """`False` : Uncertainty propagation is *not* possible for this class.
        """
        return False

    @property
    def uncertainty_type(self):
        """``"unknown"`` : `UnknownUncertainty` implements any unknown \
                           uncertainty type.
        """
        return 'unknown'

    def _convert_uncertainty(self, other_uncert):
        """Raise an Exception because unknown uncertainty types cannot
        implement propagation.
        """
        msg = "Uncertainties of unknown type cannot be propagated."
        raise IncompatibleUncertaintiesException(msg)

    def _propagate_add(self, other_uncert, result_data, correlation):
        """Not possible for unknown uncertainty types.
        """
        return None

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return None

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return None

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return None


class _VariancePropagationMixin:
    """
    Propagation of uncertainties for variances, also used to perform error
    propagation for variance-like uncertainties (standard deviation and inverse
    variance).
    """
    def _propagate_add_sub(self, other_uncert, result_data, correlation,
                           add_or_subtract,
                           to_variance=lambda x: x, from_variance=lambda x: x):

        if self.array is None:
            # Formula: sigma = dB
            if (other_uncert.unit is not None and
                    (result_data.unit**2 != to_variance(other_uncert.unit))):
                # If the other uncertainty has a unit and this unit differs
                # from the unit of the result convert it to the results unit
                return from_variance(to_variance(other_uncert.array *
                                                 other_uncert.unit).to(
                                                 result_data.unit**2).value)
            else:
                # Copy the result because _propagate will not copy it but for
                # arithmetic operations users will expect copies.
                return deepcopy(other_uncert.array)

        elif other_uncert.array is None:
            # Formula: sigma = dA

            if (self.unit is not None and
                    to_variance(self.unit) != self.parent_nddata.unit**2):
                # If the uncertainty has a different unit than the result we
                # need to convert it to the results unit.
                return from_variance(to_variance(self.array * self.unit).to(result_data.unit**2))
            else:
                # Copy the result because _propagate will not copy it but for
                # arithmetic operations users will expect copies.
                return deepcopy(self.array)

        else:
            # Formula: sigma = sqrt(dA**2 + dB**2 + 2*cor*dA*dB)

            # Calculate: dA (this) and dB (other)
            if self.unit != other_uncert.unit:
                # In case the two uncertainties (or data) have different units
                # we need to use quantity operations. The case where only one
                # has a unit and the other doesn't is not possible with
                # addition and would have raised an exception in the data
                # computation
                this = to_variance(self.array * self.unit)
                other = to_variance(other_uncert.array * other_uncert.unit)
            else:
                # Since both units are the same or None we can just use
                # numpy operations
                this = to_variance(self.array)
                other = to_variance(other_uncert.array)

            # Determine the result depending on the correlation
            if isinstance(correlation, np.ndarray) or correlation != 0:
                corr = 2 * correlation * np.sqrt(this * other)
                result = this + other + add_or_subtract * corr
            else:
                result = this + other

            if isinstance(result, Quantity):
                # In case we worked with quantities we need to return the
                # uncertainty that has the same unit as the resulting data.
                # Note that this call is fast if the units are the same.
                return from_variance(result.to_value(result_data.unit**2))
            else:
                return from_variance(result)

    def _propagate_multiply_divide(self, other_uncert, result_data,
                                   correlation,
                                   divide=False,
                                   to_variance=lambda x: x,
                                   from_variance=lambda x: x):

        # For multiplication we don't need the result as quantity
        if isinstance(result_data, Quantity):
            result_data = result_data.value

        if divide:
            correlation_sign = -1
        else:
            correlation_sign = 1

        if self.array is None:
            # We want the result to have a unit consistent with the parent, so
            # we only need to convert the unit of the other uncertainty if it
            # is different from its data's unit.
            if (other_uncert.unit and
                to_variance(1 * other_uncert.unit) !=
                    ((1 * other_uncert.parent_nddata.unit)**2).unit):
                other = to_variance((other_uncert.array * other_uncert.unit)).to(
                    (1 * other_uncert.parent_nddata.unit)**2).value
            else:
                other = to_variance(other_uncert.array)
            if divide:
                # Formula: sigma**2 = (A / B)**2 * (dB / B**2)
                return from_variance(result_data**2 * other /
                                     other_uncert.parent_nddata.data**2)
            else:
                # Formula: sigma**2 = |A|**2 * dB
                return from_variance(np.abs(self.parent_nddata.data**2 * other))

        elif other_uncert.array is None:
            # Just the reversed case
            if (self.unit and
                to_variance(1 * self.unit) !=
                    ((1 * self.parent_nddata.unit)**2).unit):
                this = to_variance(self.array * self.unit).to(
                    (1 * self.parent_nddata.unit)**2).value
            else:
                this = to_variance(self.array)
            if divide:
                # Formula: sigma**2 = dA / B**2
                return from_variance(this / other_uncert.parent_nddata.data**2)
            else:
                # Formula: sigma**2 = dA**2 * |B|**2
                return from_variance(np.abs(other_uncert.parent_nddata.data**2 * this))
        else:
            # Multiplication Formulae
            #   sigma**2 = |AB|**2*(dA/A**2+dB/B**2+2*sqrt(dA)/A*sqrt(dB)/B*cor)

            # This formula is not very handy since it generates NaNs for every
            # zero in A and B. So we rewrite it:

            # Multiplication Formula:
            #   sigma**2 = (dA*B**2 + dB*A**2 + (2 * cor * ABsqrt(dAdB)))

            # Division formula (rewritten):
            #   sigma**2 = dA/B**2 + (A/B)**2 * dB/B**2
            #                   - 2 * cor * A *sqrt(dAdB) / B**3
            #   sigma**2 = dA/B**2 + (A/B)**2 * dB/B**2
            #                   - 2*cor * sqrt(dA)/B**2  * sqrt(dB) * A / B
            #   sigma**2 = multiplication formula/B**4 (and sign change in
            #               the correlation)

            if (self.unit and
                    to_variance(self.unit) != (self.parent_nddata.unit)**2):
                # To get the unit right we need to convert the unit of
                # each uncertainty to the same unit as it's parent
                dA = (to_variance(self.array * self.unit).to(
                    self.parent_nddata.unit**2).value)
            else:
                dA = to_variance(self.array)
            # Calculate: dA * B**2 (left)
            left = dA * other_uncert.parent_nddata.data**2

            if (other_uncert.unit and
                    to_variance(other_uncert.unit) != other_uncert.parent_nddata.unit**2):
                # Make sure the unit of dB is the same as its parent
                dB = (to_variance(other_uncert.array * other_uncert.unit).to(
                         other_uncert.parent_nddata.unit**2).value)
            else:
                dB = to_variance(other_uncert.array)

            # Calculate: dB * A**2 (right)
            right = dB * self.parent_nddata.data**2

            if isinstance(correlation, np.ndarray) or correlation != 0:
                corr = (2 * correlation * np.sqrt(dA * dB) *
                        self.parent_nddata.data *
                        other_uncert.parent_nddata.data)
            else:
                corr = 0
            if divide:
                return from_variance((left + right + correlation_sign * corr) /
                                     other_uncert.parent_nddata.data**4)
            else:
                return from_variance(left + right + correlation_sign * corr)


class StdDevUncertainty(_VariancePropagationMixin, NDUncertainty):
    """Standard deviation uncertainty assuming first order gaussian error
    propagation.

    This class implements uncertainty propagation for ``addition``,
    ``subtraction``, ``multiplication`` and ``division`` with other instances
    of `StdDevUncertainty`. The class can handle if the uncertainty has a
    unit that differs from (but is convertible to) the parents `NDData` unit.
    The unit of the resulting uncertainty will have the same unit as the
    resulting data. Also support for correlation is possible but requires the
    correlation as input. It cannot handle correlation determination itself.

    Parameters
    ----------
    args, kwargs :
        see `NDUncertainty`

    Examples
    --------
    `StdDevUncertainty` should always be associated with an `NDData`-like
    instance, either by creating it during initialization::

        >>> from astropy.nddata import NDData, StdDevUncertainty
        >>> ndd = NDData([1,2,3],
        ...              uncertainty=StdDevUncertainty([0.1, 0.1, 0.1]))
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        StdDevUncertainty([0.1, 0.1, 0.1])

    or by setting it manually on the `NDData` instance::

        >>> ndd.uncertainty = StdDevUncertainty([0.2], unit='m', copy=True)
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        StdDevUncertainty([0.2])

    the uncertainty ``array`` can also be set directly::

        >>> ndd.uncertainty.array = 2
        >>> ndd.uncertainty
        StdDevUncertainty(2)

    .. note::
        The unit will not be displayed.
    """

    @property
    def supports_correlated(self):
        """`True` : `StdDevUncertainty` allows to propagate correlated \
                    uncertainties.

        ``correlation`` must be given, this class does not implement computing
        it by itself.
        """
        return True

    @property
    def uncertainty_type(self):
        """``"std"`` : `StdDevUncertainty` implements standard deviation.
        """
        return 'std'

    def _convert_uncertainty(self, other_uncert):
        if isinstance(other_uncert, StdDevUncertainty):
            return other_uncert
        else:
            raise IncompatibleUncertaintiesException

    def _propagate_add(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(other_uncert, result_data,
                                          correlation, 1,
                                          to_variance=lambda x: x**2,
                                          from_variance=np.sqrt)

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(other_uncert, result_data,
                                          correlation, -1,
                                          to_variance=lambda x: x**2,
                                          from_variance=np.sqrt)

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(other_uncert,
                                                  result_data, correlation,
                                                  divide=False,
                                                  to_variance=lambda x: x**2,
                                                  from_variance=np.sqrt)

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(other_uncert,
                                                  result_data, correlation,
                                                  divide=True,
                                                  to_variance=lambda x: x**2,
                                                  from_variance=np.sqrt)


