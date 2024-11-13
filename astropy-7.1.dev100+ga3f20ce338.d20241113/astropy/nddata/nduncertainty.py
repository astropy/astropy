# Licensed under a 3-clause BSD style license - see LICENSE.rst

import weakref
from abc import ABCMeta, abstractmethod
from copy import deepcopy

import numpy as np

# from astropy.utils.compat import ignored
from astropy import log
from astropy.units import Quantity, Unit, UnitConversionError

__all__ = [
    "MissingDataAssociationException",
    "IncompatibleUncertaintiesException",
    "NDUncertainty",
    "StdDevUncertainty",
    "UnknownUncertainty",
    "VarianceUncertainty",
    "InverseVariance",
]

# mapping from collapsing operations to the complementary methods used for `to_variance`
collapse_to_variance_mapping = {
    np.sum: np.square,
    np.mean: np.square,
}


def _move_preserved_axes_first(arr, preserve_axes):
    # When collapsing an ND array and preserving M axes, move the
    # preserved axes to the first M axes of the output. For example,
    # if arr.shape == (6, 5, 4, 3, 2) and we're preserving axes (1, 2),
    # then the output should have shape (20, 6, 3, 2). Axes 1 and 2 have
    # shape 5 and 4, so we take their product and put them both in the zeroth
    # axis.
    zeroth_axis_after_reshape = np.prod(np.array(arr.shape)[list(preserve_axes)])
    collapse_axes = [i for i in range(arr.ndim) if i not in preserve_axes]
    return arr.reshape(
        [zeroth_axis_after_reshape] + np.array(arr.shape)[collapse_axes].tolist()
    )


def _unravel_preserved_axes(arr, collapsed_arr, preserve_axes):
    # After reshaping an array with _move_preserved_axes_first and collapsing
    # the result, convert the reshaped first axis back into the shape of each
    # of the original preserved axes.
    # For example, if arr.shape == (6, 5, 4, 3, 2) and we're preserving axes (1, 2),
    # then the output of _move_preserved_axes_first should have shape (20, 6, 3, 2).
    # This method unravels the first axis in the output *after* a collapse, so the
    # output with shape (20,) becomes shape (5, 4).
    if collapsed_arr.ndim != len(preserve_axes):
        arr_shape = np.array(arr.shape)
        return collapsed_arr.reshape(arr_shape[np.asarray(preserve_axes)])
    return collapsed_arr


def from_variance_for_mean(x, axis):
    if axis is None:
        # do operation on all dimensions:
        denom = np.ma.count(x)
    else:
        denom = np.ma.count(x, axis)
    return np.sqrt(np.ma.sum(x, axis)) / denom


# mapping from collapsing operations to the complementary methods used for `from_variance`
collapse_from_variance_mapping = {
    np.sum: lambda x, axis: np.sqrt(np.ma.sum(x, axis)),
    np.mean: from_variance_for_mean,
    np.median: None,
}


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

    unit : unit-like, optional
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
            if unit is not None and unit != array._unit:
                # TODO : Clarify it (see NDData.init for same problem)?
                log.info("overwriting Uncertainty's current unit with specified unit.")
            elif array._unit is not None:
                unit = array.unit
            array = array.array

        elif isinstance(array, Quantity):
            # Check if two units are given and take the explicit one then.
            if unit is not None and array.unit is not None and unit != array.unit:
                log.info("overwriting Quantity's current unit with specified unit.")
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
        """`bool` : Supports uncertainty propagation with correlated uncertainties?

        .. versionadded:: 1.2
        """
        return False

    @property
    def array(self):
        """`numpy.ndarray` : the uncertainty's value."""
        return self._array

    @array.setter
    def array(self, value):
        if isinstance(value, (list, np.ndarray)):
            value = np.asarray(value)
        self._array = value

    @property
    def unit(self):
        """`~astropy.units.Unit` : The unit of the uncertainty, if any."""
        return self._unit

    @unit.setter
    def unit(self, value):
        """
        The unit should be set to a value consistent with the parent NDData
        unit and the uncertainty type.
        """
        if value is not None:
            # Check the hidden attribute below, not the property. The property
            # raises an exception if there is no parent_nddata.
            if self._parent_nddata is not None:
                parent_unit = self.parent_nddata.unit
                try:
                    # Check for consistency with the unit of the parent_nddata
                    self._data_unit_to_uncertainty_unit(parent_unit).to(value)
                except UnitConversionError:
                    raise UnitConversionError(
                        f"Unit {value} is incompatible with unit {parent_unit} of "
                        "parent nddata"
                    )

            self._unit = Unit(value)
        else:
            self._unit = value

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
        no_parent_message = "uncertainty is not associated with an NDData object"
        parent_lost_message = (
            "the associated NDData object was deleted and cannot be accessed "
            "anymore. You can prevent the NDData object from being deleted by "
            "assigning it to a variable. If this happened after unpickling "
            "make sure you pickle the parent not the uncertainty directly."
        )
        try:
            parent = self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException(no_parent_message)
        else:
            if parent is None:
                raise MissingDataAssociationException(no_parent_message)
            else:
                # The NDData is saved as weak reference so we must call it
                # to get the object the reference points to. However because
                # we have a weak reference here it's possible that the parent
                # was deleted because its reference count dropped to zero.
                if isinstance(self._parent_nddata, weakref.ref):
                    resolved_parent = self._parent_nddata()
                    if resolved_parent is None:
                        log.info(parent_lost_message)
                    return resolved_parent
                else:
                    log.info("parent_nddata should be a weakref to an NDData object.")
                    return self._parent_nddata

    @parent_nddata.setter
    def parent_nddata(self, value):
        if value is not None and not isinstance(value, weakref.ref):
            # Save a weak reference on the uncertainty that points to this
            # instance of NDData. Direct references should NOT be used:
            # https://github.com/astropy/astropy/pull/4799#discussion_r61236832
            value = weakref.ref(value)
        # Set _parent_nddata here and access below with the property because value
        # is a weakref
        self._parent_nddata = value
        # set uncertainty unit to that of the parent if it was not already set, unless initializing
        # with empty parent (Value=None)
        if value is not None:
            parent_unit = self.parent_nddata.unit
            # this will get the unit for masked quantity input:
            parent_data_unit = getattr(self.parent_nddata.data, "unit", None)
            if parent_unit is None and parent_data_unit is None:
                self.unit = None
            elif self.unit is None and parent_unit is not None:
                # Set the uncertainty's unit to the appropriate value
                self.unit = self._data_unit_to_uncertainty_unit(parent_unit)
            elif parent_data_unit is not None:
                # if the parent_nddata object has a unit, use it:
                self.unit = self._data_unit_to_uncertainty_unit(parent_data_unit)
            else:
                # Check that units of uncertainty are compatible with those of
                # the parent. If they are, no need to change units of the
                # uncertainty or the data. If they are not, let the user know.
                unit_from_data = self._data_unit_to_uncertainty_unit(parent_unit)
                try:
                    unit_from_data.to(self.unit)
                except UnitConversionError:
                    raise UnitConversionError(
                        f"Unit {self.unit} of uncertainty "
                        f"incompatible with unit {parent_unit} of "
                        "data"
                    )

    @abstractmethod
    def _data_unit_to_uncertainty_unit(self, value):
        """
        Subclasses must override this property. It should take in a data unit
        and return the correct unit for the uncertainty given the uncertainty
        type.
        """
        return None

    def __repr__(self):
        prefix = self.__class__.__name__ + "("
        try:
            body = np.array2string(self.array, separator=", ", prefix=prefix)
        except AttributeError:
            # In case it wasn't possible to use array2string
            body = str(self.array)
        return f"{prefix}{body})"

    def __getstate__(self):
        # Because of the weak reference the class wouldn't be picklable.
        try:
            return self._array, self._unit, self.parent_nddata
        except MissingDataAssociationException:
            # In case there's no parent
            return self._array, self._unit, None

    def __setstate__(self, state):
        if len(state) != 3:
            raise TypeError("The state should contain 3 items.")
        self._array = state[0]
        self._unit = state[1]

        parent = state[2]
        if parent is not None:
            parent = weakref.ref(parent)
        self._parent_nddata = parent

    def __getitem__(self, item):
        """Normal slicing on the array, keep the unit and return a reference."""
        return self.__class__(self.array[item], unit=self.unit, copy=False)

    def propagate(self, operation, other_nddata, result_data, correlation, axis=None):
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

        result_data : `~astropy.units.Quantity` or ndarray
            The result of the arithmetic operations on the data.

        correlation : `numpy.ndarray` or number
            The correlation (rho) is defined between the uncertainties in
            sigma_AB = sigma_A * sigma_B * rho. A value of ``0`` means
            uncorrelated operands.

        axis : int or tuple of ints, optional
            Axis over which to perform a collapsing operation.

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
                raise ValueError(
                    f"{type(self).__name__} does not support uncertainty propagation"
                    " with correlation."
                )

        if other_nddata is not None:
            # Get the other uncertainty (and convert it to a matching one)
            other_uncert = self._convert_uncertainty(other_nddata.uncertainty)

            if operation.__name__ == "add":
                result = self._propagate_add(other_uncert, result_data, correlation)
            elif operation.__name__ == "subtract":
                result = self._propagate_subtract(
                    other_uncert, result_data, correlation
                )
            elif operation.__name__ == "multiply":
                result = self._propagate_multiply(
                    other_uncert, result_data, correlation
                )
            elif operation.__name__ in ["true_divide", "divide"]:
                result = self._propagate_divide(other_uncert, result_data, correlation)
            else:
                raise ValueError(f"unsupported operation: {operation.__name__}")
        else:
            # assume this is a collapsing operation:
            result = self._propagate_collapse(operation, axis)

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

    def represent_as(self, other_uncert):
        """Convert this uncertainty to a different uncertainty type.

        Parameters
        ----------
        other_uncert : `NDUncertainty` subclass
            The `NDUncertainty` subclass to convert to.

        Returns
        -------
        resulting_uncertainty : `NDUncertainty` instance
            An instance of ``other_uncert`` subclass containing the uncertainty
            converted to the new uncertainty type.

        Raises
        ------
        TypeError
            If either the initial or final subclasses do not support
            conversion, a `TypeError` is raised.
        """
        as_variance = getattr(self, "_convert_to_variance", None)
        if as_variance is None:
            raise TypeError(
                f"{type(self)} does not support conversion to another uncertainty type."
            )
        from_variance = getattr(other_uncert, "_convert_from_variance", None)
        if from_variance is None:
            raise TypeError(
                f"{other_uncert.__name__} does not support conversion from "
                "another uncertainty type."
            )
        return from_variance(as_variance())


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
        """`False` : Uncertainty propagation is *not* possible for this class."""
        return False

    @property
    def uncertainty_type(self):
        """``"unknown"`` : `UnknownUncertainty` implements any unknown \
                           uncertainty type.
        """
        return "unknown"

    def _data_unit_to_uncertainty_unit(self, value):
        """
        No way to convert if uncertainty is unknown.
        """
        return None

    def _convert_uncertainty(self, other_uncert):
        """Raise an Exception because unknown uncertainty types cannot
        implement propagation.
        """
        msg = "Uncertainties of unknown type cannot be propagated."
        raise IncompatibleUncertaintiesException(msg)

    def _propagate_add(self, other_uncert, result_data, correlation):
        """Not possible for unknown uncertainty types."""
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

    def _propagate_collapse(self, numpy_op, axis=None):
        """
        Error propagation for collapse operations on variance or
        variance-like uncertainties. Uncertainties are calculated using the
        formulae for variance but can be used for uncertainty convertible to
        a variance.

        Parameters
        ----------
        numpy_op : function
            Numpy operation like `np.sum` or `np.max` to use in the collapse

        subtract : bool, optional
            If ``True``, propagate for subtraction, otherwise propagate for
            addition.

        axis : tuple, optional
            Axis on which to compute collapsing operations.
        """
        try:
            result_unit_sq = self.parent_nddata.unit**2
        except (AttributeError, TypeError):
            result_unit_sq = None

        if self.array is not None:
            # Formula: sigma**2 = dA

            if numpy_op in [np.min, np.max]:
                # Find the indices of the min/max in parent data along each axis,
                # return the uncertainty at the corresponding entry:
                return self._get_err_at_extremum(numpy_op, axis=axis)

            # np.sum and np.mean operations use similar pattern
            # to `_propagate_add_sub`, for example:
            else:
                # lookup the mapping for to_variance and from_variance for this
                # numpy operation:
                to_variance = collapse_to_variance_mapping[numpy_op]
                from_variance = collapse_from_variance_mapping[numpy_op]
                masked_uncertainty = np.ma.masked_array(
                    self.array, self.parent_nddata.mask
                )
                if (
                    self.unit is not None
                    and to_variance(self.unit) != self.parent_nddata.unit**2
                ):
                    # If the uncertainty has a different unit than the result we
                    # need to convert it to the results unit.
                    this = (
                        to_variance(masked_uncertainty << self.unit)
                        .to(result_unit_sq)
                        .value
                    )
                else:
                    this = to_variance(masked_uncertainty)

                return from_variance(this, axis=axis)

    def _get_err_at_extremum(self, extremum, axis):
        """
        Return the value of the ``uncertainty`` array at the indices
        which satisfy the ``extremum`` function applied to the ``measurement`` array,
        where we expect ``extremum`` to be np.argmax or np.argmin, and
        we expect a two-dimensional output.

        Assumes the ``measurement`` and ``uncertainty`` array dimensions
        are ordered such that the zeroth dimension is the one to preserve.
        For example, if you start with array with shape (a, b, c), this
        function applies the ``extremum`` function to the last two dimensions,
        with shapes b and c.

        This operation is difficult to cast in a vectorized way. Here
        we implement it with a list comprehension, which is likely not the
        most performant solution.
        """
        if axis is not None and not hasattr(axis, "__len__"):
            # this is a single axis:
            axis = [axis]

        if extremum is np.min:
            arg_extremum = np.ma.argmin
        elif extremum is np.max:
            arg_extremum = np.ma.argmax

        all_axes = np.arange(self.array.ndim)

        if axis is None:
            # collapse over all dimensions
            ind = arg_extremum(np.asanyarray(self.parent_nddata).ravel())
            return self.array.ravel()[ind]

        # collapse an ND array over arbitrary dimensions:
        preserve_axes = [ax for ax in all_axes if ax not in axis]
        meas = np.ma.masked_array(
            _move_preserved_axes_first(self.parent_nddata.data, preserve_axes),
            _move_preserved_axes_first(self.parent_nddata.mask, preserve_axes),
        )
        err = _move_preserved_axes_first(self.array, preserve_axes)

        result = np.array(
            [e[np.unravel_index(arg_extremum(m), m.shape)] for m, e in zip(meas, err)]
        )

        return _unravel_preserved_axes(
            self.parent_nddata.data,
            result,
            preserve_axes,
        )

    def _propagate_add_sub(
        self,
        other_uncert,
        result_data,
        correlation,
        subtract=False,
        to_variance=lambda x: x,
        from_variance=lambda x: x,
    ):
        """
        Error propagation for addition or subtraction of variance or
        variance-like uncertainties. Uncertainties are calculated using the
        formulae for variance but can be used for uncertainty convertible to
        a variance.

        Parameters
        ----------
        other_uncert : `~astropy.nddata.NDUncertainty` instance
            The uncertainty, if any, of the other operand.

        result_data : `~astropy.nddata.NDData` instance
            The results of the operation on the data.

        correlation : float or array-like
            Correlation of the uncertainties.

        subtract : bool, optional
            If ``True``, propagate for subtraction, otherwise propagate for
            addition.

        to_variance : function, optional
            Function that will transform the input uncertainties to variance.
            The default assumes the uncertainty is the variance.

        from_variance : function, optional
            Function that will convert from variance to the input uncertainty.
            The default assumes the uncertainty is the variance.
        """
        if subtract:
            correlation_sign = -1
        else:
            correlation_sign = 1

        try:
            result_unit_sq = result_data.unit**2
        except AttributeError:
            result_unit_sq = None

        if other_uncert.array is not None:
            # Formula: sigma**2 = dB
            if other_uncert.unit is not None and result_unit_sq != to_variance(
                other_uncert.unit
            ):
                # If the other uncertainty has a unit and this unit differs
                # from the unit of the result convert it to the results unit
                other = (
                    to_variance(other_uncert.array << other_uncert.unit)
                    .to(result_unit_sq)
                    .value
                )
            else:
                other = to_variance(other_uncert.array)
        else:
            other = 0

        if self.array is not None:
            # Formula: sigma**2 = dA

            if (
                self.unit is not None
                and to_variance(self.unit) != self.parent_nddata.unit**2
            ):
                # If the uncertainty has a different unit than the result we
                # need to convert it to the results unit.
                this = to_variance(self.array << self.unit).to(result_unit_sq).value
            else:
                this = to_variance(self.array)
        else:
            this = 0

        # Formula: sigma**2 = dA + dB +/- 2*cor*sqrt(dA*dB)
        # Formula: sigma**2 = sigma_other + sigma_self +/- 2*cor*sqrt(dA*dB)
        #     (sign depends on whether addition or subtraction)

        # Determine the result depending on the correlation
        if isinstance(correlation, np.ndarray) or correlation != 0:
            corr = 2 * correlation * np.sqrt(this * other)
            result = this + other + correlation_sign * corr
        else:
            result = this + other

        return from_variance(result)

    def _propagate_multiply_divide(
        self,
        other_uncert,
        result_data,
        correlation,
        divide=False,
        to_variance=lambda x: x,
        from_variance=lambda x: x,
    ):
        """
        Error propagation for multiplication or division of variance or
        variance-like uncertainties. Uncertainties are calculated using the
        formulae for variance but can be used for uncertainty convertible to
        a variance.

        Parameters
        ----------
        other_uncert : `~astropy.nddata.NDUncertainty` instance
            The uncertainty, if any, of the other operand.

        result_data : `~astropy.nddata.NDData` instance
            The results of the operation on the data.

        correlation : float or array-like
            Correlation of the uncertainties.

        divide : bool, optional
            If ``True``, propagate for division, otherwise propagate for
            multiplication.

        to_variance : function, optional
            Function that will transform the input uncertainties to variance.
            The default assumes the uncertainty is the variance.

        from_variance : function, optional
            Function that will convert from variance to the input uncertainty.
            The default assumes the uncertainty is the variance.
        """
        # For multiplication we don't need the result as quantity
        if isinstance(result_data, Quantity):
            result_data = result_data.value

        if divide:
            correlation_sign = -1
        else:
            correlation_sign = 1

        if other_uncert.array is not None:
            # We want the result to have a unit consistent with the parent, so
            # we only need to convert the unit of the other uncertainty if it
            # is different from its data's unit.
            if (
                other_uncert.unit
                and to_variance(1 * other_uncert.unit)
                != ((1 * other_uncert.parent_nddata.unit) ** 2).unit
            ):
                d_b = (
                    to_variance(other_uncert.array << other_uncert.unit)
                    .to((1 * other_uncert.parent_nddata.unit) ** 2)
                    .value
                )
            else:
                d_b = to_variance(other_uncert.array)
            # Formula: sigma**2 = |A|**2 * d_b
            right = np.abs(self.parent_nddata.data**2 * d_b)
        else:
            right = 0

        if self.array is not None:
            # Just the reversed case
            if (
                self.unit
                and to_variance(1 * self.unit)
                != ((1 * self.parent_nddata.unit) ** 2).unit
            ):
                d_a = (
                    to_variance(self.array << self.unit)
                    .to((1 * self.parent_nddata.unit) ** 2)
                    .value
                )
            else:
                d_a = to_variance(self.array)
            # Formula: sigma**2 = |B|**2 * d_a
            left = np.abs(other_uncert.parent_nddata.data**2 * d_a)
        else:
            left = 0

        # Multiplication
        #
        # The fundamental formula is:
        #   sigma**2 = |AB|**2*(d_a/A**2+d_b/B**2+2*sqrt(d_a)/A*sqrt(d_b)/B*cor)
        #
        # This formula is not very handy since it generates NaNs for every
        # zero in A and B. So we rewrite it:
        #
        # Multiplication Formula:
        #   sigma**2 = (d_a*B**2 + d_b*A**2 + (2 * cor * ABsqrt(dAdB)))
        #   sigma**2 = (left + right + (2 * cor * ABsqrt(dAdB)))
        #
        # Division
        #
        # The fundamental formula for division is:
        #   sigma**2 = |A/B|**2*(d_a/A**2+d_b/B**2-2*sqrt(d_a)/A*sqrt(d_b)/B*cor)
        #
        # As with multiplication, it is convenient to rewrite this to avoid
        # nans where A is zero.
        #
        # Division formula (rewritten):
        #   sigma**2 = d_a/B**2 + (A/B)**2 * d_b/B**2
        #                   - 2 * cor * A *sqrt(dAdB) / B**3
        #   sigma**2 = d_a/B**2 + (A/B)**2 * d_b/B**2
        #                   - 2*cor * sqrt(d_a)/B**2  * sqrt(d_b) * A / B
        #   sigma**2 = multiplication formula/B**4 (and sign change in
        #               the correlation)

        if isinstance(correlation, np.ndarray) or correlation != 0:
            corr = (
                2
                * correlation
                * np.sqrt(d_a * d_b)
                * self.parent_nddata.data
                * other_uncert.parent_nddata.data
            )
        else:
            corr = 0

        if divide:
            return from_variance(
                (left + right + correlation_sign * corr)
                / other_uncert.parent_nddata.data**4
            )
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
        >>> ndd = NDData([1,2,3], unit='m',
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
        """``"std"`` : `StdDevUncertainty` implements standard deviation."""
        return "std"

    def _convert_uncertainty(self, other_uncert):
        if isinstance(other_uncert, StdDevUncertainty):
            return other_uncert
        else:
            raise IncompatibleUncertaintiesException

    def _propagate_add(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert,
            result_data,
            correlation,
            subtract=False,
            to_variance=np.square,
            from_variance=np.sqrt,
        )

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert,
            result_data,
            correlation,
            subtract=True,
            to_variance=np.square,
            from_variance=np.sqrt,
        )

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert,
            result_data,
            correlation,
            divide=False,
            to_variance=np.square,
            from_variance=np.sqrt,
        )

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert,
            result_data,
            correlation,
            divide=True,
            to_variance=np.square,
            from_variance=np.sqrt,
        )

    def _propagate_collapse(self, numpy_operation, axis):
        # defer to _VariancePropagationMixin
        return super()._propagate_collapse(numpy_operation, axis=axis)

    def _data_unit_to_uncertainty_unit(self, value):
        return value

    def _convert_to_variance(self):
        new_array = None if self.array is None else self.array**2
        new_unit = None if self.unit is None else self.unit**2
        return VarianceUncertainty(new_array, unit=new_unit)

    @classmethod
    def _convert_from_variance(cls, var_uncert):
        new_array = None if var_uncert.array is None else var_uncert.array ** (1 / 2)
        new_unit = None if var_uncert.unit is None else var_uncert.unit ** (1 / 2)
        return cls(new_array, unit=new_unit)


class VarianceUncertainty(_VariancePropagationMixin, NDUncertainty):
    """
    Variance uncertainty assuming first order Gaussian error
    propagation.

    This class implements uncertainty propagation for ``addition``,
    ``subtraction``, ``multiplication`` and ``division`` with other instances
    of `VarianceUncertainty`. The class can handle if the uncertainty has a
    unit that differs from (but is convertible to) the parents `NDData` unit.
    The unit of the resulting uncertainty will be the square of the unit of the
    resulting data. Also support for correlation is possible but requires the
    correlation as input. It cannot handle correlation determination itself.

    Parameters
    ----------
    args, kwargs :
        see `NDUncertainty`

    Examples
    --------
    Compare this example to that in `StdDevUncertainty`; the uncertainties
    in the examples below are equivalent to the uncertainties in
    `StdDevUncertainty`.

    `VarianceUncertainty` should always be associated with an `NDData`-like
    instance, either by creating it during initialization::

        >>> from astropy.nddata import NDData, VarianceUncertainty
        >>> ndd = NDData([1,2,3], unit='m',
        ...              uncertainty=VarianceUncertainty([0.01, 0.01, 0.01]))
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        VarianceUncertainty([0.01, 0.01, 0.01])

    or by setting it manually on the `NDData` instance::

        >>> ndd.uncertainty = VarianceUncertainty([0.04], unit='m^2', copy=True)
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        VarianceUncertainty([0.04])

    the uncertainty ``array`` can also be set directly::

        >>> ndd.uncertainty.array = 4
        >>> ndd.uncertainty
        VarianceUncertainty(4)

    .. note::
        The unit will not be displayed.
    """

    @property
    def uncertainty_type(self):
        """``"var"`` : `VarianceUncertainty` implements variance."""
        return "var"

    @property
    def supports_correlated(self):
        """`True` : `VarianceUncertainty` allows to propagate correlated \
                    uncertainties.

        ``correlation`` must be given, this class does not implement computing
        it by itself.
        """
        return True

    def _propagate_add(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert, result_data, correlation, subtract=False
        )

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert, result_data, correlation, subtract=True
        )

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert, result_data, correlation, divide=False
        )

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert, result_data, correlation, divide=True
        )

    def _data_unit_to_uncertainty_unit(self, value):
        return value**2

    def _convert_to_variance(self):
        return self

    @classmethod
    def _convert_from_variance(cls, var_uncert):
        return var_uncert


def _inverse(x):
    """Just a simple inverse for use in the InverseVariance."""
    return 1 / x


class InverseVariance(_VariancePropagationMixin, NDUncertainty):
    """
    Inverse variance uncertainty assuming first order Gaussian error
    propagation.

    This class implements uncertainty propagation for ``addition``,
    ``subtraction``, ``multiplication`` and ``division`` with other instances
    of `InverseVariance`. The class can handle if the uncertainty has a unit
    that differs from (but is convertible to) the parents `NDData` unit. The
    unit of the resulting uncertainty will the inverse square of the unit of
    the resulting data. Also support for correlation is possible but requires
    the correlation as input. It cannot handle correlation determination
    itself.

    Parameters
    ----------
    args, kwargs :
        see `NDUncertainty`

    Examples
    --------
    Compare this example to that in `StdDevUncertainty`; the uncertainties
    in the examples below are equivalent to the uncertainties in
    `StdDevUncertainty`.

    `InverseVariance` should always be associated with an `NDData`-like
    instance, either by creating it during initialization::

        >>> from astropy.nddata import NDData, InverseVariance
        >>> ndd = NDData([1,2,3], unit='m',
        ...              uncertainty=InverseVariance([100, 100, 100]))
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        InverseVariance([100, 100, 100])

    or by setting it manually on the `NDData` instance::

        >>> ndd.uncertainty = InverseVariance([25], unit='1/m^2', copy=True)
        >>> ndd.uncertainty  # doctest: +FLOAT_CMP
        InverseVariance([25])

    the uncertainty ``array`` can also be set directly::

        >>> ndd.uncertainty.array = 0.25
        >>> ndd.uncertainty
        InverseVariance(0.25)

    .. note::
        The unit will not be displayed.
    """

    @property
    def uncertainty_type(self):
        """``"ivar"`` : `InverseVariance` implements inverse variance."""
        return "ivar"

    @property
    def supports_correlated(self):
        """`True` : `InverseVariance` allows to propagate correlated \
                    uncertainties.

        ``correlation`` must be given, this class does not implement computing
        it by itself.
        """
        return True

    def _propagate_add(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert,
            result_data,
            correlation,
            subtract=False,
            to_variance=_inverse,
            from_variance=_inverse,
        )

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return super()._propagate_add_sub(
            other_uncert,
            result_data,
            correlation,
            subtract=True,
            to_variance=_inverse,
            from_variance=_inverse,
        )

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert,
            result_data,
            correlation,
            divide=False,
            to_variance=_inverse,
            from_variance=_inverse,
        )

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return super()._propagate_multiply_divide(
            other_uncert,
            result_data,
            correlation,
            divide=True,
            to_variance=_inverse,
            from_variance=_inverse,
        )

    def _data_unit_to_uncertainty_unit(self, value):
        return 1 / value**2

    def _convert_to_variance(self):
        new_array = None if self.array is None else 1 / self.array
        new_unit = None if self.unit is None else 1 / self.unit
        return VarianceUncertainty(new_array, unit=new_unit)

    @classmethod
    def _convert_from_variance(cls, var_uncert):
        new_array = None if var_uncert.array is None else 1 / var_uncert.array
        new_unit = None if var_uncert.unit is None else 1 / var_uncert.unit
        return cls(new_array, unit=new_unit)
