# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Distribution class and associated machinery.
"""
import numpy as np

from astropy import units as u
from astropy.utils.misc import isiterable


__all__ = ['Uncertainty', 'Variable']


# Derivatives of ufuncs relative to their input(s).
UFUNC_DERIVATIVES = {
    np.add: (lambda x, y: 1.,
             lambda x, y: 1.),
    np.subtract: (lambda x, y: 1.,
                  lambda x, y: -1.),
    np.multiply: (lambda x, y: y,
                  lambda x, y: x),
    np.true_divide: (lambda x, y: 1/y,
                     lambda x, y: -x/y**2),
    np.arccos: lambda x: -1/np.sqrt(1-x**2),
    np.arccosh: lambda x: 1/np.sqrt(x**2-1),
    np.arcsin: lambda x: 1/np.sqrt(1-x**2),
    np.arcsinh: lambda x: 1/np.sqrt(1+x**2),
    np.arctan: lambda x: 1/(1+x**2),
    np.arctan2: (lambda y, x: x/(x**2+y**2),
                 lambda y, x: -y/(x**2+y**2)),
    np.arctanh: lambda x: 1/(1-x**2),
    np.copysign: (lambda x, y: np.where(x >= 0, np.copysign(1, y),
                                        -np.copysign(1, y)),
                  lambda x, y: 0),
    np.cos: lambda x: -np.sin(x),
    np.cosh: np.sinh,
    np.degrees: lambda x: np.degrees(1.),
    np.exp: np.exp,
    np.expm1: np.exp,
    np.fabs: lambda x: np.where(x >= 0., 1., -1.),
    np.hypot: (lambda x, y: x / np.hypot(x, y),
               lambda x, y: y / np.hypot(x, y)),
    np.log: lambda x: 1/x,
    np.log10: lambda x: 1/x/np.log(10),
    np.log1p: lambda x: 1/(1+x),
    np.power: (lambda x, y: x**(y-1) * y,   # explicitly zero y==0?
               lambda x, y: np.log(x)*x**y),
    # above: explicityly zero np.logical_and(x==0, y>0)?
    np.radians: lambda x: np.radians(1),
    np.sin: np.cos,
    np.sinh: np.cosh,
    np.sqrt: lambda x: 0.5/np.sqrt(x),
    np.square: lambda x: 2.0*x,
    np.tan: lambda x: 1+np.tan(x)**2,
    np.tanh: lambda x: 1-np.tanh(x)**2}

UFUNC_DERIVATIVES[np.divide] = UFUNC_DERIVATIVES[np.true_divide]
UFUNC_DERIVATIVES[np.abs] = UFUNC_DERIVATIVES[np.fabs]
UFUNC_DERIVATIVES[np.cbrt] = lambda x: 1./(3.*np.cbrt(x)**2)


class Uncertainty:
    """Measurement Uncertainty.

    Parameters
    ----------
    uncertainty : `~numpy.ndarray` or subclass
        Measurement uncertainties
    check : Bool
        Whether to convert ``uncertainty`` to an array, and check it is
        positive (default: ``True``).

    Notes
    -----
    This class is mostly for internal use by `Variable`.
    """
    def __init__(self, uncertainty, check=True):
        if check and np.any(uncertainty < 0):
                raise ValueError("The uncertainty cannot be negative.")
        self.uncertainty = uncertainty
        # Set up an ID that uses the memory location and the array strides;
        # the latter helps keep track of slices.
        self.id = (uncertainty.__array_interface__['data'][0],
                   uncertainty.__array_interface__['strides'])
        # Any derived uncertainty will use this dictionary, ensuring there is
        # a link that keeps the value in memory until it is not used any more.
        self.derivatives = {self.id: [self, 1.]}

    def __call__(self):
        """Get the uncertainty.  Returns simply the stored value."""
        return self.uncertainty

    def __getitem__(self, item):
        # Particular items are stored in a new instance, with a new ID.
        return self.__class__(uncertainty=self.uncertainty[item], check=False)

    def __repr__(self):
        return '{0}({1})'.format(type(self).__name__, self.uncertainty)


class DerivedUncertainty:
    """Derived Uncertainty.

    Parameters
    ----------
    derivatives : dict
        A dictionary keyed on uncertainty IDs, containing both the uncertainty
        and the derivative with which a variable depends on it.

    Notes
    -----
    This class is mostly for internal use by `Variable`.
    """
    def __init__(self, derivatives):
        self.derivatives = derivatives

    @property
    def uncertainty_components(self):
        """Iterator over the different components of the uncertainty."""
        for unc_id, (uncertainty, derivative) in self.derivatives.items():
            uncertainty = uncertainty.uncertainty
            component = np.abs(derivative * uncertainty)
            # The component is generally an array; if any derivatives are
            # undefined, but the uncertainty is 0, the component is zero too.
            # [Twice faster is to do: np.isnan(np.max(a))]
            if np.any(np.isnan(derivative)):
                component[np.logical_and(np.isnan(derivative),
                                         uncertainty == 0.)] = 0.
            yield component

    def __call__(self):
        """Calculate the uncertainty as the rms of the components."""
        variance = None
        for delta in self.uncertainty_components:
            variance = delta**2 if variance is None else variance + delta**2
        return np.sqrt(variance)

    def __getitem__(self, item):
        # Get the item from each of the uncertainties that contribute.
        # TODO: this needs proper support for slicing, incl. broadcasting, etc.
        derivatives = {}
        for unc_id, (uncertainty, derivative) in self.derivatives.items():
            uncertainty = uncertainty[item]
            if isiterable(derivative):
                derivative = derivative[item]
            derivatives[uncertainty.id] = [uncertainty, derivative]
        return self.__class__(derivatives)

    def __repr__(self):
        return '{0}({1})'.format(type(self).__name__, self.derivatives)


class Variable(np.ndarray):
    """Value with an uncertainty.

    Parameters
    ----------
    value : array
        The nominal value of the variable.
    uncertainty : array
        The associated measurement uncertainty.  This is assumed to have no
        dependences on any other variables.
    kwargs: dict
        Any further arguments passed on the `~numpy.array` initializer.
        By default, ``subok=True``, i.e., subclasses are passed through.
    """

    __array_priority__ = 100000

    def __new__(cls, value, uncertainty, **kwargs):
        kwargs.setdefault('subok', True)
        value = np.array(value, **kwargs)
        if uncertainty is not None:
            uncertainty = Uncertainty(np.array(uncertainty, **kwargs))
        # We keep the nominal value separately.  It points to the same data
        # as `self`, but can have a different class (e.g., Quantity).
        self = value.view(cls)
        self._nominal_value = value
        self._uncertainty = uncertainty
        return self

    @property
    def nominal_value(self):
        return self._nominal_value

    @property
    def uncertainty(self):
        """Uncertainty associated with the variable.

        Returns either the measurement uncertainty, or the uncertainty derived
        by propagating the uncertainties of the underlying measurements.
        """
        return self._uncertainty() if self._uncertainty is not None else 0

    def __getitem__(self, item):
        nominal_value = self._nominal_value[item]
        result = nominal_value.view(type(self))
        result._nominal_value = nominal_value
        result._uncertainty = self._uncertainty[item]
        return result

    def __array_finalize__(self, obj):
        if super().__array_finalize__:
            super().__array_finalize__(obj)
        self._nominal_value = getattr(obj, '_nominal_value', obj)
        self._uncertainty = getattr(obj, '_uncertainty', None)

    def __array_prepare__(self, obj, context=None):
        if context is None:
            return obj
        # The calculation will be on the data proper, so will produce the
        # nominal value.  Hence, here we prepare the object that will become
        # the nominal value using the nominal input values and the
        # __array_prepare__ of the highest ranked one.
        function = context[0]
        nominal_values = tuple(getattr(arg, 'nominal_value', arg)
                               for arg in context[1][:function.nin])
        nv_context = (context[0], nominal_values + context[1][function.nin:],
                      context[2])
        prepare = self._nominal_value.__array_prepare__
        priority = self._nominal_value.__array_priority__
        for nv in nominal_values:
            if getattr(nv, '__array_priority__', priority) > priority:
                prepare = nv.__array_prepare__
                priority = nv.__array_priority__

        return prepare(obj, nv_context)

    def __array_wrap__(self, obj, context=None):
        if context is None:
            return obj
        # First finish the calculation of the nominal value, by wrapping using
        # the nominal input context and the highest priority __array_wrap__.
        function = context[0]
        inputs = context[1][:function.nin]
        nominal_values = tuple(getattr(arg, 'nominal_value', arg)
                               for arg in inputs)
        nv_context = (context[0], nominal_values + context[1][function.nin:],
                      context[2])
        wrap = self._nominal_value.__array_wrap__
        priority = self._nominal_value.__array_priority__
        for nv in nominal_values:
            if getattr(nv, '__array_priority__', priority) > priority:
                wrap = nv.__array_wrap__
                priority = nv.__array_priority__

        nominal_value = wrap(obj, nv_context)

        obj = nominal_value.view(type(self))
        obj._nominal_value = nominal_value
        return obj._add_derivatives(function, inputs, nominal_values)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Evaluate a function on the nominal value, deriving its uncertainty.

        The derivation is done by error propagation, using the derivatives of
        the given function.
        """
        # No support yet for ufunc methods (at, reduce, reduceat, outer, etc.).
        assert method == '__call__'
        # No support yet for output arguments.
        assert 'out' not in kwargs
        # Apply the function to the nominal values.
        values = [getattr(arg, 'nominal_value', arg) for arg in inputs]
        value = ufunc(*values, **kwargs)
        # Set up the output as a Variable that contains a derived uncertainty.
        if not isinstance(value, np.ndarray):
            value = np.array(value)
        result = value.view(type(self))
        return result._add_derivatives(ufunc, inputs, values)

    def _add_derivatives(self, ufunc, inputs, nominal_values):
        # get the functions that calculate derivatives of the ufunc
        # relative to the nominal input values.
        ufunc_derivs = UFUNC_DERIVATIVES[ufunc]
        if len(inputs) == 1:
            ufunc_derivs = (ufunc_derivs,)

        # For Variable inputs, get possible derivatives to other variables.
        deriv_dicts = [getattr(getattr(arg, '_uncertainty', None),
                               'derivatives', {}) for arg in inputs]
        # Set up a new derivatives dictionary, with entries to all variables
        # the result depends on.
        derivatives = {}
        for deriv_dict in deriv_dicts:
            for k, [u, d] in deriv_dict.items():
                derivatives[k] = [u, None]
        # Add to derivatives using chain rule.
        for ufunc_deriv, deriv_dict in zip(ufunc_derivs, deriv_dicts):
            if deriv_dict:
                new_deriv = ufunc_deriv(*nominal_values)
            for unc_id, (uncertainty, derivative) in deriv_dict.items():
                if derivatives[unc_id][1] is None:
                    derivatives[unc_id][1] = new_deriv * derivative
                else:
                    derivatives[unc_id][1] += new_deriv * derivative
        self._uncertainty = DerivedUncertainty(derivatives)
        return self

    def __str__(self):
        return '{0}±{1}'.format(self.nominal_value, self.uncertainty)

    def __repr__(self):
        return '{0}(value={1}, uncertainty={2})'.format(
            type(self).__name__, self.nominal_value, self.uncertainty)
