# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module provides utility functions for the models package.
"""

import warnings

# pylint: disable=invalid-name
from collections import UserDict
from inspect import signature

import numpy as np

from astropy import units as u
from astropy.utils.decorators import deprecated

__all__ = ["poly_map_domain", "comb", "ellipse_extent"]


def make_binary_operator_eval(oper, f, g):
    """
    Given a binary operator (as a callable of two arguments) ``oper`` and
    two callables ``f`` and ``g`` which accept the same arguments,
    returns a *new* function that takes the same arguments as ``f`` and ``g``,
    but passes the outputs of ``f`` and ``g`` in the given ``oper``.

    ``f`` and ``g`` are assumed to return tuples (which may be 1-tuples).  The
    given operator is applied element-wise to tuple outputs).

    Example
    -------
    >>> from operator import add
    >>> def prod(x, y):
    ...     return (x * y,)
    ...
    >>> sum_of_prod = make_binary_operator_eval(add, prod, prod)
    >>> sum_of_prod(3, 5)
    (30,)
    """
    return lambda inputs, params: tuple(
        oper(x, y) for x, y in zip(f(inputs, params), g(inputs, params))
    )


def poly_map_domain(oldx, domain, window):
    """
    Map domain into window by shifting and scaling.

    Parameters
    ----------
    oldx : array
          original coordinates
    domain : list or tuple of length 2
          function domain
    window : list or tuple of length 2
          range into which to map the domain
    """
    domain = np.array(domain, dtype=np.float64)
    window = np.array(window, dtype=np.float64)
    if domain.shape != (2,) or window.shape != (2,):
        raise ValueError('Expected "domain" and "window" to be a tuple of size 2.')
    scl = (window[1] - window[0]) / (domain[1] - domain[0])
    off = (window[0] * domain[1] - window[1] * domain[0]) / (domain[1] - domain[0])
    return off + scl * oldx


def _validate_domain_window(value):
    if value is not None:
        if np.asanyarray(value).shape != (2,):
            raise ValueError("domain and window should be tuples of size 2.")
        return tuple(value)
    return value


@deprecated("5.3", alternative="math.comb")
def comb(N, k):
    """
    The number of combinations of N things taken k at a time.

    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.

    """
    if (k > N) or (N < 0) or (k < 0):
        return 0
    val = 1
    for j in range(min(k, N - k)):
        val = (val * (N - j)) / (j + 1)
    return val


def array_repr_oneline(array):
    """
    Represents a multi-dimensional Numpy array flattened onto a single line.
    """
    r = np.array2string(array, separator=", ", suppress_small=True)
    return " ".join(line.strip() for line in r.splitlines())


def combine_labels(left, right):
    """
    For use with the join operator &: Combine left input/output labels with
    right input/output labels.

    If none of the labels conflict then this just returns a sum of tuples.
    However if *any* of the labels conflict, this appends '0' to the left-hand
    labels and '1' to the right-hand labels so there is no ambiguity).
    """
    if set(left).intersection(right):
        left = tuple(label + "0" for label in left)
        right = tuple(label + "1" for label in right)

    return left + right


def ellipse_extent(a, b, theta):
    """
    Calculates the half size of a box encapsulating a rotated 2D
    ellipse.

    Parameters
    ----------
    a : float or `~astropy.units.Quantity`
        The ellipse semimajor axis.
    b : float or `~astropy.units.Quantity`
        The ellipse semiminor axis.
    theta : float or `~astropy.units.Quantity` ['angle']
        The rotation angle as an angular quantity
        (`~astropy.units.Quantity` or `~astropy.coordinates.Angle`) or
        a value in radians (as a float). The rotation angle increases
        counterclockwise.

    Returns
    -------
    offsets : tuple
        The absolute value of the offset distances from the ellipse center that
        define its bounding box region, ``(dx, dy)``.

    Examples
    --------
    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.modeling.models import Ellipse2D
        from astropy.modeling.utils import ellipse_extent, render_model

        amplitude = 1
        x0 = 50
        y0 = 50
        a = 30
        b = 10
        theta = np.pi / 4

        model = Ellipse2D(amplitude, x0, y0, a, b, theta)
        dx, dy = ellipse_extent(a, b, theta)
        limits = [x0 - dx, x0 + dx, y0 - dy, y0 + dy]
        model.bounding_box = limits

        image = render_model(model)

        plt.imshow(image, cmap='binary', interpolation='nearest', alpha=.5,
                  extent = limits)
        plt.show()
    """
    from .parameters import Parameter  # prevent circular import

    if isinstance(theta, Parameter):
        if theta.quantity is None:
            theta = theta.value
        else:
            theta = theta.quantity

    t = np.arctan2(-b * np.tan(theta), a)
    dx = a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)

    t = np.arctan2(b, a * np.tan(theta))
    dy = b * np.sin(t) * np.cos(theta) + a * np.cos(t) * np.sin(theta)

    if isinstance(dx, u.Quantity) or isinstance(dy, u.Quantity):
        return np.abs(u.Quantity([dx, dy], subok=True))
    return np.abs([dx, dy])


def get_inputs_and_params(func):
    """
    Given a callable, determine the input variables and the
    parameters.

    Parameters
    ----------
    func : callable

    Returns
    -------
    inputs, params : tuple
        Each entry is a list of inspect.Parameter objects
    """
    sig = signature(func)

    inputs = []
    params = []
    for param in sig.parameters.values():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            raise ValueError("Signature must not have *args or **kwargs")
        if param.default == param.empty:
            inputs.append(param)
        else:
            params.append(param)

    return inputs, params


def _combine_equivalency_dict(keys, eq1=None, eq2=None):
    # Given two dictionaries that give equivalencies for a set of keys, for
    # example input value names, return a dictionary that includes all the
    # equivalencies
    eq = {}
    for key in keys:
        eq[key] = []
        if eq1 is not None and key in eq1:
            eq[key].extend(eq1[key])
        if eq2 is not None and key in eq2:
            eq[key].extend(eq2[key])
    return eq


def _to_radian(value):
    """Convert ``value`` to radian."""
    if isinstance(value, u.Quantity):
        return value.to(u.rad)
    return np.deg2rad(value)


def _to_orig_unit(value, raw_unit=None, orig_unit=None):
    """Convert value with ``raw_unit`` to ``orig_unit``."""
    if raw_unit is not None:
        return (value * raw_unit).to(orig_unit)
    return np.rad2deg(value)


class _ConstraintsDict(UserDict):
    """
    Wrapper around UserDict to allow updating the constraints
    on a Parameter when the dictionary is updated.
    """

    def __init__(self, model, constraint_type):
        self._model = model
        self.constraint_type = constraint_type
        c = {}
        for name in model.param_names:
            param = getattr(model, name)
            c[name] = getattr(param, constraint_type)
        super().__init__(c)

    def __setitem__(self, key, val):
        super().__setitem__(key, val)
        param = getattr(self._model, key)
        setattr(param, self.constraint_type, val)


class _SpecialOperatorsDict(UserDict):
    """
    Wrapper around UserDict to allow for better tracking of the Special
    Operators for CompoundModels. This dictionary is structured so that
    one cannot inadvertently overwrite an existing special operator.

    Parameters
    ----------
    unique_id: int
        the last used unique_id for a SPECIAL OPERATOR
    special_operators: dict
        a dictionary containing the special_operators

    Notes
    -----
    Direct setting of operators (`dict[key] = value`) into the
    dictionary has been deprecated in favor of the `.add(name, value)`
    method, so that unique dictionary keys can be generated and tracked
    consistently.
    """

    def __init__(self, unique_id=0, special_operators={}):
        super().__init__(special_operators)
        self._unique_id = unique_id

    def _set_value(self, key, val):
        if key in self:
            raise ValueError(f'Special operator "{key}" already exists')
        else:
            super().__setitem__(key, val)

    def __setitem__(self, key, val):
        self._set_value(key, val)
        warnings.warn(
            DeprecationWarning(
                """
            Special operator dictionary assignment has been deprecated.
            Please use `.add` instead, so that you can capture a unique
            key for your operator.
            """
            )
        )

    def _get_unique_id(self):
        self._unique_id += 1

        return self._unique_id

    def add(self, operator_name, operator):
        """
        Adds a special operator to the dictionary, and then returns the
        unique key that the operator is stored under for later reference.

        Parameters
        ----------
        operator_name: str
            the name for the operator
        operator: function
            the actual operator function which will be used

        Returns
        -------
        the unique operator key for the dictionary
            `(operator_name, unique_id)`
        """
        key = (operator_name, self._get_unique_id())

        self._set_value(key, operator)

        return key
