# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Comparison functions for `astropy.cosmology.Cosmology`.

This module is **NOT** public API. To use these functions, import them from
the top-level namespace -- :mod:`astropy.cosmology`.
This module will be moved.
"""

from __future__ import annotations

import functools
import inspect
from typing import Any, Callable, Tuple, Union

import numpy as np
from numpy import False_, True_, ndarray

from astropy.cosmology.core import Cosmology
from astropy import table

__all__ = []  # Nothing is scoped here


##############################################################################
# PARAMETERS

_FormatType = Union[bool, None, str]
_FormatsType = Union[_FormatType, Tuple[_FormatType, ...]]


##############################################################################
# UTILITIES


class _CosmologyWrapper:
    """
    A private wrapper class to hide things from :mod:`numpy`.
    This should never be exposed to the user.
    """

    _cantbroadcast: Tuple[type, ...] = (table.Row, table.Table)
    """
    Have to deal with things that do not broadcast well.
    e.g. astropy.row cannot be used in an array, even if dtype=object
    and will raise a segfault when used in a ufunc.
    """

    value: Any

    def __init__(self, value: Any) -> None:
        self.value = value


# https://github.com/numpy/numpy/issues/9477
# TODO segfault on astropy.row and np.vectorize can't coerce table to dtypes
def _wrap_to_ufunc(pyfunc: Callable[[Any, _FormatType], Cosmology]) -> np.ufunc:
    ufunc = np.frompyfunc(pyfunc, 2, 1)
    return ufunc


@_wrap_to_ufunc  # TODO! use when py3.9+: @functools.partial(np.frompyfunc, nin=2, nout=1)
def _parse_format(cosmo: Any, format: _FormatType, /,) -> Cosmology:
    """Parse Cosmology-like input into Cosmologies, given a format hint.

    Parameters
    ----------
    cosmo : |Cosmology|-like, positional-only
        Cosmology to parse.
    format : bool or None or str or tuple thereof, positional-only
        Whether to allow, before equivalence is checked, the object to be
        converted to a |Cosmology|. This allows, e.g. a |Table| to be
        equivalent to a Cosmology.
        `False` (default) will not allow conversion. `True` or `None` will,
        and will use the auto-identification to try to infer the correct
        format. A `str` is assumed to be the correct format to use when
        converting.

    Returns
    -------
    |Cosmology| or generator thereof

    Raises
    ------
    TypeError
        If 'cosmo' is not a |Cosmology| and 'format' equals `numpy.False_`.
        If 'cosmo' is a |Cosmology| and 'format' is not `None` or equal to `numpy.True_`.
    """
    # Deal with private wrapper
    if isinstance(cosmo, _CosmologyWrapper):
        cosmo = cosmo.value

    # Shortcut if already a cosmology
    if isinstance(cosmo, Cosmology):
        if format not in (None, True_, False_):  # also catches real bool
            raise ValueError(f"for a Cosmology, 'format' must be None or True, not {format}")
        return cosmo

    if format != False_:  # catches False and False_
        format = None if format == True_ else format  # str->str, None/True->None
        out = Cosmology.from_format(cosmo, format=format)  # this can error!
    elif not isinstance(cosmo, Cosmology):
        raise TypeError(f"if 'format' is False, arguments must be a Cosmology, not {cosmo}")
    else:
        out = cosmo

    return out


def _parse_formats(*cosmos: object, format: _FormatsType) -> ndarray:
    """Parse Cosmology-like to |Cosmology|, using provided formats.

    `format` is broadcast to match the shape of the cosmology arguments.
    Note that the cosmology arguments are not broadcast against ``format``,
    so it cannot determine the output shape.

    Parameters
    ----------
    *cosmos : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as
        specified by the corresponding `format`.

    Raises
    ------
    TypeError
        If any in 'cosmos' is not a |Cosmology| and the corresponding 'format'
        equals `numpy.False_`.
    """
    formats = np.broadcast_to(format, len(cosmos))
    # parse each cosmo & format

    # Have to deal with things that do not broadcast well.
    # astropy.row cannot be used in an array, even if dtype=object
    # and will raise a segfault when used in a ufunc.
    towrap = (isinstance(cosmo, _CosmologyWrapper._cantbroadcast) for cosmo in cosmos)
    wcosmos = [(c if not wrap else _CosmologyWrapper(c)) for c, wrap in zip(cosmos, towrap)]

    return _parse_format(wcosmos, formats)


def _comparison_decorator(pyfunc: Callable[..., Any]) -> Callable[..., Any]:
    """
    Decorator to make wrapper function that parses cosmology-like inputs.

    Parameters
    ----------
    pyfunc : Python function object
        An arbitrary Python function.

    Returns
    -------
    callable[..., Any]
        Wrapped `pyfunc`, as described above.

    Notes
    -----
    All decorated functions should add the following to 'Parameters'.

    format : bool or None or str or tuple thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|.
        This allows, e.g. a |Table| to be given instead a Cosmology.
        `False` (default) will not allow conversion. `True` or `None` will,
        and will use the auto-identification to try to infer the correct
        format. A `str` is assumed to be the correct format to use when
        converting.
        `format` is broadcast to match the shape of the cosmology arguments.
        Note that the cosmology arguments are not broadcast against ``format``,
        so it cannot determine the output shape.
    """
    sig = inspect.signature(pyfunc)
    nin = len([p.kind == 0 for p in sig.parameters.values()])

    # Make wrapper function that parses cosmology-like inputs
    @functools.wraps(pyfunc)
    def wrapper(*cosmos: Any, format: _FormatsType=False, **kwargs: Any) -> bool:
        if len(cosmos) > nin:
            raise TypeError
        # Parse cosmologies to format. Only do specified number.
        cosmos = _parse_formats(*cosmos, format=format)
        # Evaluate pyfunc, erroring if didn't match specified number.
        result = wrapper.__wrapped__(*cosmos, **kwargs)
        # Return, casting to correct type casting is possible.
        return result

    return wrapper


##############################################################################
# COMPARISON FUNCTIONS


@_comparison_decorator
def cosmology_equal(cosmo1: Any, cosmo2: Any, /, *, allow_equivalent: bool=False) -> bool:
    r"""Return element-wise equality check on the cosmologies.

    .. note::

        Cosmologies are currently scalar in their parameters.

    Parameters
    ----------
    cosmo1, cosmo2 : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as
        specified by `format`.

    format : bool or None or str or tuple thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|.
        This allows, e.g. a |Table| to be given instead a Cosmology.
        `False` (default) will not allow conversion. `True` or `None` will,
        and will use the auto-identification to try to infer the correct
        format. A `str` is assumed to be the correct format to use when
        converting.
        `format` is broadcast to match the shape of the cosmology arguments.
        Note that the cosmology arguments are not broadcast against ``format``,
        so it cannot determine the output shape.

    allow_equivalent : bool, optional keyword-only
        Whether to allow cosmologies to be equal even if not of the same class.
        For example, an instance of ``LambdaCDM`` might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like ``FlatLambdaCDM``.

    See Also
    --------
    astropy.cosmology.funcs.cosmology_not_equal
        Element-wise non-equality check, with argument conversion to Cosmology.

    Examples
    --------
    Assuming the following imports

        >>> import astropy.units as u
        >>> from astropy.cosmology import FlatLambdaCDM

    Two identical cosmologies are equal.

        >>> cosmo1 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
        >>> cosmo2 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
        >>> cosmology_equal(cosmo1, cosmo2)
        True

    And cosmologies with different parameters are not.

        >>> cosmo3 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.4)
        >>> cosmology_equal(cosmo1, cosmo3)
        False

    Two cosmologies may be equivalent even if not of the same class.
    In this examples the ``LambdaCDM`` has ``Ode0`` set to the same value
    calculated in ``FlatLambdaCDM``.

        >>> from astropy.cosmology import LambdaCDM
        >>> cosmo3 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
        >>> cosmology_equal(cosmo1, cosmo3)
        False
        >>> cosmology_equal(cosmo1, cosmo3, allow_equivalent=True)
        True

    While in this example, the cosmologies are not equivalent.

        >>> cosmo4 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, Tcmb0=3 * u.K)
        >>> cosmology_equal(cosmo3, cosmo4, allow_equivalent=True)
        False

    Also, using the keyword argument, the notion of equality is extended
    to any Python object that can be converted to a |Cosmology|.

        >>> mapping = cosmo2.to_format("mapping")
        >>> cosmology_equal(cosmo1, mapping, format=True)
        True

    Either (or both) arguments can be |Cosmology|-like.

        >>> cosmology_equal(mapping, cosmo2, format=True)
        True

    The list of valid formats, e.g. the |Table| in this example, may be
    checked with ``Cosmology.from_format.list_formats()``.

    As can be seen in the list of formats, not all formats can be
    auto-identified by ``Cosmology.from_format.registry``. Objects of
    these kinds can still be checked for equality, but the correct
    format string must be used.

        >>> yml = cosmo2.to_format("yaml")
        >>> cosmology_equal(cosmo1, yml, format="yaml")
        True

    This also works with an array of 'format' matching the number of cosmologies.

        >>> cosmology_equal(mapping, yml, format=[True, "yaml"])
        True
    """
    # Check parameter equality
    if not allow_equivalent:
        eq = (cosmo1 == cosmo2)

    else:
        # Check parameter equivalence
        # The options are: 1) same class & parameters; 2) same class, different
        # parameters; 3) different classes, equivalent parameters; 4) different
        # classes, different parameters. (1) & (3) => True, (2) & (4) => False.
        eq = cosmo1.__equiv__(cosmo2)
        if eq is NotImplemented:
            eq = cosmo2.__equiv__(cosmo1)  # that failed, try from 'other'

        eq = eq if eq is not NotImplemented else False

    # TODO! include equality check of metadata

    return eq


@_comparison_decorator
def cosmology_not_equal(cosmo1: Any, cosmo2: Any, /, *, allow_equivalent: bool=False) -> bool:
    r"""Return element-wise cosmology non-equality check.

    .. note::

        Cosmologies are currently scalar in their parameters.

    Parameters
    ----------
    cosmo1, cosmo2 : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as
        specified by `format`.

    out : ndarray, None, optional
        A location into which the result is stored. If provided, it must have
        a shape that the inputs broadcast to. If not provided or None,
        a freshly-allocated array is returned.

    format : bool or None or str or tuple thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|.
        This allows, e.g. a |Table| to be given instead a Cosmology.
        `False` (default) will not allow conversion. `True` or `None` will,
        and will use the auto-identification to try to infer the correct
        format. A `str` is assumed to be the correct format to use when
        converting.
        `format` is broadcast to match the shape of the cosmology arguments.
        Note that the cosmology arguments are not broadcast against ``format``,
        so it cannot determine the output shape.

    allow_equivalent : bool, optional keyword-only
        Whether to allow cosmologies to be equal even if not of the same class.
        For example, an instance of ``LambdaCDM`` might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like ``FlatLambdaCDM``.

    See Also
    --------
    astropy.cosmology.funcs.cosmology_equal
        Element-wise equality check, with argument conversion to Cosmology.

    Examples
    --------
    Assuming the following imports

        >>> import astropy.units as u
        >>> from astropy.cosmology import FlatLambdaCDM

    Two identical cosmologies are equal.

        >>> cosmo1 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
        >>> cosmo2 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3)
        >>> cosmology_not_equal(cosmo1, cosmo2)
        False

    And cosmologies with different parameters are not.

        >>> cosmo3 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.4)
        >>> cosmology_not_equal(cosmo1, cosmo3)
        True

    Two cosmologies may be equivalent even if not of the same class.
    In this examples the ``LambdaCDM`` has ``Ode0`` set to the same value
    calculated in ``FlatLambdaCDM``.

        >>> from astropy.cosmology import LambdaCDM
        >>> cosmo4 = LambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, 0.7)
        >>> cosmology_not_equal(cosmo1, cosmo4)
        True
        >>> cosmology_not_equal(cosmo1, cosmo4, allow_equivalent=True)
        False

    While in this example, the cosmologies are not equivalent.

        >>> cosmo5 = FlatLambdaCDM(70 * (u.km/u.s/u.Mpc), 0.3, Tcmb0=3 * u.K)
        >>> cosmology_not_equal(cosmo4, cosmo5, allow_equivalent=True)
        True

    Also, using the keyword argument, the notion of equality is extended
    to any Python object that can be converted to a |Cosmology|.

        >>> mapping = cosmo2.to_format("mapping")
        >>> cosmology_not_equal(cosmo1, mapping, format=True)
        False

    Either (or both) arguments can be |Cosmology|-like.

        >>> cosmology_not_equal(mapping, cosmo2, format=True)
        False

    The list of valid formats, e.g. the |Table| in this example, may be
    checked with ``Cosmology.from_format.list_formats()``.

    As can be seen in the list of formats, not all formats can be
    auto-identified by ``Cosmology.from_format.registry``. Objects of
    these kinds can still be checked for equality, but the correct
    format string must be used.

        >>> yml = cosmo3.to_format("yaml")
        >>> cosmology_not_equal(cosmo1, yml, format="yaml")
        True

    This also works with an array of 'format' matching the number of cosmologies.

        >>> cosmology_not_equal(mapping, yml, format=[True, "yaml"])
        True
    """
    neq = not cosmology_equal(cosmo1, cosmo2, allow_equivalent=allow_equivalent)
    # TODO! it might eventually be worth the speed boost to implement some of
    #       the internals of cosmology_equal here, but for now it's a hassle.

    return neq
