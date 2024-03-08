# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Comparison functions for `astropy.cosmology.Cosmology`.

This module is **NOT** public API. To use these functions, import them from
the top-level namespace -- :mod:`astropy.cosmology`. This module will be
moved.
"""

from __future__ import annotations

import functools
import inspect
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy import False_, True_, ndarray

from astropy import table
from astropy.cosmology.core import Cosmology

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, TypeAlias

__all__: list[str] = []  # Nothing is scoped here


##############################################################################
# PARAMETERS

_FormatType: TypeAlias = bool | None | str
_FormatsType: TypeAlias = _FormatType | tuple[_FormatType, ...]

_COSMO_AOK: set[Any] = {None, True_, False_, "astropy.cosmology"}
# The numpy bool also catches real bool for ops "==" and "in"


##############################################################################
# UTILITIES


_CANT_BROADCAST: tuple[type, ...] = (table.Row, table.Table)
"""Things that cannot broadcast.

Have to deal with things that do not broadcast well. e.g.
`~astropy.table.Row` cannot be used in an array, even if ``dtype=object``
and will raise a segfault when used in a `numpy.ufunc`.
"""


@dataclass(frozen=True)
class _CosmologyWrapper:
    """A private wrapper class to hide things from :mod:`numpy`.

    This should never be exposed to the user.
    """

    __slots__ = ("wrapped",)

    wrapped: Any


@functools.partial(np.frompyfunc, nin=2, nout=1)
def _parse_format(cosmo: Any, format: _FormatType, /) -> Cosmology:
    """Parse Cosmology-like input into Cosmologies, given a format hint.

    Parameters
    ----------
    cosmo : |Cosmology|-like, positional-only
        |Cosmology| to parse.
    format : bool or None or str, positional-only
        Whether to allow, before equivalence is checked, the object to be
        converted to a |Cosmology|. This allows, e.g. a |Table| to be equivalent
        to a |Cosmology|. `False` (default) will not allow conversion. `True` or
        `None` will, and will use the auto-identification to try to infer the
        correct format. A `str` is assumed to be the correct format to use when
        converting.

    Returns
    -------
    |Cosmology| or generator thereof

    Raises
    ------
    TypeError
        If ``cosmo`` is not a |Cosmology| and ``format`` equals `False`.
    TypeError
        If ``cosmo`` is a |Cosmology| and ``format`` is not `None` or equal to
        `True`.
    """
    # Deal with private wrapper
    if isinstance(cosmo, _CosmologyWrapper):
        cosmo = cosmo.wrapped

    # Shortcut if already a cosmology
    if isinstance(cosmo, Cosmology):
        if format not in _COSMO_AOK:
            allowed = "/".join(map(str, _COSMO_AOK))
            raise ValueError(
                f"for parsing a Cosmology, 'format' must be {allowed}, not {format}"
            )
        return cosmo
    # Convert, if allowed.
    elif format == False_:  # catches False and False_
        raise TypeError(
            f"if 'format' is False, arguments must be a Cosmology, not {cosmo}"
        )
    else:
        format = None if format == True_ else format  # str->str, None/True/True_->None
        out = Cosmology.from_format(cosmo, format=format)  # this can error!

    return out


def _parse_formats(*cosmos: object, format: _FormatsType) -> ndarray:
    """Parse Cosmology-like to |Cosmology|, using provided formats.

    ``format`` is broadcast to match the shape of the cosmology arguments. Note
    that the cosmology arguments are not broadcast against ``format``, so it
    cannot determine the output shape.

    Parameters
    ----------
    *cosmos : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as specified
        by the corresponding ``format``.

    format : bool or None or str or array-like thereof, positional-only
        Whether to allow, before equivalence is checked, the object to be
        converted to a |Cosmology|. This allows, e.g. a |Table| to be equivalent
        to a |Cosmology|. `False` (default) will not allow conversion. `True` or
        `None` will, and will use the auto-identification to try to infer the
        correct format. A `str` is assumed to be the correct format to use when
        converting. Note ``format`` is broadcast as an object array to match the
        shape of ``cosmos`` so ``format`` cannot determine the output shape.

    Raises
    ------
    TypeError
        If any in ``cosmos`` is not a |Cosmology| and the corresponding
        ``format`` equals `False`.
    """
    formats = np.broadcast_to(np.array(format, dtype=object), len(cosmos))
    # parse each cosmo & format

    # Have to deal with things that do not broadcast well.
    # astropy.row cannot be used in an array, even if dtype=object
    # and will raise a segfault when used in a ufunc.
    wcosmos = [
        c if not isinstance(c, _CANT_BROADCAST) else _CosmologyWrapper(c)
        for c in cosmos
    ]

    return _parse_format(wcosmos, formats)


def _comparison_decorator(pyfunc: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to make wrapper function that parses |Cosmology|-like inputs.

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

    format : bool or None or str or array-like thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|. This
        allows, e.g. a |Table| to be given instead a |Cosmology|. `False`
        (default) will not allow conversion. `True` or `None` will, and will use
        the auto-identification to try to infer the correct format. A `str` is
        assumed to be the correct format to use when converting. Note ``format``
        is broadcast as an object array to match the shape of ``cosmos`` so
        ``format`` cannot determine the output shape.
    """
    sig = inspect.signature(pyfunc)
    nin = sum(p.kind == 0 for p in sig.parameters.values())

    # Make wrapper function that parses cosmology-like inputs
    @functools.wraps(pyfunc)
    def wrapper(*cosmos: Any, format: _FormatsType = False, **kwargs: Any) -> bool:
        if len(cosmos) > nin:
            raise TypeError(
                f"{wrapper.__wrapped__.__name__} takes {nin} positional"
                f" arguments but {len(cosmos)} were given"
            )
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
def cosmology_equal(
    cosmo1: Any, cosmo2: Any, /, *, allow_equivalent: bool = False
) -> bool:
    r"""Return element-wise equality check on the cosmologies.

    .. note::

        Cosmologies are currently scalar in their parameters.

    Parameters
    ----------
    cosmo1, cosmo2 : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as specified
        by ``format``.

    format : bool or None or str or tuple thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|. This
        allows, e.g. a |Table| to be given instead a |Cosmology|. `False`
        (default) will not allow conversion. `True` or `None` will, and will use
        the auto-identification to try to infer the correct format. A `str` is
        assumed to be the correct format to use when converting. Note ``format``
        is broadcast as an object array to match the shape of ``cosmos`` so
        ``format`` cannot determine the output shape.

    allow_equivalent : bool, optional keyword-only
        Whether to allow cosmologies to be equal even if not of the same class.
        For example, an instance of |LambdaCDM| might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like |FlatLambdaCDM|.

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

    Two cosmologies may be equivalent even if not of the same class. In these
    examples the |LambdaCDM| has :attr:`~astropy.cosmology.LambdaCDM.Ode0` set
    to the same value calculated in |FlatLambdaCDM|.

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

    Also, using the keyword argument, the notion of equality is extended to any
    Python object that can be converted to a |Cosmology|.

        >>> mapping = cosmo2.to_format("mapping")
        >>> cosmology_equal(cosmo1, mapping, format=True)
        True

    Either (or both) arguments can be |Cosmology|-like.

        >>> cosmology_equal(mapping, cosmo2, format=True)
        True

    The list of valid formats, e.g. the |Table| in this example, may be checked
    with ``Cosmology.from_format.list_formats()``.

    As can be seen in the list of formats, not all formats can be
    auto-identified by ``Cosmology.from_format.registry``. Objects of these
    kinds can still be checked for equality, but the correct format string must
    be used.

        >>> yml = cosmo2.to_format("yaml")
        >>> cosmology_equal(cosmo1, yml, format=(None, "yaml"))
        True

    This also works with an array of ``format`` matching the number of
    cosmologies.

        >>> cosmology_equal(mapping, yml, format=[True, "yaml"])
        True
    """
    # Check parameter equality
    if not allow_equivalent:
        eq = cosmo1 == cosmo2

    else:
        # Check parameter equivalence
        # The options are: 1) same class & parameters; 2) same class, different
        # parameters; 3) different classes, equivalent parameters; 4) different
        # classes, different parameters. (1) & (3) => True, (2) & (4) => False.
        eq = cosmo1.__equiv__(cosmo2)
        if eq is NotImplemented:
            eq = cosmo2.__equiv__(cosmo1)  # that failed, try from 'other'

        eq = False if eq is NotImplemented else eq

    # TODO! include equality check of metadata

    return eq


@_comparison_decorator
def _cosmology_not_equal(
    cosmo1: Any, cosmo2: Any, /, *, allow_equivalent: bool = False
) -> bool:
    r"""Return element-wise cosmology non-equality check.

    .. note::

        Cosmologies are currently scalar in their parameters.

    Parameters
    ----------
    cosmo1, cosmo2 : |Cosmology|-like
        The objects to compare. Must be convertible to |Cosmology|, as specified
        by ``format``.

    out : ndarray, None, optional
        A location into which the result is stored. If provided, it must have a
        shape that the inputs broadcast to. If not provided or None, a
        freshly-allocated array is returned.

    format : bool or None or str or tuple thereof, optional keyword-only
        Whether to allow the arguments to be converted to a |Cosmology|. This
        allows, e.g. a |Table| to be given instead a Cosmology. `False`
        (default) will not allow conversion. `True` or `None` will, and will use
        the auto-identification to try to infer the correct format. A `str` is
        assumed to be the correct format to use when converting. ``format`` is
        broadcast to match the shape of the cosmology arguments. Note that the
        cosmology arguments are not broadcast against ``format``, so it cannot
        determine the output shape.

    allow_equivalent : bool, optional keyword-only
        Whether to allow cosmologies to be equal even if not of the same class.
        For example, an instance of |LambdaCDM| might have :math:`\Omega_0=1`
        and :math:`\Omega_k=0` and therefore be flat, like |FlatLambdaCDM|.

    See Also
    --------
    astropy.cosmology.cosmology_equal
        Element-wise equality check, with argument conversion to Cosmology.
    """
    neq = not cosmology_equal(cosmo1, cosmo2, allow_equivalent=allow_equivalent)
    # TODO! it might eventually be worth the speed boost to implement some of
    #       the internals of cosmology_equal here, but for now it's a hassle.

    return neq
