# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Miscellaneous utilities for `astropy.units`.

None of the functions in the module are meant for use outside of the
package.
"""

import inspect
import io
import re
from collections.abc import Sequence
from fractions import Fraction

import numpy as np
from numpy import finfo


_float_finfo = finfo(float)
# take float here to ensure comparison with another float is fast
# give a little margin since often multiple calculations happened
_JUST_BELOW_UNITY = float(1.-4.*_float_finfo.epsneg)
_JUST_ABOVE_UNITY = float(1.+4.*_float_finfo.eps)


def _get_first_sentence(s):
    """
    Get the first sentence from a string and remove any carriage
    returns.
    """

    x = re.match(r".*?\S\.\s", s)
    if x is not None:
        s = x.group(0)
    return s.replace('\n', ' ')


def _iter_unit_summary(namespace):
    """
    Generates the ``(unit, doc, represents, aliases, prefixes)``
    tuple used to format the unit summary docs in `generate_unit_summary`.
    """

    from . import core

    # Get all of the units, and keep track of which ones have SI
    # prefixes
    units = []
    has_prefixes = set()
    for key, val in namespace.items():
        # Skip non-unit items
        if not isinstance(val, core.UnitBase):
            continue

        # Skip aliases
        if key != val.name:
            continue

        if isinstance(val, core.PrefixUnit):
            # This will return the root unit that is scaled by the prefix
            # attached to it
            has_prefixes.add(val._represents.bases[0].name)
        else:
            units.append(val)

    # Sort alphabetically, case insensitive
    units.sort(key=lambda x: x.name.lower())

    for unit in units:
        doc = _get_first_sentence(unit.__doc__).strip()
        represents = ''
        if isinstance(unit, core.Unit):
            represents = f":math:`{unit._represents.to_string('latex')[1:-1]}`"
        aliases = ', '.join(f'``{x}``' for x in unit.aliases)

        yield (unit, doc, represents, aliases, 'Yes' if unit.name in has_prefixes else 'No')


def generate_unit_summary(namespace):
    """
    Generates a summary of units from a given namespace.  This is used
    to generate the docstring for the modules that define the actual
    units.

    Parameters
    ----------
    namespace : dict
        A namespace containing units.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """

    docstring = io.StringIO()

    docstring.write("""
.. list-table:: Available Units
   :header-rows: 1
   :widths: 10 20 20 20 1

   * - Unit
     - Description
     - Represents
     - Aliases
     - SI Prefixes
""")

    for unit_summary in _iter_unit_summary(namespace):
        docstring.write("""
   * - ``{}``
     - {}
     - {}
     - {}
     - {}
""".format(*unit_summary))

    return docstring.getvalue()


def generate_prefixonly_unit_summary(namespace):
    """
    Generates table entries for units in a namespace that are just prefixes
    without the base unit.  Note that this is intended to be used *after*
    `generate_unit_summary` and therefore does not include the table header.

    Parameters
    ----------
    namespace : dict
        A namespace containing units that are prefixes but do *not* have the
        base unit in their namespace.

    Returns
    -------
    docstring : str
        A docstring containing a summary table of the units.
    """
    from . import PrefixUnit

    faux_namespace = {}
    for nm, unit in namespace.items():
        if isinstance(unit, PrefixUnit):
            base_unit = unit.represents.bases[0]
            faux_namespace[base_unit.name] = base_unit

    docstring = io.StringIO()

    for unit_summary in _iter_unit_summary(faux_namespace):
        docstring.write("""
   * - Prefixes for ``{}``
     - {} prefixes
     - {}
     - {}
     - Only
""".format(*unit_summary))

    return docstring.getvalue()


def is_effectively_unity(value):
    # value is *almost* always real, except, e.g., for u.mag**0.5, when
    # it will be complex.  Use try/except to ensure normal case is fast
    try:
        return _JUST_BELOW_UNITY <= value <= _JUST_ABOVE_UNITY
    except TypeError:  # value is complex
        return (_JUST_BELOW_UNITY <= value.real <= _JUST_ABOVE_UNITY and
                _JUST_BELOW_UNITY <= value.imag + 1 <= _JUST_ABOVE_UNITY)


def sanitize_scale(scale):
    if is_effectively_unity(scale):
        return 1.0

    # Maximum speed for regular case where scale is a float.
    if scale.__class__ is float:
        return scale

    # We cannot have numpy scalars, since they don't autoconvert to
    # complex if necessary.  They are also slower.
    if hasattr(scale, 'dtype'):
        scale = scale.item()

    # All classes that scale can be (int, float, complex, Fraction)
    # have an "imag" attribute.
    if scale.imag:
        if abs(scale.real) > abs(scale.imag):
            if is_effectively_unity(scale.imag/scale.real + 1):
                return scale.real

        elif is_effectively_unity(scale.real/scale.imag + 1):
            return complex(0., scale.imag)

        return scale

    else:
        return scale.real


def maybe_simple_fraction(p, max_denominator=100):
    """Fraction very close to x with denominator at most max_denominator.

    The fraction has to be such that fraction/x is unity to within 4 ulp.
    If such a fraction does not exist, returns the float number.

    The algorithm is that of `fractions.Fraction.limit_denominator`, but
    sped up by not creating a fraction to start with.
    """
    if p == 0 or p.__class__ is int:
        return p
    n, d = p.as_integer_ratio()
    a = n // d
    # Normally, start with 0,1 and 1,0; here we have applied first iteration.
    n0, d0 = 1, 0
    n1, d1 = a, 1
    while d1 <= max_denominator:
        if _JUST_BELOW_UNITY <= n1/(d1*p) <= _JUST_ABOVE_UNITY:
            return Fraction(n1, d1)
        n, d = d, n-a*d
        a = n // d
        n0, n1 = n1, n0+a*n1
        d0, d1 = d1, d0+a*d1

    return p


def validate_power(p):
    """Convert a power to a floating point value, an integer, or a Fraction.

    If a fractional power can be represented exactly as a floating point
    number, convert it to a float, to make the math much faster; otherwise,
    retain it as a `fractions.Fraction` object to avoid losing precision.
    Conversely, if the value is indistinguishable from a rational number with a
    low-numbered denominator, convert to a Fraction object.

    Parameters
    ----------
    p : float, int, Rational, Fraction
        Power to be converted
    """
    denom = getattr(p, 'denominator', None)
    if denom is None:
        try:
            p = float(p)
        except Exception:
            if not np.isscalar(p):
                raise ValueError("Quantities and Units may only be raised "
                                 "to a scalar power")
            else:
                raise

        # This returns either a (simple) Fraction or the same float.
        p = maybe_simple_fraction(p)
        # If still a float, nothing more to be done.
        if isinstance(p, float):
            return p

        # Otherwise, check for simplifications.
        denom = p.denominator

    if denom == 1:
        p = p.numerator

    elif (denom & (denom - 1)) == 0:
        # Above is a bit-twiddling hack to see if denom is a power of two.
        # If so, float does not lose precision and will speed things up.
        p = float(p)

    return p


def resolve_fractions(a, b):
    """
    If either input is a Fraction, convert the other to a Fraction
    (at least if it does not have a ridiculous denominator).
    This ensures that any operation involving a Fraction will use
    rational arithmetic and preserve precision.
    """
    # We short-circuit on the most common cases of int and float, since
    # isinstance(a, Fraction) is very slow for any non-Fraction instances.
    a_is_fraction = (a.__class__ is not int and a.__class__ is not float and
                     isinstance(a, Fraction))
    b_is_fraction = (b.__class__ is not int and b.__class__ is not float and
                     isinstance(b, Fraction))
    if a_is_fraction and not b_is_fraction:
        b = maybe_simple_fraction(b)
    elif not a_is_fraction and b_is_fraction:
        a = maybe_simple_fraction(a)
    return a, b


def quantity_asanyarray(a, dtype=None):
    from .quantity import Quantity
    if not isinstance(a, np.ndarray) and not np.isscalar(a) and any(isinstance(x, Quantity) for x in a):
        return Quantity(a, dtype=dtype)
    else:
        return np.asanyarray(a, dtype=dtype)


# ------------------------------------------------------------------------------


def quantity_frompyfunc(func, nin, nout, inunits=None, ounits=None,
                        *, identity=None, assume_correct_units=False):
    """Quantity-aware `~numpy.frompyfunc`.

    `~numpy.ufunc`s operate on only recognized `~numpy.dtype`s (e.g. float32),
    so units must be removed beforehand and replaced afterwards. Therefore units
    MUST BE KNOWN a priori, as they will not be propagated by the astropy
    machinery.

    Parameters
    ----------
    func : callable
    nin, nout : int
        Number of ufunc's inputs and outputs.
    inunits, ounits : unit-like or sequence thereof (optional)
        Sequence of the input and output units, respectively.

        .. warning::
            Inputs will be converted to these units before being passed to the
            returned `~numpy.ufunc`. Outputs will be assigned these units.
            Make sure these are the correct units.

    identity : object (optional, keyword-only)
        The value to use for the `~numpy.ufunc.identity` attribute of the
        resulting object. If specified, this is equivalent to setting the
        underlying C ``identity`` field to ``PyUFunc_IdentityValue``.
        If omitted, the identity is set to ``PyUFunc_None``. Note that this is
        _not_ equivalent to setting the identity to ``None``, which implies the
        operation is reorderable.

    assume_correct_units : bool (optional, keyword-only)
        When input arrays are given without units, but the ufunc has 'inunits',
        whether the array is assumed to have dimensionless units (default) or
        have 'inunits'.

    Returns
    -------
    `~numpy.ufunc`
        Registered into `astropy.units.Quantity` `numpy.ufunc` registry.

    See Also
    --------
    `astropy.units.quantity_helper.helpers.register_ufunc`

    Examples
    --------
    We first need to import relevant packages:

        >>> import numpy as np
        >>> import astropy.units as u
        >>> from astropy.units.imperial import Fahrenheit

    Now we can define a python function. For this example we will define the
    conversion between Celsius and Fahrenheit.

        >>> def c2f(x): return 9./5. * x + 32

    With numpy this function can be turned into a `numpy.ufunc`. This is useful
    if the python function works only on scalars, but we want to be able to
    pass in arrays. One of the limitations of a `numpy.ufunc` is that it cannot
    work with `~astropy.units.Quantity`. This is a partially solved problem as
    numpy allows for `numpy.ufunc` evaluation to be overridden. We register
    this ``ufunc`` and provide the input and output units. The internal
    calculation will be done on the unitless arrays (by converting to the input
    units) and then the output units will be assigned.
    ``c2f`` will work on Quantities, but pretending it didn't...

        >>> ufunc = quantity_frompyfunc(c2f, nin=1, nout=1,
        ...                             inunits=u.Celsius, ounits=Fahrenheit)
        >>> ufunc
        <ufunc 'c2f (vectorized)'>

        >>> ufunc(36 * u.Celsius)
        <Quantity 96.8 deg_F>
        >>> ufunc(np.array([0, 10, 20]) * u.Celsius)
        <Quantity [32.0, 50.0, 68.0] deg_F>


    **There are two caveats to note**:

    1. The `numpy.ufunc` overrides only work when at least one argument
       is a `~astropy.units.Quantity`. In the above example ``c2f`` takes only
       one argument, so if a scalar or `~numpy.ndarray` were passed instead of
       a Quantity, the output will also be an ndarray.

        >>> ufunc(36)
        96.8
        >>> ufunc(np.array([0, 10, 20]))  # note dtype is an object
        array([32.0, 50.0, 68.0], dtype=object)

    2. The function cannot return a Quantity with units. If so, an object array
       of Quantity will be returned instead of a Quantity array.

       >>> def badc2f(x): return (9./5. * x + 32) << Fahrenheit
       >>> badufunc = quantity_frompyfunc(badc2f, 1, 1, u.Celsius, Fahrenheit)
       >>> badufunc(np.array([0, 10, 20]) * u.Celsius)
       <Quantity [<Quantity 32. deg_F>, <Quantity 50. deg_F>,
                  <Quantity 68. deg_F>] deg_F>


    **Extra features**:

    As a convenience, ``quantity_frompyfunc`` can also introspect function
    annotations and use these to determine the input and output units,
    obviating the need for arguments ``inunits`` and ``ounits``.

        >>> def c2f(x: u.Celsius) -> Fahrenheit: return 9./5. * x + 32
        >>> ufunc = quantity_frompyfunc(c2f, 1, 1)

        >>> ufunc(-40 * u.Celsius)
        <Quantity -40. deg_F>

    When a ufunc has at least 2 inputs, if one of the arguments does not have
    units it is assumed to be `~astropy.units.dimensionless_unscaled`. However,
    ``quantity_frompyfunc`` takes the keyword argument "assume_correct_units",
    in which case the ufunc will instead interpret a unitless argument as
    having units 'inunits' -- i.e. the correct units.

        >>> def exf(x: u.km, y: u.s) -> u.km**2/u.s: return x ** 2 / y
        >>> exufunc = quantity_frompyfunc(exf, 2, 1, assume_correct_units=True)
        >>> exufunc(3 * u.km, 2)
        <Quantity 4.5 km2 / s>

    """
    from astropy.units import Unit
    from astropy.units.quantity_helper.helpers import (
        _is_ulike, _is_seq_ulike, register_ufunc)

    # -------------------------
    # determine units by introspection

    sig = inspect.signature(func)

    # input units
    if inunits is None:
        svals = tuple(sig.parameters.values())  # sequence[Parameter]
        inunits = [Unit(p.annotation) if _is_seq_ulike(p.annotation) else None
                   for p in svals]

    # output units
    if ounits is None:
        ra = sig.return_annotation
        ra = [ra] if _is_ulike(ra) else ra  # now a sequence, if unit-like
        ounits = ra if _is_seq_ulike(ra) else [None]

        if ounits != [None] and len(ounits) != nout:
            raise ValueError(
                "function annotation is a sequence of unit-like, but "
                f"its length ({len(ra)}) does not equal `nout` ({nout})")

    # -------------------------
    # make and register ufunc

    ufunc = np.frompyfunc(func, nin, nout, identity=identity)
    register_ufunc(ufunc, inunits=inunits, ounits=ounits,
                   assume_correct_units=assume_correct_units)

    return ufunc
