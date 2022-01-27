# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Helper functions for testing :mod:`astropy.cosmology`."""

##############################################################################
# IMPORTS

# STDLIB
from collections.abc import Mapping

# THIRD-PARTY
import numpy as np

# LOCAL
from astropy.cosmology import Cosmology
from astropy.cosmology import io  # make sure everything is registered


###############################################################################
# FUNCTIONS


def _recursive_dict_eq(map1, map2):
    if not isinstance(map1, Mapping) or not isinstance(map2, Mapping):
        return False

    elif set(map1.keys()) != set(map2.keys()):
        return False

    for k, v in map1.items():
        try:  # Attempt normal comparison
            eq = (map2[k] == v)
        except ValueError:  # Maybe it's a NumPy array
            try:
                eq = np.array_equal(map2[k], v)
            except DeprecationWarning:  # Some element-wise failures. Maybe mappings?
                if isinstance(v, Mapping) and isinstance(map2[k], Mapping):
                    eq = _recursive_dict_eq(v, map2[k])
                else:
                    eq = False
        if not np.all(eq):
            return False

    return True


def cosmology_equal(cosmo1, cosmo2, *, check_meta=True, format=False):
    r"""Check equality between Cosmologies.

    Parameters
    ----------
    cosmo1, cosmo2 : `~astropy.cosmology.Cosmology` subclass instance or Any
        The objects to check for equality .
    check_meta : bool, optional keyword-only
        Whether to also check the metadata when determining equality.
    format : bool or None or str, optional keyword-only
        Whether to allow, before equality is checked, ``cosmo1`` and ``cosmo2``
        to be converted to a |Cosmology|. This allows, e.g. a |Table| to be
        equal to a Cosmology.
        `False` (default) will not allow conversion. `True` or `None` will,
        and will use the auto-identification to try to infer the correct
        format. A `str` is assumed to be the correct format to use when
        converting.

    Returns
    -------
    bool
        `True` if cosmologies are equal, `False` cosmo2wise.

    See Also
    --------
    astropy.cosmology.Cosmology.is_equivalent
        Cosmologies may be equivalent, even if not the same class or name
        or metadata.

    Examples
    --------
    For a simple check that the two objects are |Cosmology|, with the same
    names and parameter values, use the standard python equality operator.

        >>> from astropy.cosmology import Planck13, Planck18
        >>> Planck18 == Planck18
        True

        >>> Planck13 != Planck18
        True

    This can also be checked with ``cosmology_equal``:

        >>> cosmology_equal(Planck18, Planck18)
        True

        >>> cosmology_equal(Planck18, Planck13)
        False

    By default ``cosmology_equal`` also check that the metadata
    (:attr:`astropy.Cosmology.meta`) are equal, with the keyword argument
    ``check_meta``.

        >>> cosmo = Planck18.clone(name="Planck18", meta=dict(info="new"))
        >>> cosmology_equal(Planck18, cosmo, check_meta=True)  # <- the default
        False

        >>> cosmology_equal(Planck18, cosmo, check_meta=False)
        True

    When the cosmologies have different names or parameter values they
    are never equal.

        >>> cosmology_equal(Planck18, Planck13)
        False

    Using the keyword argument ``format``, the notion of equality is
    extended to any Python object that can be converted to a |Cosmology|.

        >>> tbl = Planck18.to_format("astropy.table")
        >>> cosmology_equal(Planck18, tbl, format=True)
        True

        >>> cosmology_equal(tbl, Planck18, format=True)
        True

        >>> cosmology_equal(tbl, tbl, format=True)
        True

    The list of valid formats, e.g. the |Table| in this example, may be
    checked with ``Cosmology.from_format.list_formats()``.

    As can be seen in the list of formats, not all formats can be
    auto-identified by ``Cosmology.from_format.registry``. Objects of
    these kinds can still be checked for equality, but the correct
    format string must be used.

        >>> tbl = Planck18.to_format("yaml")
        >>> cosmology_equal(Planck18, tbl, format="yaml")
        True

    Only one ``format`` may be specified, so if both objects are different
    types and neither can be auto-identified, at least one of the objects must
    be converted to a |Cosmology|.
    """
    # Allow for different formats to be considered equal.
    if format is not False:
        format = None if format is True else format  # str->str, None/True->None
        try:
            cosmo1 = Cosmology.from_format(cosmo1, format=format)
        except Exception:  # TODO! should enforce only TypeError
            pass
        try:
            cosmo2 = Cosmology.from_format(cosmo2, format=format)
        except Exception:  # TODO! should enforce only TypeError
            pass

    # Parameter equality
    eq = cosmo1.__eq__(cosmo2) if hasattr(cosmo1, "__eq__") else NotImplemented
    if eq is NotImplemented and hasattr(cosmo2, "__eq__"):
        eq = cosmo2.__eq__(cosmo1)  # that failed, try from 'cosmo2'

    # Metadata
    if check_meta and eq is True:  # Only check if required.
        if not hasattr(cosmo1, "meta") or not hasattr(cosmo2, "meta"):
            eq = False
        else:
            eq &= _recursive_dict_eq(cosmo1.meta, cosmo2.meta)

    return eq if isinstance(eq, bool) else False  # Ensure boolean
