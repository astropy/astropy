# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Configure the tests for :mod:`astropy.cosmology`."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from astropy.cosmology.tests.helper import clean_registry  # noqa: F401
from astropy.tests.helper import pickle_protocol  # noqa: F401

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping, Sequence

K = TypeVar("K")
V = TypeVar("V")


def filter_keys_from_items(
    m: Mapping[K, V], /, filter_out: Sequence[K]
) -> Iterable[K, V]:
    """Filter ``m``, returning key-value pairs not including keys in ``filter``.

    Parameters
    ----------
    m : mapping[K, V]
        A mapping from which to remove keys in ``filter_out``.
    filter_out : sequence[K]
        Sequence of keys to filter out from ``m``.

    Returns
    -------
    iterable[K, V]
        Iterable of ``(key, value)`` pairs with the ``filter_out`` keys removed.
    """
    return ((k, v) for k, v in m.items() if k not in filter_out)
