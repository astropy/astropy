# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from importlib.metadata import entry_points

from astropy.utils.compat import PY_LT_310

__all__ = ['load_all_entry_points']


def load_all_entry_points(group):
    """Load all entry points in a given group.

    Parameters
    ----------
    group : str
        Name of the entry points section, e.g. 'astropy_io_registry_table'.

    Raises
    ------
    Exception
        If ``group`` not found in the set of entry points.
    """
    if eps := (entry_points().get(group) if PY_LT_310 else entry_points(group=group)):
        for entry_point in (eps if PY_LT_310 else (eps[n] for n in eps.names)):
            try:
                registration_func = entry_point.load()
            except Exception as e:  # TODO! which exceptio type
                warnings.warn(e)
            else:
                registration_func()
    else:
        raise Exception(f'No group {group} found in entry points')
