from importlib.metadata import entry_points

__all__ = ['load_all_entry_points']


def load_all_entry_points(group):
    """
    Load all entry points in a given group for the I/O registry.

    Paramters
    ---------
    name : str
        Name of the entry points section, e.g. 'astropy_io_registry_table'.
    """
    eps = entry_points()
    if group in eps:
        for entry_point in eps[group]:
            func = entry_point.load()
            func()
    else:
        raise Exception(f'No group {group} found in entry points')