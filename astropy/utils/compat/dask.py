from .optional_deps import HAS_DASK


def _is_dask_array(data):
    """Check whether data is a dask array."""
    if not HAS_DASK or not hasattr(data, "compute"):
        return False

    from dask.array import Array

    return isinstance(data, Array)
