import inspect


def is_positional_only(func, param_name="z"):
    """Return True if ``param_name`` is a positional-only parameter.

    Parameters
    ----------
    func : callable
        Function to inspect.
    param_name : str
        Parameter name to check (default: ``'z'``).
    """
    sig = inspect.signature(func)
    p = sig.parameters.get(param_name)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY


__all__ = ["is_positional_only"]
