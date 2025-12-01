import inspect
from collections.abc import Callable


def is_positional_only(func: Callable, /, param: str) -> bool:
    """Return True if ``param:str`` is a positional-only parameter.

    Parameters
    ----------
    param : str
        Parameter name to check
    """
    sig = inspect.signature(func)
    p = sig.parameters.get(param)
    return p is not None and p.kind == inspect.Parameter.POSITIONAL_ONLY
