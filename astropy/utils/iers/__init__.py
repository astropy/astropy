from .iers import *


def __getattr__(name):
    # Deprecated leap-second constants are not in iers.__all__ and are
    # only available through iers.__getattr__, which warns.
    from . import iers

    return iers.__getattr__(name)
