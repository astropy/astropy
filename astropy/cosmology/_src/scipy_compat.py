"""Scipy compatibility."""

__all__ = ("ellipkinc", "hyp2f1", "quad")

from typing import Never

from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.integrate import quad
    from scipy.special import ellipkinc, hyp2f1

else:

    def quad(*args, **kwargs) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.integrate'")

    def ellipkinc(*args, **kwargs) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def hyp2f1(*args, **kwargs) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")
