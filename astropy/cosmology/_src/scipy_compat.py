"""Scipy compatibility."""

__all__ = ("ellipkinc", "hyp2f1", "minimize_scalar", "quad")

from typing import Any, Never

from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.integrate import quad
    from scipy.optimize import minimize_scalar
    from scipy.special import ellipkinc, hyp2f1

else:

    def quad(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.integrate'")

    def ellipkinc(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def hyp2f1(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def minimize_scalar(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.optimize'")
