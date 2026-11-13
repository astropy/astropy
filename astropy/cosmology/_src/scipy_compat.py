"""Scipy compatibility."""

__all__ = ("ellipkinc", "hyp2f1", "quad", "solve_ivp")

from typing import Any, Never

from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.integrate import quad, solve_ivp
    from scipy.special import ellipkinc, hyp2f1

else:

    def quad(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.integrate'")

    def solve_ivp(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.integrate'")

    def ellipkinc(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")

    def hyp2f1(*args: Any, **kwargs: Any) -> Never:
        raise ModuleNotFoundError("No module named 'scipy.special'")
