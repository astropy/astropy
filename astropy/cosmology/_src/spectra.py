from __future__ import annotations

from abc import ABC, abstractmethod
from numbers import Integral, Real
from typing import Any

import numpy as np


def _window_tophat(x: np.ndarray) -> np.ndarray:
    """Fourier-space spherical top-hat window.

    W(x) = 3 (sin x - x cos x) / x^3
    """
    out = np.ones_like(x, dtype=float)
    m = x != 0.0
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / (xm * xm * xm)
    return out


def _as_float_array(x: Any, *, name: str) -> np.ndarray:
    """Convert input to a float ndarray and require positivity.

    Args:
        x: Input value(s).
        name: Name used in error messages.

    Returns
    -------
        Float ndarray view/copy of `x`.

    Raises
    ------
        ValueError: If any element is <= 0.
    """
    arr = np.asarray(x, dtype=float)
    if np.any(arr <= 0.0):
        msg = f"All {name} values must be positive."
        raise ValueError(msg)
    return arr


class PowerSpectrum(ABC):
    """
    Abstract base class for power spectrum models.

    Subclasses implement `_compute_spectrum(k)` and `_compute_integral(x)`.

    The public call `ps(k)` evaluates the spectrum with:
        - an overall prefactor `prefactor = 1 / (2*pi)^ndims`
        - an optional multiplicative `scale`
        - an optional Gaussian cutoff exp(-(k/k_cutoff)^2)

    Notes
    -----
        The `integrate(x)` method is intentionally *not* multiplied by the
        prefactor again. Subclasses should define integrals in terms of the
        evaluated spectrum `self(k)` where appropriate.
    """

    def __init__(
        self,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        """Initialize the spectrum model.

        Args:
            ndims: Dimensionality used to define the prefactor
                `prefactor = 1 / (2*pi)^ndims`.
            scale: Optional multiplicative scale applied to the spectrum.
                Defaults to 1.0.
            k_cutoff: Optional Gaussian cutoff scale. If provided, the spectrum
                is multiplied by exp(-(k/k_cutoff)^2).
        """
        self.ndims = int(ndims)
        self.scale = 1.0 if scale is None else float(scale)
        self.k_cutoff = None if k_cutoff is None else float(k_cutoff)

        self.prefactor = 1.0 / (2.0 * np.pi) ** self.ndims
        self.spectrum_type: str | None = None

        # integration defaults
        self.k_min = 1.0e-5
        self.k_max = 1.0e2
        self.nk = 4096

    @abstractmethod
    def _compute_spectrum(self, k: Any) -> Any:
        """Compute the spectrum P(k) without prefactor/scale/cutoff.

        Args:
            k: Wavenumber(s).

        Returns
        -------
            Spectrum values with the same shape as `k`.
        """
        raise NotImplementedError

    def _cutoff_multiplier(self, k: Any) -> Any:
        """
        Compute Gaussian cutoff multiplier for the given k.

        Args:
            k: Wavenumber(s).

        Returns
        -------
            exp(-(k/k_cutoff)^2) evaluated elementwise, or 1 if no cutoff.
        """
        if self.k_cutoff is None:
            return 1.0

        k_arr = np.asarray(k, dtype=float)
        return np.exp(-(k_arr * k_arr) / (self.k_cutoff * self.k_cutoff))

    def __call__(self, k: Any) -> Any:
        """
        Evaluate the spectrum at wavenumber(s) k.

        Args:
            k: Wavenumber(s).

        Returns
        -------
            Evaluated spectrum including prefactor, scale, and optional cutoff.
        """
        base = self._compute_spectrum(k)
        cutoff = self._cutoff_multiplier(k)
        return self.scale * self.prefactor * base * cutoff

    def integrate(self, R: Any) -> Any:
        """Compute sigma(R) using a spherical top-hat window."""
        R_arr = np.asarray(R, dtype=float)

        if np.any(R_arr <= 0):
            msg = "R must be positive."
            raise ValueError(msg)

        k = np.logspace(np.log10(self.k_min), np.log10(self.k_max), self.nk)
        Pk = np.asarray(self(k), dtype=float)

        k_col = k[:, None]
        R_row = np.atleast_1d(R_arr)[None, :]

        W = _window_tophat(k_col * R_row)

        integrand = (k_col * k_col) * Pk[:, None] * (W * W)

        sigma2 = (1.0 / (2.0 * np.pi * np.pi)) * np.trapezoid(integrand, k, axis=0)
        sigma = np.sqrt(sigma2)

        if np.ndim(R) == 0:
            return float(sigma[0])

        return sigma


class AnalyticalPowerSpectrum(PowerSpectrum):
    """Base class for analytical (closed-form) power spectra."""

    def __init__(
        self,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        """Initialize an analytical spectrum model.

        Args:
            ndims: Dimensionality used to define the prefactor.
            scale: Optional multiplicative scale applied to the spectrum.
            k_cutoff: Optional Gaussian cutoff scale.
        """
        super().__init__(ndims=ndims, scale=scale, k_cutoff=k_cutoff)
        self.spectrum_type = "analytical"


class ScaleInvariantSpectrum(AnalyticalPowerSpectrum):
    """
    Harrison-Zel'dovich (scale-invariant) primordial spectrum.

    This spectrum is the strictly scale-invariant primordial curvature spectrum:

        P(k) = A_s.

    Args:
        A_s: Amplitude (dimensionless).
        ndims: Dimensionality used by the base prefactor.
        scale: Optional multiplicative scale applied to the evaluated spectrum.
        k_cutoff: Optional Gaussian cutoff scale applied to the evaluated spectrum.

    References
    ----------
        Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*,
        arXiv:1807.06209 (2018). Parameter conventions and typical amplitudes.
    """

    def __init__(
        self,
        A_s: Real = 2.1e-9,
        *,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        super().__init__(
            ndims=ndims,
            scale=scale,
            k_cutoff=k_cutoff,
        )
        self.A_s = float(A_s)

    def _compute_spectrum(self, k: Any) -> Any:
        k_arr = _as_float_array(k, name="k")
        out = np.full_like(k_arr, self.A_s, dtype=float)

        if np.ndim(k) == 0:
            return float(out.item())

        return out


class PowerLawSpectrum(AnalyticalPowerSpectrum):
    r"""Primordial power-law curvature spectrum.

    Implements the standard primordial (inflationary) power-law form:

        P(k) = A_s * (k / k_pivot)^(n_s - 1).

    Args:
        A_s: Primordial curvature amplitude at the pivot scale (dimensionless).
        n_s: Scalar spectral index.
        k_pivot: Pivot scale in the same units as k. Must be positive.
        ndims: Dimensionality used by the base prefactor `1 / (2*pi)^ndims`.
        scale: Optional multiplicative scale applied to the evaluated spectrum.
        k_cutoff: Optional Gaussian cutoff scale applied to the evaluated spectrum.

    Notes
    -----
        The common defaults `A_s=2.1e-9`, `n_s=0.965`, and `k_pivot=0.05` are
        consistent with Planck 2018 base-ΛCDM conventions (A_s defined at
        k = 0.05 Mpc^-1) and typical best-fit values.

    References
    ----------
        Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*,
        arXiv:1807.06209 (2018). (See parameter definitions: A_s at k=0.05 Mpc^-1,
        and reported n_s ≈ 0.965.)
    """

    def __init__(
        self,
        A_s: Real = 2.1e-9,
        n_s: Real = 0.965,
        k_pivot: Real = 0.05,
        *,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        if k_pivot <= 0:
            msg = "k_pivot must be positive."
            raise ValueError(msg)

        super().__init__(ndims=ndims, scale=scale, k_cutoff=k_cutoff)
        self.A_s = float(A_s)
        self.n_s = float(n_s)
        self.k_pivot = float(k_pivot)

    def _compute_spectrum(self, k: Any) -> Any:
        k_arr = _as_float_array(k, name="k")
        out = self.A_s * (k_arr / self.k_pivot) ** (self.n_s - 1.0)

        if np.ndim(k) == 0:
            return float(out.item())

        return out


class RunningPowerLawSpectrum(AnalyticalPowerSpectrum):
    r"""
    Primordial power spectrum with running of the tilt.

    Uses the common expansion:
        ln P(k) = ln A_s + (n_s - 1) ln(k/k_pivot) + 0.5 * alpha_s [ln(k/k_pivot)]^2.

    Equivalently:
        P(k) = A_s * exp((n_s - 1)L + 0.5*alpha_s*L^2),  L = ln(k/k_pivot).

    Args:
        A_s: Amplitude at k_pivot.
        n_s: Scalar tilt at k_pivot.
        alpha_s: Running dn_s/dlnk at k_pivot.
        k_pivot: Pivot scale (same units as k). Must be positive.
        ndims: Dimensionality used by the base prefactor.
        scale: Optional multiplicative scale applied to the evaluated spectrum.
        k_cutoff: Optional Gaussian cutoff scale.

    References
    ----------
        Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*,
        arXiv:1807.06209 (2018). (Standard parameterization and constraints on
        running in common extensions.)
    """

    def __init__(
        self,
        A_s: Real = 2.1e-9,
        n_s: Real = 0.965,
        alpha_s: Real = 0.0,
        k_pivot: Real = 0.05,
        *,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        if k_pivot <= 0:
            msg = "k_pivot must be positive."
            raise ValueError(msg)

        super().__init__(ndims=ndims, scale=scale, k_cutoff=k_cutoff)
        self.A_s = float(A_s)
        self.n_s = float(n_s)
        self.alpha_s = float(alpha_s)
        self.k_pivot = float(k_pivot)

    def _compute_spectrum(self, k: Any) -> Any:
        k_arr = _as_float_array(k, name="k")
        L = np.log(k_arr / self.k_pivot)
        out = np.exp(
            np.log(self.A_s) + (self.n_s - 1.0) * L + 0.5 * self.alpha_s * (L * L)
        )

        if np.ndim(k) == 0:
            return float(out.item())

        return out


class LogOscillationSpectrum(AnalyticalPowerSpectrum):
    r"""Power-law spectrum with a log-oscillatory feature.

    P(k) = P_pl(k) * [1 + A_osc * sin(omega * ln(k/k_pivot) + phi)],

    where P_pl(k) is the underlying power-law.

    Args:
        A_s: Power-law amplitude at k_pivot.
        n_s: Power-law tilt.
        k_pivot: Pivot scale.
        A_osc: Oscillation amplitude (non-negative).
        omega: Oscillation frequency in log-k.
        phi: Phase in radians.
        ndims: Dimensionality used by the base prefactor.
        scale: Optional multiplicative scale applied to the evaluated spectrum.
        k_cutoff: Optional Gaussian cutoff scale.

    Notes
    -----
        For physical positivity of P(k), typically require A_osc << 1
        (since 1 + A_osc*sin(...) can become negative if A_osc > 1).

    References
    ----------
        Planck Collaboration, *Planck 2018 results. VI. Cosmological parameters*,
        arXiv:1807.06209 (2018). (Underlying power-law primordial convention.)
    """

    def __init__(
        self,
        A_s: Real = 2.1e-9,
        n_s: Real = 0.965,
        k_pivot: Real = 0.05,
        *,
        A_osc: Real = 0.0,
        omega: Real = 0.0,
        phi: Real = 0.0,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        if k_pivot <= 0:
            msg = "k_pivot must be positive."
            raise ValueError(msg)
        if A_osc < 0:
            msg = "A_osc must be non-negative."
            raise ValueError(msg)

        super().__init__(
            ndims=ndims,
            scale=scale,
            k_cutoff=k_cutoff,
        )

        self.A_s = float(A_s)
        self.n_s = float(n_s)
        self.k_pivot = float(k_pivot)
        self.A_osc = float(A_osc)
        self.omega = float(omega)
        self.phi = float(phi)

    def _compute_spectrum(self, k: Any) -> Any:
        k_arr = _as_float_array(k, name="k")
        pl = self.A_s * (k_arr / self.k_pivot) ** (self.n_s - 1.0)

        if self.A_osc == 0.0 or self.omega == 0.0:
            out = pl
        else:
            L = np.log(k_arr / self.k_pivot)
            out = pl * (1.0 + self.A_osc * np.sin(self.omega * L + self.phi))

        if np.ndim(k) == 0:
            return float(out.item())

        return out


class BrokenPowerLawSpectrum(AnalyticalPowerSpectrum):
    r"""
    Broken power-law spectrum with a smooth transition.

    This is a phenomenological form often used for “knee” spectra:

        P(k) = A * (k/k_break)^n1 * [1 + (k/k_break)^(1/s)]^{s(n2-n1)},

    where `s = smooth` controls the sharpness of the transition.

    Args:
        A: Amplitude scale.
        k_break: Break scale (same units as k). Must be positive.
        n1: Low-k slope exponent.
        n2: High-k slope exponent.
        smooth: Smoothness parameter (>0). Smaller -> sharper.
        ndims: Dimensionality used by the base prefactor.
        scale: Optional multiplicative scale applied to the evaluated spectrum.
        k_cutoff: Optional Gaussian cutoff scale.

    Notes
    -----
        The implementation uses a log-form evaluation for stability.

    References
    ----------
        (Phenomenological; used broadly across cosmology/SGWB fitting. No single
        canonical source. If you have a target application (PTA, strings,
        phase transitions), we should cite that specific literature.)
    """

    def __init__(
        self,
        A: Real = 1.0,
        k_break: Real = 1.0,
        n1: Real = 0.0,
        n2: Real = -3.0,
        *,
        smooth: Real = 0.1,
        ndims: Integral = 3,
        scale: Real | None = None,
        k_cutoff: Real | None = None,
    ) -> None:
        if k_break <= 0:
            msg = "k_break must be positive."
            raise ValueError(msg)
        if smooth <= 0:
            msg = "smooth must be positive."
            raise ValueError(msg)

        super().__init__(ndims=ndims, scale=scale, k_cutoff=k_cutoff)
        self.A = float(A)
        self.k_break = float(k_break)
        self.n1 = float(n1)
        self.n2 = float(n2)
        self.smooth = float(smooth)

    def _compute_spectrum(self, k: Any) -> Any:
        k_arr = _as_float_array(k, name="k")
        x = k_arr / self.k_break

        lnP = np.log(self.A) + self.n1 * np.log(x)
        ln_fac = self.smooth * (self.n2 - self.n1) * np.log1p(x ** (1.0 / self.smooth))
        out = np.exp(lnP + ln_fac)

        if np.ndim(k) == 0:
            return float(out.item())

        return out
