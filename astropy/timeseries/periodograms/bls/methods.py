# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["bls_fast", "bls_slow"]

import numpy as np
from functools import partial

from ._impl import bls_impl


def bls_slow(t, y, ivar, period, duration, oversample, use_likelihood):
    """Compute the periodogram using a brute force reference method

    t : array_like
        Sequence of observation times.
    y : array_like
        Sequence of observations associated with times t.
    ivar : array_like
        The inverse variance of ``y``.
    period : array_like
        The trial periods where the periodogram should be computed.
    duration : array_like
        The durations that should be tested.
    oversample :
        The resolution of the phase grid in units of durations.
    use_likeliood : bool
        If true, maximize the log likelihood over phase, duration, and depth.

    Returns
    -------
    power : array_like
        The periodogram evaluated at the periods in ``period``.
    depth : array_like
        The estimated depth of the maximum power model at each period.
    depth_err : array_like
        The 1-sigma uncertainty on ``depth``.
    duration : array_like
        The maximum power duration at each period.
    transit_time : array_like
        The maximum power phase of the transit in units of time. This
        indicates the mid-transit time and it will always be in the range
        (0, period).
    depth_snr : array_like
        The signal-to-noise with which the depth is measured at maximum power.
    log_likelihood : array_like
        The log likelihood of the maximum power model.

    """
    f = partial(_bls_slow_one, t, y, ivar, duration,
                oversample, use_likelihood)
    return _apply(f, period)


def bls_fast(t, y, ivar, period, duration, oversample, use_likelihood):
    """Compute the periodogram using an optimized Cython implementation

    t : array_like
        Sequence of observation times.
    y : array_like
        Sequence of observations associated with times t.
    ivar : array_like
        The inverse variance of ``y``.
    period : array_like
        The trial periods where the periodogram should be computed.
    duration : array_like
        The durations that should be tested.
    oversample :
        The resolution of the phase grid in units of durations.
    use_likeliood : bool
        If true, maximize the log likelihood over phase, duration, and depth.

    Returns
    -------
    power : array_like
        The periodogram evaluated at the periods in ``period``.
    depth : array_like
        The estimated depth of the maximum power model at each period.
    depth_err : array_like
        The 1-sigma uncertainty on ``depth``.
    duration : array_like
        The maximum power duration at each period.
    transit_time : array_like
        The maximum power phase of the transit in units of time. This
        indicates the mid-transit time and it will always be in the range
        (0, period).
    depth_snr : array_like
        The signal-to-noise with which the depth is measured at maximum power.
    log_likelihood : array_like
        The log likelihood of the maximum power model.

    """
    return bls_impl(
        t, y, ivar, period, duration, oversample, use_likelihood
    )


def _bls_slow_one(t, y, ivar, duration, oversample, use_likelihood, period):
    """A private function to compute the brute force periodogram result"""
    best = (-np.inf, None)
    hp = 0.5*period
    min_t = np.min(t)
    for dur in duration:

        # Compute the phase grid (this is set by the duration and oversample).
        d_phase = dur / oversample
        phase = np.arange(0, period+d_phase, d_phase)

        for t0 in phase:
            # Figure out which data points are in and out of transit.
            m_in = np.abs((t-min_t-t0+hp) % period - hp) < 0.5*dur
            m_out = ~m_in

            # Compute the estimates of the in and out-of-transit flux.
            ivar_in = np.sum(ivar[m_in])
            ivar_out = np.sum(ivar[m_out])
            y_in = np.sum(y[m_in] * ivar[m_in]) / ivar_in
            y_out = np.sum(y[m_out] * ivar[m_out]) / ivar_out

            # Use this to compute the best fit depth and uncertainty.
            depth = y_out - y_in
            depth_err = np.sqrt(1.0 / ivar_in + 1.0 / ivar_out)
            snr = depth / depth_err

            # Compute the log likelihood of this model.
            loglike = -0.5*np.sum((y_in - y[m_in])**2 * ivar[m_in])
            loglike += 0.5*np.sum((y_out - y[m_in])**2 * ivar[m_in])

            # Choose which objective should be used for the optimization.
            if use_likelihood:
                objective = loglike
            else:
                objective = snr

            # If this model is better than any before, keep it.
            if depth > 0 and objective > best[0]:
                best = (
                    objective,
                    (objective, depth, depth_err, dur, (t0+min_t) % period,
                     snr, loglike)
                )

    return best[1]


def _apply(f, period):
    return tuple(map(np.array, zip(*map(f, period))))
