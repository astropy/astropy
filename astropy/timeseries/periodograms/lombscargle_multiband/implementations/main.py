"""
Main Lomb-Scargle Multiband Implementation

The ``lombscarglemultiband`` function here is essentially a switch
statement for the various implementations available in this submodule
"""

__all__ = ["lombscargle_multiband", "available_methods"]

from .mbfast_impl import lombscargle_mbfast
from .mbflex_impl import lombscargle_mbflex


def available_methods():
    methods = ["fast", "flexible"]
    return methods


def lombscargle_multiband(
    t,
    y,
    bands,
    dy=None,
    frequency=None,
    method="flexible",
    sb_method="auto",
    assume_regular_frequency=False,
    normalization="standard",
    fit_mean=True,
    center_data=True,
    method_kwds=None,
    nterms_base=1,
    nterms_band=1,
    reg_base=None,
    reg_band=1e-6,
    regularize_by_trace=True,
    fit_period=False,
):
    methods = available_methods()

    if method == "flexible":
        power = lombscargle_mbflex(
            t,
            y,
            bands,
            frequency,
            dy=dy,
            nterms_base=nterms_base,
            nterms_band=nterms_band,
            reg_base=reg_base,
            reg_band=reg_band,
            regularize_by_trace=regularize_by_trace,
            center_data=center_data,
        )
    elif method == "fast":
        power = lombscargle_mbfast(
            t,
            y,
            bands,
            dy,
            frequency=frequency,
            sb_method=sb_method,
            assume_regular_frequency=assume_regular_frequency,
            normalization=normalization,
            fit_mean=fit_mean,
            center_data=center_data,
            method_kwds=method_kwds,
            nterms=nterms_base,
        )  # nterms_base used for single_band nterms
    if method not in methods:
        raise ValueError(f"Invalid Method: {method}")

    return power
