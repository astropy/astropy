import numpy as np

from astropy.timeseries.periodograms.lombscargle.implementations import lombscargle

__all__ = ["lombscargle_mbfast"]


def lombscargle_mbfast(
    t,
    y,
    bands,
    dy=None,
    frequency=None,
    sb_method="auto",
    assume_regular_frequency=False,
    normalization="standard",
    fit_mean=True,
    center_data=True,
    method_kwds=None,
    nterms=1,
):
    # create masks for each filter/bandpass
    unique_bands = np.unique(bands)
    masks = [(bands == band) for band in unique_bands]

    # calculate singleband powers for each filter/bandpass
    if dy is not None:
        models = [
            lombscargle(
                t[mask],
                y[mask],
                dy=dy[mask],
                frequency=frequency,
                method=sb_method,
                assume_regular_frequency=assume_regular_frequency,
                normalization=normalization,
                fit_mean=fit_mean,
                center_data=center_data,
                method_kwds=method_kwds,
                nterms=nterms,
            )
            for mask in masks
        ]
    else:
        models = [
            lombscargle(
                t[mask],
                y[mask],
                dy=None,
                frequency=frequency,
                method=sb_method,
                assume_regular_frequency=assume_regular_frequency,
                normalization=normalization,
                fit_mean=fit_mean,
                center_data=center_data,
                method_kwds=method_kwds,
                nterms=nterms,
            )
            for mask in masks
        ]

    # Total score is the sum of powers weighted by chi2-normalization
    powers = np.array(models)
    chi2_0 = np.array([np.sum(model**2) for model in models])
    return np.dot(chi2_0 / chi2_0.sum(), powers)
