import numpy as np

__all__ = ["construct_regularization", "design_matrix", "periodic_fit"]


def design_matrix(t, bands, frequency, dy=None, nterms_base=1, nterms_band=1):
    t = np.asarray(t)
    omega = np.asarray(2 * np.pi * frequency)
    unique_bands = np.unique(bands)

    # Construct X - design matrix of the stationary sinusoid model
    cols = [np.ones(len(t))]
    cols = sum(
        (
            [np.sin((i + 1) * omega * t), np.cos((i + 1) * omega * t)]
            for i in range(nterms_base)
        ),
        cols,
    )

    # Add band columns for the multiband model, binary flag indicating the band of a given observation
    for band in unique_bands:
        cols.append(np.ones(len(t)))  # A bias column is added by default
        cols = sum(
            (
                [np.sin((i + 1) * omega * t), np.cos((i + 1) * omega * t)]
                for i in range(nterms_band)
            ),
            cols,
        )
        mask = bands == band
        for i in range(-1 - 2 * nterms_band, 0):
            cols[i][~mask] = 0
    if dy is not None:
        XT = np.transpose(np.vstack(cols) / dy)  # weighted
    else:
        XT = np.transpose(np.vstack(cols))

    return XT


def construct_regularization(
    bands, nterms_base=1, nterms_band=1, reg_base=None, reg_band=1e-6
):
    unique_bands = np.unique(bands)
    # Construct Regularization
    if reg_base is None and reg_band is None:
        regularization = 0
    else:
        Nbase = 1 + 2 * nterms_base
        Nband = 1 + 2 * nterms_band
        regularization = np.zeros(Nbase + len(unique_bands) * Nband)
        if reg_base is not None:
            regularization[:Nbase] = reg_base
        if reg_band is not None:
            regularization[Nbase:] = reg_band

    return regularization


def periodic_fit(
    t,
    y,
    dy,
    bands,
    frequency,
    t_fit,
    bands_fit,
    center_data=True,
    nterms_base=1,
    nterms_band=1,
    reg_base=None,
    reg_band=1e-6,
    regularize_by_trace=True,
):
    """Compute the Lomb-Scargle model fit at a given frequency

    Parameters
    ----------
    t, y, dy : float or array-like
        The times, observations, and uncertainties to fit
    bands : str, or array-like
        The bands of each observation
    frequency : float
        The frequency at which to compute the model
    t_fit : float or array-like
        The times at which the fit should be computed
    center_data : bool (default=True)
        If True, center the input data before applying the fit
    nterms : int (default=1)
        The number of Fourier terms to include in the fit

    Returns
    -------
    y_fit : ndarray
        The model fit evaluated at each value of t_fit
    """
    t, y, bands, frequency = map(np.asarray, (t, y, bands, frequency))
    bands_fit = bands_fit[:, np.newaxis]
    unique_bands = np.unique(bands)
    unique_bands_fit = np.unique(bands_fit)

    if not set(unique_bands_fit).issubset(set(unique_bands)):
        raise ValueError(
            "bands_fit does not match training data: "
            f"input: {set(unique_bands_fit)} output: {set(unique_bands)}"
        )

    t_fit, bands_fit = np.broadcast_arrays(t_fit, bands_fit)

    if dy is None:
        dy = np.ones_like(y)
    else:
        dy = np.asarray(dy)

    if t.ndim != 1 or y.ndim != 1 or dy.ndim != 1:
        raise ValueError("t, y, dy should be one dimensional")
    if frequency.ndim != 0:
        raise ValueError("frequency should be a scalar")

    # need to make sure all unique filters are represented
    u, i = np.unique(
        np.concatenate([bands_fit.ravel(), unique_bands]), return_inverse=True
    )

    # Calculate ymeans
    ymeans = np.zeros(
        y.shape
    )  # An array of shape y, with each index given a filter specific mean
    if center_data:
        for band in unique_bands:
            mask = bands == band
            ymeans[mask] = np.average(y[mask], weights=1 / dy[mask] ** 2)
        y = y - ymeans
    ymeans_fit = ymeans[i[: -len(unique_bands)]]

    # Theta -- Construct X and M from t and bands, using weighting
    X = design_matrix(
        t, bands, frequency, dy=dy, nterms_base=nterms_base, nterms_band=nterms_band
    )
    M = np.dot(X.T, X)
    regularization = construct_regularization(
        bands,
        nterms_base=nterms_base,
        nterms_band=nterms_band,
        reg_base=reg_base,
        reg_band=reg_band,
    )

    if regularization is not None:
        diag = M.ravel(order="K")[
            :: M.shape[0] + 1
        ]  # M is being affected by operations on diag
        if regularize_by_trace:
            diag += diag.sum() * np.asarray(regularization)
        else:
            diag += np.asarray(regularization)

    theta_MLE = np.linalg.solve(M, np.dot(X.T, y / dy))

    # Fit to t_fit and bands_fit
    X_fit = design_matrix(
        t_fit.ravel(),
        bands_fit.ravel(),
        frequency,
        dy=None,
        nterms_base=nterms_base,
        nterms_band=nterms_band,
    )

    if center_data:
        y_fit = ymeans_fit + np.dot(X_fit, theta_MLE)
    else:
        y_fit = np.dot(X_fit, theta_MLE)

    return y_fit.reshape(np.shape(t_fit))
