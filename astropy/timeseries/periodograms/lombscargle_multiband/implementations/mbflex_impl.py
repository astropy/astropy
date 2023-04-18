import numpy as np

__all__ = ["lombscargle_mbflex"]


def lombscargle_mbflex(
    t,
    y,
    bands,
    frequency,
    dy=None,
    nterms_base=1,
    nterms_band=1,
    reg_base=None,
    reg_band=1e-6,
    regularize_by_trace=True,
    center_data=True,
):
    if (nterms_base == 0) and (nterms_band == 0):
        raise ValueError(
            "At least one of nterms_base and nterms_band must be greater than 0."
        )
    # Inputs
    unique_bands = np.unique(bands)
    t = np.asarray(t)
    y = np.asarray(y)
    bands = np.asarray(bands)
    frequency = np.asarray(frequency)

    # Create a ones array for dy (errors) if not provided
    if dy is not None:
        dy = np.asarray(dy)
    else:
        dy = np.ones(y.shape)

    # Calculate ymeans -- per band/filter
    ymeans = np.zeros(
        y.shape
    )  # An array of shape y, with each index given a filter specific mean
    for band in unique_bands:
        mask = bands == band
        ymeans[mask] = np.average(y[mask], weights=1 / dy[mask] ** 2)
    # Construct weighted y matrix
    if center_data:
        y = y - ymeans

    yw = y / dy  # weighted by dy, as above; one's array if not provided

    # Construct Regularization
    if reg_base is None and reg_band is None:
        regularization = 0
    else:
        n_base = 1 + 2 * nterms_base
        n_band = 1 + 2 * nterms_band
        regularization = np.zeros(n_base + len(unique_bands) * n_band)
        if reg_base is not None:
            regularization[:n_base] = reg_base
        if reg_band is not None:
            regularization[n_base:] = reg_band

    # Scoring
    omegas = 2 * np.pi * frequency

    # Calculate chi-squared
    chi2_0 = np.dot(yw.T, yw)
    chi2_ref = np.copy(chi2_0)  # reference chi2 for later comparison

    # Iterate through the omegas and compute the power for each
    chi2_0_minus_chi2 = []
    for i, omega in enumerate(omegas.flat):
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
            cols.append(np.ones(len(t)))
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

        X = np.transpose(np.vstack(cols) / dy)  # weighted
        M = np.dot(X.T, X)

        if regularization is not None:
            diag = M.ravel(order="K")[
                :: M.shape[0] + 1
            ]  # M is being affected by operations on diag
            if regularize_by_trace:
                diag += diag.sum() * np.asarray(regularization)
            else:
                diag += np.asarray(regularization)

        # Construct Xw, XTX, XTy
        Xw = X
        XTX = M
        XTy = np.dot(Xw.T, yw)

        # Matrix Algebra to calculate the Lomb-Scargle power at each omega step
        try:
            chi2_0_minus_chi2.append(np.dot(XTy.T, np.linalg.solve(XTX, XTy)))
        # If X'X is not invertible, use pseudoinverse instead
        except np.linalg.LinAlgError:
            chi2_0_minus_chi2.append(
                np.dot(XTy.T, np.linalg.lstsq(XTX, XTy, rcond=None)[0])
            )

    # Construct and return the power from the chi2 difference
    if center_data:
        P = chi2_0_minus_chi2 / chi2_ref
    else:
        P = 1 + (chi2_0_minus_chi2 - chi2_0) / chi2_ref

    return P
