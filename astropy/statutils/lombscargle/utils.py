import numpy as np


NORMALIZATIONS = ['standard', 'psd', 'model', 'log']


def compute_chi2_ref(y, dy=None, center_data=True, fit_mean=True):
    """Compute the reference chi-square for a particular dataset.

    Note: this is not valid center_data=False and fit_mean=False.

    Parameters
    ----------
    y : array_like
        data values
    dy : float, array, or None (optional)
        data uncertainties
    center_data : boolean
        specify whether data should be pre-centered
    fit_mean : boolean
        specify whether model should fit the mean of the data

    Returns
    -------
    chi2_ref : float
        The reference chi-square for the periodogram of this data
    """
    if dy is None:
        dy = 1
    y, dy = np.broadcast_arrays(y, dy)
    w = dy ** -2.0
    if center_data or fit_mean:
        mu = np.dot(w, y) / w.sum()
    else:
        mu = 0
    yw = (y - mu) / dy
    return np.dot(yw, yw)


def convert_normalization(Z, N, from_normalization, to_normalization,
                          chi2_ref=None):
    """Convert power from one normalization to another.

    This currently only works for standard & floating-mean models.

    Parameters
    ----------
    Z : array_like
        the periodogram output
    N : integer
        the number of data points
    from_normalization, to_normalization : strings
        the normalization to convert from and to. Options are
        ['standard', 'model', 'log', 'psd']
    chi2_ref : float
        The reference chi-square, required for converting to or from the
        psd normalization.

    Returns
    -------
    Z_out : ndarray
        The periodogram in the new normalization
    """
    Z = np.asarray(Z)
    from_to = (from_normalization, to_normalization)

    for norm in from_to:
        if norm not in NORMALIZATIONS:
            raise ValueError("{0} is not a valid normalization"
                             "".format(from_normalization))

    if from_normalization == to_normalization:
        return Z

    if "psd" in from_to and chi2_ref is None:
        raise ValueError("must supply reference chi^2 when converting "
                         "to or from psd normalization")

    if from_to == ('log', 'standard'):
        return 1 - np.exp(-Z)
    elif from_to == ('standard', 'log'):
        return -np.log(1 - Z)
    elif from_to == ('log', 'model'):
        return np.exp(Z) - 1
    elif from_to == ('model', 'log'):
        return np.log(Z + 1)
    elif from_to == ('model', 'standard'):
        return Z / (1 + Z)
    elif from_to == ('standard', 'model'):
        return Z / (1 - Z)
    elif from_normalization == "psd":
        return convert_normalization(2 / chi2_ref * Z, N,
                                     from_normalization='standard',
                                     to_normalization=to_normalization)
    elif to_normalization == "psd":
        Z_standard = convert_normalization(Z, N,
                                           from_normalization=from_normalization,
                                           to_normalization='standard')
        return 0.5 * chi2_ref * Z_standard
    else:
        raise NotImplementedError("conversion from '{0}' to '{1}'"
                                  "".format(from_normalization,
                                            to_normalization))
