import numpy as np


NORMALIZATIONS = ['standard', 'psd', 'model', 'log']


def compute_chi2_ref(y, dy, center_data=True, fit_mean=True):
    """Compute the reference chi-square for a particular dataset.

    Note: this is not valid center_data=False and fit_mean=False.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO
    """
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
    TODO

    Returns
    -------
    TODO
    """
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
