# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Internal helpers for uncertainty propagation in coordinate transformations.

This module is intentionally minimal and experimental.
"""


def propagate_uncertainty(coord_in, coord_out, uncertainty):
    """
    Placeholder for uncertainty propagation.

    Parameters
    ----------
    coord_in : SkyCoord
        Input coordinate.
    coord_out : SkyCoord
        Output coordinate after transformation.
    uncertainty : dict
        Mapping of component name to Quantity uncertainty.

    Returns
    -------
    dict
        Propagated uncertainties (currently unmodified).
    """
    # Phase-0: no propagation yet
    return uncertainty.copy()
