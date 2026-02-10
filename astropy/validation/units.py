import numpy as np

from astropy import units as u

from .errors import PhysicalInconsistencyError


def validate_units_consistency(flux, distance, strict=True):
    """
    Validate whether flux and distance units are physically
    compatible for luminosity computation.

    Parameters
    ----------
    flux : astropy.units.Quantity
        Observed flux (e.g. erg / s / cm^2)

    distance : astropy.units.Quantity
        Distance to the source (e.g. pc, m)

    strict : bool, optional
        If True, raise PhysicalInconsistencyError on failure.
        If False, return False instead.

    Returns
    -------
    bool
        True if units are physically consistent, False otherwise.

    Raises
    ------
    PhysicalInconsistencyError
        If units are incompatible and strict=True.
    """
    try:
        # Luminosity relation: L = 4 * pi * d^2 * F
        luminosity = flux * (4 * np.pi * distance**2)

        # Force unit resolution (must reduce to power)
        luminosity.to(u.erg / u.s)

    except Exception:
        if strict:
            raise PhysicalInconsistencyError(
                "Flux unit is incompatible with distance unit "
                "for luminosity computation."
            )
        return False

    return True
