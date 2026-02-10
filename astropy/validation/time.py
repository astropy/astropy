from .errors import PhysicalInconsistencyError


def validate_time_scale_consistency(times, strict=True):
    """
    Validate that all provided Time objects use the same time scale.

    Parameters
    ----------
    times : astropy.time.Time or iterable of Time
        One or more Time objects to validate.

    strict : bool, optional
        If True, raise PhysicalInconsistencyError on failure.
        If False, return False instead.

    Returns
    -------
    bool
        True if time scales are consistent, False otherwise.
    """
    try:
        # Normalize input to an iterable of Time-like objects
        iterable = times if hasattr(times, "__iter__") else [times]
        scales = {t.scale for t in iterable}
    except Exception:
        if strict:
            raise PhysicalInconsistencyError("Invalid time object(s) provided.")
        return False

    if len(scales) > 1:
        if strict:
            raise PhysicalInconsistencyError(
                f"Inconsistent time scales detected: {scales}"
            )
        return False

    return True
