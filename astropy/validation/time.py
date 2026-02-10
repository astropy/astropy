from .errors import PhysicalInconsistencyError


def validate_time_scale_consistency(times, strict=True):
    scales = {t.scale for t in times}

    if len(scales) > 1:
        if strict:
            raise PhysicalInconsistencyError(
                f"Inconsistent time scales detected: {scales}"
            )
        return False
    return True
