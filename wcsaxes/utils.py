import numpy as np
from astropy import units as u

# Modified from axis_artist, supports astropy.units


def select_step_degree(dv):

    # Modified from axis_artist, supports astropy.units

    degree_limits_ = [1.5, 3, 7, 13, 20, 40, 70, 120, 270, 520]
    degree_steps_ = [1, 2, 5, 10, 15, 30, 45, 90, 180, 360]
    degree_units = [u.degree] * len(degree_steps_)

    minsec_limits_ = [1.5, 2.5, 3.5, 8, 11, 18, 25, 45]
    minsec_steps_ = [1, 2, 3, 5, 10, 15, 20, 30]

    minute_limits_ = np.array(minsec_limits_) / 60.
    minute_units = [u.arcmin] * len(minute_limits_)

    second_limits_ = np.array(minsec_limits_) / 3600.
    second_units = [u.arcsec] * len(second_limits_)

    degree_limits = np.concatenate([second_limits_,
                                    minute_limits_,
                                    degree_limits_])

    degree_steps = minsec_steps_ + minsec_steps_ + degree_steps_
    degree_units = second_units + minute_units + degree_units

    n = degree_limits.searchsorted(dv)
    step = degree_steps[n]
    unit = degree_units[n]

    return step * unit


def select_step_hour(dv):

    hour_limits_ = [1.5, 2.5, 3.5, 5, 7, 10, 15, 21, 36]
    hour_steps_ = [1, 2, 3, 4, 6, 8, 12, 18, 24]
    hour_units = [u.hour] * len(hour_steps_)

    minsec_limits_ = [1.5, 2.5, 3.5, 4.5, 5.5, 8, 11, 14, 18, 25, 45]
    minsec_steps_ = [1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30]

    minute_limits_ = np.array(minsec_limits_) / 60.
    minute_units = [15. * u.arcmin] * len(minute_limits_)

    second_limits_ = np.array(minsec_limits_) / 3600.
    second_units = [15. * u.arcsec] * len(second_limits_)

    hour_limits = np.concatenate([second_limits_,
                                  minute_limits_,
                                  hour_limits_])

    hour_steps = minsec_steps_ + minsec_steps_ + hour_steps_
    hour_units = second_units + minute_units + hour_units

    n = hour_limits.searchsorted(dv)
    step = hour_steps[n]
    unit = hour_units[n]

    return step * unit


def get_coordinate_system(wcs):

    xcoord = wcs.wcs.ctype[0][0:4]
    ycoord = wcs.wcs.ctype[1][0:4]

    from astropy.coordinates import FK5Coordinates, GalacticCoordinates

    if xcoord == 'RA--' and ycoord == 'DEC-':
        coordinate_class = FK5Coordinates
    elif xcoord == 'GLON' and ycoord == 'GLAT':
        coordinate_class = GalacticCoordinates
    else:
        raise ValueError("System not supported (yet): {0}/{1}".format(xcoord, ycoord))

    return coordinate_class
