# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy import units as u
from astropy.extern import six
from astropy.coordinates import BaseCoordinateFrame

# Modified from axis_artist, supports astropy.units


def select_step_degree(dv):

    # Modified from axis_artist, supports astropy.units

    if dv > 1. * u.arcsec:

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

        n = degree_limits.searchsorted(dv.to(u.degree))
        step = degree_steps[n]
        unit = degree_units[n]

        return step * unit

    else:

        return select_step_scalar(dv.to(u.arcsec).value) * u.arcsec


def select_step_hour(dv):

    if dv > 15. * u.arcsec:

        hour_limits_ = [1.5, 2.5, 3.5, 5, 7, 10, 15, 21, 36]
        hour_steps_ = [1, 2, 3, 4, 6, 8, 12, 18, 24]
        hour_units = [u.hourangle] * len(hour_steps_)

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

        n = hour_limits.searchsorted(dv.to(u.hourangle))
        step = hour_steps[n]
        unit = hour_units[n]

        return step * unit

    else:

        return select_step_scalar(dv.to(15. * u.arcsec).value) * (15. * u.arcsec)


def select_step_scalar(dv):

    log10_dv = np.log10(dv)

    base = np.floor(log10_dv)
    frac = log10_dv - base

    steps = np.log10([1, 2, 5, 10])

    imin = np.argmin(np.abs(frac - steps))

    return 10. ** (base + steps[imin])


def get_coord_meta(frame):

    coord_meta = {}
    coord_meta['type'] = ('longitude', 'latitude')
    coord_meta['wrap'] = (None, None)
    coord_meta['unit'] = (u.deg, u.deg)

    try:

        from astropy.coordinates import frame_transform_graph

        if isinstance(frame, six.string_types):
            initial_frame = frame
            frame = frame_transform_graph.lookup_name(frame)
            if frame is None:
                raise ValueError("Unknown frame: {0}".format(initial_frame))

        if not isinstance(frame, BaseCoordinateFrame):
            frame = frame()

        names = list(frame.representation_component_names.keys())
        coord_meta['name'] = names[:2]

    except ImportError:

        if isinstance(frame, six.string_types):
            if frame in ('fk4', 'fk5', 'icrs'):
                coord_meta['name'] = ('ra', 'dec')
            elif frame == 'galactic':
                coord_meta['name'] = ('l', 'b')
            else:
                raise ValueError("Unknown frame: {0}".format(frame))

    return coord_meta


def coord_type_from_ctype(ctype):
    """
    Determine whether a particular WCS ctype corresponds to an angle or scalar
    coordinate.
    """
    if ctype[:4] in ['RA--'] or ctype[1:4] == 'LON':
        return 'longitude', None
    elif ctype[:4] in ['HPLN']:
        return 'longitude', 180.
    elif ctype[:4] in ['DEC-', 'HPLT'] or ctype[1:4] == 'LAT':
        return 'latitude', None
    else:
        return 'scalar', None
