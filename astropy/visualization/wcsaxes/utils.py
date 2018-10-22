# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from ... import units as u
from ...coordinates import BaseCoordinateFrame

__all__ = ['select_step_degree', 'select_step_hour', 'select_step_scalar',
           'coord_type_from_ctype', 'transform_contour_set_inplace']

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

        return select_step_scalar(dv.to_value(u.arcsec)) * u.arcsec


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

        return select_step_scalar(dv.to_value(15. * u.arcsec)) * (15. * u.arcsec)


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

    from astropy.coordinates import frame_transform_graph

    if isinstance(frame, str):
        initial_frame = frame
        frame = frame_transform_graph.lookup_name(frame)
        if frame is None:
            raise ValueError("Unknown frame: {0}".format(initial_frame))

    if not isinstance(frame, BaseCoordinateFrame):
        frame = frame()

    names = list(frame.representation_component_names.keys())
    coord_meta['name'] = names[:2]

    return coord_meta


def coord_type_from_ctype(ctype):
    """
    Determine whether a particular WCS ctype corresponds to an angle or scalar
    coordinate.
    """
    if ctype[:4] == 'RA--':
        return 'longitude', u.hourangle, None
    elif ctype[:4] == 'HPLN':
        return 'longitude', u.arcsec, 180.
    elif ctype[:4] == 'HPLT':
        return 'latitude', u.arcsec, None
    elif ctype[:4] == 'HGLN':
        return 'longitude', None, 180.
    elif ctype[1:4] == 'LON' or ctype[2:4] == 'LN':
        return 'longitude', None, None
    elif ctype[:4] == 'DEC-' or ctype[1:4] == 'LAT' or ctype[2:4] == 'LT':
        return 'latitude', None, None
    else:
        return 'scalar', None, None


def transform_contour_set_inplace(cset, transform):
    """
    Transform a contour set in-place using a specified
    :class:`matplotlib.transform.Transform`

    Using transforms with the native Matplotlib contour/contourf can be slow if
    the transforms have a non-negligible overhead (which is the case for
    WCS/SkyCoord transforms) since the transform is called for each individual
    contour line. It is more efficient to stack all the contour lines together
    temporarily and transform them in one go.
    """

    # The contours are represented as paths grouped into levels. Each can have
    # one or more paths. The approach we take here is to stack the vertices of
    # all paths and transform them in one go. The pos_level list helps us keep
    # track of where the set of segments for each overall contour level ends.
    # The pos_segments list helps us keep track of where each segmnt ends for
    # each contour level.
    all_paths = []
    pos_level = []
    pos_segments = []

    for collection in cset.collections:
        paths = collection.get_paths()
        all_paths.append(paths)
        # The last item in pos isn't needed for np.split and in fact causes
        # issues if we keep it because it will cause an extra empty array to be
        # returned.
        pos = np.cumsum([len(x) for x in paths])
        pos_segments.append(pos[:-1])
        pos_level.append(pos[-1])

    # As above the last item isn't needed
    pos_level = np.cumsum(pos_level)[:-1]

    # Stack all the segments into a single (n, 2) array
    vertices = [path.vertices for paths in all_paths for path in paths]
    if len(vertices) > 0:
        vertices = np.concatenate(vertices)
    else:
        return

    # Transform all coordinates in one go
    vertices = transform.transform(vertices)

    # Split up into levels again
    vertices = np.split(vertices, pos_level)

    # Now re-populate the segments in the line collections
    for ilevel, vert in enumerate(vertices):
        vert = np.split(vert, pos_segments[ilevel])
        for iseg, ivert in enumerate(vert):
            all_paths[ilevel][iseg].vertices = ivert
