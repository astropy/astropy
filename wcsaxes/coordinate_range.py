import numpy as np

from . import settings

# Algorithm inspired by PGSBOX from WCSLIB by M. Calabretta


def wrap_180(values):
    values_new = values % 360.
    values_new[values_new > 180.] -= 360
    return values_new


def find_coordinate_range(transform, extent, coord_types):
    '''
    Find the range of coordinates to use for ticks/grids

    Parameters
    ----------
    pix2world : func
        Function to transform pixel to world coordinates. Should take two
        values (the pixel coordinates) and return two values (the world
        coordinates).
    extent : iterable
        The range of the image viewport in pixel coordinates, given as [xmin,
        xmax, ymin, ymax].
    coord_types : list of str
        Whether each coordinate is a ``'longitude'``, ``'latitude'``, or
        ``'scalar'`` value.
    '''

    # Sample coordinates on a NX x NY grid.
    NX = NY = settings.COORDINATE_RANGE_SAMPLES
    x = np.linspace(extent[0], extent[1], NX + 1)
    y = np.linspace(extent[2], extent[3], NY + 1)
    xp, yp = np.meshgrid(x, y)
    world = transform.transform(np.vstack([xp.ravel(), yp.ravel()]).transpose())

    ranges = []

    for coord_index, coord_type in enumerate(coord_types):

        xw = world[:, coord_index].reshape(xp.shape)

        if coord_type in ['longitude', 'latitude']:

            # Iron out coordinates along first row
            wjump = xw[0, 1:] - xw[0, :-1]
            reset = np.abs(wjump) > 180.
            if np.any(reset):
                wjump = wjump + np.sign(wjump) * 180.
                wjump = 360. * (wjump / 360.).astype(int)
                xw[0, 1:][reset] -= wjump[reset]

            # Now iron out coordinates along all columns, starting with first row.
            wjump = xw[1:] - xw[:1]
            reset = np.abs(wjump) > 180.
            if np.any(reset):
                wjump = wjump + np.sign(wjump) * 180.
                wjump = 360. * (wjump / 360.).astype(int)
                xw[1:][reset] -= wjump[reset]

        xw_min = np.nanmin(xw)
        xw_max = np.nanmax(xw)

        # Check if range is smaller when normalizing to the range 0 to 360

        if coord_type in ['longitude', 'latitude']:

            xw_min_check = np.min(xw % 360.)
            xw_max_check = np.max(xw % 360.)

            if xw_max - xw_min < 360. and xw_max - xw_min >= xw_max_check - xw_min_check:
                xw_min = xw_min_check
                xw_max = xw_max_check

        # Check if range is smaller when normalizing to the range -180 to 180

        if coord_type in ['longitude', 'latitude']:

            xw_min_check = np.min(wrap_180(xw))
            xw_max_check = np.max(wrap_180(xw))

            if xw_max_check - xw_min_check < 360. and xw_max - xw_min >= xw_max_check - xw_min_check:
                xw_min = xw_min_check
                xw_max = xw_max_check

        x_range = xw_max - xw_min
        if coord_type == 'longitude':
            if x_range > 300.:
                xw_min = 0.
                xw_max = 360.
            elif xw_min < 0.:
                xw_min = max(-180., xw_min - 0.1 * x_range)
                xw_max = min(+180., xw_max + 0.1 * x_range)
            else:
                xw_min = max(0., xw_min - 0.1 * x_range)
                xw_max = min(360., xw_max + 0.1 * x_range)
        elif coord_type == 'latitude':
            xw_min = max(-90., xw_min - 0.1 * x_range)
            xw_max = min(+90., xw_max + 0.1 * x_range)

        ranges.append((xw_min, xw_max))

    return ranges
