import numpy as np

# Algorithm inspired by PGSBOX from WCSLIB by M. Calabretta


def wrap_180(values):
    values_new = values % 360.
    values_new[values_new > 180.] -= 360
    return values_new


def find_coordinate_range(transform, extent, x_type='scalar', y_type='scalar'):
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
    x_type : bool
        Whether the x coordinate is a ``'longitude'``, ``'latitude'``, or
        ``'scalar'`` value.
    y_type : bool
        Whether the y coordinate is a ``'longitude'``, ``'latitude'``, or
        ``'scalar'`` value.
    '''

    # Initialize the ranges
    wmin = np.repeat(np.inf, 2)
    wmax = np.repeat(-np.inf, 2)
    vmin = np.repeat(np.inf, 4).reshape(2, 2)
    vmax = np.repeat(-np.inf, 4).reshape(2, 2)

    # Sample coordinates on a 50 x 50 grid.
    NX = NY = 50
    x = np.linspace(extent[0], extent[1], NX + 1)
    y = np.linspace(extent[2], extent[3], NY + 1)
    xp, yp = np.meshgrid(x, y)
    world = transform.transform(np.vstack([xp.ravel(), yp.ravel()]).transpose())
    xw, yw = world[:, 0].reshape(xp.shape), world[:, 1].reshape(yp.shape)

    if x_type in ['longitude', 'latitude']:

        # Iron out coordinates along first row

        for ix in range(1, NX + 1):

            wjump = xw[0, ix] - xw[0, ix - 1]
            if np.abs(wjump) > 180.:
                wjump = wjump + np.sign(wjump) * 180.
                wjump = 360. * (wjump / 360.).astype(int)
                xw[0, ix] -= wjump

        # Now iron out coordinates along all columns, starting with first row.

        for iy in range(1, NY + 1):

            wjump = xw[iy,:] - xw[iy - 1,:]
            reset = np.abs(wjump) > 180.
            if np.any(reset):
                wjump = wjump + np.sign(wjump) * 180.
                wjump = 360. * (wjump / 360.).astype(int)
                xw[iy,:][reset] -= wjump[reset]

    if y_type in ['longitude', 'latitude']:

        # Iron out coordinates along first row

        wjump = yw[0, ix] - yw[0, ix - 1]
        if np.abs(wjump) > 180.:
            wjump = wjump + np.sign(wjump) * 180.
            wjump = 360. * (wjump / 360.).astype(int)
            yw[0, ix] -= wjump

        # Now iron out coordinates along all columns, starting with first row.

        for iy in range(1, NY + 1):

            wjump = yw[iy,:] - yw[iy - 1,:]
            reset = np.abs(wjump) > 180.
            if np.any(reset):
                wjump = wjump + np.sign(wjump) * 180.
                wjump = 360. * (wjump / 360.).astype(int)
                yw[iy,:][reset] -= wjump[reset]

    xw_min = np.nanmin(xw)
    xw_max = np.nanmax(xw)
    yw_min = np.nanmin(yw)
    yw_max = np.nanmax(yw)

    # Check if range is smaller when normalizing to the range 0 to 360

    if x_type in ['longitude', 'latitude']:

        xw_min_check = np.min(xw % 360.)
        xw_max_check = np.max(xw % 360.)

        if xw_max - xw_min < 360. and xw_max - xw_min >= xw_max_check - xw_min_check:
            xw_min = xw_min_check
            xw_max = xw_max_check

    if y_type in ['longitude', 'latitude']:

        yw_min_check = np.min(yw % 360.)
        yw_max_check = np.max(yw % 360.)

        if yw_max - yw_min < 360. and yw_max - yw_min >= yw_max_check - yw_min_check:
            yw_min = yw_min_check
            yw_max = yw_max_check

    # Check if range is smaller when normalizing to the range -180 to 180

    if x_type in ['longitude', 'latitude']:

        xw_min_check = np.min(wrap_180(xw))
        xw_max_check = np.max(wrap_180(xw))


        if xw_max_check - xw_min_check < 360. and xw_max - xw_min >= xw_max_check - xw_min_check:
            xw_min = xw_min_check
            xw_max = xw_max_check

    if y_type in ['longitude', 'latitude']:

        yw_min_check = np.min(wrap_180(yw))
        yw_max_check = np.max(wrap_180(yw))

        if yw_max_check - yw_min_check < 360. and yw_max - yw_min >= yw_max_check - yw_min_check:
            yw_min = yw_min_check
            yw_max = yw_max_check

    x_range = xw_max - xw_min
    if x_type == 'longitude':
        if x_range > 300.:
            xw_min = 0.
            xw_max = 360.
        elif xw_min < 0.:
            xw_min = max(-180., xw_min - 0.1 * x_range)
            xw_max = min(+180., xw_max + 0.1 * x_range)
        else:
            xw_min = max(0., xw_min - 0.1 * x_range)
            xw_max = min(360., xw_max + 0.1 * x_range)
    elif y_type == 'latitude':
        xw_min = max(-90., xw_min - 0.1 * x_range)
        xw_max = min(+90., xw_max + 0.1 * x_range)

    y_range = yw_max - yw_min
    if y_type == 'longitude':
        if y_range > 360.:
            yw_min = -180.
            yw_max = 180.
        elif yw_min < 0.:
            yw_min = max(-180., yw_min - 0.1 * y_range)
            yw_max = min(+180., yw_max + 0.1 * y_range)
        else:
            yw_min = max(0., yw_min - 0.1 * y_range)
            yw_max = min(360., yw_max + 0.1 * y_range)
    elif y_type == 'latitude':
        yw_min = max(-90., yw_min - 0.1 * y_range)
        yw_max = min(+90., yw_max + 0.1 * y_range)

    return (xw_min, xw_max), (yw_min, yw_max)
