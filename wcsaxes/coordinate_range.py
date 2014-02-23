import numpy as np

# Algorithm inspired by PGSBOX from WCSLIB by M. Calabretta


def wrap_180(values):
    values_new = values % 360.
    values_new[values_new > 180.] -= 360
    return values_new


def find_coordinate_range(transform, extent, x_angle=False, y_angle=False):
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
    x_angle : bool
        Whether the x coordinate is an angle
    y_angle : bool
        Whether the y coordinate is an angle
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

    if x_angle:

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

    if y_angle:

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

    if x_angle:

        xw_min_check = np.min(xw % 360.)
        xw_max_check = np.max(xw % 360.)

        if xw_max - xw_min < 360. and xw_max - xw_min > xw_max_check - xw_min_check:
            if xw_max > 0.:
                xw_min = xw_min_check
                xw_max = xw_max_check
            else:
                xw_min = xw_min_check - 360.
                xw_max = xw_max_check - 360.

    if y_angle:

        yw_min_check = np.min(yw % 360.)
        yw_max_check = np.max(yw % 360.)

        if yw_max - yw_min < 360. and yw_max - yw_min > yw_max_check - yw_min_check:
            if yw_max > 0.:
                yw_min = yw_min_check
                yw_max = yw_max_check
            else:
                yw_min = yw_min_check - 360.
                yw_max = yw_max_check - 360.

    # Check if range is smaller when normalizing to the range -180 to 180

    if x_angle:

        xw_min_check = np.min(wrap_180(xw))
        xw_max_check = np.max(wrap_180(xw))

        if xw_max - xw_min < 360. and xw_max - xw_min > xw_max_check - xw_min_check:
            if xw_max > 0.:
                if xw_max_check > 0:
                    xw_min = xw_min_check
                    xw_max = xw_max_check
                else:
                    xw_min = xw_min_check + 360.
                    xw_max = xw_max_check + 360.
            else:
                if xw_max_check < 0:
                    xw_min = xw_min_check
                    xw_max = xw_max_check
                else:
                    xw_min = xw_min_check - 360.
                    xw_max = xw_max_check - 360.

    if y_angle:

        yw_min_check = np.min(wrap_180(yw))
        yw_max_check = np.max(wrap_180(yw))

        if yw_max - yw_min < 360. and yw_max - yw_min > yw_max_check - yw_min_check:
            if yw_max > 0.:
                if yw_max_check > 0:
                    yw_min = yw_min_check
                    yw_max = yw_max_check
                else:
                    yw_min = yw_min_check + 360.
                    yw_max = yw_max_check + 360.
            else:
                if yw_max_check < 0:
                    yw_min = yw_min_check
                    yw_max = yw_max_check
                else:
                    yw_min = yw_min_check - 360.
                    yw_max = yw_max_check - 360.

    return (xw_min, xw_max), (yw_min, yw_max)
