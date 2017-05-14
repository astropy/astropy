# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, unicode_literals

import numpy as np

def scatter_contour(x, y,
                    levels=10,
                    threshold=100,
                    log_counts=False,
                    histogram2d_kwargs=None,
                    plot_kwargs=None,
                    contour_kwargs=None,
                    filled_contour=True,
                    ax=None):
    x = np.asarray(x)
    y = np.asarray(y)

    # plot contours on top
    default_contour_kwargs = dict(zorder=2)
    default_plot_kwargs = dict(marker='.', linestyle='none', zorder=1)

    if plot_kwargs is not None:
        default_plot_kwargs.update(plot_kwargs)
    plot_kwargs = default_plot_kwargs

    if contour_kwargs is not None:
        default_contour_kwargs.update(contour_kwargs)
    contour_kwargs = default_contour_kwargs

    if histogram2d_kwargs is None:
        histogram2d_kwargs = {}

    if contour_kwargs is None:
        contour_kwargs = {}

    if ax is None:
        # Import here so that testing with Agg will work
        from matplotlib import pyplot as plt
        ax = plt.gca()

    H, xbins, ybins = np.histogram2d(x, y, **histogram2d_kwargs)

    if log_counts:
        H = np.log10(1 + H)
        threshold = np.log10(1 + threshold)

    levels = np.asarray(levels)

    if levels.size == 1: # if an integer is passed

        if threshold > H.max():
            raise ValueError('Counts threshold is larger than the maximum '
                             'number of points per bin. Either decrease the '
                             'threshold or use larger bins (use the '
                             'histogram2d_kwargs argument to pass kwargs to '
                             'numpy.histogram2d).')

        levels = np.linspace(threshold, H.max(), levels)

    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]

    if filled_contour:
        contours = ax.contourf(H.T, levels, extent=extent, **contour_kwargs)
    else:
        contours = ax.contour(H.T, levels, extent=extent, **contour_kwargs)

    X = np.hstack([x[:, None], y[:, None]])

    # This grabs the outermost / lowest-level contour path vertices so we
    # can reduce the number of points to draw
    paths = contours.allsegs[0]
    if len(paths) > 0:
        if filled_contour:
            # We take only the first half because the second half is the
            # inner edge of the outer contour shape.
            n_vertices = paths[0].shape[0]
            outer_poly = paths[0][:n_vertices//2]
        else:
            outer_poly = paths[0]

        try:
            # this works in newer matplotlib versions
            from matplotlib.path import Path
            points_inside = Path(outer_poly).contains_points(X)
        except:
            # this works in older matplotlib versions
            import matplotlib.nxutils as nx
            points_inside = nx.points_inside_poly(X, outer_poly)

        Xplot = X[~points_inside]

    else:
        Xplot = X

    points = ax.plot(Xplot[:, 0], Xplot[:, 1], **plot_kwargs)

    return points, contours
