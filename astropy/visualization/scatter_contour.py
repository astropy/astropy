# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, unicode_literals

import numpy as np

__all__ = ['scatter_contour']

def scatter_contour(x, y,
                    levels=10,
                    threshold=100,
                    log_counts=False,
                    histogram2d_kwargs=None,
                    plot_kwargs=None,
                    contour_kwargs=None,
                    filled_contour=True,
                    ax=None):
    """Scatter plot with contours over dense regions

    The contours are either drawn with :meth:`~matplotlib.axes.Axes.contour` or
    :meth:`~matplotlib.axes.Axes.contourf`. The points in sparse regions are
    plotted using :meth:`~matplotlib.axes.Axes.plot`.

    This function was ported from astroML: http://astroML.org/

    Original code:
    Copyright (c) 2012-2013, Jacob Vanderplas. All rights reserved.

    Parameters
    ----------
    x : array_like
        x data for the contour plot.
    y : array_like
        y data for the contour plot.
    levels : integer or array_like (optional, default=10)
        Number of contour levels, or array of contour levels.
    threshold : float (optional)
        Number of points per 2D bin at which to begin drawing contours.
    log_counts : boolean (optional)
        If ``True``, contour levels are the base-10 logarithm of the bin counts.
    histogram2d_kwargs : dict
        Keyword arguments passed to :func:`numpy.histogram2d`, which is used to
        generate the input for the contour plotting function used. See docstring
        of :func:`numpy.histogram2d` for more information.
    plot_kwargs : dict
        Keyword arguments passed to :meth:`~matplotlib.axes.Axes.plot`. By
        default it will use ``marker='.'`` and ``linestyle='none'``. See
        docstring of :meth:`matplotlib.axes.Axes.plot` for more information.
    contour_kwargs : dict
        Keyword arguments passed to :meth:`~matplotlib.axes.Axes.contour` or
        :meth:`matplotlib.axes.Axes.contourf`, depending on the value of
        ``filled_contour``. See docstrings of these functions for more
        information.
    filled_contour : bool
        If ``True`` (default) use filled contours,
        :meth:`matplotlib.axes.Axes.contourf`. Otherwise, use contour outlines,
        :meth:`matplotlib.axes.Axes.contour`.
    ax : `matplotlib.axes.Axes`
        The axes on which to plot. If not specified, the current axes are used,
        from :func:`~matplotlib.pyplot.gca()`.

    Returns
    -------
    points, contours :
       ``points`` is the return value of :meth:`~matplotlib.axes.Axes.plot`.
       ``contours`` is the return value of the contour function that is called.

    """
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
