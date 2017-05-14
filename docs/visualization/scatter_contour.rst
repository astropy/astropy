.. _astropy-visualization-scatter-contour:

***********************************************
Scatter plots with large numbers of data points
***********************************************

This page demonstrates how to use the
:func:`~astropy.visualization.scatter_contour` function.

When creating a scatter plot with a large number of data points, dense regions
are often hard to interpret because of overlapping markers. For example:

.. plot::
    :align: center
    :context: close-figs
    :include-source: True

    import numpy as np
    import matplotlib.pyplot as plt

    figsize = (6, 6)
    lim = (-8, 8)

    np.random.seed(42)
    x, y = np.random.normal([0., 0.], [1., 1.], size=(16384, 2)).T

    fig,ax = plt.subplots(1, 1, figsize=figsize)
    ax.plot(x, y, marker='.', linestyle='none')
    ax.set_xlim(lim)
    ax.set_ylim(lim)

One solution that works for some situations is to use marker transparency.
However, when there is a large dynamic range between the most dense and most
sparse regions, even this can lead to confusion:

.. plot::
    :align: center
    :context: close-figs
    :include-source: True

    fig,ax = plt.subplots(1, 1, figsize=figsize)
    ax.plot(x, y, marker='.', linestyle='none', alpha=0.2)
    ax.set_xlim(lim)
    ax.set_ylim(lim)

The :mod:`astropy.visualization` module provides an alternate solution with the
function :func:`~astropy.visualization.scatter_contour`. The function switches
from drawing markers to drawing contours given some density threshold:

.. plot::
    :align: center
    :context: close-figs
    :include-source: True

    from astropy.visualization import scatter_contour
    fig,ax = plt.subplots(1, 1, figsize=figsize)
    scatter_contour(x, y, ax=ax)
    ax.set_xlim(lim)
    ax.set_ylim(lim)

This function can be customized to change the number and placement of the bins
used to compute the contours (using the argument ``histogram2d_kwargs``), the
marker plotting style (using the argument ``plot_kwargs``), and the contour
plotting style (using the arguments ``filled_contour`` and ``contour_kwargs``).
Here are some examples of other visualizations using
:func:`~astropy.visualization.scatter_contour` with the same data:

.. plot::
    :align: center
    :context: close-figs
    :include-source: True

    fig,ax = plt.subplots(1, 1, figsize=figsize)
    scatter_contour(x, y, ax=ax, contour_kwargs=dict(cmap='Greys_r'))
    ax.set_title('alternate colormap')
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    fig,ax = plt.subplots(1, 1, figsize=figsize)
    scatter_contour(x, y, ax=ax, threshold=10,
                    log_counts=True, levels=8,
                    histogram2d_kwargs=dict(bins=(16,16)),
                    filled_contour=False,
                    contour_kwargs=dict(colors='k'))
    ax.set_title('log-spaced contour lines')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
