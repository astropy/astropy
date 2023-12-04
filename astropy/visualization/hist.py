# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.stats.histogram import calculate_bin_edges

__all__ = ["hist"]


def hist(x, bins=10, ax=None, max_bins=1e5, **kwargs):
    """Enhanced histogram function.

    This is a histogram function that enables the use of more sophisticated
    algorithms for determining bins.  Aside from the ``bins`` argument allowing
    a string specified how bins are computed, the parameters are the same
    as matplotlib.pyplot.hist().

    This function was ported from astroML: https://www.astroml.org/

    Parameters
    ----------
    x : array-like
        array of data to be histogrammed

    bins : int, list, or str, optional
        If bins is a string, then it must be one of:

        - 'blocks' : use bayesian blocks for dynamic bin widths

        - 'knuth' : use Knuth's rule to determine bins

        - 'scott' : use Scott's rule to determine bins

        - 'freedman' : use the Freedman-Diaconis rule to determine bins

    ax : `~matplotlib.axes.Axes` instance, optional
        Specify the Axes on which to draw the histogram. If not specified,
        then the current active axes will be used.

    max_bins : int, optional
        Maximum number of bins allowed. With more than a few thousand bins
        the performance of matplotlib will not be great. If the number of
        bins is large *and* the number of input data points is large then
        the it will take a very long time to compute the histogram.

    **kwargs :
        other keyword arguments are described in ``plt.hist()``.

    Notes
    -----
    Return values are the same as for ``plt.hist()``

    See Also
    --------
    astropy.stats.histogram
    """
    # Note that we only calculate the bin edges...matplotlib will calculate
    # the actual histogram.
    range = kwargs.get("range", None)
    weights = kwargs.get("weights", None)
    bins = calculate_bin_edges(x, bins, range=range, weights=weights)

    if len(bins) > max_bins:
        raise ValueError(
            "Histogram has too many bins: "
            f"{len(bins)}. Use max_bins to increase the number "
            "of allowed bins or range to restrict "
            "the histogram range."
        )

    if ax is None:
        # optional dependency; only import if strictly needed.
        import matplotlib.pyplot as plt

        ax = plt.gca()

    return ax.hist(x, bins, **kwargs)
