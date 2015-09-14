# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import absolute_import, division, unicode_literals

import numpy as np
from ..stats import histogram
from ..utils.compat.funcsigs import signature

__all__ = ['hist']


def hist(x, bins=10, ax=None, **kwargs):
    """Enhanced histogram function

    This is a histogram function that enables the use of more sophisticated
    algorithms for determining bins.  Aside from the ``bins`` argument allowing
    a string specified how bins are computed, the parameters are the same
    as pylab.hist().

    This function was ported from astroML: http://astroML.org/

    Parameters
    ----------
    x : array_like
        array of data to be histogrammed

    bins : int or list or str (optional)
        If bins is a string, then it must be one of:

        - 'blocks' : use bayesian blocks for dynamic bin widths

        - 'knuth' : use Knuth's rule to determine bins

        - 'scott' : use Scott's rule to determine bins

        - 'freedman' : use the Freedman-diaconis rule to determine bins

    ax : Axes instance (optional)
        specify the Axes on which to draw the histogram.  If not specified,
        then the current active axes will be used.

    **kwargs :
        other keyword arguments are described in ``plt.hist()``.

    Notes
    -----
    Return values are the same as for ``plt.hist()``

    See Also
    --------
    astropy.stats.histogram
    """
    # arguments of np.histogram should be passed to astropy.stats.histogram
    arglist = list(signature(np.histogram).parameters.keys())[1:]
    np_hist_kwds = dict((key, kwargs[key]) for key in arglist if key in kwargs)
    hist, bins = histogram(x, bins, **np_hist_kwds)

    if ax is None:
        # optional dependency; only import if strictly needed.
        import matplotlib.pyplot as plt
        ax = plt.gca()

    return ax.hist(x, bins, **kwargs)
