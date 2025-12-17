# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import contextmanager
from .units import quantity_support
from .time import time_support

__all__ = ['astropy_support']

@contextmanager
def astropy_support():
    """
    Enable support for plotting both `astropy.units.Quantity` and
    `astropy.time.Time` instances in matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u, time
      >>> from astropy import visualization
      >>> with visualization.astropy_support():
      ...     plt.figure()
      ...     plt.plot(time.Time(['2016-03-22T12:30:31', '2016-03-22T12:30:38']), [1, 2] * u.m)
      ...     plt.draw()
    """
    with quantity_support(), time_support():
        yield
