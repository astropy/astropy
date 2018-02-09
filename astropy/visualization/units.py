# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from .. import units as u

try:
    from matplotlib.units import (ConversionInterface, AxisInfo,
                                  registry as units_registry)
except ImportError:  # no matplotlib
    class ConversionInterface(object):
        def __init__(self, *args, **kwargs):
            raise ImportError('matplotlib is required in order to use this '
                              'class.')


__doctest_skip__ = ['quantity_support']


def rad_fn(x, pos=None):
    """Represent a value in radians as a multiple of pi.
    """
    n = int((x / np.pi) * 2.0 + 0.25)
    if n == 0:
        return '0'
    elif n == 1:
        return 'π/2'
    elif n == 2:
        return 'π'
    elif n % 2 == 0:
        return '{0}π'.format(n / 2)
    else:
        return '{0}π/2'.format(n)


class MplQuantityConverter(ConversionInterface):
    """`ConversionInterface` with support for Astropy units
    """
    def __init__(self, format='latex_inline'):
        if u.Quantity not in units_registry:
            units_registry[u.Quantity] = self
            self._remove = True
        else:
            self._remove = False
        self.format = format

    def axisinfo(self, unit, axis):
        from matplotlib import ticker

        if unit == u.radian:
            return AxisInfo(
                majloc=ticker.MultipleLocator(base=np.pi/2),
                majfmt=ticker.FuncFormatter(rad_fn),
                label=unit.to_string(),
            )
        elif unit == u.degree:
            return AxisInfo(
                majloc=ticker.AutoLocator(),
                majfmt=ticker.FormatStrFormatter('%i°'),
                label=unit.to_string(),
            )
        elif unit is not None:
            return AxisInfo(label=unit.to_string(self.format))
        return None

    @staticmethod
    def convert(val, unit, axis):
        if isinstance(val, u.Quantity):
            return val.to_value(unit)
        else:
            return val

    @staticmethod
    def default_units(x, axis):
        if hasattr(x, 'unit'):
            return x.unit
        return None

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        if self._remove:
            del units_registry[u.Quantity]


def quantity_support(format='latex_inline'):
    """
    Enable support for plotting `astropy.units.Quantity` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> with visualization.quantity_support():
      ...     plt.figure()
      ...     plt.plot([1, 2, 3] * u.m)
      [...]
      ...     plt.plot([101, 125, 150] * u.cm)
      [...]
      ...     plt.draw()

    Parameters
    ----------
    format : `astropy.units.format.Base` instance or str
        The name of a format or a formatter object.  If not
        provided, defaults to ``latex_inline``.

    """
    return MplQuantityConverter()
