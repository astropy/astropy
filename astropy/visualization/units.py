# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)


import numpy as np


def quantity_support():
    """
    Enable support for plotting `astropy.units.Quantity` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> visualization.quantity_support()
      <...>
      >>> plt.plot([1, 2, 3] * u.m)
      [...]
      >>> plt.plot([101, 125, 150] * u.cm)
      [...]
      >>> plt.draw()
    """
    from .. import units as u

    from matplotlib import units
    from matplotlib import ticker

    def rad_fn(x, pos=None):
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

    class MplQuantityConverter(units.ConversionInterface):
        def __init__(self):
            if u.Quantity not in units.registry[u.Quantity]:
                units.registry[u.Quantity] = mpl_quantity_converter
                self._remove = True
            else:
                self._remove = False

        @staticmethod
        def axisinfo(unit, axis):
            if unit == u.radian:
                return units.AxisInfo(
                    majloc=ticker.MultipleLocator(base=np.pi/2),
                    majfmt=ticker.FuncFormatter(rad_fn),
                    label=unit.to_string(),
                )
            elif unit == u.degree:
                return units.AxisInfo(
                    majloc=ticker.AutoLocator(),
                    majfmt=ticker.FormatStrFormatter('%i°'),
                    label=unit.to_string(),
                )
            elif unit is not None:
                return units.AxisInfo(label=unit.to_string('unicode'))
            return None

        @staticmethod
        def convert(val, unit, axis):
            if isinstance(val, u.Quantity):
                return val.to(unit)
            else:
                return val

        @staticmethod
        def default_units(x, axis):
            if hasattr(x, 'unit'):
                if x.unit.is_equivalent(u.radian):
                    return u.radian
                return x.unit
            return None

        def __enter__(self):
            return self

        def __exit__(self, type, value, tb):
            if self._remove:
                del units.registry[u.Quantity]

    return MplQuantityConverter()
