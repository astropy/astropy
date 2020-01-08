# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np


__doctest_skip__ = ['quantity_support']


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
    from astropy import units as u
    # import Angle just so we have a more or less complete list of Quantity
    # subclasses loaded - matplotlib needs them all separately!
    # NOTE: in matplotlib >=3.2, subclasses will be recognized automatically,
    # and once that becomes our minimum version, we can remove this,
    # adding just u.Quantity itself to the registry.
    from astropy.coordinates import Angle  # noqa

    from matplotlib import units
    from matplotlib import ticker

    # Get all subclass for Quantity, since matplotlib checks on class,
    # not subclass.
    def all_issubclass(cls):
        return {cls}.union(
            [s for c in cls.__subclasses__() for s in all_issubclass(c)])

    def rad_fn(x, pos=None):
        n = int((x / np.pi) * 2.0 + 0.25)
        if n == 0:
            return '0'
        elif n == 1:
            return 'π/2'
        elif n == 2:
            return 'π'
        elif n % 2 == 0:
            return '{}π'.format(n / 2)
        else:
            return f'{n}π/2'

    class MplQuantityConverter(units.ConversionInterface):

        _all_issubclass_quantity = all_issubclass(u.Quantity)

        def __init__(self):

            # Keep track of original converter in case the context manager is
            # used in a nested way.
            self._original_converter = {}

            for cls in self._all_issubclass_quantity:
                self._original_converter[cls] = units.registry.get(cls)
                units.registry[cls] = self

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
                return units.AxisInfo(label=unit.to_string(format))
            return None

        @staticmethod
        def convert(val, unit, axis):
            if isinstance(val, u.Quantity):
                return val.to_value(unit)
            elif isinstance(val, list) and val and isinstance(val[0], u.Quantity):
                return [v.to_value(unit) for v in val]
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
            for cls in self._all_issubclass_quantity:
                if self._original_converter[cls] is None:
                    del units.registry[cls]
                else:
                    units.registry[cls] = self._original_converter[cls]

    return MplQuantityConverter()
