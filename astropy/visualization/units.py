# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import ContextDecorator

import numpy as np

__all__ = ["quantity_support"]
__doctest_skip__ = ["quantity_support"]


def quantity_support(format="latex_inline"):
    """
    Enable support for plotting `astropy.units.Quantity` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> with visualization.quantity_support():
      ...     fig, ax = plt.subplots()
      ...     ax.plot([1, 2, 3] * u.m)
      [...]
      ...     ax.plot([101, 125, 150] * u.cm)
      [...]
      ...     ax.yaxis.set_units(u.km)
      ...     plt.draw()

    The default axis unit is inferred from the first plot using a Quantity.
    To override it, you can explicitly set the axis unit using
    :meth:`matplotlib.axis.Axis.set_units`, for example,
    ``ax.yaxis.set_units(u.km)``.

    Parameters
    ----------
    format : `astropy.units.format.Base` subclass or str
        The name of a format or a formatter class.  If not
        provided, defaults to ``latex_inline``.

    """
    from matplotlib import ticker, units

    from astropy import units as u

    def rad_fn(x, pos=None):
        n = int((x / np.pi) * 2.0 + 0.25)
        if n == 0:
            return "0"
        elif n == 1:
            return "π/2"
        elif n == 2:
            return "π"
        elif n % 2 == 0:
            return f"{n // 2}π"
        else:
            return f"{n}π/2"

    class MplQuantityConverter(units.ConversionInterface, ContextDecorator):
        def __init__(self):
            # Keep track of original converter in case the context manager is
            # used in a nested way.
            self._original_converter = {u.Quantity: units.registry.get(u.Quantity)}
            units.registry[u.Quantity] = self

        @staticmethod
        def axisinfo(unit, axis):
            if unit == u.radian:
                return units.AxisInfo(
                    majloc=ticker.MultipleLocator(base=np.pi / 2),
                    majfmt=ticker.FuncFormatter(rad_fn),
                    label=unit.to_string(),
                )
            elif unit == u.degree:
                return units.AxisInfo(
                    majloc=ticker.AutoLocator(),
                    majfmt=ticker.FormatStrFormatter("%g°"),
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
            if hasattr(x, "unit"):
                return x.unit
            return None

        def __enter__(self):
            return self

        def __exit__(self, type, value, tb):
            if self._original_converter[u.Quantity] is None:
                del units.registry[u.Quantity]
            else:
                units.registry[u.Quantity] = self._original_converter[u.Quantity]

    return MplQuantityConverter()
from contextlib import contextmanager, ExitStack
from .units import quantity_support   # this is already in the same file
from .time import time_support        # import time_support from time.py

@contextmanager
def astro_support(quantity_format="latex_inline", time_scale=None, time_format=None, simplify=True):
    """
    Enable both Quantity and Time support for matplotlib in one context.

    This context manager wraps `quantity_support()` and `time_support()` so that
    matplotlib can plot `astropy.units.Quantity` and `astropy.time.Time` objects
    together seamlessly.

    Parameters
    ----------
    quantity_format : str, optional
        Format for Quantity axis labels (default is "latex_inline").
    time_scale : str, optional
        Time scale to use for Time objects (default None, uses first object's scale).
    time_format : str, optional
        Time format to use for Time objects (default None, uses first object's format).
    simplify : bool, optional
        Simplify Time labels if possible (default True).

    Usage
    -----
    >>> from astropy.visualization import astro_support
    >>> with astro_support():
    ...     # plot Quantity and Time objects
    """
    with ExitStack() as stack:
        stack.enter_context(quantity_support(format=quantity_format))
        stack.enter_context(time_support(scale=time_scale, format=time_format, simplify=simplify))
        yield