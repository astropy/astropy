# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ["MplQuantityConverter", "quantity_support"]

import types
from contextlib import ContextDecorator

import numpy as np

from astropy import units as u
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB

if HAS_MATPLOTLIB:
    from matplotlib.ticker import (
        AutoLocator,
        FormatStrFormatter,
        FuncFormatter,
        MultipleLocator,
    )
    from matplotlib.units import AxisInfo, ConversionInterface, registry
else:  # Create mock-up classes to avoid import errors when matplotlib is not available.

    class MplMockUp:
        def __init__(self, *args, **kwargs):
            raise ImportError("matplotlib is required in order to use this class.")

    ConversionInterface = types.new_class("ConversionInterface", (MplMockUp,))
    AxisInfo = types.new_class("AxisInfo", (MplMockUp,))
    MultipleLocator = types.new_class("MultipleLocator", (MplMockUp,))
    FuncFormatter = types.new_class("FuncFormatter", (MplMockUp,))
    AutoLocator = types.new_class("AutoLocator", (MplMockUp,))
    FormatStrFormatter = types.new_class("FormatStrFormatter", (MplMockUp,))

    registry = {}


__doctest_skip__ = ["quantity_support"]

_default_format = "latex_inline"


def quantity_support(format=None):
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
    format : `astropy.units.format.Base` subclass or str or None, optional
        The name of a format or a formatter class.  If not
        provided, defaults to ``latex_inline``.

    """

    class MplQuantityConverterFormatted(MplQuantityConverter):
        @staticmethod
        def axisinfo(unit, axis):
            return MplQuantityConverter.axisinfo(unit, axis, format)

    return MplQuantityConverterFormatted()


class MplQuantityConverter(ConversionInterface, ContextDecorator):
    """Matplotlib converter for ``astropy.units.Quantity``.

    Registers itself to matplotlib as the converter for
    `astropy.units.Quantity` when initialized. If used as a context manager,
    it will restore the original converter upon exit. Also see
    ``quantity_support`` for a convenient way to use this converter
    with an optional format for the ``axisinfo``.
    """

    def __init__(self):
        # Keep track of original converter in case the context manager is
        # used in a nested way.
        self._original_converter = {u.Quantity: registry.get(u.Quantity)}
        registry[u.Quantity] = self

    @staticmethod
    def axisinfo(unit, axis, format=None):
        if not format:
            format = _default_format
        if unit == u.radian:

            def rad_fn(x, pos=None):
                n = int((x / np.pi) * 2.0 + 0.25)
                if n < 3:
                    return ("0", "π/2", "π")[n]
                elif n % 2 == 0:
                    return f"{n // 2}π"
                else:
                    return f"{n}π/2"

            return AxisInfo(
                majloc=MultipleLocator(base=np.pi / 2),
                majfmt=FuncFormatter(rad_fn),
                label=unit.to_string(),
            )
        elif unit == u.degree:
            return AxisInfo(
                majloc=AutoLocator(),
                majfmt=FormatStrFormatter("%g°"),
                label=unit.to_string(),
            )
        elif unit is not None:
            return AxisInfo(label=unit.to_string(format))
        return None

    @staticmethod
    def convert(val, unit, axis):
        if isinstance(val, u.Quantity):
            return val.to_value(unit)
        elif (
            all(hasattr(val, attr) for attr in ["__getitem__", "__iter__", "__len__"])
            and len(val) > 0
            and isinstance(val[0], u.Quantity)
        ):
            return [v.to_value(unit) for v in val]
        else:
            return val

    @staticmethod
    def default_units(x, axis):
        if hasattr(x, "unit"):
            return x.unit
        elif (
            all(hasattr(x, attr) for attr in ["__getitem__", "__iter__", "__len__"])
            and len(x) > 0
            and hasattr(x[0], "unit")
        ):
            return x[0].unit
        return None

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        if self._original_converter[u.Quantity] is None:
            del registry[u.Quantity]
        else:
            registry[u.Quantity] = self._original_converter[u.Quantity]
