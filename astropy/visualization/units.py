# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import ContextDecorator

import numpy as np

from astropy import units as u
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB

__all__ = ["MplQuantityConverter", "quantity_support"]

__doctest_skip__ = ["quantity_support"]

_default_format = "latex_inline"

if HAS_MATPLOTLIB:
    from matplotlib import ticker, units

    class MplQuantityConverter(units.ConversionInterface, ContextDecorator):
        """Matplotlib converter for ``astropy.units.Quantity``.

        Registers itself to matplotlib as the converter for
        ``astropy.units.Quantity`` when initialized. If used as a context manager,
        it will restore the original converter upon exit. Also see
        :meth:`quantity_support` for a convenient way to use this converter
        with an optional format for the ``axisinfo``.
        """

        def __init__(self):
            # Keep track of original converter in case the context manager is
            # used in a nested way.
            self._original_converter = {u.Quantity: units.registry.get(u.Quantity)}
            units.registry[u.Quantity] = self

        @staticmethod
        def axisinfo(unit, axis, format=None):
            """Return a :class:`matplotlib.units.AxisInfo` for *unit* and *axis*.

            Parameters
            ----------
            unit : `~astropy.units.UnitBase`
                The unit to format the axis for.
            axis : `matplotlib.axis.Axis`
                The matplotlib axis being formatted.
            format : `astropy.units.format.Base` subclass or str or None, optional
                The name of a format or a formatter class used to render the
                axis label.  If `None`, the module-level default
                (``"latex_inline"``) is used.
            """
            if format is None:
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
            elif (
                all(
                    hasattr(val, attr)
                    for attr in ["__getitem__", "__iter__", "__len__"]
                )
                and len(val) > 0
                and isinstance(val[0], u.Quantity)
            ):
                return np.array([v.to_value(unit) for v in val])
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
                del units.registry[u.Quantity]
            else:
                units.registry[u.Quantity] = self._original_converter[u.Quantity]

else:
    # Create mock-up class to avoid import errors when matplotlib is not available.
    class MplQuantityConverter:
        def __init__(self, *args, **kwargs):
            raise ImportError("matplotlib is required in order to use this class.")


def quantity_support(format=None):
    """
    Enable support for plotting :class:`astropy.units.Quantity` instances in
    matplotlib.

    May be (optionally) used with a ``with`` statement.

      >>> import matplotlib.pyplot as plt
      >>> from astropy import units as u
      >>> from astropy import visualization
      >>> with visualization.quantity_support():
      ...     fig, ax = plt.subplots()
      ...     ax.plot([1, 2, 3] * u.m)
      ...  # doctest: +ELLIPSIS
      ...     ax.plot([101, 125, 150] * u.cm)
      ...  # doctest: +ELLIPSIS
      ...     ax.yaxis.set_units(u.km)
      ...     plt.draw()

    The default axis unit is inferred from the first plot using a Quantity.
    To override it, you can explicitly set the axis unit using
    :meth:`matplotlib.axis.Axis.set_units`, for example,
    ``ax.yaxis.set_units(u.km)``.

    Parameters
    ----------
    format : :class:`astropy.units.format.Base` subclass or str or None, optional
        The name of a format or a formatter class.  If not
        provided, defaults to ``latex_inline``.

    """

    class MplQuantityConverterFormatted(MplQuantityConverter):
        # Override to pass on `format`
        @staticmethod
        def axisinfo(unit, axis):
            return MplQuantityConverter.axisinfo(unit, axis, format)

    return MplQuantityConverterFormatted()
