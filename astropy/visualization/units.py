# Licensed under a 3-clause BSD style license - see LICENSE.rst

from contextlib import ContextDecorator

import numpy as np

__all__ = ["quantity_support"]
__doctest_skip__ = ["quantity_support"]


def quantity_support(*args,format="latex_inline"):
    outer_args = args
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
        axis_labels_val = outer_args
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
                    label=MplQuantityConverter.decide_label(axis,unit),
                )
            elif unit == u.degree:
                return units.AxisInfo(
                    majloc=ticker.AutoLocator(),
                    majfmt=ticker.FormatStrFormatter("%g°"),
                    label=MplQuantityConverter.decide_label(axis,unit),
                )
            elif unit is not None:
                return units.AxisInfo(label=MplQuantityConverter.decide_label(axis,unit))
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
        
        @staticmethod
        def decide_label(axis, unit):
            axis_labels = MplQuantityConverter.axis_labels_val 
            
            default_label = f"{unit.physical_type}({unit.to_string()})"
            exceptional_units = {u.radian, u.degree}
            axis_name_ent = axis.axes
            
            has_label_0 = len(axis_labels) > 0 and bool(axis_labels[0])
            has_label_1 = len(axis_labels) > 1 and bool(axis_labels[1])
            
            # Get information about which axis we're dealing with
            is_x_axis = isinstance(axis, type(axis_name_ent.xaxis))
            is_y_axis = isinstance(axis, type(axis_name_ent.yaxis))
            
            # Handle exceptional units (radians, degrees)
            if unit in exceptional_units:
                if is_x_axis and has_label_0:
                    return f"{axis_labels[0]}({unit.to_string()})"
                elif is_y_axis and has_label_1:
                    return f"{axis_labels[1]}({unit.to_string()})"
                else:
                    return default_label
                
            # Handle normal units
            else:
                if is_x_axis and has_label_0:
                    return f"{axis_labels[0]}({unit.to_string()})"
                elif is_y_axis and has_label_1:
                    return f"{axis_labels[1]}({unit.to_string()})"
                else:
                    return default_label
                
        def __exit__(self, type, value, tb):
            if self._original_converter[u.Quantity] is None:
                del units.registry[u.Quantity]
            else:
                units.registry[u.Quantity] = self._original_converter[u.Quantity]

    return MplQuantityConverter()
