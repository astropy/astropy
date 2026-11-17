# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
import inspect

from astropy import units as u

# Matplotlib discards the units of quantities before they reach the transform
# of an artist, so the plotting methods of WCSAxes convert quantities to the
# units that the transform expects before calling Matplotlib. For each of the
# x and y coordinates these units are either a unit, None if values are not
# coordinates (in which case quantities are passed on unchanged), or the
# following if the transform is not one for which the units are known.
UNKNOWN_UNIT = object()


def to_value(value, unit):
    # Convert a Quantity, or a list or tuple that contains quantities, to plain
    # values in the given unit, and return anything else unchanged
    if unit is None or not (
        isinstance(value, u.Quantity)
        or (
            isinstance(value, (list, tuple))
            and any(isinstance(item, u.Quantity) for item in value)
        )
    ):
        return value
    if unit is UNKNOWN_UNIT:
        raise TypeError(
            "Quantities can only be plotted on WCSAxes without a transform, in "
            "which case they have to be in pixels, or with a transform returned "
            "by WCSAxes.get_transform, since the units that other transforms "
            "expect are not known"
        )
    try:
        return u.Quantity(value).to_value(unit)
    except u.UnitConversionError as exc:
        if unit == u.pix:
            expected = (
                "pixel coordinates are expected unless a transform is given, for "
                "example transform=ax.get_transform('world')"
            )
        else:
            expected = f"the transform expects units of '{unit}'"
        raise u.UnitConversionError(f"{exc}: {expected}") from None


def convert_plot_args(args, x_unit, y_unit):
    # Convert quantities in arguments given as ([x], y, [fmt]) groups
    args = list(args)
    converted = []
    while args:
        group = [args.pop(0)]
        if args and not isinstance(args[0], str):
            group.append(args.pop(0))
        units = (x_unit, y_unit)[-len(group) :]
        converted += [to_value(value, unit) for value, unit in zip(group, units)]
        if args and isinstance(args[0], str):
            converted.append(args.pop(0))
    return converted


def convert_grid_args(args, x_unit, y_unit):
    # Convert quantities in arguments given as ([X, Y], Z, ...), where X and Y
    # are optional and present if there are at least three arguments
    if len(args) < 3:
        return args
    return [to_value(args[0], x_unit), to_value(args[1], y_unit), *args[2:]]


def convert_quantities(method, x=(), y=(), convert_args=None):
    # Wrap a Matplotlib plotting method so that quantities are converted to
    # the units expected by the transform, for the x and y coordinates given
    # either by the names of arguments or, for methods with a variable number
    # of positional arguments, found by one of the functions above
    names = list(inspect.signature(method).parameters)[1:]

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        x_unit, y_unit = self._get_transform_input_units(kwargs.get("transform"))
        if convert_args is not None:
            args = convert_args(args, x_unit, y_unit)
        args = list(args)
        for name, unit in [(name, x_unit) for name in x] + [
            (name, y_unit) for name in y
        ]:
            if (index := names.index(name)) < len(args):
                args[index] = to_value(args[index], unit)
            elif name in kwargs:
                kwargs[name] = to_value(kwargs[name], unit)
        return method(self, *args, **kwargs)

    wrapper.__doc__ = f"""
        Call :meth:`matplotlib.axes.Axes.{method.__name__}`.

        This is the same as the Matplotlib method, except that coordinates
        given as `~astropy.units.Quantity` objects are converted to the units
        that the transform expects (see
        :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform`), or to
        pixels if no transform is used.
        """
    return wrapper
