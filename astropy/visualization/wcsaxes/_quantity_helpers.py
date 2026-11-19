# Licensed under a 3-clause BSD style license - see LICENSE.rst

import functools
import inspect

from astropy import units as u
from astropy.coordinates import BaseCoordinateFrame
from astropy.wcs.wcsapi import BaseHighLevelWCS, BaseLowLevelWCS

from .frame import RectangularFrame1D

# Matplotlib discards the units of quantities before they reach the transform
# of an artist, so the plotting methods of WCSAxes convert quantities to the
# units that the transform expects before calling Matplotlib. For each of the
# x and y coordinates these units are either a unit, None if values are not
# coordinates (in which case quantities are passed on unchanged), or the
# following if the transform is not one for which the units are known.
UNKNOWN_UNIT = object()


def get_frame_input_units(ax, frame):
    # The units in which the transform returned by ax.get_transform(frame)
    # expects its input x and y values to be (see UNKNOWN_UNIT)

    if isinstance(frame, (BaseLowLevelWCS, BaseHighLevelWCS)) or (
        isinstance(frame, str) and frame == "pixel"
    ):
        # For a WCS, the values are pixel coordinates in that WCS
        return get_transform_input_units(ax, None)

    if isinstance(frame, str) and frame == "world":
        # The world units are not known if the axes were created with a
        # transform rather than a WCS, and there is no single unit for
        # each of x and y if there are more than two world coordinates
        units = getattr(ax._transform_pixel2world, "units_out", None)
        if units is not None and len(units) <= 2:
            # For a 1D WCS, y is not a world coordinate
            return (*(u.Unit(unit) for unit in units), None)[:2]

    elif isinstance(frame, (str, BaseCoordinateFrame)):
        # The name of a coordinate frame or a coordinate frame, for which the
        # values are the longitude and latitude in degrees
        return (u.deg, u.deg)

    return (UNKNOWN_UNIT, UNKNOWN_UNIT)


def get_transform_input_units(ax, transform):
    # The units in which the transform given to a plotting method of the axes
    # expects its input x and y values to be (see UNKNOWN_UNIT). These are
    # recorded by WCSAxes.get_transform, and are not known for any other
    # transform.

    if transform is None or transform is ax.transData:
        # Pixel coordinates, except for the y axis of 1D plots, which
        # shows data values
        if ax.frame_class is RectangularFrame1D:
            return (u.pix, None)
        return (u.pix, u.pix)

    return getattr(transform, "_wcsaxes_input_units", (UNKNOWN_UNIT, UNKNOWN_UNIT))


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
        if not isinstance(value, u.Quantity):
            value = u.Quantity(value)
        return value.to_value(unit)
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
        x_unit, y_unit = get_transform_input_units(self, kwargs.get("transform"))
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
