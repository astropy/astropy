# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
from collections.abc import Sequence

import numpy as np

from astropy import units as u
from astropy.units import IrreducibleUnit, Unit
from astropy.utils.compat import COPY_IF_NEEDED

from .baseframe import (
    BaseCoordinateFrame,
    _get_diff_cls,
    _get_repr_cls,
    frame_transform_graph,
)
from .representation import (
    BaseRepresentation,
    SphericalRepresentation,
    UnitSphericalRepresentation,
)

"""
This module contains utility functions to make the SkyCoord initializer more modular
and maintainable. No functionality here should be in the public API, but rather used as
part of creating SkyCoord objects.
"""

PLUS_MINUS_RE = re.compile(r"(\+|\-)")
J_PREFIXED_RA_DEC_RE = re.compile(
    r"""J                              # J prefix
    ([0-9]{6,7}\.?[0-9]{0,2})          # RA as HHMMSS.ss or DDDMMSS.ss, optional decimal digits
    ([\+\-][0-9]{6}\.?[0-9]{0,2})\s*$  # Dec as DDMMSS.ss, optional decimal digits
    """,
    re.VERBOSE,
)


def _get_frame_class(frame):
    """
    Get a frame class from the input `frame`, which could be a frame name
    string, or frame class.
    """
    if isinstance(frame, str):
        frame_names = frame_transform_graph.get_names()
        if frame not in frame_names:
            raise ValueError(
                f'Coordinate frame name "{frame}" is not a known '
                f"coordinate frame ({sorted(frame_names)})"
            )
        frame_cls = frame_transform_graph.lookup_name(frame)

    elif isinstance(frame, type) and issubclass(frame, BaseCoordinateFrame):
        frame_cls = frame

    else:
        raise ValueError(
            "Coordinate frame must be a frame name or frame class, not a"
            f" '{frame.__class__.__name__}'"
        )

    return frame_cls


_conflict_err_msg = (
    "Coordinate attribute '{0}'={1!r} conflicts with keyword argument '{0}'={2!r}. This"
    " usually means an attribute was set on one of the input objects and also in the "
    "keyword arguments to {3}"
)


def _get_frame_without_data(args, kwargs):
    """
    Determines the coordinate frame from input SkyCoord args and kwargs.

    This function extracts (removes) all frame attributes from the kwargs and
    determines the frame class either using the kwargs, or using the first
    element in the args (if a single frame object is passed in, for example).
    This function allows a frame to be specified as a string like 'icrs' or a
    frame class like ICRS, or an instance ICRS(), as long as the instance frame
    attributes don't conflict with kwargs passed in (which could require a
    three-way merge with the coordinate data possibly specified via the args).
    """
    from .sky_coordinate import SkyCoord

    # We eventually (hopefully) fill and return these by extracting the frame
    # and frame attributes from the input:
    frame_cls = None
    frame_cls_kwargs = {}

    # The first place to check: the frame could be specified explicitly
    frame = kwargs.pop("frame", None)

    if isinstance(frame, SkyCoord):
        # If the frame was passed as a SkyCoord, we also want to preserve
        # any extra attributes (e.g., obstime) if they are not already
        # specified in the kwargs. We preserve these extra attributes by
        # adding them to the kwargs dict:
        for attr in frame._extra_frameattr_names:
            value = getattr(frame, attr)
            if attr not in kwargs:
                kwargs[attr] = value
            elif np.any(value != kwargs[attr]):
                # This SkyCoord attribute passed in with the frame= object
                # conflicts with an attribute passed in directly to the
                # SkyCoord initializer as a kwarg:
                raise ValueError(
                    _conflict_err_msg.format(attr, value, kwargs[attr], "SkyCoord")
                )
        frame = frame.frame

    if isinstance(frame, BaseCoordinateFrame):
        # Extract any frame attributes
        for attr in frame.frame_attributes:
            # If the frame was specified as an instance, we have to make
            # sure that no frame attributes were specified as kwargs - this
            # would require a potential three-way merge:
            if attr in kwargs:
                raise ValueError(
                    f"Cannot specify frame attribute '{attr}' directly as an"
                    " argument to SkyCoord because a frame instance was passed in."
                    " Either pass a frame class, or modify the frame attributes of"
                    " the input frame instance."
                )
            if not frame.is_frame_attr_default(attr):
                kwargs[attr] = getattr(frame, attr)

        frame_cls = frame.__class__

        # Make sure we propagate representation/differential _type choices,
        # unless these are specified directly in the kwargs:
        kwargs.setdefault("representation_type", frame.representation_type)
        kwargs.setdefault("differential_type", frame.differential_type)

    elif frame_cls is None and frame is not None:  # frame probably a string
        frame_cls = _get_frame_class(frame)

    # Check that the new frame doesn't conflict with existing coordinate frame
    # if a coordinate is supplied in the args list.  If the frame still had not
    # been set by this point and a coordinate was supplied, then use that frame.
    for arg in args:
        # this catches the "single list passed in" case.  For that case we want
        # to allow the first argument to set the class.  That's OK because
        # _parse_coordinate_arg goes and checks that the frames match between
        # the first and all the others
        if isinstance(arg, (Sequence, np.ndarray)) and len(args) == 1 and len(arg) > 0:
            arg = arg[0]

        if isinstance(arg, BaseCoordinateFrame):
            coord_frame_obj = arg
        elif isinstance(arg, SkyCoord):
            coord_frame_obj = arg.frame
        else:
            continue

        coord_frame_cls = type(coord_frame_obj)
        if frame_cls is None:
            frame_cls = coord_frame_cls
        elif frame_cls is not coord_frame_cls:
            raise ValueError(
                f"Cannot override frame='{coord_frame_cls.__name__}' of input "
                f"coordinate with new frame='{frame_cls.__name__}'. Instead, "
                "transform the coordinate."
            )

        for attr in coord_frame_obj.frame_attributes:
            if not coord_frame_obj.is_frame_attr_default(attr):
                value = getattr(coord_frame_obj, attr)
                if attr not in kwargs:
                    kwargs[attr] = value
                elif np.any(value != kwargs[attr]):
                    raise ValueError(
                        f"Frame attribute '{attr}' has conflicting values between the"
                        " input coordinate data and either keyword arguments or the "
                        f"frame specification (frame=...): {value} =/= {kwargs[attr]}"
                    )

    if frame_cls is None:
        from .builtin_frames import ICRS

        frame_cls = ICRS

    # By now, frame_cls should be set - if it's not, something went wrong
    if not issubclass(frame_cls, BaseCoordinateFrame):
        # We should hopefully never get here...
        raise ValueError(f"Frame class has unexpected type: {frame_cls.__name__}")

    for attr in frame_cls.frame_attributes:
        if attr in kwargs:
            frame_cls_kwargs[attr] = kwargs.pop(attr)

    if (representation_type := kwargs.pop("representation_type", None)) is not None:
        frame_cls_kwargs["representation_type"] = _get_repr_cls(representation_type)

    if (differential_type := kwargs.pop("differential_type", None)) is not None:
        frame_cls_kwargs["differential_type"] = _get_diff_cls(differential_type)

    return frame_cls, frame_cls_kwargs


def _parse_coordinate_data(frame, args, kwargs):
    """
    Extract coordinate data from the args and kwargs passed to SkyCoord.

    By this point, we assume that all of the frame attributes have been
    extracted from kwargs (see _get_frame_without_data()), so all that are left
    are (1) extra SkyCoord attributes, and (2) the coordinate data, specified in
    any of the valid ways.
    """
    valid_skycoord_kwargs = {}
    info = None

    # Look through the remaining kwargs to see if any are valid attribute names
    # by asking the frame transform graph:
    for attr in tuple(kwargs):
        if attr in frame_transform_graph.frame_attributes:
            valid_skycoord_kwargs[attr] = kwargs.pop(attr)

    # By this point in parsing the arguments, anything left in the args and
    # kwargs should be data. Either as individual components, or a list of
    # objects, or a representation, etc.

    # Get units of components
    units = _get_representation_component_units(args, kwargs)

    # Grab any frame-specific attr names like `ra` or `l` or `distance` from
    # kwargs and move them to valid_components.
    valid_components = _get_representation_attrs(frame, units, kwargs)

    # Error if anything is still left in kwargs
    if kwargs:
        # The next few lines add a more user-friendly error message to a
        # common and confusing situation when the user specifies, e.g.,
        # `pm_ra` when they really should be passing `pm_ra_cosdec`. The
        # extra error should only turn on when the positional representation
        # is spherical, and when the component 'pm_<lon>' is passed.
        pm_message = ""
        if frame.representation_type == SphericalRepresentation:
            lon_name, lat_name, _ = frame.get_representation_component_names()
            if f"pm_{lon_name}" in kwargs:
                pm_message = (
                    "\n\n By default, most frame classes expect the longitudinal proper"
                    " motion to include the cos(latitude) term, named"
                    f" `pm_{lon_name}_cos{lat_name}`. Did you mean to pass in this"
                    " component?"
                )

        raise ValueError(
            "Unrecognized keyword argument(s) {}{}".format(
                ", ".join(f"'{key}'" for key in kwargs), pm_message
            )
        )

    if len(args) > 3:
        raise ValueError(
            f"Must supply no more than three positional arguments, got {len(args)}"
        )

    # Finally deal with the unnamed args.  This figures out what the arg[0]
    # is and returns a dict with appropriate key/values for initializing
    # frame class. Note that differentials are *never* valid args, only
    # kwargs.  So they are not accounted for here (unless they're in a frame
    # or SkyCoord object)
    _skycoord_kwargs = {}
    _components = {}
    if len(args) > 1:
        for arg, (frame_attr_name, repr_attr_name), unit in zip(
            args, frame.representation_component_names.items(), units
        ):
            attr_class = frame.representation_type.attr_classes[repr_attr_name]
            _components[frame_attr_name] = attr_class(arg, unit=unit)

    elif len(args) == 1:
        # One arg which must be a coordinate.  In this case coord_kwargs
        # will contain keys like 'ra', 'dec', 'distance' along with any
        # frame attributes like equinox or obstime which were explicitly
        # specified in the coordinate object (i.e. non-default).
        _skycoord_kwargs, _components = _parse_coordinate_arg(args[0], frame, units)

        # Copy other 'info' attr only if it has actually been defined.
        if "info" in getattr(args[0], "__dict__", ()):
            info = args[0].info

    # The next two loops copy the component and skycoord attribute data into
    # their final, respective "valid_" dictionaries. For each, we check that
    # there are no relevant conflicts with values specified by the user
    # through other means:

    # First validate the component data
    for attr, coord_value in _components.items():
        if attr in valid_components:
            raise ValueError(
                _conflict_err_msg.format(
                    attr, coord_value, valid_components[attr], "SkyCoord"
                )
            )
        valid_components[attr] = coord_value

    # Now validate the custom SkyCoord attributes
    for attr, value in _skycoord_kwargs.items():
        if attr in valid_skycoord_kwargs and np.any(
            valid_skycoord_kwargs[attr] != value
        ):
            raise ValueError(
                _conflict_err_msg.format(
                    attr, value, valid_skycoord_kwargs[attr], "SkyCoord"
                )
            )
        valid_skycoord_kwargs[attr] = value

    return valid_skycoord_kwargs, valid_components, info


def _get_representation_component_units(args, kwargs):
    """
    Get the unit from kwargs for the *representation* components (not the
    differentials).
    """
    if "unit" not in kwargs:
        units = [None, None, None]

    else:
        units = kwargs.pop("unit")

        if isinstance(units, str):
            units = [x.strip() for x in units.split(",")]
            # Allow for input like unit='deg' or unit='m'
            if len(units) == 1:
                units = [units[0], units[0], units[0]]
        elif isinstance(units, (Unit, IrreducibleUnit)):
            units = [units, units, units]

        try:
            units = [(Unit(x) if x else None) for x in units]
            units.extend(None for x in range(3 - len(units)))
            if len(units) > 3:
                raise ValueError()
        except Exception as err:
            raise ValueError(
                "Unit keyword must have one to three unit values as "
                "tuple or comma-separated string."
            ) from err

    return units


def _parse_coordinate_arg(coords, frame, units):
    from .sky_coordinate import SkyCoord

    is_scalar = False  # Differentiate between scalar and list input
    # valid_kwargs = {}  # Returned dict of lon, lat, and distance (optional)
    components = {}
    skycoord_kwargs = {}

    frame_attr_names = list(frame.representation_component_names.keys())
    repr_attr_names = list(frame.representation_component_names.values())
    repr_attr_classes = list(frame.representation_type.attr_classes.values())
    n_attr_names = len(repr_attr_names)

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, str):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `frame` it is checked that any coordinate
        # args have the same frame as explicitly supplied, so don't worry here.

        if not coords.has_data:
            raise ValueError("Cannot initialize from a frame without coordinate data")

        data = coords.data.represent_as(frame.representation_type)

        # If coords did not have an explicit distance then don't include in initializers.
        if isinstance(coords.data, UnitSphericalRepresentation):
            try:
                index = repr_attr_names.index("distance")
            except ValueError:
                pass
            else:
                del repr_attr_names[index]
                del units[index]
                del frame_attr_names[index]
                del repr_attr_classes[index]

        # List of values corresponding to representation attrs
        values = [getattr(data, name) for name in repr_attr_names]

        if coords.data.differentials and "s" in coords.data.differentials:
            orig_vel = coords.data.differentials["s"]
            vel = coords.data.represent_as(
                frame.representation_type, frame.get_representation_cls("s")
            ).differentials["s"]
            for frname, reprname in frame.get_representation_component_names(
                "s"
            ).items():
                if (
                    reprname == "d_distance"
                    and not hasattr(orig_vel, reprname)
                    and "unit" in orig_vel.name
                ):
                    continue
                values.append(getattr(vel, reprname))
                units.append(None)
                frame_attr_names.append(frname)
                repr_attr_names.append(reprname)
                repr_attr_classes.append(vel.attr_classes[reprname])

        is_skycoord = isinstance(coords, SkyCoord)
        for attr in frame_transform_graph.frame_attributes:
            if (value := getattr(coords, attr, None)) is not None and (
                is_skycoord or attr not in coords.frame_attributes
            ):
                skycoord_kwargs[attr] = value

    elif isinstance(coords, BaseRepresentation):
        if coords.differentials and "s" in coords.differentials:
            diffs = frame.get_representation_cls("s")
            data = coords.represent_as(frame.representation_type, diffs)
            values = [getattr(data, name) for name in repr_attr_names]
            for frname, reprname in frame.get_representation_component_names(
                "s"
            ).items():
                values.append(getattr(data.differentials["s"], reprname))
                units.append(None)
                frame_attr_names.append(frname)
                repr_attr_names.append(reprname)
                repr_attr_classes.append(data.differentials["s"].attr_classes[reprname])

        else:
            data = coords.represent_as(frame.representation_type)
            values = [getattr(data, name) for name in repr_attr_names]

    elif (
        isinstance(coords, np.ndarray)
        and coords.dtype.kind in "if"
        and coords.ndim == 2
        and coords.shape[1] <= 3
    ):
        # 2-d array of coordinate values.  Handle specially for efficiency.
        values = coords.transpose()  # Iterates over repr attrs

    elif isinstance(coords, (Sequence, np.ndarray)):
        # Handles list-like input.
        coord_types = (SkyCoord, BaseCoordinateFrame, BaseRepresentation)
        if any(isinstance(coord, coord_types) for coord in coords):
            # this parsing path is used when there are coordinate-like objects
            # in the list - instead of creating lists of values, we create
            # SkyCoords from the list elements and then combine them.
            scs = list(map(SkyCoord, coords))

            # Check that all frames are equivalent
            for sc in scs[1:]:
                if not sc.is_equivalent_frame(scs[0]):
                    raise ValueError(
                        f"List of inputs don't have equivalent frames: {sc} != {scs[0]}"
                    )

            # Now use the first to determine if they are all UnitSpherical
            not_unit_sphere = not isinstance(scs[0].data, UnitSphericalRepresentation)

            # get the frame attributes from the first coord in the list, because
            # from the above we know it matches all the others.  First copy over
            # the attributes that are in the frame itself, then copy over any
            # extras in the SkyCoord
            for fattrnm in scs[0].frame.frame_attributes:
                skycoord_kwargs[fattrnm] = getattr(scs[0].frame, fattrnm)
            for fattrnm in scs[0]._extra_frameattr_names:
                skycoord_kwargs[fattrnm] = getattr(scs[0], fattrnm)

            # Now combine the values, to be used below
            values = [
                np.concatenate([np.atleast_1d(getattr(sc, data_attr)) for sc in scs])
                for data_attr, repr_attr in zip(frame_attr_names, repr_attr_names)
                if not_unit_sphere or repr_attr != "distance"
            ]
        else:
            is_radec = (
                "ra" in frame.representation_component_names
                and "dec" in frame.representation_component_names
            )
            # none of the elements are "frame-like", create a list of sequences like
            # [[v1_0, v2_0, v3_0], ... [v1_N, v2_N, v3_N]]
            vals = [
                _parse_one_coord_str(c, is_radec=is_radec) if isinstance(c, str) else c
                for c in coords
            ]

            # Do some basic validation of the list elements: all have a length and all
            # lengths the same
            try:
                n_coords = {len(x) for x in vals}
            except Exception as err:
                raise ValueError(
                    "One or more elements of input sequence does not have a length."
                ) from err

            if len(n_coords) > 1:
                raise ValueError(
                    "Input coordinate values must have same number of elements, found"
                    f" {sorted(n_coords)}"
                )
            n_coords = n_coords.pop()

            # Must have no more coord inputs than representation attributes
            if n_coords > n_attr_names:
                raise ValueError(
                    f"Input coordinates have {n_coords} values but representation"
                    f" {frame.representation_type.name} only accepts"
                    f" {n_attr_names}"
                )

            # Now transpose vals to get [(v1_0 .. v1_N), (v2_0 .. v2_N), (v3_0 .. v3_N)]
            # (ok since we know it is exactly rectangular).  (Note: can't just use zip(*values)
            # because Longitude et al distinguishes list from tuple so [a1, a2, ..] is needed
            # while (a1, a2, ..) doesn't work.
            values = [list(x) for x in zip(*vals)]

            if is_scalar:
                values = [x[0] for x in values]
    else:
        raise ValueError("Cannot parse coordinates from first argument")

    # Finally we have a list of values from which to create the keyword args
    # for the frame initialization.  Validate by running through the appropriate
    # class initializer and supply units (which might be None).
    try:
        for frame_attr_name, repr_attr_class, value, unit in zip(
            frame_attr_names, repr_attr_classes, values, units
        ):
            components[frame_attr_name] = repr_attr_class(
                value, unit=unit, copy=COPY_IF_NEEDED
            )
    except Exception as err:
        raise ValueError(
            f'Cannot parse first argument data "{value}" for attribute'
            f" {frame_attr_name}"
        ) from err
    return skycoord_kwargs, components


def _get_representation_attrs(frame, units, kwargs):
    """
    Find instances of the "representation attributes" for specifying data
    for this frame.  Pop them off of kwargs, run through the appropriate class
    constructor (to validate and apply unit), and put into the output
    valid_kwargs.  "Representation attributes" are the frame-specific aliases
    for the underlying data values in the representation, e.g. "ra" for "lon"
    for many equatorial spherical representations, or "w" for "x" in the
    cartesian representation of Galactic.

    This also gets any *differential* kwargs, because they go into the same
    frame initializer later on.
    """
    frame_attr_names = frame.representation_component_names.keys()
    repr_attr_classes = frame.representation_type.attr_classes.values()

    valid_kwargs = {}
    for frame_attr_name, repr_attr_class, unit in zip(
        frame_attr_names, repr_attr_classes, units
    ):
        value = kwargs.pop(frame_attr_name, None)
        if value is not None:
            try:
                valid_kwargs[frame_attr_name] = repr_attr_class(value, unit=unit)
            except u.UnitConversionError as err:
                error_message = (
                    f"Unit '{unit}' ({unit.physical_type}) could not be applied to"
                    f" '{frame_attr_name}'. This can occur when passing units for some"
                    " coordinate components when other components are specified as"
                    " Quantity objects. Either pass a list of units for all components"
                    " (and unit-less coordinate data), or pass Quantities for all"
                    " components."
                )
                raise u.UnitConversionError(error_message) from err

    # also check the differentials.  They aren't included in the units keyword,
    # so we only look for the names.

    differential_type = frame.differential_type
    if differential_type is not None:
        for frame_name, repr_name in frame.get_representation_component_names(
            "s"
        ).items():
            diff_attr_class = differential_type.attr_classes[repr_name]
            value = kwargs.pop(frame_name, None)
            if value is not None:
                valid_kwargs[frame_name] = diff_attr_class(value)

    return valid_kwargs


def _parse_one_coord_str(coord_str: str, *, is_radec: bool = True) -> tuple[str, str]:
    """Parse longitude-like and latitude-like values from a string.

    Currently the following formats are always supported:

     * space separated 2-value or 6-value format

    If the input can be assumed to represent an RA and Dec then the
    following are additionally supported:

     * space separated <6-value format, this requires a plus or minus sign
       separation between RA and Dec
     * sign separated format
     * JHHMMSS.ss+DDMMSS.ss format, with up to two optional decimal digits
     * JDDDMMSS.ss+DDMMSS.ss format, with up to two optional decimal digits

    Parameters
    ----------
    coord_str : str
        Coordinate string to parse.
    is_radec : bool, keyword-only
        Whether the coordinates represent an RA and Dec.

    Returns
    -------
    longitude-like, latitude-like : str
        Parsed coordinate values. If ``is_radec`` is `True` then they are
        RA and Dec.
    """
    if not isinstance(coord_str, str):
        # This exception should never be raised from SkyCoord
        raise TypeError("coord_str must be a single str")
    split_coord = coord_str.split()
    if len(split_coord) == 6:
        return " ".join(split_coord[:3]), " ".join(split_coord[3:])
    if len(split_coord) == 2 or not is_radec:
        return tuple(split_coord)
    if len(split_coord) == 1 and (match_j := J_PREFIXED_RA_DEC_RE.match(coord_str)):
        ra, dec = match_j.groups()
        if len(ra.split(".", 1)[0]) == 7:
            ra = f"{ra[0:3]} {ra[3:5]} {ra[5:]}"
        else:
            ra = f"{ra[0:2]} {ra[2:4]} {ra[4:]}"
        return ra, f"{dec[0:3]} {dec[3:5]} {dec[5:]}"
    ra, *dec = PLUS_MINUS_RE.split(coord_str)
    return ra, " ".join(dec)
