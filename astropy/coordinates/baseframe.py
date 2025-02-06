# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Framework and base classes for coordinate frames/"low-level" coordinate
classes.
"""

from __future__ import annotations

__all__ = [
    "BaseCoordinateFrame",
    "CoordinateFrameInfo",
    "GenericFrame",
    "RepresentationMapping",
    "frame_transform_graph",
]

import copy
import operator
import warnings
from collections import defaultdict
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from astropy import units as u
from astropy.table import QTable
from astropy.utils import ShapedLikeNDArray
from astropy.utils.data_info import MixinInfo
from astropy.utils.decorators import format_doc, lazyproperty
from astropy.utils.exceptions import AstropyWarning
from astropy.utils.masked import MaskableShapedLikeNDArray, combine_masks

from . import representation as r
from .angles import Angle, position_angle
from .attributes import Attribute
from .errors import NonRotationTransformationError, NonRotationTransformationWarning
from .transformations import (
    DynamicMatrixTransform,
    StaticMatrixTransform,
    TransformGraph,
)

if TYPE_CHECKING:
    from typing import Literal

    from astropy.coordinates import Latitude, Longitude, SkyCoord
    from astropy.units import Unit

# the graph used for all transformations between frames
frame_transform_graph = TransformGraph()


def _get_repr_cls(value):
    """
    Return a valid representation class from ``value`` or raise exception.
    """
    if value in r.REPRESENTATION_CLASSES:
        value = r.REPRESENTATION_CLASSES[value]
    elif not isinstance(value, type) or not issubclass(value, r.BaseRepresentation):
        raise ValueError(
            f"Representation is {value!r} but must be a BaseRepresentation class "
            f"or one of the string aliases {list(r.REPRESENTATION_CLASSES)}"
        )
    return value


def _get_diff_cls(value):
    """
    Return a valid differential class from ``value`` or raise exception.

    As originally created, this is only used in the SkyCoord initializer, so if
    that is refactored, this function my no longer be necessary.
    """
    if value in r.DIFFERENTIAL_CLASSES:
        value = r.DIFFERENTIAL_CLASSES[value]
    elif not isinstance(value, type) or not issubclass(value, r.BaseDifferential):
        raise ValueError(
            f"Differential is {value!r} but must be a BaseDifferential class "
            f"or one of the string aliases {list(r.DIFFERENTIAL_CLASSES)}"
        )
    return value


def _get_repr_classes(base, **differentials):
    """Get valid representation and differential classes.

    Parameters
    ----------
    base : str or `~astropy.coordinates.BaseRepresentation` subclass
        class for the representation of the base coordinates.  If a string,
        it is looked up among the known representation classes.
    **differentials : dict of str or `~astropy.coordinates.BaseDifferentials`
        Keys are like for normal differentials, i.e., 's' for a first
        derivative in time, etc.  If an item is set to `None`, it will be
        guessed from the base class.

    Returns
    -------
    repr_classes : dict of subclasses
        The base class is keyed by 'base'; the others by the keys of
        ``diffferentials``.
    """
    base = _get_repr_cls(base)
    repr_classes = {"base": base}

    for name, differential_type in differentials.items():
        if differential_type == "base":
            # We don't want to fail for this case.
            differential_type = r.DIFFERENTIAL_CLASSES.get(base.name, None)

        elif differential_type in r.DIFFERENTIAL_CLASSES:
            differential_type = r.DIFFERENTIAL_CLASSES[differential_type]

        elif differential_type is not None and (
            not isinstance(differential_type, type)
            or not issubclass(differential_type, r.BaseDifferential)
        ):
            raise ValueError(
                "Differential is {differential_type!r} but must be a BaseDifferential"
                f" class or one of the string aliases {list(r.DIFFERENTIAL_CLASSES)}"
            )
        repr_classes[name] = differential_type
    return repr_classes


class RepresentationMapping(NamedTuple):
    """
    This :class:`~typing.NamedTuple` is used with the
    ``frame_specific_representation_info`` attribute to tell frames what
    attribute names (and default units) to use for a particular representation.
    ``reprname`` and ``framename`` should be strings, while ``defaultunit`` can
    be either an astropy unit, the string ``'recommended'`` (which is degrees
    for Angles, nothing otherwise), or None (to indicate that no unit mapping
    should be done).
    """

    reprname: str
    framename: str
    defaultunit: str | Unit = "recommended"


class CoordinateFrameInfo(MixinInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """

    attrs_from_parent = {"unit"}  # Unit is read-only
    _supports_indexing = False
    mask_val = np.ma.masked

    @staticmethod
    def default_format(val):
        repr_data = val.info._repr_data
        formats = ["{0." + compname + ".value:}" for compname in repr_data.components]
        return ",".join(formats).format(repr_data)

    @property
    def unit(self):
        repr_data = self._repr_data
        return (
            ",".join(
                str(getattr(repr_data, comp).unit) or "None"
                for comp in repr_data.components
            )
            if repr_data is not None
            else None
        )

    @property
    def _repr_data(self):
        coord = self._parent
        if coord is None or not coord.has_data:
            return None

        if issubclass(
            coord.representation_type, r.SphericalRepresentation
        ) and isinstance(coord.data, r.UnitSphericalRepresentation):
            repr_data = coord.represent_as(coord.data.__class__, in_frame_units=True)
        else:
            repr_data = coord.represent_as(
                coord.representation_type, in_frame_units=True
            )
        return repr_data

    def _represent_as_dict(self):
        coord = self._parent
        attrs = []
        representation_type = coord.representation_type
        differential_type = None
        if coord.has_data:
            attrs.extend(coord.representation_component_names)
            # Don't output distance unless it's actually distance.
            if issubclass(
                representation_type, r.SphericalRepresentation
            ) and isinstance(coord.data, r.UnitSphericalRepresentation):
                attrs = attrs[:-1]

            diff = coord.data.differentials.get("s")
            if diff is not None:
                differential_type = coord.differential_type
                diff_attrs = list(coord.get_representation_component_names("s"))
                if issubclass(
                    differential_type,
                    (r.SphericalDifferential, r.SphericalCosLatDifferential),
                ):
                    # Don't output proper motions or radial velocities unless
                    # they've been specifically passed in.
                    if isinstance(diff, r.RadialDifferential):
                        diff_attrs = diff_attrs[2:]
                    elif isinstance(
                        diff,
                        (
                            r.UnitSphericalDifferential,
                            r.UnitSphericalCosLatDifferential,
                        ),
                    ):
                        diff_attrs = diff_attrs[:-1]
                attrs.extend(diff_attrs)

        attrs.extend(frame_transform_graph.frame_attributes.keys())

        out = super()._represent_as_dict(attrs)

        out["representation_type"] = representation_type.name
        if differential_type is not None:
            out["differential_type"] = differential_type.name

        # Note that coord.info.unit is a fake composite unit (e.g. 'deg,deg,None'
        # or None,None,m) and is not stored. The individual attributes have
        # units.

        return out

    def new_like(self, coords, length, metadata_conflicts="warn", name=None):
        """A new consistent coordinate instance with the given length.

        Return a new SkyCoord or BaseCoordinateFrame instance which is
        consistent with the input coordinate objects ``coords`` and has
        ``length`` rows.  Being "consistent" is defined as being able to set an
        item from one to each of the rest without any exception being raised.

        This is intended for creating a new coordinate instance whose elements
        can be set in-place for table operations like join or vstack.  This is
        used when a coordinate object is used as a mixin column in an astropy
        Table.

        The data values are not predictable and it is expected that the consumer
        of the object will fill in all values.

        Parameters
        ----------
        coords : list
            List of input SkyCoord or BaseCoordinateFrame objects
        length : int
            Length of the output SkyCoord or BaseCoordinateFrame object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output name (sets output coord.info.name)

        Returns
        -------
        coord : |SkyCoord|, |BaseFrame|
            Instance of this class consistent with ``coords``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(
            coords, metadata_conflicts, name, ("meta", "description")
        )

        # Make a new coordinate with the desired length.
        coord0 = coords[0]
        out = coord0._apply(np.zeros_like, shape=(length,) + coord0.shape[1:])

        # Use __setitem__ machinery to check for consistency of all coords.
        # We use :0 to ensure we do not break on empty coordinates (with the
        # side benefit that we do not actually set anything).
        for coord in coords[1:]:
            try:
                out[:0] = coord[:0]
            except Exception as err:
                raise ValueError("Input coords are inconsistent.") from err

        # Set (merged) info attributes
        for attr in ("name", "meta", "description"):
            if attr in attrs:
                setattr(out.info, attr, attrs[attr])

        return out

    def _insert(self, obj, values, axis=0):
        """
        Make a copy with coordinate values inserted before the given indices.

        The values to be inserted must conform to the rules for in-place setting
        of the object.

        The API signature matches the ``np.insert`` API, but is more limited.
        The specification of insert index ``obj`` must be a single integer,
        and the ``axis`` must be ``0`` for simple insertion before the index.

        Parameters
        ----------
        obj : int
            Integer index before which ``values`` is inserted.
        values : array-like
            Value(s) to insert.  If the type of ``values`` is different
            from that of quantity, ``values`` is converted to the matching type.
        axis : int, optional
            Axis along which to insert ``values``.  Default is 0, which is the
            only allowed value and will insert a row.

        Returns
        -------
        coord : |SkyCoord|, |BaseFrame|
            Copy of instance with new values inserted.
        """
        # TODO: move this up to MixinInfo, so it can be used for Time too?
        coord = self._parent
        # Validate inputs: obj arg is integer, axis=0, coord is not a scalar, and
        # input index is in bounds.
        try:
            idx0 = operator.index(obj)
        except TypeError:
            raise TypeError("obj arg must be an integer")

        if axis != 0:
            raise ValueError("axis must be 0")

        if not coord.shape:
            raise TypeError(
                f"cannot insert into scalar {coord.__class__.__name__} object"
            )

        if abs(idx0) > len(coord):
            raise IndexError(
                f"index {idx0} is out of bounds for axis 0 with size {len(coord)}"
            )

        # Turn negative index into positive
        if idx0 < 0:
            idx0 = len(coord) + idx0

        n_values = len(values) if values.shape else 1

        # Finally make the new object with the correct length and set values for the
        # three sections, before insert, the insert, and after the insert.
        out = self.new_like([coord], len(coord) + n_values, name=self.name)

        # Set the output values. This is where validation of `values` takes place to ensure
        # that it can indeed be inserted.
        out[:idx0] = coord[:idx0]
        out[idx0 : idx0 + n_values] = values
        out[idx0 + n_values :] = coord[idx0:]

        return out


base_doc = """{__doc__}
    Parameters
    ----------
    data : `~astropy.coordinates.BaseRepresentation` subclass instance
        A representation object or ``None`` to have no data (or use the
        coordinate component arguments, see below).
    {components}
    representation_type : `~astropy.coordinates.BaseRepresentation` subclass, str, optional
        A representation class or string name of a representation class. This
        sets the expected input representation class, thereby changing the
        expected keyword arguments for the data passed in. For example, passing
        ``representation_type='cartesian'`` will make the classes expect
        position data with cartesian names, i.e. ``x, y, z`` in most cases
        unless overridden via ``frame_specific_representation_info``. To see this
        frame's names, check out ``<this frame>().representation_info``.
    differential_type : `~astropy.coordinates.BaseDifferential` subclass, str, dict, optional
        A differential class or dictionary of differential classes (currently
        only a velocity differential with key 's' is supported). This sets the
        expected input differential class, thereby changing the expected keyword
        arguments of the data passed in. For example, passing
        ``differential_type='cartesian'`` will make the classes expect velocity
        data with the argument names ``v_x, v_y, v_z`` unless overridden via
        ``frame_specific_representation_info``. To see this frame's names,
        check out ``<this frame>().representation_info``.
    copy : bool, optional
        If `True` (default), make copies of the input coordinate arrays.
        Can only be passed in as a keyword argument.
    {footer}
"""

_components = """
    *args, **kwargs
        Coordinate components, with names that depend on the subclass.
"""


@format_doc(base_doc, components=_components, footer="")
class BaseCoordinateFrame(MaskableShapedLikeNDArray):
    """
    The base class for coordinate frames.

    This class is intended to be subclassed to create instances of specific
    systems.  Subclasses can implement the following attributes:

    * `default_representation`
        A subclass of `~astropy.coordinates.BaseRepresentation` that will be
        treated as the default representation of this frame.  This is the
        representation assumed by default when the frame is created.

    * `default_differential`
        A subclass of `~astropy.coordinates.BaseDifferential` that will be
        treated as the default differential class of this frame.  This is the
        differential class assumed by default when the frame is created.

    * `~astropy.coordinates.Attribute` class attributes
       Frame attributes such as ``FK4.equinox`` or ``FK4.obstime`` are defined
       using a descriptor class.  See the narrative documentation or
       built-in classes code for details.

    * `frame_specific_representation_info`
        A dictionary mapping the name or class of a representation to a list of
        `~astropy.coordinates.RepresentationMapping` objects that tell what
        names and default units should be used on this frame for the components
        of that representation.

    Unless overridden via `frame_specific_representation_info`, velocity name
    defaults are:

      * ``pm_{lon}_cos{lat}``, ``pm_{lat}`` for `~astropy.coordinates.SphericalCosLatDifferential` velocity components
      * ``pm_{lon}``, ``pm_{lat}`` for `~astropy.coordinates.SphericalDifferential` velocity components
      * ``radial_velocity`` for any ``d_distance`` component
      * ``v_{x,y,z}`` for `~astropy.coordinates.CartesianDifferential` velocity components

    where ``{lon}`` and ``{lat}`` are the frame names of the angular components.
    """

    default_representation = None
    default_differential = None

    # Specifies special names and units for representation and differential
    # attributes.
    frame_specific_representation_info = {}

    frame_attributes = {}
    # Default empty frame_attributes dict

    # Declare that BaseCoordinateFrame can be used as a Table column by defining
    # the info property.
    info = CoordinateFrameInfo()

    def __init_subclass__(cls, **kwargs):
        # We first check for explicitly set values for these:
        default_repr = getattr(cls, "default_representation", None)
        default_diff = getattr(cls, "default_differential", None)
        repr_info = getattr(cls, "frame_specific_representation_info", None)
        # Then, to make sure this works for subclasses-of-subclasses, we also
        # have to check for cases where the attribute names have already been
        # replaced by underscore-prefaced equivalents by the logic below:
        if default_repr is None or isinstance(default_repr, property):
            default_repr = getattr(cls, "_default_representation", None)

        if default_diff is None or isinstance(default_diff, property):
            default_diff = getattr(cls, "_default_differential", None)

        if repr_info is None or isinstance(repr_info, property):
            repr_info = getattr(cls, "_frame_specific_representation_info", None)

        repr_info = cls._infer_repr_info(repr_info)

        # Make read-only properties for the frame class attributes that should
        # be read-only to make them immutable after creation.
        # We copy attributes instead of linking to make sure there's no
        # accidental cross-talk between classes
        cls._create_readonly_property(
            "default_representation",
            default_repr,
            "Default representation for position data",
        )
        cls._create_readonly_property(
            "default_differential",
            default_diff,
            "Default representation for differential data (e.g., velocity)",
        )
        cls._create_readonly_property(
            "frame_specific_representation_info",
            copy.deepcopy(repr_info),
            "Mapping for frame-specific component names",
        )

        # Set the frame attributes. We first construct the attributes from
        # superclasses, going in reverse order to keep insertion order,
        # and then add any attributes from the frame now being defined
        # (if any old definitions are overridden, this keeps the order).
        # Note that we cannot simply start with the inherited frame_attributes
        # since we could be a mixin between multiple coordinate frames.
        # TODO: Should this be made to use readonly_prop_factory as well or
        # would it be inconvenient for getting the frame_attributes from
        # classes?
        frame_attrs = {}
        for basecls in reversed(cls.__bases__):
            if issubclass(basecls, BaseCoordinateFrame):
                frame_attrs.update(basecls.frame_attributes)

        frame_attrs |= {k: v for k, v in vars(cls).items() if isinstance(v, Attribute)}

        cls.frame_attributes = frame_attrs

        # Deal with setting the name of the frame:
        if not hasattr(cls, "name"):
            cls.name = cls.__name__.lower()
        elif BaseCoordinateFrame not in cls.__bases__ and cls.name in [
            getattr(base, "name", None) for base in cls.__bases__
        ]:
            # This may be a subclass of a subclass of BaseCoordinateFrame,
            # like ICRS(BaseRADecFrame). In this case, cls.name will have been
            # set by init_subclass
            cls.name = cls.__name__.lower()

        # A cache that *must be unique to each frame class* - it is
        # insufficient to share them with superclasses, hence the need to put
        # them in the meta
        cls._frame_class_cache = {}

        super().__init_subclass__(**kwargs)

        # call this once here to initialize defaults
        # (via FrameAttribute.__get__/convert_input)
        cls.get_frame_attr_defaults()

    def __init__(
        self,
        *args,
        copy=True,
        representation_type=None,
        differential_type=None,
        **kwargs,
    ):
        self._attr_names_with_defaults = []

        self._representation = self._infer_representation(
            representation_type, differential_type
        )
        data = self._infer_data(args, copy, kwargs)  # possibly None.

        shapes = [] if data is None else [data.shape]

        # Set frame attributes, if any.
        # Keep track of their shapes, but do not broadcast them yet.
        values = {}
        for fnm, fdefault in self.get_frame_attr_defaults().items():
            # Read-only frame attributes are defined as FrameAttribute
            # descriptors which are not settable, so set 'real' attributes as
            # the name prefaced with an underscore.

            if fnm in kwargs:
                value = kwargs.pop(fnm)
                setattr(self, "_" + fnm, value)
                # Validate attribute by getting it.
                values[fnm] = value = getattr(self, fnm)
                shapes.append(getattr(value, "shape", ()))
            else:
                setattr(self, "_" + fnm, fdefault)
                self._attr_names_with_defaults.append(fnm)

        if kwargs:
            raise TypeError(
                f"Coordinate frame {self.__class__.__name__} got unexpected "
                f"keywords: {list(kwargs)}"
            )

        # Determine the overall shape of the frame.
        try:
            self._shape = np.broadcast_shapes(*shapes)
        except ValueError as exc:
            exc.add_note(
                f"non-scalar data and/or attributes with inconsistent shapes: {shapes}",
            )
            raise exc

        # Broadcast the data if necessary and set it
        if data is not None and data.shape != self._shape:
            try:
                # if broadcasting isn't strictly needed, avoid it
                # see https://github.com/astropy/astropy/issues/16219
                data = data.reshape(self._shape)
            except Exception:
                data = data._apply(np.broadcast_to, shape=self._shape, subok=True)
                if copy:
                    data = data.copy()
        self._data = data
        # Broadcast the attributes if necessary by getting them again
        # (we now know the shapes will be OK).
        for key in values:
            getattr(self, key)

        # The logic of this block is not related to the previous one
        if self.has_data:
            # This makes the cache keys backwards-compatible, but also adds
            # support for having differentials attached to the frame data
            # representation object.
            if "s" in self._data.differentials:
                # TODO: assumes a velocity unit differential
                key = (
                    self._data.__class__.__name__,
                    self._data.differentials["s"].__class__.__name__,
                    False,
                )
            else:
                key = (self._data.__class__.__name__, False)

            # Set up representation cache.
            self.cache["representation"][key] = self._data

    def _infer_representation(self, representation_type, differential_type):
        if representation_type is None and differential_type is None:
            return {"base": self.default_representation, "s": self.default_differential}

        if representation_type is None:
            representation_type = self.default_representation

        if isinstance(differential_type, type) and issubclass(
            differential_type, r.BaseDifferential
        ):
            # TODO: assumes the differential class is for the velocity
            # differential
            differential_type = {"s": differential_type}

        elif isinstance(differential_type, str):
            # TODO: assumes the differential class is for the velocity
            # differential
            diff_cls = r.DIFFERENTIAL_CLASSES[differential_type]
            differential_type = {"s": diff_cls}

        elif differential_type is None:
            if representation_type == self.default_representation:
                differential_type = {"s": self.default_differential}
            else:
                differential_type = {"s": "base"}  # see set_representation_cls()

        return _get_repr_classes(representation_type, **differential_type)

    def _infer_data(self, args, copy, kwargs):
        # if not set below, this is a frame with no data
        representation_data = None
        differential_data = None

        args = list(args)  # need to be able to pop them
        if args and (isinstance(args[0], r.BaseRepresentation) or args[0] is None):
            representation_data = args.pop(0)  # This can still be None
            if len(args) > 0:
                raise TypeError(
                    "Cannot create a frame with both a representation object "
                    "and other positional arguments"
                )

            if representation_data is not None:
                diffs = representation_data.differentials
                differential_data = diffs.get("s", None)
                if (differential_data is None and len(diffs) > 0) or (
                    differential_data is not None and len(diffs) > 1
                ):
                    raise ValueError(
                        "Multiple differentials are associated with the representation"
                        " object passed in to the frame initializer. Only a single"
                        f" velocity differential is supported. Got: {diffs}"
                    )

        else:
            representation_cls = self.get_representation_cls()
            # Get any representation data passed in to the frame initializer
            # using keyword or positional arguments for the component names
            repr_kwargs = {}
            for nmkw, nmrep in self.representation_component_names.items():
                if len(args) > 0:
                    # first gather up positional args
                    repr_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    repr_kwargs[nmrep] = kwargs.pop(nmkw)

            # special-case the Spherical->UnitSpherical if no `distance`

            if repr_kwargs:
                # TODO: determine how to get rid of the part before the "try" -
                # currently removing it has a performance regression for
                # unitspherical because of the try-related overhead.
                # Also frames have no way to indicate what the "distance" is
                if repr_kwargs.get("distance", True) is None:
                    del repr_kwargs["distance"]

                if (
                    issubclass(representation_cls, r.SphericalRepresentation)
                    and "distance" not in repr_kwargs
                ):
                    representation_cls = representation_cls._unit_representation

                try:
                    representation_data = representation_cls(copy=copy, **repr_kwargs)
                except TypeError as e:
                    # this except clause is here to make the names of the
                    # attributes more human-readable.  Without this the names
                    # come from the representation instead of the frame's
                    # attribute names.
                    try:
                        representation_data = representation_cls._unit_representation(
                            copy=copy, **repr_kwargs
                        )
                    except Exception:
                        msg = str(e)
                        names = self.get_representation_component_names()
                        for frame_name, repr_name in names.items():
                            msg = msg.replace(repr_name, frame_name)
                        msg = msg.replace("__init__()", f"{self.__class__.__name__}()")
                        e.args = (msg,)
                        raise e

            # Now we handle the Differential data:
            # Get any differential data passed in to the frame initializer
            # using keyword or positional arguments for the component names
            differential_cls = self.get_representation_cls("s")
            diff_component_names = self.get_representation_component_names("s")
            diff_kwargs = {}
            for nmkw, nmrep in diff_component_names.items():
                if len(args) > 0:
                    # first gather up positional args
                    diff_kwargs[nmrep] = args.pop(0)
                elif nmkw in kwargs:
                    diff_kwargs[nmrep] = kwargs.pop(nmkw)

            if diff_kwargs:
                if (
                    hasattr(differential_cls, "_unit_differential")
                    and "d_distance" not in diff_kwargs
                ):
                    differential_cls = differential_cls._unit_differential

                elif len(diff_kwargs) == 1 and "d_distance" in diff_kwargs:
                    differential_cls = r.RadialDifferential

                try:
                    differential_data = differential_cls(copy=copy, **diff_kwargs)
                except TypeError as e:
                    # this except clause is here to make the names of the
                    # attributes more human-readable.  Without this the names
                    # come from the representation instead of the frame's
                    # attribute names.
                    msg = str(e)
                    names = self.get_representation_component_names("s")
                    for frame_name, repr_name in names.items():
                        msg = msg.replace(repr_name, frame_name)
                    msg = msg.replace("__init__()", f"{self.__class__.__name__}()")
                    e.args = (msg,)
                    raise

        if len(args) > 0:
            raise TypeError(
                f"{type(self).__name__}.__init__ had {len(args)} remaining "
                "unhandled arguments"
            )

        if representation_data is None and differential_data is not None:
            raise ValueError(
                "Cannot pass in differential component data "
                "without positional (representation) data."
            )

        if differential_data:
            # Check that differential data provided has units compatible
            # with time-derivative of representation data.
            # NOTE: there is no dimensionless time while lengths can be
            # dimensionless (u.dimensionless_unscaled).
            for comp in representation_data.components:
                if (diff_comp := f"d_{comp}") in differential_data.components:
                    current_repr_unit = representation_data._units[comp]
                    current_diff_unit = differential_data._units[diff_comp]
                    expected_unit = current_repr_unit / u.s
                    if not current_diff_unit.is_equivalent(expected_unit):
                        for (
                            key,
                            val,
                        ) in self.get_representation_component_names().items():
                            if val == comp:
                                current_repr_name = key
                                break
                        for key, val in self.get_representation_component_names(
                            "s"
                        ).items():
                            if val == diff_comp:
                                current_diff_name = key
                                break
                        raise ValueError(
                            f'{current_repr_name} has unit "{current_repr_unit}" with'
                            f' physical type "{current_repr_unit.physical_type}", but'
                            f" {current_diff_name} has incompatible unit"
                            f' "{current_diff_unit}" with physical type'
                            f' "{current_diff_unit.physical_type}" instead of the'
                            f' expected "{(expected_unit).physical_type}".'
                        )

            representation_data = representation_data.with_differentials(
                {"s": differential_data}
            )

        return representation_data

    @classmethod
    def _infer_repr_info(cls, repr_info):
        # Unless overridden via `frame_specific_representation_info`, velocity
        # name defaults are (see also docstring for BaseCoordinateFrame):
        #   * ``pm_{lon}_cos{lat}``, ``pm_{lat}`` for
        #     `SphericalCosLatDifferential` proper motion components
        #   * ``pm_{lon}``, ``pm_{lat}`` for `SphericalDifferential` proper
        #     motion components
        #   * ``radial_velocity`` for any `d_distance` component
        #   * ``v_{x,y,z}`` for `CartesianDifferential` velocity components
        # where `{lon}` and `{lat}` are the frame names of the angular
        # components.
        if repr_info is None:
            repr_info = {}

        # the tuple() call below is necessary because if it is not there,
        # the iteration proceeds in a difficult-to-predict manner in the
        # case that one of the class objects hash is such that it gets
        # revisited by the iteration.  The tuple() call prevents this by
        # making the items iterated over fixed regardless of how the dict
        # changes
        for cls_or_name in tuple(repr_info.keys()):
            if isinstance(cls_or_name, str):
                # TODO: this provides a layer of backwards compatibility in
                # case the key is a string, but now we want explicit classes.
                repr_info[_get_repr_cls(cls_or_name)] = repr_info.pop(cls_or_name)

        # The default spherical names are 'lon' and 'lat'
        sph_repr = repr_info.setdefault(
            r.SphericalRepresentation,
            [RepresentationMapping("lon", "lon"), RepresentationMapping("lat", "lat")],
        )

        sph_component_map = {m.reprname: m.framename for m in sph_repr}
        lon = sph_component_map["lon"]
        lat = sph_component_map["lat"]

        ang_v_unit = u.mas / u.yr
        lin_v_unit = u.km / u.s

        sph_coslat_diff = repr_info.setdefault(
            r.SphericalCosLatDifferential,
            [
                RepresentationMapping("d_lon_coslat", f"pm_{lon}_cos{lat}", ang_v_unit),
                RepresentationMapping("d_lat", f"pm_{lat}", ang_v_unit),
                RepresentationMapping("d_distance", "radial_velocity", lin_v_unit),
            ],
        )
        sph_diff = repr_info.setdefault(
            r.SphericalDifferential,
            [
                RepresentationMapping("d_lon", f"pm_{lon}", ang_v_unit),
                RepresentationMapping("d_lat", f"pm_{lat}", ang_v_unit),
                RepresentationMapping("d_distance", "radial_velocity", lin_v_unit),
            ],
        )
        repr_info.setdefault(
            r.RadialDifferential,
            [RepresentationMapping("d_distance", "radial_velocity", lin_v_unit)],
        )
        repr_info.setdefault(
            r.CartesianDifferential,
            [RepresentationMapping(f"d_{c}", f"v_{c}", lin_v_unit) for c in "xyz"],
        )

        # Unit* classes should follow the same naming conventions
        # TODO: this adds some unnecessary mappings for the Unit classes, so
        # this could be cleaned up, but in practice doesn't seem to have any
        # negative side effects
        repr_info.setdefault(r.UnitSphericalRepresentation, sph_repr)
        repr_info.setdefault(r.UnitSphericalCosLatDifferential, sph_coslat_diff)
        repr_info.setdefault(r.UnitSphericalDifferential, sph_diff)

        return repr_info

    @classmethod
    def _create_readonly_property(cls, attr_name, value, doc=None):
        private_attr = "_" + attr_name

        def getter(self):
            return getattr(self, private_attr)

        setattr(cls, private_attr, value)
        setattr(cls, attr_name, property(getter, doc=doc))

    @lazyproperty
    def cache(self):
        """Cache for this frame, a dict.

        It stores anything that should be computed from the coordinate data (*not* from
        the frame attributes). This can be used in functions to store anything that
        might be expensive to compute but might be reused by some other function.
        E.g.::

            if 'user_data' in myframe.cache:
                data = myframe.cache['user_data']
            else:
                myframe.cache['user_data'] = data = expensive_func(myframe.lat)

        If in-place modifications are made to the frame data, the cache should
        be cleared::

            myframe.cache.clear()

        """
        return defaultdict(dict)

    @property
    def data(self):
        """
        The coordinate data for this object.  If this frame has no data, an
        `ValueError` will be raised.  Use `has_data` to
        check if data is present on this frame object.
        """
        if self._data is None:
            raise ValueError(
                f'The frame object "{self!r}" does not have associated data'
            )
        return self._data

    @property
    def has_data(self):
        """
        True if this frame has `data`, False otherwise.
        """
        return self._data is not None

    @property
    def shape(self):
        return self._shape

    def __bool__(self):
        return self.has_data and self.size > 0

    @property
    def size(self):
        return self.data.size

    @property
    def masked(self):
        """Whether the underlying data is masked.

        Raises
        ------
        ValueError
            If the frame has no associated data.
        """
        return self.data.masked

    def get_mask(self, *attrs):
        """Get the mask associated with these coordinates.

        Parameters
        ----------
        *attrs : str
            Attributes from which to get the masks to combine. Items can be
            dotted, like ``"data.lon", "data.lat"``. By default, get the
            combined mask of all components (including from differentials),
            ignoring possible masks of attributes.

        Returns
        -------
        mask : ~numpy.ndarray of bool
            The combined, read-only mask. If the instance is not masked, it
            is an array of `False` with the correct shape.

        Raises
        ------
        ValueError
            If the coordinate frame has no associated data.

        """
        if attrs:
            values = operator.attrgetter(*attrs)(self)
            if not isinstance(values, tuple):
                values = (values,)
            masks = [getattr(v, "mask", None) for v in values]
        elif self.data.masked:
            masks = [diff.mask for diff in self.data.differentials.values()]
            masks.append(self.data.mask)
        else:
            # Short-cut if the data is not masked.
            masks = []

        # Broadcast makes it readonly too.
        return np.broadcast_to(combine_masks(masks), self.shape)

    mask = property(
        get_mask,
        doc="""The mask associated with these coordinates.

    Combines the masks of all components of the underlying representation,
    including possible differentials.
    """,
    )

    @classmethod
    def get_frame_attr_defaults(cls):
        """Return a dict with the defaults for each frame attribute."""
        return {name: getattr(cls, name).default for name in cls.frame_attributes}

    def get_representation_cls(self, which="base"):
        """The class used for part of this frame's data.

        Parameters
        ----------
        which : ('base', 's', `None`)
            The class of which part to return.  'base' means the class used to
            represent the coordinates; 's' the first derivative to time, i.e.,
            the class representing the proper motion and/or radial velocity.
            If `None`, return a dict with both.

        Returns
        -------
        representation : `~astropy.coordinates.BaseRepresentation` or `~astropy.coordinates.BaseDifferential`.
        """
        return self._representation if which is None else self._representation[which]

    def set_representation_cls(self, base=None, s="base"):
        """Set representation and/or differential class for this frame's data.

        Parameters
        ----------
        base : str, `~astropy.coordinates.BaseRepresentation` subclass, optional
            The name or subclass to use to represent the coordinate data.
        s : `~astropy.coordinates.BaseDifferential` subclass, optional
            The differential subclass to use to represent any velocities,
            such as proper motion and radial velocity.  If equal to 'base',
            which is the default, it will be inferred from the representation.
            If `None`, the representation will drop any differentials.
        """
        if base is None:
            base = self._representation["base"]
        self._representation = _get_repr_classes(base=base, s=s)

    representation_type = property(
        fget=get_representation_cls,
        fset=set_representation_cls,
        doc="""The representation class used for this frame's data.

        This will be a subclass from `~astropy.coordinates.BaseRepresentation`.
        Can also be *set* using the string name of the representation. If you
        wish to set an explicit differential class (rather than have it be
        inferred), use the ``set_representation_cls`` method.
        """,
    )

    @property
    def differential_type(self):
        """
        The differential used for this frame's data.

        This will be a subclass from `~astropy.coordinates.BaseDifferential`.
        For simultaneous setting of representation and differentials, see the
        ``set_representation_cls`` method.
        """
        return self.get_representation_cls("s")

    @differential_type.setter
    def differential_type(self, value):
        self.set_representation_cls(s=value)

    @classmethod
    def _get_representation_info(cls):
        # This exists as a class method only to support handling frame inputs
        # without units, which are deprecated and will be removed.  This can be
        # moved into the representation_info property at that time.
        # note that if so moved, the cache should be acceessed as
        # self.__class__._frame_class_cache

        if (
            cls._frame_class_cache.get("last_reprdiff_hash", None)
            != r.get_reprdiff_cls_hash()
        ):
            repr_attrs = {}
            for repr_diff_cls in list(r.REPRESENTATION_CLASSES.values()) + list(
                r.DIFFERENTIAL_CLASSES.values()
            ):
                repr_attrs[repr_diff_cls] = {"names": [], "units": []}
                for c, c_cls in repr_diff_cls.attr_classes.items():
                    repr_attrs[repr_diff_cls]["names"].append(c)
                    rec_unit = u.deg if issubclass(c_cls, Angle) else None
                    repr_attrs[repr_diff_cls]["units"].append(rec_unit)

            for (
                repr_diff_cls,
                mappings,
            ) in cls._frame_specific_representation_info.items():
                # take the 'names' and 'units' tuples from repr_attrs,
                # and then use the RepresentationMapping objects
                # to update as needed for this frame.
                nms = repr_attrs[repr_diff_cls]["names"]
                uns = repr_attrs[repr_diff_cls]["units"]
                comptomap = {m.reprname: m for m in mappings}
                for i, c in enumerate(repr_diff_cls.attr_classes.keys()):
                    if (mapping := comptomap.get(c)) is not None:
                        nms[i] = mapping.framename
                        defaultunit = mapping.defaultunit

                        # need the isinstance because otherwise if it's a unit it
                        # will try to compare to the unit string representation
                        if not (
                            isinstance(defaultunit, str)
                            and defaultunit == "recommended"
                        ):
                            uns[i] = defaultunit
                            # else we just leave it as recommended_units says above

                # Convert to tuples so that this can't mess with frame internals
                repr_attrs[repr_diff_cls]["names"] = tuple(nms)
                repr_attrs[repr_diff_cls]["units"] = tuple(uns)

            cls._frame_class_cache["representation_info"] = repr_attrs
            cls._frame_class_cache["last_reprdiff_hash"] = r.get_reprdiff_cls_hash()
        return cls._frame_class_cache["representation_info"]

    @lazyproperty
    def representation_info(self):
        """
        A dictionary with the information of what attribute names for this frame
        apply to particular representations.
        """
        return self._get_representation_info()

    def get_representation_component_names(self, which="base"):
        cls = self.get_representation_cls(which)
        if cls is None:
            return {}
        return dict(zip(self.representation_info[cls]["names"], cls.attr_classes))

    def get_representation_component_units(self, which="base"):
        repr_or_diff_cls = self.get_representation_cls(which)
        if repr_or_diff_cls is None:
            return {}
        repr_attrs = self.representation_info[repr_or_diff_cls]
        return {k: v for k, v in zip(repr_attrs["names"], repr_attrs["units"]) if v}

    representation_component_names = property(get_representation_component_names)

    representation_component_units = property(get_representation_component_units)

    def _replicate(self, data, copy=False, **kwargs):
        """Base for replicating a frame, with possibly different attributes.

        Produces a new instance of the frame using the attributes of the old
        frame (unless overridden) and with the data given.

        Parameters
        ----------
        data : `~astropy.coordinates.BaseRepresentation` or None
            Data to use in the new frame instance.  If `None`, it will be
            a data-less frame.
        copy : bool, optional
            Whether data and the attributes on the old frame should be copied
            (default), or passed on by reference.
        **kwargs
            Any attributes that should be overridden.
        """
        # This is to provide a slightly nicer error message if the user tries
        # to use frame_obj.representation instead of frame_obj.data to get the
        # underlying representation object [e.g., #2890]
        if isinstance(data, type):
            raise TypeError(
                "Class passed as data instead of a representation instance. If you"
                " called frame.representation, this returns the representation class."
                " frame.data returns the instantiated object - you may want to  use"
                " this instead."
            )
        if copy and data is not None:
            data = data.copy()

        for attr in self.frame_attributes:
            if attr not in self._attr_names_with_defaults and attr not in kwargs:
                value = getattr(self, attr)
                kwargs[attr] = value.copy() if copy else value
        return self.__class__(data, copy=False, **kwargs)

    def replicate(self, copy=False, **kwargs):
        """
        Return a replica of the frame, optionally with new frame attributes.

        The replica is a new frame object that has the same data as this frame
        object and with frame attributes overridden if they are provided as extra
        keyword arguments to this method. If ``copy`` is set to `True` then a
        copy of the internal arrays will be made.  Otherwise the replica will
        use a reference to the original arrays when possible to save memory. The
        internal arrays are normally not changeable by the user so in most cases
        it should not be necessary to set ``copy`` to `True`.

        Parameters
        ----------
        copy : bool, optional
            If True, the resulting object is a copy of the data.  When False,
            references are used where  possible. This rule also applies to the
            frame attributes.
        **kwargs
            Any additional keywords are treated as frame attributes to be set on the
            new frame object.

        Returns
        -------
        frameobj : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            Replica of this object, but possibly with new frame attributes.
        """
        return self._replicate(self.data, copy=copy, **kwargs)

    def replicate_without_data(self, copy=False, **kwargs):
        """
        Return a replica without data, optionally with new frame attributes.

        The replica is a new frame object without data but with the same frame
        attributes as this object, except where overridden by extra keyword
        arguments to this method.  The ``copy`` keyword determines if the frame
        attributes are truly copied vs being references (which saves memory for
        cases where frame attributes are large).

        This method is essentially the converse of `realize_frame`.

        Parameters
        ----------
        copy : bool, optional
            If True, the resulting object has copies of the frame attributes.
            When False, references are used where  possible.
        **kwargs
            Any additional keywords are treated as frame attributes to be set on the
            new frame object.

        Returns
        -------
        frameobj : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            Replica of this object, but without data and possibly with new frame
            attributes.
        """
        return self._replicate(None, copy=copy, **kwargs)

    def realize_frame(self, data, **kwargs):
        """
        Generates a new frame with new data from another frame (which may or
        may not have data). Roughly speaking, the converse of
        `replicate_without_data`.

        Parameters
        ----------
        data : `~astropy.coordinates.BaseRepresentation`
            The representation to use as the data for the new frame.
        **kwargs
            Any additional keywords are treated as frame attributes to be set on the
            new frame object. In particular, `representation_type` can be specified.

        Returns
        -------
        frameobj : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
            A new object in *this* frame, with the same frame attributes as
            this one, but with the ``data`` as the coordinate data.

        """
        return self._replicate(data, **kwargs)

    def represent_as(self, base, s="base", in_frame_units=False):
        """
        Generate and return a new representation of this frame's `data`
        as a Representation object.

        Note: In order to make an in-place change of the representation
        of a Frame or SkyCoord object, set the ``representation``
        attribute of that object to the desired new representation, or
        use the ``set_representation_cls`` method to also set the differential.

        Parameters
        ----------
        base : subclass of BaseRepresentation or string
            The type of representation to generate.  Must be a *class*
            (not an instance), or the string name of the representation
            class.
        s : subclass of `~astropy.coordinates.BaseDifferential`, str, optional
            Class in which any velocities should be represented. Must be
            a *class* (not an instance), or the string name of the
            differential class.  If equal to 'base' (default), inferred from
            the base class.  If `None`, all velocity information is dropped.
        in_frame_units : bool, keyword-only
            Force the representation units to match the specified units
            particular to this frame

        Returns
        -------
        newrep : BaseRepresentation-derived object
            A new representation object of this frame's `data`.

        Raises
        ------
        AttributeError
            If this object had no `data`

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.coordinates import SkyCoord, CartesianRepresentation
        >>> coord = SkyCoord(0*u.deg, 0*u.deg)
        >>> coord.represent_as(CartesianRepresentation)  # doctest: +FLOAT_CMP
        <CartesianRepresentation (x, y, z) [dimensionless]
                (1., 0., 0.)>

        >>> coord.representation_type = CartesianRepresentation
        >>> coord  # doctest: +FLOAT_CMP
        <SkyCoord (ICRS): (x, y, z) [dimensionless]
            (1., 0., 0.)>
        """
        # In the future, we may want to support more differentials, in which
        # case one probably needs to define **kwargs above and use it here.
        # But for now, we only care about the velocity.
        repr_classes = _get_repr_classes(base=base, s=s)
        representation_cls = repr_classes["base"]
        # We only keep velocity information
        if "s" in self.data.differentials:
            # For the default 'base' option in which _get_repr_classes has
            # given us a best guess based on the representation class, we only
            # use it if the class we had already is incompatible.
            if s == "base" and (
                self.data.differentials["s"].__class__
                in representation_cls._compatible_differentials
            ):
                differential_cls = self.data.differentials["s"].__class__
            else:
                differential_cls = repr_classes["s"]
        elif s is None or s == "base":
            differential_cls = None
        else:
            raise TypeError(
                "Frame data has no associated differentials (i.e. the frame has no"
                " velocity data) - represent_as() only accepts a new representation."
            )

        if differential_cls:
            cache_key = (
                representation_cls.__name__,
                differential_cls.__name__,
                in_frame_units,
            )
        else:
            cache_key = (representation_cls.__name__, in_frame_units)

        if cached_repr := self.cache["representation"].get(cache_key):
            return cached_repr

        if differential_cls:
            # Sanity check to ensure we do not just drop radial
            # velocity.  TODO: should Representation.represent_as
            # allow this transformation in the first place?
            if (
                isinstance(self.data, r.UnitSphericalRepresentation)
                and issubclass(representation_cls, r.CartesianRepresentation)
                and not isinstance(
                    self.data.differentials["s"],
                    (
                        r.UnitSphericalDifferential,
                        r.UnitSphericalCosLatDifferential,
                        r.RadialDifferential,
                    ),
                )
            ):
                raise u.UnitConversionError(
                    "need a distance to retrieve a cartesian representation "
                    "when both radial velocity and proper motion are present, "
                    "since otherwise the units cannot match."
                )

            # TODO NOTE: only supports a single differential
            data = self.data.represent_as(representation_cls, differential_cls)
            diff = data.differentials["s"]  # TODO: assumes velocity
        else:
            data = self.data.represent_as(representation_cls)

        # If the new representation is known to this frame and has a defined
        # set of names and units, then use that.
        if in_frame_units and (
            new_attrs := self.representation_info.get(representation_cls)
        ):
            datakwargs = {comp: getattr(data, comp) for comp in data.components}
            for comp, new_attr_unit in zip(data.components, new_attrs["units"]):
                if new_attr_unit:
                    datakwargs[comp] = datakwargs[comp].to(new_attr_unit)
            data = data.__class__(copy=False, **datakwargs)

        if differential_cls:
            # the original differential
            data_diff = self.data.differentials["s"]

            # If the new differential is known to this frame and has a
            # defined set of names and units, then use that.
            if in_frame_units and (
                new_attrs := self.representation_info.get(differential_cls)
            ):
                diffkwargs = {comp: getattr(diff, comp) for comp in diff.components}
                for comp, new_attr_unit in zip(diff.components, new_attrs["units"]):
                    # Some special-casing to treat a situation where the
                    # input data has a UnitSphericalDifferential or a
                    # RadialDifferential. It is re-represented to the
                    # frame's differential class (which might be, e.g., a
                    # dimensional Differential), so we don't want to try to
                    # convert the empty component units
                    if (
                        isinstance(
                            data_diff,
                            (
                                r.UnitSphericalDifferential,
                                r.UnitSphericalCosLatDifferential,
                                r.RadialDifferential,
                            ),
                        )
                        and comp not in data_diff.__class__.attr_classes
                    ):
                        continue

                    # Try to convert to requested units. Since that might
                    # not be possible (e.g., for a coordinate with proper
                    # motion but without distance, one cannot convert to a
                    # cartesian differential in km/s), we allow the unit
                    # conversion to fail.  See gh-7028 for discussion.
                    if new_attr_unit and hasattr(diff, comp):
                        try:
                            diffkwargs[comp] = diffkwargs[comp].to(new_attr_unit)
                        except Exception:
                            pass

                diff = diff.__class__(copy=False, **diffkwargs)

                # Here we have to bypass using with_differentials() because
                # it has a validation check. But because
                # .representation_type and .differential_type don't point to
                # the original classes, if the input differential is a
                # RadialDifferential, it usually gets turned into a
                # SphericalCosLatDifferential (or whatever the default is)
                # with strange units for the d_lon and d_lat attributes.
                # This then causes the dictionary key check to fail (i.e.
                # comparison against `diff._get_deriv_key()`)
                data._differentials.update({"s": diff})

        self.cache["representation"][cache_key] = data
        return data

    def transform_to(self, new_frame):
        """
        Transform this object's coordinate data to a new frame.

        Parameters
        ----------
        new_frame : coordinate-like
            The frame to transform this coordinate frame into.

        Returns
        -------
        transframe : coordinate-like
            A new object with the coordinate data represented in the
            ``newframe`` system.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from .errors import ConvertError

        if self._data is None:
            raise ValueError("Cannot transform a frame with no data")

        if (
            getattr(self.data, "differentials", None)
            and hasattr(self, "obstime")
            and hasattr(new_frame, "obstime")
            and np.any(self.obstime != new_frame.obstime)
        ):
            raise NotImplementedError(
                "You cannot transform a frame that has velocities to another frame at a"
                " different obstime. If you think this should (or should not) be"
                " possible, please comment at"
                " https://github.com/astropy/astropy/issues/6280"
            )

        if hasattr(new_frame, "_sky_coord_frame"):
            # Input new_frame is not a frame instance or class and is most
            # likely a SkyCoord object.
            new_frame = new_frame._sky_coord_frame

        trans = frame_transform_graph.get_transform(self.__class__, new_frame.__class__)
        if trans is None:
            if new_frame is self.__class__:
                # no special transform needed, but should update frame info
                return new_frame.realize_frame(self.data)
            msg = "Cannot transform from {0} to {1}"
            raise ConvertError(msg.format(self.__class__, new_frame.__class__))
        return trans(self, new_frame)

    def is_transformable_to(self, new_frame):
        """
        Determines if this coordinate frame can be transformed to another
        given frame.

        Parameters
        ----------
        new_frame : `~astropy.coordinates.BaseCoordinateFrame` subclass or instance
            The proposed frame to transform into.

        Returns
        -------
        transformable : bool or str
            `True` if this can be transformed to ``new_frame``, `False` if
            not, or the string 'same' if ``new_frame`` is the same system as
            this object but no transformation is defined.

        Notes
        -----
        A return value of 'same' means the transformation will work, but it will
        just give back a copy of this object.  The intended usage is::

            if coord.is_transformable_to(some_unknown_frame):
                coord2 = coord.transform_to(some_unknown_frame)

        This will work even if ``some_unknown_frame``  turns out to be the same
        frame class as ``coord``.  This is intended for cases where the frame
        is the same regardless of the frame attributes (e.g. ICRS), but be
        aware that it *might* also indicate that someone forgot to define the
        transformation between two objects of the same frame class but with
        different attributes.
        """
        new_frame_cls = new_frame if isinstance(new_frame, type) else type(new_frame)
        if frame_transform_graph.get_transform(type(self), new_frame_cls):
            return True
        return "same" if new_frame_cls is type(self) else False

    def is_frame_attr_default(self, attrnm):
        """
        Determine whether or not a frame attribute has its value because it's
        the default value, or because this frame was created with that value
        explicitly requested.

        Parameters
        ----------
        attrnm : str
            The name of the attribute to check.

        Returns
        -------
        isdefault : bool
            True if the attribute ``attrnm`` has its value by default, False if
            it was specified at creation of this frame.
        """
        return attrnm in self._attr_names_with_defaults

    @staticmethod
    def _frameattr_equiv(left_fattr, right_fattr):  # noqa: PLR0911
        """
        Determine if two frame attributes are equivalent.  Implemented as a
        staticmethod mainly as a convenient location, although conceivable it
        might be desirable for subclasses to override this behavior.

        Primary purpose is to check for equality of representations.
        Secondary purpose is to check for equality of coordinate attributes,
        which first checks whether they themselves are in equivalent frames
        before checking for equality in the normal fashion.  This is because
        checking for equality with non-equivalent frames raises an error.
        """
        if left_fattr is right_fattr:
            # shortcut if it's exactly the same object
            return True
        elif left_fattr is None or right_fattr is None:
            # shortcut if one attribute is unspecified and the other isn't
            return False

        left_is_repr = isinstance(left_fattr, r.BaseRepresentationOrDifferential)
        if left_is_repr ^ isinstance(right_fattr, r.BaseRepresentationOrDifferential):
            return False
        if left_is_repr:
            # both are representations.
            if getattr(left_fattr, "differentials", False) or getattr(
                right_fattr, "differentials", False
            ):
                warnings.warn(
                    "Two representation frame attributes were checked for equivalence"
                    " when at least one of them has differentials.  This yields False"
                    " even if the underlying representations are equivalent (although"
                    " this may change in future versions of Astropy)",
                    AstropyWarning,
                )
                return False
            return np.all(
                left_fattr == right_fattr
                if type(left_fattr) is type(right_fattr)
                else left_fattr.to_cartesian() == right_fattr.to_cartesian()
            )

        left_is_coord = isinstance(left_fattr, BaseCoordinateFrame)
        if left_is_coord ^ isinstance(right_fattr, BaseCoordinateFrame):
            return False
        if left_is_coord:
            # both are coordinates
            return left_fattr.is_equivalent_frame(right_fattr) and np.all(
                left_fattr == right_fattr
            )

        return np.all(left_fattr == right_fattr)

    def is_equivalent_frame(self, other):
        """
        Checks if this object is the same frame as the ``other`` object.

        To be the same frame, two objects must be the same frame class and have
        the same frame attributes.  Note that it does *not* matter what, if any,
        data either object has.

        Parameters
        ----------
        other : :class:`~astropy.coordinates.BaseCoordinateFrame`
            the other frame to check

        Returns
        -------
        isequiv : bool
            True if the frames are the same, False if not.

        Raises
        ------
        TypeError
            If ``other`` isn't a `~astropy.coordinates.BaseCoordinateFrame` or subclass.
        """
        if self.__class__ == other.__class__:
            return all(
                self._frameattr_equiv(getattr(self, attr), getattr(other, attr))
                for attr in self.frame_attributes
            )
        elif not isinstance(other, BaseCoordinateFrame):
            raise TypeError(
                "Tried to do is_equivalent_frame on something that isn't a frame"
            )
        else:
            return False

    def __repr__(self):
        frameattrs = self._frame_attrs_repr()
        data_repr = self._data_repr()

        if frameattrs:
            frameattrs = f" ({frameattrs})"

        if data_repr:
            return f"<{self.__class__.__name__} Coordinate{frameattrs}: {data_repr}>"
        else:
            return f"<{self.__class__.__name__} Frame{frameattrs}>"

    def _data_repr(self):
        """Returns a string representation of the coordinate data."""
        if not self.has_data:
            return ""

        if rep_cls := self.representation_type:
            if isinstance(self.data, getattr(rep_cls, "_unit_representation", ())):
                rep_cls = self.data.__class__

            dif_cls = None
            if "s" in self.data.differentials:
                dif_cls = self.get_representation_cls("s")
                if isinstance(
                    dif_data := self.data.differentials["s"],
                    (
                        r.UnitSphericalDifferential,
                        r.UnitSphericalCosLatDifferential,
                        r.RadialDifferential,
                    ),
                ):
                    dif_cls = dif_data.__class__

            data = self.represent_as(rep_cls, dif_cls, in_frame_units=True)

            data_repr = repr(data)
            # Generate the list of component names out of the repr string
            part1, _, remainder = data_repr.partition("(")
            if remainder:
                comp_str, part2 = remainder.split(")", 1)
                # Swap in frame-specific component names
                invnames = {
                    nmrepr: nmpref
                    for nmpref, nmrepr in self.representation_component_names.items()
                }
                comp_names = (invnames.get(name, name) for name in comp_str.split(", "))
                # Reassemble the repr string
                data_repr = f"{part1}({', '.join(comp_names)}){part2}"

        else:
            data = self.data
            data_repr = repr(self.data)

        if data_repr.startswith(class_prefix := f"<{type(data).__name__} "):
            data_repr = data_repr.removeprefix(class_prefix).removesuffix(">")
        else:
            data_repr = "Data:\n" + data_repr

        if "s" not in self.data.differentials:
            return data_repr

        data_repr_spl = data_repr.split("\n")
        first, *middle, last = repr(data.differentials["s"]).split("\n")
        if first.startswith("<"):
            first = " " + first.split(" ", 1)[1]
        for frm_nm, rep_nm in self.get_representation_component_names("s").items():
            first = first.replace(rep_nm, frm_nm)
        data_repr_spl[-1] = "\n".join((first, *middle, last.removesuffix(">")))
        return "\n".join(data_repr_spl)

    def _frame_attrs_repr(self):
        """
        Returns a string representation of the frame's attributes, if any.
        """
        attr_strs = []
        for attribute_name in self.frame_attributes:
            attr = getattr(self, attribute_name)
            # Check to see if this object has a way of representing itself
            # specific to being an attribute of a frame. (Note, this is not the
            # Attribute class, it's the actual object).
            if hasattr(attr, "_astropy_repr_in_frame"):
                attrstr = attr._astropy_repr_in_frame()
            else:
                attrstr = str(attr)
            attr_strs.append(f"{attribute_name}={attrstr}")

        return ", ".join(attr_strs)

    def _apply(self, method, *args, **kwargs):
        """Create a new instance, applying a method to the underlying data.

        In typical usage, the method is any of the shape-changing methods for
        `~numpy.ndarray` (``reshape``, ``swapaxes``, etc.), as well as those
        picking particular elements (``__getitem__``, ``take``, etc.), which
        are all defined in `~astropy.utils.shapes.ShapedLikeNDArray`. It will be
        applied to the underlying arrays in the representation (e.g., ``x``,
        ``y``, and ``z`` for `~astropy.coordinates.CartesianRepresentation`),
        as well as to any frame attributes that have a shape, with the results
        used to create a new instance.

        Internally, it is also used to apply functions to the above parts
        (in particular, `~numpy.broadcast_to`).

        Parameters
        ----------
        method : str or callable
            If str, it is the name of a method that is applied to the internal
            ``components``. If callable, the function is applied.
        *args : tuple
            Any positional arguments for ``method``.
        **kwargs : dict
            Any keyword arguments for ``method``.
        """

        def apply_method(value):
            if isinstance(value, ShapedLikeNDArray):
                return value._apply(method, *args, **kwargs)
            else:
                if callable(method):
                    return method(value, *args, **kwargs)
                else:
                    return getattr(value, method)(*args, **kwargs)

        new = super().__new__(self.__class__)
        if hasattr(self, "_representation"):
            new._representation = self._representation.copy()
        new._attr_names_with_defaults = self._attr_names_with_defaults.copy()

        new_shape = ()
        for attr in self.frame_attributes:
            _attr = "_" + attr
            if attr in self._attr_names_with_defaults:
                setattr(new, _attr, getattr(self, _attr))
            else:
                value = getattr(self, _attr)
                if getattr(value, "shape", ()):
                    value = apply_method(value)
                    new_shape = new_shape or value.shape
                elif method == "copy" or method == "flatten":
                    # flatten should copy also for a single element array, but
                    # we cannot use it directly for array scalars, since it
                    # always returns a one-dimensional array. So, just copy.
                    value = copy.copy(value)

                setattr(new, _attr, value)

        if self.has_data:
            new._data = apply_method(self.data)
            new_shape = new_shape or new._data.shape
        else:
            new._data = None

        new._shape = new_shape

        # Copy other 'info' attr only if it has actually been defined.
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if "info" in self.__dict__:
            new.info = self.info

        return new

    def __setitem__(self, item, value):
        if value is np.ma.masked or value is np.ma.nomask:
            self.data.__setitem__(item, value)
            self.cache.clear()
            return

        if self.__class__ is not value.__class__:
            raise TypeError(
                f"can only set from object of same class: {self.__class__.__name__} vs."
                f" {value.__class__.__name__}"
            )

        if not self.is_equivalent_frame(value):
            raise ValueError("can only set frame item from an equivalent frame")

        if value._data is None:
            raise ValueError("can only set frame with value that has data")

        if self._data is None:
            raise ValueError("cannot set frame which has no data")

        if self.shape == ():
            raise TypeError(
                f"scalar '{self.__class__.__name__}' frame object "
                "does not support item assignment"
            )

        if self._data is None:
            raise ValueError("can only set frame if it has data")

        if self._data.__class__ is not value._data.__class__:
            raise TypeError(
                "can only set from object of same class: "
                f"{self._data.__class__.__name__} vs. {value._data.__class__.__name__}"
            )

        if self._data._differentials:
            # Can this ever occur? (Same class but different differential keys).
            # This exception is not tested since it is not clear how to generate it.
            if self._data._differentials.keys() != value._data._differentials.keys():
                raise ValueError("setitem value must have same differentials")

            for key, self_diff in self._data._differentials.items():
                if self_diff.__class__ is not value._data._differentials[key].__class__:
                    raise TypeError(
                        "can only set from object of same class: "
                        f"{self_diff.__class__.__name__} vs. "
                        f"{value._data._differentials[key].__class__.__name__}"
                    )

        # Set representation data
        self._data[item] = value._data

        # Frame attributes required to be identical by is_equivalent_frame,
        # no need to set them here.

        self.cache.clear()

    def insert(self, obj, values, axis=0):
        return self.info._insert(obj, values, axis)

    insert.__doc__ = CoordinateFrameInfo._insert.__doc__

    def __dir__(self):
        """
        Override the builtin `dir` behavior to include representation
        names.

        TODO: dynamic representation transforms (i.e. include cylindrical et al.).
        """
        return sorted(
            set(super().__dir__())
            | set(self.representation_component_names)
            | set(self.get_representation_component_names("s"))
        )

    def __getattr__(self, attr):
        """
        Allow access to attributes on the representation and differential as
        found via ``self.get_representation_component_names``.

        TODO: We should handle dynamic representation transforms here (e.g.,
        `.cylindrical`) instead of defining properties as below.
        """
        # attr == '_representation' is likely from the hasattr() test in the
        # representation property which is used for
        # self.representation_component_names.
        #
        # Prevent infinite recursion here.
        if attr.startswith("_"):
            return self.__getattribute__(attr)  # Raise AttributeError.

        repr_names = self.representation_component_names
        if attr in repr_names:
            if self._data is None:
                # this raises the "no data" error by design - doing it this way means we
                # don't have to replicate the error message here.
                self.data  # noqa: B018

            rep = self.represent_as(self.representation_type, in_frame_units=True)
            val = getattr(rep, repr_names[attr])
            return val

        diff_names = self.get_representation_component_names("s")
        if attr in diff_names:
            if self._data is None:
                self.data  # noqa: B018  # see above.
            # TODO: this doesn't work for the case when there is only
            # unitspherical information. The differential_type gets set to the
            # default_differential, which expects full information, so the
            # units don't work out
            rep = self.represent_as(
                in_frame_units=True, **self.get_representation_cls(None)
            )
            val = getattr(rep.differentials["s"], diff_names[attr])
            return val

        return self.__getattribute__(attr)  # Raise AttributeError.

    def __setattr__(self, attr, value):
        # Don't slow down access of private attributes!
        if not attr.startswith("_"):
            if hasattr(self, "representation_info"):
                repr_attr_names = set()
                for representation_attr in self.representation_info.values():
                    repr_attr_names.update(representation_attr["names"])

                if attr in repr_attr_names:
                    raise AttributeError(f"Cannot set any frame attribute {attr}")

        super().__setattr__(attr, value)

    def __eq__(self, value):
        """Equality operator for frame.

        This implements strict equality and requires that the frames are
        equivalent and that the representation data are exactly equal.
        """
        if not isinstance(value, BaseCoordinateFrame):
            return NotImplemented

        is_equiv = self.is_equivalent_frame(value)

        if self._data is None and value._data is None:
            # For Frame with no data, == compare is same as is_equivalent_frame()
            return is_equiv

        if not is_equiv:
            raise TypeError(
                "cannot compare: objects must have equivalent frames: "
                f"{self.replicate_without_data()} vs. {value.replicate_without_data()}"
            )

        if (value._data is None) != (self._data is None):
            raise ValueError(
                "cannot compare: one frame has data and the other does not"
            )

        return self._data == value._data

    def __ne__(self, value):
        return np.logical_not(self == value)

    def _prepare_unit_sphere_coords(
        self,
        other: BaseCoordinateFrame | SkyCoord,
        origin_mismatch: Literal["ignore", "warn", "error"],
    ) -> tuple[Longitude, Latitude, Longitude, Latitude]:
        other_frame = getattr(other, "frame", other)
        if not (
            origin_mismatch == "ignore"
            or self.is_equivalent_frame(other_frame)
            or all(
                isinstance(comp, (StaticMatrixTransform, DynamicMatrixTransform))
                for comp in frame_transform_graph.get_transform(
                    type(self), type(other_frame)
                ).transforms
            )
        ):
            if origin_mismatch == "warn":
                warnings.warn(NonRotationTransformationWarning(self, other_frame))
            elif origin_mismatch == "error":
                raise NonRotationTransformationError(self, other_frame)
            else:
                raise ValueError(
                    f"{origin_mismatch=} is invalid. Allowed values are 'ignore', "
                    "'warn' or 'error'."
                )
        self_sph = self.represent_as(r.UnitSphericalRepresentation)
        other_sph = other_frame.transform_to(self).represent_as(
            r.UnitSphericalRepresentation
        )
        return self_sph.lon, self_sph.lat, other_sph.lon, other_sph.lat

    def position_angle(self, other: BaseCoordinateFrame | SkyCoord) -> Angle:
        """Compute the on-sky position angle to another coordinate.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The other coordinate to compute the position angle to.  It is
            treated as the "head" of the vector of the position angle.

        Returns
        -------
        `~astropy.coordinates.Angle`
            The (positive) position angle of the vector pointing from ``self``
            to ``other``, measured East from North.  If either ``self`` or
            ``other`` contain arrays, this will be an array following the
            appropriate `numpy` broadcasting rules.

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.coordinates import ICRS, SkyCoord
        >>> c1 = SkyCoord(0*u.deg, 0*u.deg)
        >>> c2 = ICRS(1*u.deg, 0*u.deg)
        >>> c1.position_angle(c2).to(u.deg)
        <Angle 90. deg>
        >>> c2.position_angle(c1).to(u.deg)
        <Angle 270. deg>
        >>> c3 = SkyCoord(1*u.deg, 1*u.deg)
        >>> c1.position_angle(c3).to(u.deg)  # doctest: +FLOAT_CMP
        <Angle 44.995636455344844 deg>
        """
        return position_angle(*self._prepare_unit_sphere_coords(other, "ignore"))

    def separation(
        self,
        other: BaseCoordinateFrame | SkyCoord,
        *,
        origin_mismatch: Literal["ignore", "warn", "error"] = "warn",
    ) -> Angle:
        """
        Computes on-sky separation between this coordinate and another.

        For more on how to use this (and related) functionality, see the
        examples in :ref:`astropy-coordinates-separations-matching`.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The coordinate to get the separation to.
        origin_mismatch : {"warn", "ignore", "error"}, keyword-only
            If the ``other`` coordinates are in a different frame then they
            will have to be transformed, and if the transformation is not a
            pure rotation then ``self.separation(other)`` can be
            different from ``other.separation(self)``. With
            ``origin_mismatch="warn"`` (default) the transformation is
            always performed, but a warning is emitted if it is not a
            pure rotation. If ``origin_mismatch="ignore"`` then the
            required transformation is always performed without warnings.
            If ``origin_mismatch="error"`` then only transformations
            that are pure rotations are allowed.

        Returns
        -------
        sep : `~astropy.coordinates.Angle`
            The on-sky separation between this and the ``other`` coordinate.

        Notes
        -----
        The separation is calculated using the Vincenty formula, which
        is stable at all locations, including poles and antipodes [1]_.

        .. [1] https://en.wikipedia.org/wiki/Great-circle_distance

        """
        from .angles import Angle, angular_separation

        return Angle(
            angular_separation(
                *self._prepare_unit_sphere_coords(other, origin_mismatch)
            ),
            unit=u.degree,
        )

    def separation_3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        For more on how to use this (and related) functionality, see the
        examples in :ref:`astropy-coordinates-separations-matching`.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseCoordinateFrame` or `~astropy.coordinates.SkyCoord`
            The coordinate system to get the distance to.

        Returns
        -------
        sep : `~astropy.coordinates.Distance`
            The real-space distance between these two coordinates.

        Raises
        ------
        ValueError
            If this or the other coordinate do not have distances.
        """
        from .distances import Distance

        if isinstance(self.data, r.UnitSphericalRepresentation):
            raise ValueError(
                "This object does not have a distance; cannot compute 3d separation."
            )

        # do this first just in case the conversion somehow creates a distance
        other = getattr(other, "frame", other).transform_to(self)

        if isinstance(other, r.UnitSphericalRepresentation):
            raise ValueError(
                "The other object does not have a distance; "
                "cannot compute 3d separation."
            )

        # drop the differentials to ensure they don't do anything odd in the
        # subtraction
        dist = (
            self.data.without_differentials().represent_as(r.CartesianRepresentation)
            - other.data.without_differentials().represent_as(r.CartesianRepresentation)
        ).norm()
        return dist if dist.unit == u.one else Distance(dist)

    @property
    def cartesian(self):
        """
        Shorthand for a cartesian representation of the coordinates in this
        object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as("cartesian", in_frame_units=True)

    @property
    def cylindrical(self):
        """
        Shorthand for a cylindrical representation of the coordinates in this
        object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as("cylindrical", in_frame_units=True)

    @property
    def spherical(self):
        """
        Shorthand for a spherical representation of the coordinates in this
        object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as("spherical", in_frame_units=True)

    @property
    def sphericalcoslat(self):
        """
        Shorthand for a spherical representation of the positional data and a
        `~astropy.coordinates.SphericalCosLatDifferential` for the velocity
        data in this object.
        """
        # TODO: if representations are updated to use a full transform graph,
        #       the representation aliases should not be hard-coded like this
        return self.represent_as("spherical", "sphericalcoslat", in_frame_units=True)

    @property
    def velocity(self):
        """
        Shorthand for retrieving the Cartesian space-motion as a
        `~astropy.coordinates.CartesianDifferential` object.

        This is equivalent to calling ``self.cartesian.differentials['s']``.
        """
        if "s" not in self.data.differentials:
            raise ValueError(
                "Frame has no associated velocity (Differential) data information."
            )

        return self.cartesian.differentials["s"]

    @property
    def proper_motion(self):
        """
        Shorthand for the two-dimensional proper motion as a
        `~astropy.units.Quantity` object with angular velocity units. In the
        returned `~astropy.units.Quantity`, ``axis=0`` is the longitude/latitude
        dimension so that ``.proper_motion[0]`` is the longitudinal proper
        motion and ``.proper_motion[1]`` is latitudinal. The longitudinal proper
        motion already includes the cos(latitude) term.
        """
        if "s" not in self.data.differentials:
            raise ValueError(
                "Frame has no associated velocity (Differential) data information."
            )

        sph = self.represent_as("spherical", "sphericalcoslat", in_frame_units=True)
        pm_lon = sph.differentials["s"].d_lon_coslat
        pm_lat = sph.differentials["s"].d_lat
        return (
            np.stack((pm_lon.value, pm_lat.to(pm_lon.unit).value), axis=0) * pm_lon.unit
        )

    @property
    def radial_velocity(self):
        """
        Shorthand for the radial or line-of-sight velocity as a
        `~astropy.units.Quantity` object.
        """
        if "s" not in self.data.differentials:
            raise ValueError(
                "Frame has no associated velocity (Differential) data information."
            )

        sph = self.represent_as("spherical", in_frame_units=True)
        return sph.differentials["s"].d_distance

    def to_table(self) -> QTable:
        """
        Convert this |BaseFrame| to a |QTable|.

        Any attributes that have the same length as the |BaseFrame| will be
        converted to columns of the |QTable|. All other attributes will be
        recorded as metadata.

        Returns
        -------
        `~astropy.table.QTable`
            A |QTable| containing the data of this |BaseFrame|.

        Examples
        --------
        >>> from astropy.coordinates import ICRS
        >>> coord = ICRS(ra=[40, 70]*u.deg, dec=[0, -20]*u.deg)
        >>> t =  coord.to_table()
        >>> t
        <QTable length=2>
           ra     dec
          deg     deg
        float64 float64
        ------- -------
           40.0     0.0
           70.0   -20.0
        >>> t.meta
        {'representation_type': 'spherical'}
        """
        columns = {}
        metadata = {}
        # Record attributes that have the same length as self as columns in the
        # table, and the other attributes as table metadata.  This matches
        # table.serialize._represent_mixin_as_column().
        for key, value in self.info._represent_as_dict().items():
            if getattr(value, "shape", ())[:1] == (len(self),):
                columns[key] = value
            else:
                metadata[key] = value
        return QTable(columns, meta=metadata)


class GenericFrame(BaseCoordinateFrame):
    """
    A frame object that can't store data but can hold any arbitrary frame
    attributes. Mostly useful as a utility for the high-level class to store
    intermediate frame attributes.

    Parameters
    ----------
    frame_attrs : dict
        A dictionary of attributes to be used as the frame attributes for this
        frame.
    """

    name = None  # it's not a "real" frame so it doesn't have a name

    def __init__(self, frame_attrs):
        self.frame_attributes = {}
        for name, default in frame_attrs.items():
            self.frame_attributes[name] = Attribute(default)
            setattr(self, "_" + name, default)

        super().__init__(None)

    def __getattr__(self, name):
        if "_" + name in self.__dict__:
            return getattr(self, "_" + name)
        else:
            raise AttributeError(f"no {name}")

    def __setattr__(self, name, value):
        if name in self.frame_attributes:
            raise AttributeError(f"can't set frame attribute '{name}'")
        super().__setattr__(name, value)
