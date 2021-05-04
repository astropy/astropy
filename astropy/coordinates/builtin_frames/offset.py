"""Offset Frames in 6D."""

from __future__ import annotations

# BUILT-IN
import warnings
from typing import Any, Callable, Dict, Type, TypeVar

# THIRD-PARTY
import numpy as np

# LOCAL
from astropy import units as u
from astropy.coordinates import representation as r
from astropy.coordinates.attributes.core import CoordinateAttribute
from astropy.coordinates.attributes.optional import RotationAttribute
from astropy.coordinates.baseframe import BaseCoordinateFrame, frame_transform_graph
from astropy.coordinates.errors import ConvertError
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix
from astropy.coordinates.transformations import AffineTransform, FunctionTransform
from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    # THIRD-PARTY
    from scipy.spatial.transform import Rotation as ScipyRotation
else:  # make no-op version
    warnings.warn("`scipy` must be installed for the Rotation attribute")


__all__ = ["OffsetFrame"]


Self = TypeVar("Self", bound="OffsetFrame")
RepT = TypeVar("RepT", bound=r.BaseRepresentation)


_offset_cache: dict[
    type[BaseCoordinateFrame], type[OffsetFrame]
] = {}  # Cache of offset frames for fast retrieval

zero_v = r.CartesianDifferential(  # zero velocity base
    u.Quantity(0, u.km / u.s), u.Quantity(0, u.km / u.s), u.Quantity(0, u.km / u.s)
)


def _check_coord_repr_diff_types(c: BaseCoordinateFrame) -> None:
    """Check Representation and Differential type of the Frame's data.

    Parameters
    ----------
    c : BaseCoordinateFrame
        Frame to check.

    Raises
    ------
    ConvertError
        If the data is 2D in either the Representation or Differential.
    """
    # TODO! need a more robust check for 2D- data.
    if isinstance(c.data, (r.UnitSphericalRepresentation, r.RadialRepresentation)):
        raise ConvertError(
            "Transforming to/from an OffsetFrame frame "
            "requires a 3D coordinate, e.g. (angle, angle, "
            "distance) or (x, y, z)."
        )

    # TODO! need a more robust check for 2D- data.
    need3d = (
        r.UnitSphericalDifferential,
        r.UnitSphericalCosLatDifferential,
        r.RadialDifferential,
    )
    if isinstance(c.data.differentials.get("s", None), need3d):
        raise ConvertError(
            "Transforming to/from a OffsetFrame frame "
            "requires a 3D velocity, e.g., proper motion "
            "components and radial velocity."
        )


def get_matrix_vectors(
    offset_frame: OffsetFrame, inverse: bool = False
) -> tuple[np.ndarray, r.CartesianRepresentation]:
    """Get matrix vectors of offset.

    Use the ``inverse`` argument to get the inverse transformation, matrix
    and offsets to go from OffsetFrame to from_frame.

    .. todo::  allow for rotations about the axes

    """
    of = offset_frame  # shorthand

    # Make translation in 6D.
    if "s" in of.origin.data.differentials:
        translation = of.origin.represent_as(
            r.CartesianRepresentation, s=r.CartesianDifferential
        )
    else:
        translation = of.origin.represent_as(r.CartesianRepresentation)
        translation = translation.with_differentials(zero_v)

    R = np.identity(3)
    offset = -translation

    return R, offset


def offset_to_offset(
    from_offset_coord: OffsetFrame, to_offset_frame: OffsetFrame
) -> OffsetFrame:
    """Transform between two offset frames."""
    # This transform goes through the parent frames on each side.
    # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
    intermediate_from = from_offset_coord.transform_to(from_offset_coord.origin)
    intermediate_to = intermediate_from.transform_to(to_offset_frame.origin)
    return intermediate_to.transform_to(to_offset_frame)


def reference_to_offset(reference_coord: BaseCoordinateFrame, offset_frame: OffsetFrame) -> tuple[np.ndarray, r.CartesianRepresentation]:
    """Convert a reference coordinate to an offset frame."""
    _check_coord_repr_diff_types(reference_coord)
    return get_matrix_vectors(offset_frame)


def offset_to_reference(offset_coord: OffsetFrame, _: BaseCoordinateFrame) -> tuple[np.ndarray, r.CartesianRepresentation]:
    """Convert an offset frame coordinate to the reference frame."""
    _check_coord_repr_diff_types(offset_coord)
    return get_matrix_vectors(offset_coord, inverse=True)


def make_offset_cls(framecls: type[BaseCoordinateFrame]) -> type[OffsetFrame]:
    """Make frame-offset class.

    Create a new class that is the sky offset frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.

    Parameters
    ----------
    framecls : `~astropy.coordinates.BaseCoordinateFrame` class
        The class to create the OffsetFrame of.

    Returns
    -------
    offsetframecls : class
        The class for the new offset frame.

    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct offset frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given python session.
    """
    # Check if already made
    if framecls in _offset_cache:
        return _offset_cache[framecls]

    # Create a new OffsetFrame subclass for this frame class.
    _OffsetFrameCls = type(
        "SkyOffset" + framecls.__name__,
        (OffsetFrame, framecls),
        {
            "origin": CoordinateAttribute(frame=framecls, default=None),
            # The following two have to be done because otherwise we use the
            # defaults of OffsetFrame set by BaseCoordinateFrame.
            "_default_representation": framecls._default_representation,
            "_default_differential": framecls._default_differential,
            "__doc__": OffsetFrame.__doc__,
        },
    )

    # Register transformations
    FunctionTransform(offset_to_offset, _OffsetFrameCls, _OffsetFrameCls,
                      register_graph=frame_transform_graph, priority=1)
    AffineTransform(reference_to_offset, framecls, _OffsetFrameCls,
                    register_graph=frame_transform_graph, priority=1)
    AffineTransform(offset_to_reference, _OffsetFrameCls, framecls,
                    register_graph=frame_transform_graph, priority=1)

    # Cache and return
    _offset_cache[framecls] = _OffsetFrameCls
    return _OffsetFrameCls


class OffsetFrame(BaseCoordinateFrame):
    """Offset Frame.

    A frame which is relative to some specific position and oriented to match
    its frame. OffsetFrames always have component names for spherical
    coordinates of ``lon``/``lat``, *not* the component names for the frame of
    ``origin``. This is useful for calculating offsets and dithers in the
    frame of the sky relative to an arbitrary position. Coordinates in this
    frame are both centered on the position specified by the ``origin``
    coordinate, *and* they are oriented in the same manner as the ``origin``
    frame.  E.g., if ``origin`` is `~astropy.coordinates.ICRS`, this object's
    ``lat`` will be pointed in the direction of Dec, while ``lon`` will point
    in the direction of RA. For more on offset frames, see
    :ref:`astropy-offset-frames`.

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data (or use the other
        keywords)
    origin : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies the origin of this frame.
    rotation : (3, 3) ndarray[float], optional
        A rotation matrix satisfying
        :func:`~astropy.coordinates.matrix_utilities.is_rotation`.
        If not specified, defaults to the identity (no-rotation) matrix.
        Can also be

    Notes
    -----
    ``OffsetFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``OffsetFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``origin``.
    """

    origin = CoordinateAttribute(default=None, frame=None)
    Rotation = RotationAttribute(default=None)

    def __new__(cls: type[Self], *args: Any, **kw: Any) -> Self:
        # We don't want to call this method if we've already set up
        # an offset frame for this class.
        if not (issubclass(cls, OffsetFrame) and cls is not OffsetFrame):
            # We get the origin argument, and handle it here.
            # this ensures that all subclasses define an origin.
            try:
                origin_frame = kw["origin"]
            except KeyError:
                raise TypeError(
                    "Can't initialize an OffsetFrame without origin= keyword."
                )
            # handle SkyCoord
            if hasattr(origin_frame, "frame"):
                origin_frame = origin_frame.frame

            newcls = make_offset_cls(origin_frame.__class__)
            return newcls.__new__(newcls, *args, **kw)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kw)

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        # deal with possible values of Rotation
        Rotation = kwargs.pop("Rotation", None)
        rotation = kwargs.pop("rotation", None)

        if rotation is not None and Rotation is not None:
            raise TypeError("can only specify one of 'Rotation' or 'rotation'.")
        elif rotation is None and Rotation is None:
            R = ScipyRotation.from_matrix(np.identity(3))
        else:
            rot = rotation if Rotation is None else Rotation
            if isinstance(rot, ScipyRotation):
                R = rot
            elif isinstance(rot, np.ndarray):
                R = ScipyRotation.from_matrix(rot.T)
            else:
                raise ValueError
        kwargs["Rotation"] = R

        super().__init__(*args, **kwargs)

        if self.origin is not None and not self.origin.has_data:
            raise ValueError("the origin supplied to OffsetFrame has no data.")
        elif self.has_data:
            self._set_offset_data_lon_wrap_angle(self.data)

    @property
    def rotation(self) -> np.ndarray:
        """Rotation matrix, from `~scipy.spatial.transform.Rotation` attribute.

        Raises
        ------
        ModuleNotFoundError
            if :mod:`scipy` is not installed.
        """
        return self.Rotation.as_matrix().T  # scipy uses active transform formalism

    @staticmethod
    def _set_offset_data_lon_wrap_angle(data: RepT) -> RepT:
        if hasattr(data, "lon"):
            data.lon.wrap_angle = u.Quantity(180.0, u.deg)
        return data

    def represent_as(
        self,
        base: type[RepT],
        s: type[r.BaseDifferential] | str = "base",
        in_frame_units: bool = False,
    ) -> RepT:
        """Represent frame data.

        Parameters
        ----------
        base : `astropy.coordinates.BaseRepresentation` class
            The type of the repreentation in which to represent the data.
        s : `astropy.coordinates.BaseDifferential` or str, optional
            The type of the differential in which to represent the data, or the
            string name of the differential type, or 'base'. Default is "base".
        in_frame_units : bool, optional
            Whether to keep the representation in the frame's units (`True`) or
            to convert to the representation's defaault units (`True`, default)

        Returns
        -------
        `astropy.coordinates.BaseRepresentation`
            The frame data, transformed to the specified representation (and
            possibly differential).
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        self._set_offset_data_lon_wrap_angle(data)
        return data

    def __reduce__(
        self,
    ) -> tuple[
        Callable[[BaseCoordinateFrame, np.ndarray], OffsetFrame],
        tuple[BaseCoordinateFrame, np.ndarray],
        dict[str, Any],
    ]:
        return (_offset_reducer, (self.origin, self.rotation), self.__dict__)


def _offset_reducer(origin: BaseCoordinateFrame, rotation: np.ndarray) -> OffsetFrame:
    return OffsetFrame.__new__(OffsetFrame, origin=origin, rotation=rotation)
