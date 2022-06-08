# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import inspect
from astropy import units as u
from astropy.coordinates.transformations import DynamicMatrixTransform, FunctionTransform
from astropy.coordinates.baseframe import (frame_transform_graph,
                                           BaseCoordinateFrame)
from astropy.coordinates.attributes import CoordinateAttribute, QuantityAttribute
from astropy.coordinates.matrix_utilities import (rotation_matrix,
                                                  matrix_product,
                                                  matrix_transpose)
from astropy.utils.metaclasses import InheritanceInMixMeta

_skyoffset_cache = {}


def make_skyoffset_cls(framecls):
    """
    Create a new class that is the sky offset frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.

    The new class will always have component names for spherical coordinates of
    ``lon``/``lat``.

    Parameters
    ----------
    framecls : `~astropy.coordinates.BaseCoordinateFrame` subclass
        The class to create the SkyOffsetFrame of.

    Returns
    -------
    skyoffsetframecls : class
        The class for the new skyoffset frame.

    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct skyoffset frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given python session.
    """
    return SkyOffsetFrame(framecls)


def _skyoffset_to_skyoffset(from_skyoffset_coord, to_skyoffset_frame):
    """Transform between two skyoffset frames."""

    # This transform goes through the parent frames on each side.
    # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
    intermediate_from = from_skyoffset_coord.transform_to(from_skyoffset_coord.origin)
    intermediate_to = intermediate_from.transform_to(to_skyoffset_frame.origin)
    return intermediate_to.transform_to(to_skyoffset_frame)


def _reference_to_skyoffset(reference_frame, skyoffset_frame):
    """Convert a reference coordinate to an sky offset frame."""

    # Define rotation matrices along the position angle vector, and
    # relative to the origin.
    origin = skyoffset_frame.origin.spherical
    mat1 = rotation_matrix(-skyoffset_frame.rotation, 'x')
    mat2 = rotation_matrix(-origin.lat, 'y')
    mat3 = rotation_matrix(origin.lon, 'z')
    return matrix_product(mat1, mat2, mat3)


def _skyoffset_to_reference(skyoffset_coord, reference_frame):
    """Convert an sky offset frame coordinate to the reference frame"""

    # use the forward transform, but just invert it
    R = _reference_to_skyoffset(reference_frame, skyoffset_coord)
    # transpose is the inverse because R is a rotation matrix
    return matrix_transpose(R)


class SkyOffsetFrame(BaseCoordinateFrame, metaclass=InheritanceInMixMeta):
    """
    A frame which is relative to some specific position and oriented to match
    its frame.

    SkyOffsetFrames always have component names for spherical coordinates
    of ``lon``/``lat``, *not* the component names for the frame of ``origin``.

    This is useful for calculating offsets and dithers in the frame of the sky
    relative to an arbitrary position. Coordinates in this frame are both centered on the position specified by the
    ``origin`` coordinate, *and* they are oriented in the same manner as the
    ``origin`` frame.  E.g., if ``origin`` is `~astropy.coordinates.ICRS`, this
    object's ``lat`` will be pointed in the direction of Dec, while ``lon``
    will point in the direction of RA.

    For more on skyoffset frames, see :ref:`astropy:astropy-skyoffset-frames`.

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    origin : coordinate-like
        The coordinate which specifies the origin of this frame. Note that this
        origin is used purely for on-sky location/rotation.  It can have a
        ``distance`` but it will not be used by this ``SkyOffsetFrame``.
    rotation : angle-like
        The final rotation of the frame about the ``origin``. The sign of
        the rotation is the left-hand rule.  That is, an object at a
        particular position angle in the un-rotated system will be sent to
        the positive latitude (z) direction in the final frame.


    Notes
    -----
    ``SkyOffsetFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``SkyOffsetFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``origin``.
    """

    rotation = QuantityAttribute(default=0, unit=u.deg)
    origin = CoordinateAttribute(default=None, frame=None)

    # ---------------------------------------------------------------
    # `astropy.utils.metaclasses.InheritanceInMixMeta` customizations

    @classmethod
    def _inmix_prepare_type(cls, framecls, base_cls):
        name = 'SkyOffset' + framecls.__name__
        bases = (base_cls, framecls)
        namespace = {
            'origin': CoordinateAttribute(frame=framecls, default=None),
            # The following two have to be done because otherwise we use the
            # defaults of SkyOffsetFrame set by BaseCoordinateFrame.
            '_default_representation': framecls._default_representation,
            '_default_differential': framecls._default_differential,
            '__doc__': SkyOffsetFrame.__doc__,
         }

        return name, bases, namespace

    @classmethod
    def _inmix_make_class(cls, framecls):
        if not issubclass(framecls, BaseCoordinateFrame):
            raise TypeError

        # Call super
        skyoffset_framecls = type(cls)._inmix_make_class(cls, framecls)

        # Register transformations
        frame_transform_graph.transform(FunctionTransform, skyoffset_framecls, skyoffset_framecls)(_skyoffset_to_skyoffset)
        frame_transform_graph.transform(DynamicMatrixTransform, framecls, skyoffset_framecls)(_reference_to_skyoffset)
        frame_transform_graph.transform(DynamicMatrixTransform, skyoffset_framecls, framecls)(_skyoffset_to_reference)

        # TODO! deprecate _skyoffset_cache
        _skyoffset_cache[framecls] = skyoffset_framecls

        return skyoffset_framecls

    @classmethod
    def _inmix_make_instance(cls, *args, **kwargs):
        try:
            origin_frame = kwargs['origin']
        except KeyError:
            raise TypeError("Can't initialize an SkyOffsetFrame without origin= keyword.") from None
        if hasattr(origin_frame, 'frame'):
            origin_frame = origin_frame.frame

        inmixcls = cls._inmix_make_class(type(origin_frame))
        return inmixcls(*args, **kwargs)

    # ---------------------------------------------------------------

    def __new__(cls, *args, **kwargs):
        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.origin is not None and not self.origin.has_data:
            raise ValueError('The origin supplied to SkyOffsetFrame has no '
                             'data.')
        if self.has_data:
            self._set_skyoffset_data_lon_wrap_angle(self.data)

    @staticmethod
    def _set_skyoffset_data_lon_wrap_angle(data):
        if hasattr(data, 'lon'):
            data.lon.wrap_angle = 180.0 * u.deg
        return data

    def represent_as(self, base, s='base', in_frame_units=False):
        """
        Ensure the wrap angle for any spherical
        representations.
        """
        data = super().represent_as(base, s, in_frame_units=in_frame_units)
        self._set_skyoffset_data_lon_wrap_angle(data)
        return data

    def __reduce__(self):
        return (_skyoffset_reducer, (self.origin,), self.__dict__)


def _skyoffset_reducer(origin):
    return SkyOffsetFrame.__new__(SkyOffsetFrame, origin=origin)
