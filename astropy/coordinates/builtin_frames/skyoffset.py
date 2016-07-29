# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Note: `from __future__ import unicode_literals` is omitted here on purpose.
# Adding it leads to str / unicode errors on Python 2
from __future__ import (absolute_import, division, print_function)

from ... import units as u
from ..transformations import DynamicMatrixTransform, FunctionTransform
from ..baseframe import (CoordinateAttribute, QuantityFrameAttribute,
                         frame_transform_graph, RepresentationMapping,
                         BaseCoordinateFrame)
from ..matrix_utilities import (rotation_matrix,
                                matrix_product, matrix_transpose)
from ...utils.compat import namedtuple_asdict

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
    framecls : coordinate frame class (i.e., subclass of `~astropy.coordinates.BaseCoordinateFrame`)
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

    if framecls in _skyoffset_cache:
        return _skyoffset_cache[framecls]

    # the class of a class object is the metaclass
    framemeta = framecls.__class__

    class SkyOffsetMeta(framemeta):
        """
        This metaclass renames the class to be "SkyOffset<framecls>" and also
        adjusts the frame specific representation info so that spherical names
        are always "lon" and "lat" (instead of e.g. "ra" and "dec").
        """

        def __new__(cls, name, bases, members):
            # Only 'origin' is needed here, to set the origin frame properly.
            members['origin'] = CoordinateAttribute(frame=framecls, default=None)

            # This has to be done because FrameMeta will set these attributes
            # to the defaults from BaseCoordinateFrame when it creates the base
            # SkyOffsetFrame class initially.
            members['_frame_specific_representation_info'] = framecls._frame_specific_representation_info
            members['_default_representation'] = framecls._default_representation

            newname = name[:-5] if name.endswith('Frame') else name
            newname += framecls.__name__

            res = super(SkyOffsetMeta, cls).__new__(cls, newname, bases, members)

            # now go through all the component names and make any spherical names be "lon" and "lat"
            # instead of e.g. "ra" and "dec"

            lists_done = []
            for nm, component_list in res._frame_specific_representation_info.items():
                if nm in ('spherical', 'unitspherical'):
                    gotlatlon = []
                    for i, comp in enumerate(component_list):
                        if component_list in lists_done:
                            # we need this because sometimes the component_
                            # list's are the exact *same* object for both
                            # spherical and unitspherical.  So looping then makes
                            # the change *twice*.  This hack bypasses that.
                            continue

                        if comp.reprname in ('lon', 'lat'):
                            dct = namedtuple_asdict(comp)
                            # this forces the component names to be 'lat' and
                            # 'lon' regardless of what the actual base frame
                            # might use
                            dct['framename'] = comp.reprname
                            component_list[i] = type(comp)(**dct)
                            gotlatlon.append(comp.reprname)
                    if 'lon' not in gotlatlon:
                        rmlon = RepresentationMapping('lon', 'lon', 'recommended')
                        component_list.insert(0, rmlon)
                    if 'lat' not in gotlatlon:
                        rmlat = RepresentationMapping('lat', 'lat', 'recommended')
                        component_list.insert(0, rmlat)
                    lists_done.append(component_list)

            return res

    # We need this to handle the intermediate metaclass correctly, otherwise we could
    # just subclass SkyOffsetFrame.
    _SkyOffsetFramecls = SkyOffsetMeta('SkyOffsetFrame', (SkyOffsetFrame, framecls),
                                 {'__doc__': SkyOffsetFrame.__doc__})

    @frame_transform_graph.transform(FunctionTransform, _SkyOffsetFramecls, _SkyOffsetFramecls)
    def skyoffset_to_skyoffset(from_skyoffset_coord, to_skyoffset_frame):
        """Transform between two skyoffset frames."""

        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
        intermediate_from = from_skyoffset_coord.transform_to(from_skyoffset_coord.origin)
        intermediate_to = intermediate_from.transform_to(to_skyoffset_frame.origin)
        return intermediate_to.transform_to(to_skyoffset_frame)

    @frame_transform_graph.transform(DynamicMatrixTransform, framecls, _SkyOffsetFramecls)
    def reference_to_skyoffset(reference_frame, skyoffset_frame):
        """Convert a reference coordinate to an sky offset frame."""

        # Define rotation matrices along the position angle vector, and
        # relative to the origin.
        origin = skyoffset_frame.origin.spherical
        mat1 = rotation_matrix(-skyoffset_frame.rotation, 'x')
        mat2 = rotation_matrix(-origin.lat, 'y')
        mat3 = rotation_matrix(origin.lon, 'z')
        return matrix_product(mat1, mat2, mat3)

    @frame_transform_graph.transform(DynamicMatrixTransform, _SkyOffsetFramecls, framecls)
    def skyoffset_to_reference(skyoffset_coord, reference_frame):
        """Convert an sky offset frame coordinate to the reference frame"""

        # use the forward transform, but just invert it
        R = reference_to_skyoffset(reference_frame, skyoffset_coord)
        # transpose is the inverse because R is a rotation matrix
        return matrix_transpose(R)

    _skyoffset_cache[framecls] = _SkyOffsetFramecls
    return _SkyOffsetFramecls


class SkyOffsetFrame(BaseCoordinateFrame):
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

    For more on skyoffset frames, see :ref:`astropy-skyoffset-frames`.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    origin : `SkyCoord` or low-level coordinate object.
        the coordinate which specifies the origin of this frame.
    rotation : `~astropy.coordinates.Angle` or `~astropy.units.Quantity` with angle units
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

    rotation = QuantityFrameAttribute(default=0, unit=u.deg)
    origin = CoordinateAttribute(default=None, frame=None)

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an skyoffset frame for this class.
        if not (issubclass(cls, SkyOffsetFrame) and cls is not SkyOffsetFrame):
            # We get the origin argument, and handle it here.
            try:
                origin_frame = kwargs['origin']
            except KeyError:
                raise TypeError("Can't initialize an SkyOffsetFrame without origin= keyword.")
            if hasattr(origin_frame, 'frame'):
                origin_frame = origin_frame.frame
            newcls = make_skyoffset_cls(origin_frame.__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super(SkyOffsetFrame, cls).__new__ is object.__new__:
            return super(SkyOffsetFrame, cls).__new__(cls)
        return super(SkyOffsetFrame, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        super(SkyOffsetFrame, self).__init__(*args, **kwargs)
        if self.origin is not None and not self.origin.has_data:
            raise ValueError('The origin supplied to SkyOffsetFrame has no '
                             'data.')
        if self.has_data and hasattr(self.data, 'lon'):
            self.data.lon.wrap_angle = 180*u.deg
