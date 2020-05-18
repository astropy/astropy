# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np

from astropy import coordinates
from astropy import log
from astropy import units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.baseframe import (frame_transform_graph,
                                           BaseCoordinateFrame)
from astropy.coordinates.attributes import CoordinateAttribute, TimeAttribute

from .icrs import ICRS
from .utils import DEFAULT_OBSTIME

_astrometric_cache = {}


def make_astrometric_cls(framecls):
    """
    Create a new class that is the astrometric frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.

    Parameters
    ----------
    framecls : coordinate frame class (i.e., subclass of `~astropy.coordinates.BaseCoordinateFrame`)
        The class to create the AstrometricFrame of.

    Returns
    -------
    astrometricframecls : class
        The class for the new astrometric frame.

    Notes
    -----
    This function is necessary because Astropy's frame transformations depend
    on connection between specific frame *classes*.  So each type of frame
    needs its own distinct astrometric frame class.  This function generates
    just that class, as well as ensuring that only one example of such a class
    actually gets created in any given python session.
    """

    if framecls in _astrometric_cache:
        return _astrometric_cache[framecls]

    # the class of a class object is the metaclass
    framemeta = framecls.__class__

    class AstrometricMeta(framemeta):
        """
        This metaclass renames the class to be "Astrometric<framecls>".
        """

        def __new__(cls, name, bases, members):
            # Define the frame attributes that are not copied from `origin`
            members['origin'] = CoordinateAttribute(frame=framecls, default=None)
            if 'obstime' not in framecls.get_frame_attr_names():
                members['obstime'] = TimeAttribute(default=DEFAULT_OBSTIME)
            members['ref_epoch'] = TimeAttribute(default=DEFAULT_OBSTIME)

            # This has to be done because FrameMeta will set these attributes
            # to the defaults from BaseCoordinateFrame when it creates the base
            # AstrometricFrame class initially.
            members['_default_representation'] = framecls._default_representation
            members['_default_differential'] = framecls._default_differential

            members['_frame_specific_representation_info'] = framecls._frame_specific_representation_info

            newname = name[:-5] if name.endswith('Frame') else name
            newname += framecls.__name__

            return super().__new__(cls, newname, bases, members)

    # We need this to handle the intermediate metaclass correctly, otherwise we could
    # just subclass AstrometricFrame.
    _AstrometricFramecls = AstrometricMeta('AstrometricFrame',
                                           (AstrometricFrame, framecls),
                                           {'__doc__': AstrometricFrame.__doc__})

    @frame_transform_graph.transform(FunctionTransform, _AstrometricFramecls, _AstrometricFramecls)
    def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
        """Transform between two astrometric frames."""

        # This transform goes through the parent frames on each side.
        # from_frame -> from_frame.origin -> to_frame.origin -> to_frame
        intermediate_from = from_astrometric_coord.transform_to(from_astrometric_coord.origin)
        intermediate_to = intermediate_from.transform_to(to_astrometric_frame.origin)
        return intermediate_to.transform_to(to_astrometric_frame)

    @frame_transform_graph.transform(FunctionTransform, framecls, _AstrometricFramecls)
    def reference_to_astrometric(reference_coord, astrometric_frame):
        """Add motion to the reference coordinate."""
        reference_coord = reference_coord.transform_to(astrometric_frame.origin)
        icrs_coord = reference_coord.transform_to(ICRS)

        posvel = icrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < astrometric_frame._distance_threshold):  # solar-system body
                pos = posvel.without_differentials()
                vel = posvel.differentials['s'].to_cartesian()

                # The light travel time (dt) for a body in linear motion satisfies the equation:
                #   norm(D - V * dt) = c * dt
                # where D is the observer-body vector, V is the body velocity vector,
                # and c is the speed of light.  This turns into a quadratic equation in dt, and
                # only the positive solution is valid:
                #   dt = (sqrt((D dot V)^2 + (D dot D) * c2_V2) - D dot V) / c2_V2
                # where c2_V2 = c^2 - V dot V
                D = pos - astrometric_frame._observer_icrs_cartesian
                c2_V2 = speed_of_light ** 2 - vel.dot(vel)
                DdotV = D.dot(vel)
                dt = (np.sqrt(DdotV ** 2 + D.dot(D) * c2_V2) - DdotV) / c2_V2
                log.debug(f"Retarding for {dt.to(u.s)} of light travel time")

                pos -= vel * dt
                posvel = pos.with_differentials(posvel.differentials)
                icrs_coord = ICRS(posvel)
            elif np.all(distance >= astrometric_frame._distance_threshold):  # cosmic object
                dt = (astrometric_frame.obstime - astrometric_frame.ref_epoch).to(u.yr)
                log.debug(f"Propagating forward by {dt} after {astrometric_frame.ref_epoch}")
                icrs_coord = coordinates.SkyCoord(icrs_coord).apply_space_motion(dt=dt).frame
            else:  # a mix
                raise ConvertError("The transformation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")

        reference_coord = icrs_coord.transform_to(reference_coord)
        return astrometric_frame.realize_frame(reference_coord.data)

    @frame_transform_graph.transform(FunctionTransform, _AstrometricFramecls, framecls)
    def astrometric_to_reference(astrometric_coord, reference_frame):
        """Remove motion from the astrometric coordinate."""

        base_coord = astrometric_coord.as_base()
        icrs_coord = base_coord.transform_to(ICRS)

        posvel = icrs_coord.cartesian
        if 's' in posvel.differentials:
            distance = posvel.norm()
            if np.all(distance < astrometric_coord._distance_threshold):  # solar-system body
                pos = posvel.without_differentials()

                # The light travel time is straightforward to compute for an astrometric location
                dt = (pos - astrometric_coord._observer_icrs_cartesian).norm() / speed_of_light
                log.debug(f"Advancing for {dt.to(u.s)} of light travel time")

                pos += posvel.differentials['s'] * dt
                posvel = pos.with_differentials(posvel.differentials)
                icrs_coord = ICRS(posvel)
            elif np.all(distance >= astrometric_coord._distance_threshold):  # cosmic object
                dt = (astrometric_coord.obstime - astrometric_coord.ref_epoch).to(u.yr)
                log.debug(f"Propagating backward by {dt} before {astrometric_coord.ref_epoch}")
                icrs_coord = coordinates.SkyCoord(icrs_coord).apply_space_motion(dt=-dt).frame
            else:  # a mix
                raise ConvertError("The transformation cannot handle a mix of solar-system bodies "
                                   "and cosmic objects.")

        base_coord = icrs_coord.transform_to(base_coord)
        return base_coord.transform_to(reference_frame)

    _astrometric_cache[framecls] = _AstrometricFramecls
    return _AstrometricFramecls


class AstrometricFrame(BaseCoordinateFrame):
    """
    An astrometric version of a base coordinate frame, relative to a specified
    observer location.

    .. note::

        For more on astrometric frames, see
        :ref:`astropy-coordinates-astrometricframe`.

    The astrometric location of a body is the location of the body as measured
    by an observer, which depends on the motion of that body:

    * For a body in the solar system, its coordinate normally represents the
      instantaneous location of that body.  However, due to the finite speed of
      light, the light that reaches an observer at a certain observation time
      was emitted at an earlier time.  If the body is moving, then the observer
      will measure the body to be where it was at that earlier time.  This
      effect is also known as "planetary aberration".

    * For a cosmic object, its coordinate normally represents its catalog
      location at a reference epoch (e.g., J2000.0).  If the body is moving
      (i.e., has non-zero "proper motion" and/or "radial velocity"), its
      measured location will evolve over time.  In this case, light travel time
      is already included in the catalog location.

    .. warning::
        The astrometric calculation assumes that the body is moving in a
        straight line with its specified velocity.  The accuracy of this
        assumption depends on the body's true motion during the time it takes
        for light to travel the observer-body distance (for solar-system
        bodies) or the time since the reference epoch (for cosmic objects).

    Parameters
    ----------
    representation : `~astropy.coordinates.BaseRepresentation` or None
        A representation object or ``None`` to have no data.  Alternatively,
        use coordinate component keyword arguments, which depend on the base
        coordinate frame.
    origin : `~astropy.coordinates.SkyCoord` or low-level coordinate object.
        The coordinate which specifies both the base coordinate frame and the
        3D location of the observer.
    ref_epoch : `~astropy.time.Time`
        The reference epoch for the location of cosmic objects.  Defaults to
        J2000.0.
    obstime : `~astropy.time.Time`
        The observation time to use when the base coordinate frame does not
        have such a frame attribute (e.g., `~astropy.coordinates.ICRS`).
        Defaults to J2000.0

    Notes
    -----
    The astrometric calculation distinguishes between bodies in the solar
    system and cosmic objects by the distance from the solar-system barycenter
    (SSB).  The coordinates of bodies less than 1 light-year from the SSB are
    treated as instantaneous locations.  The coordinates of bodies greater than
    1 light-year from the SSB are treated as locations at the reference epoch.

    As these are "astrometric" coordinates rather than "apparent" coordinates,
    effects such as stellar aberration and gravitational deflection are not
    included, at least not explicitly.  Such effects may be implicitly
    included if they are accounted for in the base coordinate frame.

    ``AstrometricFrame`` is a factory class.  That is, the objects that it
    yields are *not* actually objects of class ``AstrometricFrame``.  Instead,
    distinct classes are created on-the-fly for whatever the frame class is
    of ``origin``.
    """
    # The distance threshold used to distinguish between solar-system bodies and cosmic objects
    _distance_threshold = 1*u.lyr

    def __new__(cls, *args, **kwargs):
        # We don't want to call this method if we've already set up
        # an astrometric frame for this class.
        if not (issubclass(cls, AstrometricFrame) and cls is not AstrometricFrame):
            # We get the origin argument, and handle it here.
            origin_frame = kwargs.get('origin', None)
            if origin_frame is None:
                raise TypeError("Can't initialize an AstrometricFrame without origin= keyword.")
            if hasattr(origin_frame, 'frame'):
                origin_frame = origin_frame.frame
            newcls = make_astrometric_cls(origin_frame.__class__)
            return newcls.__new__(newcls, *args, **kwargs)

        # http://stackoverflow.com/questions/19277399/why-does-object-new-work-differently-in-these-three-cases
        # See above for why this is necessary. Basically, because some child
        # may override __new__, we must override it here to never pass
        # arguments to the object.__new__ method.
        if super().__new__ is object.__new__:
            return super().__new__(cls)
        return super().__new__(cls, *args, **kwargs)

    def _frame_attrs_repr(self):
        # Overrides BaseCoordinateFrame._frame_attrs_repr() to hide the copied frame attributes
        attr_names = self.get_frame_attr_names()
        for attribute_name in self.origin.get_frame_attr_names():
            del attr_names[attribute_name]

        attr_strs = []
        for attribute_name in attr_names:
            attr = getattr(self, attribute_name)
            if hasattr(attr, "_astropy_repr_in_frame"):
                attrstr = attr._astropy_repr_in_frame()
            else:
                attrstr = str(attr)
            attr_strs.append(f"{attribute_name}={attrstr}")

        return ', '.join(attr_strs)

    def __setattr__(self, attr, value):
        if attr == "_origin":
            if not value.has_data:
                raise ValueError('The origin supplied to AstrometricFrame has no data.')

            if value.data.norm().unit is u.one:
                raise ValueError("The data for `origin` must have distance units.")

            # Pre-compute the Cartesian representation of the observer location in ICRS
            self._observer_icrs_cartesian = value.transform_to(ICRS).cartesian

            # Copy out all frame attributes from the `origin` frame attribute
            for attr2 in value.get_frame_attr_names():
                super().__setattr__("_" + attr2, getattr(value, attr2))

        super().__setattr__(attr, value)

    def as_base(self):
        """
        Returns a coordinate with the current representation and in the base
        coordinate frame.

        This method can be thought of as "removing" the
        `~astropy.coordinates.AstrometricFrame` layer.  Be aware that this
        method is not merely a coordinate transformation, because this method
        changes the location in inertial space that is being pointed to.
        """
        return self.origin.realize_frame(self.data)
