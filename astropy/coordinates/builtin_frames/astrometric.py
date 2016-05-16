# -*- coding: utf-8 -*-

import numpy as np

from ... import units as u
from .icrs import ICRS
from .altaz import AltAz
from ..transformations import DynamicMatrixTransform, FunctionTransform
from ..baseframe import (CoordinateAttribute, QuantityFrameAttribute,
                         frame_transform_graph, RepresentationMapping)
from ..angles import rotation_matrix
from ...utils.compat import namedtuple_asdict
from ...extern import six


_astrometric_cache = {}


def make_astrometric_cls(framecls):
    """
    Create a new class that is the Astrometric frame for a specific class of
    origin frame. If such a class has already been created for this frame, the
    same class will be returned.

    The resulting frame class will be subtly different from the base class in
    that its spherical component names will be d<lat> and d<lon>.  E.g., for
    ICRS the astrometric frame had components ``dra`` and ``ddec`` instead of
    ``ra`` and ``dec``.

    Parameters
    ----------
    framecls : coordinate frame class (i.e., subclass of `~astropy.coordinates.BaseCoordinateFrame`)
        The class to create the Astrometric frame of.

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
        This metaclass renames the class to be "Astrometric<framecls>" and also
        adjusts the frame specific representation info to have the "d" in front
        of spherical names
        """
        def __new__(cls, name, bases, members):
            newname = name + framecls.__name__
            res = super(AstrometricMeta, cls).__new__(cls, newname, bases, members)

            # now go through all the component names and make any spherical
            # lat/lon names be "d<lon>"/"d<lat>"

            lists_done = []
            for nm, component_list in res._frame_specific_representation_info.items():
                if nm in ('spherical', 'unitspherical'):
                    gotlatlon = []
                    for i, comp in enumerate(component_list):
                        if component_list in lists_done:
                            # we need this because sometimes the component_
                            # list's are the exact *same* object for both
                            # spherical and unitspherical.  So looping then adds
                            # the 'd' *twice*.  This hack bypasses that.
                            continue

                        if comp.reprname in ('lon', 'lat'):
                            dct = namedtuple_asdict(comp)
                            dct['framename'] = 'd' + dct['framename']
                            component_list[i] = type(comp)(**dct)
                            gotlatlon.append(comp.reprname)
                    if 'lon' not in gotlatlon:
                        rmlon = RepresentationMapping('lon', 'dlon', 'recommended')
                        component_list.insert(0, rmlon)
                    if 'lat' not in gotlatlon:
                        rmlat = RepresentationMapping('lat', 'dlat', 'recommended')
                        component_list.insert(0, rmlat)
                    lists_done.append(component_list)

            return res

    class Astrometric(six.with_metaclass(AstrometricMeta, framecls)):
        """
        A frame which is relative to some position on the sky. Useful for
        calculating offsets and dithers in the frame of the sky.

        Parameters
        ----------
        representation : `BaseRepresentation` or None
            A representation object or None to have no data (or use the other keywords)
        origin : `SkyCoord` or low-level coordinate object.
            the coordinate which specifiy the origin of this frame.
        rotation : Quantity with angle units
            The final rotation of the frame about the ``origin``. The sign of
            the rotation is the left-hand rule.  That is, an object at a
            particular position angle in the un-rotated system will be sent to
            the positive latitude (z) direction in the final frame.

        """
        origin = CoordinateAttribute(default=None, frame=framecls)
        rotation = QuantityFrameAttribute(default=0, unit=u.deg)

    @frame_transform_graph.transform(FunctionTransform, Astrometric, Astrometric)
    def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
        """Transform between two astrometric frames."""

        # If both frames have on-sky positions, then the transform should happen relative to both origins.
        if (from_astrometric_coord.origin is not None) and (to_astrometric_frame.origin is not None):
            return from_astrometric_coord.transform_to(framecls).transform_to(to_astrometric_frame)

        # Otherwise, the transform occurs just by setting the new origin.
        return to_astrometric_frame.realize_frame(from_astrometric_coord.cartesian)

    @frame_transform_graph.transform(DynamicMatrixTransform, framecls, Astrometric)
    def icrs_to_astrometric(reference_frame, astrometric_frame):
        """Convert an ICRS coordinate to an Astrometric frame."""

        # Define rotation matricies along the position angle vector, and
        # relative to the origin.
        origin = astrometric_frame.origin.spherical
        mat1 = rotation_matrix(-astrometric_frame.rotation, 'x')
        mat2 = rotation_matrix(-origin.lat, 'y')
        mat3 = rotation_matrix(origin.lon, 'z')
        R = mat1 * mat2 * mat3
        return R

    @frame_transform_graph.transform(DynamicMatrixTransform, Astrometric, framecls)
    def astrometric_to_icrs(astrometric_coord, reference_frame):
        """Convert an Astrometric frame coordinate to an ICRS"""

        # use the forward transform, but just invert it
        R = icrs_to_astrometric(reference_frame, astrometric_coord)
        return R.T  # this is the inverse because R is a rotation matrix

    _astrometric_cache[framecls] = Astrometric
    return Astrometric

AstrometricICRS = make_astrometric_cls(ICRS)
AstrometricAltAz = make_astrometric_cls(AltAz)
