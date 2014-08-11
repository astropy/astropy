# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the coordinate frames actually implemented by astropy.
"""
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Standard library
import inspect

# Dependencies
import numpy as np

# Project
from ..extern import six
from ..utils.compat.odict import OrderedDict
from .. import units as u
from ..time import Time
from .angles import Angle
from .representation import (SphericalRepresentation, CartesianRepresentation,
                             UnitSphericalRepresentation)
from .baseframe import (BaseCoordinateFrame, frame_transform_graph, GenericFrame,
                        FrameAttribute, TimeFrameAttribute,
                        RepresentationMapping)
from .transformations import FunctionTransform, DynamicMatrixTransform


__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'AltAz']

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')


class ICRS(BaseCoordinateFrame):
    """
    A coordinate or frame in the ICRS system.

    If you're looking for "J2000" coordinates, and aren't sure if you want to
    use this or `FK5`, you probably want to use ICRS. It's more well-defined as
    a catalog coordinate and is an inertial system, and is very close (within
    tens of milliarcseconds) to J2000 equatorial.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    @staticmethod
    def _icrs_to_fk5_matrix():
        """
        B-matrix from USNO circular 179.  Used by the ICRS->FK5 transformation
        functions.
        """
        from .angles import rotation_matrix

        eta0 = -19.9 / 3600000.
        xi0 = 9.1 / 3600000.
        da0 = -22.9 / 3600000.

        m1 = rotation_matrix(-eta0, 'x')
        m2 = rotation_matrix(xi0, 'y')
        m3 = rotation_matrix(da0, 'z')

        return m1 * m2 * m3

# define this because it only needs to be computed once
ICRS._ICRS_TO_FK5_J2000_MAT = ICRS._icrs_to_fk5_matrix()


class FK5(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK5 system.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    equinox : `~astropy.time.Time`, optional, must be keyword
        The equinox of this frame.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_J2000)

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK5 based on Capitaine et
        al. 2003/IAU2006.  Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        from .earth_orientation import precession_matrix_Capitaine

        return precession_matrix_Capitaine(oldequinox, newequinox)


# Has to be defined at module level, because `transform` needs an FK5 reference
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK5)
def fk5_to_fk5(fk5coord1, fk5frame2):
    return fk5coord1._precession_matrix(fk5coord1.equinox, fk5frame2.equinox)


class FK4(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK4 system.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    equinox : astropy.time.Time, optional, must be keyword
        The equinox of this frame.
    obstime : astropy.time.Time, optional, must be keyword
        The time this frame was observed.  If None, will be the same as
        ``equinox``.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_B1950)
    obstime = TimeFrameAttribute(default=None, secondary_attribute='equinox')


@frame_transform_graph.transform(FunctionTransform, FK4, FK4)
def fk4_to_fk4(fk4coord1, fk4frame2):
    # deceptively complicated: need to transform to No E-terms FK4, precess, and
    # then come back, because precession is non-trivial with E-terms
    fnoe_w_eqx1 = fk4coord1.transform_to(FK4NoETerms(equinox=fk4coord1.equinox))
    fnoe_w_eqx2 = fnoe_w_eqx1.transform_to(FK4NoETerms(equinox=fk4frame2.equinox))
    return fnoe_w_eqx2.transform_to(fk4frame2)


class FK4NoETerms(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK4 system, but with the E-terms of aberration
    removed.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    obstime : astropy.time.Time, optional, must be keyword
        The time this frame was observed.  If None, will be the same as
        ``equinox``.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_B1950)
    obstime = TimeFrameAttribute(default=None, secondary_attribute='equinox')

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK4 using Newcomb's method.
        Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        from .earth_orientation import _precession_matrix_besselian

        return _precession_matrix_besselian(oldequinox.byear, newequinox.byear)

    @staticmethod
    def _fk4_B_matrix(obstime):
        """
        This is a correction term in the FK4 transformations because FK4 is a
        rotating system - see Murray 89 eqn 29
        """
        # Note this is *julian century*, not besselian
        T = (obstime.jyear - 1950.) / 100.
        return _B1950_TO_J2000_M + _FK4_CORR * T


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK4NoETerms)
def fk4noe_to_fk4noe(fk4necoord1, fk4neframe2):
    return fk4necoord1._precession_matrix(fk4necoord1.equinox, fk4neframe2.equinox)


class Galactic(BaseCoordinateFrame):
    """
    Galactic Coordinates.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    l : `Angle`, optional, must be keyword
        The Galactic longitude for this object (``b`` must also be given and
        ``representation`` must be None).
    b : `Angle`, optional, must be keyword
        The Galactic latitude for this object (``l`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'l'),
                      RepresentationMapping('lat', 'b')],
        'cartesian': [RepresentationMapping('x', 'w'),
                      RepresentationMapping('y', 'u'),
                      RepresentationMapping('z', 'v')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5
    # These are from Reid & Brunthaler 2004
    _ngp_J2000 = FK5(ra=192.859508*u.degree, dec=27.128336*u.degree)
    _lon0_J2000 = Angle(122.932, u.degree)
    # These are from the IAU's definition of galactic coordinates
    _ngp_B1950 = FK4(ra=192.25*u.degree, dec=27.4*u.degree)
    _lon0_B1950 = Angle(123, u.degree)


class AltAz(BaseCoordinateFrame):
    """
    A coordinate or frame in the Altitude-Azimuth system (i.e., Horizontal
    coordinates).

    .. warning::
        The AltAz class currently does not support any transformations. In a
        future version, it will support the standard IAU2000 AltAz<->ICRS
        transformations.  It is provided right now as a placeholder for storing
        as-observed horizontal coordinates.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    az : `Angle`, optional, must be keyword
        The Azimuth for this object (``alt`` must also be given and
        ``representation`` must be None).
    alt : `Angle`, optional, must be keyword
        The Altitude for this object (``az`` must also be given and
        ``representation`` must be None).
    distance : :class:`~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
    """

    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'az'),
                      RepresentationMapping('lat', 'alt')],
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_B1950)
    location = FrameAttribute(default=None)
    obstime = TimeFrameAttribute(default=None)

    def __init__(self, *args, **kwargs):
        from warnings import warn
        from astropy.utils.exceptions import AstropyWarning

        warn(AstropyWarning('The AltAz class currently does not support any '
                            'transformations.  In a future version, it will '
                            'support the standard IAU2000 AltAz<->ICRS '
                            'transformations.'))
        super(AltAz, self).__init__(*args, **kwargs)


# <--------------------------------transformations------------------------------>
# Transformations are defined here instead of in the classes themselves, because
# we need references to the various objects to give to the decorators.

# ICRS to/from FK5 -------------------------------->
@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
def icrs_to_fk5(icrscoord, fk5frame):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5frame._precession_matrix(_EQUINOX_J2000, fk5frame.equinox)
    return pmat * icrscoord._ICRS_TO_FK5_J2000_MAT


# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
def fk5_to_icrs(fk5coord, icrsframe):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    return icrsframe._ICRS_TO_FK5_J2000_MAT.T * pmat


# FK4-NO-E to/from FK4 ----------------------------->

# In the present framework, we include two coordinate classes for FK4
# coordinates - one including the E-terms of aberration (FK4), and
# one not including them (FK4NoETerms). The following functions
# implement the transformation between these two.
def fk4_e_terms(equinox):
    """
    Return the e-terms of aberation vector

    Parameters
    ----------
    equinox : Time object
        The equinox for which to compute the e-terms
    """

    from . import earth_orientation as earth

    # Constant of aberration at J2000
    k = 0.0056932

    # Eccentricity of the Earth's orbit
    e = earth.eccentricity(equinox.jd)
    e = np.radians(e)

    # Mean longitude of perigee of the solar orbit
    g = earth.mean_lon_of_perigee(equinox.jd)
    g = np.radians(g)

    # Obliquity of the ecliptic
    o = earth.obliquity(equinox.jd, algorithm=1980)
    o = np.radians(o)

    return e * k * np.sin(g), \
           -e * k * np.cos(g) * np.cos(o), \
           -e * k * np.cos(g) * np.sin(o)


@frame_transform_graph.transform(FunctionTransform, FK4, FK4NoETerms)
def fk4_to_fk4_no_e(fk4coord, fk4noeframe):
    from .representation import CartesianRepresentation, UnitSphericalRepresentation

    # Extract cartesian vector
    c = fk4coord.cartesian.xyz
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = np.asarray(fk4_e_terms(fk4coord.equinox))
    r = r - eterms_a.reshape(3, 1) + np.dot(eterms_a, r) * r

    # Find new distance (for re-normalization)
    d_new = np.sqrt(np.sum(r ** 2))

    # Renormalize
    r = r * d_orig / d_new

    subshape = c.shape[1:]
    x = r[0].reshape(subshape)
    y = r[1].reshape(subshape)
    z = r[2].reshape(subshape)

    #now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4coord.data, UnitSphericalRepresentation):
        representation = CartesianRepresentation(x=x*u.one, y=y*u.one, z=z*u.one)
        representation = representation.represent_as(UnitSphericalRepresentation)
    else:
        representation = CartesianRepresentation(x=x*c.unit, y=y*c.unit, z=z*c.unit)

    # if no obstime was given in the new frame, use the old one for consistency
    newobstime = fk4coord._obstime if fk4noeframe._obstime is None else fk4noeframe._obstime

    fk4noe = FK4NoETerms(representation, equinox=fk4coord.equinox,
                         obstime=newobstime)
    if fk4coord.equinox != fk4noeframe.equinox:
        #precession
        fk4noe = fk4noe.transform_to(fk4noeframe)
    return fk4noe


@frame_transform_graph.transform(FunctionTransform, FK4NoETerms, FK4)
def fk4_no_e_to_fk4(fk4noecoord, fk4frame):
    from .representation import CartesianRepresentation, UnitSphericalRepresentation

    #first precess, if necessary
    if fk4noecoord.equinox != fk4frame.equinox:
        fk4noe_w_fk4equinox = FK4NoETerms(equinox=fk4frame.equinox,
                                          obstime=fk4noecoord.obstime)
        fk4noecoord = fk4noecoord.transform_to(fk4noe_w_fk4equinox)

    # Extract cartesian vector
    c = fk4noecoord.cartesian.xyz
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration. Note that this depends on the equinox (not
    # the observing time/epoch) of the coordinates. See issue #1496 for a
    # discussion of this.
    eterms_a = np.asarray(fk4_e_terms(fk4noecoord.equinox))
    r0 = r.copy()
    for _ in range(10):
        r = (eterms_a.reshape(3, 1) + r0) / (1. + np.dot(eterms_a, r))

    # Find new distance (for re-normalization)
    d_new = np.sqrt(np.sum(r ** 2))

    # Renormalize
    r = r * d_orig / d_new

    subshape = c.shape[1:]
    x = r[0].reshape(subshape)
    y = r[1].reshape(subshape)
    z = r[2].reshape(subshape)

    #now re-cast into an appropriate Representation, and precess if need be
    if isinstance(fk4noecoord.data, UnitSphericalRepresentation):
        representation = CartesianRepresentation(x=x*u.one, y=y*u.one, z=z*u.one)
        representation = representation.represent_as(UnitSphericalRepresentation)
    else:
        representation = CartesianRepresentation(x=x*c.unit, y=y*c.unit, z=z*c.unit)

    return fk4frame.realize_frame(representation)

# FK5 to/from FK4 ------------------->

# B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
_B1950_TO_J2000_M = \
    np.mat([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
            [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
            [0.0048590037723143, -0.0000271702937440,  0.9999881946023742]])

_FK4_CORR = \
    np.mat([[-0.0026455262, -1.1539918689, +2.1111346190],
            [+1.1540628161, -0.0129042997, +0.0236021478],
            [-2.1112979048, -0.0056024448, +0.0102587734]]) * 1.e-6


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, FK5)
def fk4_no_e_to_fk5(fk4noecoord, fk5frame):
    # Correction terms for FK4 being a rotating system
    B = FK4NoETerms._fk4_B_matrix(fk4noecoord.obstime)

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk4noecoord._precession_matrix(fk4noecoord.equinox, _EQUINOX_B1950)
    pmat2 = fk5frame._precession_matrix(_EQUINOX_J2000, fk5frame.equinox)

    return pmat2 * B * pmat1


# This transformation can't be static because the observation date is needed.
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK4NoETerms)
def fk5_to_fk4_no_e(fk5coord, fk4noeframe):
    # Get transposed version of the rotating correction terms... so with the
    # transpose this takes us from FK5/J200 to FK4/B1950
    B = FK4NoETerms._fk4_B_matrix(fk4noeframe.obstime).T

    # construct both precession matricies - if the equinoxes are B1950 and
    # J2000, these are just identity matricies
    pmat1 = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    pmat2 = fk4noeframe._precession_matrix(_EQUINOX_B1950, fk4noeframe.equinox)

    return pmat2 * B * pmat1


# Galactic to/from FK4/FK5 ----------------------->
# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, Galactic)
def fk5_to_gal(fk5coord, galframe):
    from .angles import rotation_matrix

    #need precess to J2000 first
    pmat = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    mat1 = rotation_matrix(180 - Galactic._lon0_J2000.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_J2000.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_J2000.ra.degree, 'z')

    return mat1 * mat2 * mat3 * pmat


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK5)
def _gal_to_fk5(galcoord, fk5frame):
    return fk5_to_gal(fk5frame, galcoord).T


@frame_transform_graph.transform(DynamicMatrixTransform, FK4NoETerms, Galactic)
def fk4_to_gal(fk4coords, galframe):
    from .angles import rotation_matrix

    mat1 = rotation_matrix(180 - Galactic._lon0_B1950.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_B1950.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_B1950.ra.degree, 'z')
    matprec = fk4coords._precession_matrix(fk4coords.equinox, _EQUINOX_B1950)

    return mat1 * mat2 * mat3 * matprec


@frame_transform_graph.transform(DynamicMatrixTransform, Galactic, FK4NoETerms)
def gal_to_fk4(galcoords, fk4frame):
    return fk4_to_gal(fk4frame, galcoords).T


def _make_transform_graph_docs():
    """
    Generates a string for use with the coordinate package's docstring
    to show the available transforms and coordinate systems
    """
    from textwrap import dedent

    isclass = inspect.isclass
    coosys = [item for item in list(six.itervalues(globals()))
              if isclass(item) and issubclass(item, BaseCoordinateFrame)]
    coosys.remove(BaseCoordinateFrame)
    coosys.remove(GenericFrame)
    graphstr = frame_transform_graph.to_dot_graph(addnodes=coosys)

    docstr = """
    The diagram below shows all of the coordinate systems built into the
    `~astropy.coordinates` package, their aliases (useful for converting
    other coordinates to them using attribute-style access) and the
    pre-defined transformations between them.  The user is free to
    override any of these transformations by defining new transformations
    between these systems, but the pre-defined transformations should be
    sufficient for typical usage.

    The graph also indicates the priority for each transformation as a
    number next to the arrow.  These priorities are used to decide the
    preferred order when two transformation paths have the same number
    of steps.  These priorities are defined such that the path with a
    *smaller* total priority is favored.


    .. graphviz::

    """

    return dedent(docstr) + '    ' + graphstr.replace('\n', '\n    ')
_transform_graph_docs = _make_transform_graph_docs()
