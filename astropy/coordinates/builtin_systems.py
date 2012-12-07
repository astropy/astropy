# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the implementations of specific coordinate systems
and the conversions between them.
"""

import numpy as np

from .angles import Angle
from .coordsystems import SphericalCoordinatesBase
from ..time import Time
from . import transformations
from .. import units as u

__all__ = ['ICRSCoordinates', 'FK5Coordinates', 'FK4Coordinates',
           'GalacticCoordinates', 'HorizontalCoordinates'
          ]


_EQUINOX_J2000 = Time('J2000', scale='utc')


#<--------------Coordinate definitions; transformations are below-------------->
@transformations.coordinate_alias('icrs')
class ICRSCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the ICRS.

    If you're looking for "J2000" coordinates, and aren't sure if you
    want to use this or `FK5Coordinates`, you probably want to use ICRS.
    It's more well-defined as a catalog coordinate and is an inertial
    system.


    Parameters
    ----------
    {params}
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to ICRSCoordinates and used as this
    coordinate.

    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()

        self._obstime = kwargs.pop('obstime', None)

        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(ICRSCoordinates, self)._initialize_latlon('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def lonangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def equinox(self):
        return _EQUINOX_J2000

    @property
    def obstime(self):
        if self._obstime is None:
            return self.equinox
        else:
            return self._obstime


@transformations.coordinate_alias('fk5')
class FK5Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK5 system.

    Parameters
    ----------
    {params}
    equinox : `~astropy.time.Time`, optional
        The equinox for these coordinates.  Defaults to J2000.
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK5Coordinates` and used as this
    coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK5Coordinates, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_J2000)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK5Coordinates, self)._initialize_latlon('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def lonangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def equinox(self):
        return self._equinox

    @property
    def obstime(self):
        if self._obstime is None:
            return self.equinox
        else:
            return self._obstime

    def precess_to(self, newequinox):
        """
        Precesses the coordinates from their current `equinox` to a new equinox and
        returns the resulting coordinate.

        Parameters
        ----------
        newequinox : `~astropy.time.Time`
            The equinox to precess these coordinates to.

        Returns
        -------
        newcoord : FK5Coordinates
            The new coordinate
        """
        from .earth_orientation import precession_matrix_Capitaine

        pmat = precession_matrix_Capitaine(self._equinox, newequinox)

        v = [self.x, self.y, self.z]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance._unit, equinox=newequinox)
        else:
            return self.__class__(x=x, y=y, z=z, equinox=newequinox)


@transformations.coordinate_alias('fk4')
class FK4Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK4 system.


    Parameters
    ----------
    {params}
    equinox : `~astropy.time.Time`, optional
        The equinox for these coordinates.  Defaults to B1950.
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK4Coordinates` and used as this
    coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK4Coordinates, self).__init__()

        self._equinox = kwargs.pop('equinox', Time('B1950', scale='utc'))
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK4Coordinates, self)._initialize_latlon('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def lonangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def equinox(self):
        return self._equinox

    @property
    def obstime(self):
        if self._obstime is None:
            return self.equinox
        else:
            return self._obstime

    def precess_to(self, newequinox):
        """
        Precesses the coordinates from their current `equinox` to a new equinox.

        Parameters
        ----------
        newequinox : `~astropy.time.Time`
            The equinox to precess these coordinates to.

        Returns
        -------
        newcoord : FK4Coordinates
            The new coordinate
        """
        from .earth_orientation import _precession_matrix_besselian

        pmat = _precession_matrix_besselian(self._equinox.byear, newequinox.byear)

        v = [self.x, self.y, self.z]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance._unit, equinox=newequinox)
        else:
            return self.__class__(x=x, y=y, z=z, equinox=newequinox)


@transformations.coordinate_alias('galactic')
class GalacticCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates.

    Parameters
    ----------
    {params}
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `GalacticCoordinates` and
    used as this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='l', latnm='b'))

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5
    _ngp_J2000 = FK5Coordinates(192.859508, 27.128336, unit=(u.degree, u.degree))
    _lon0_J2000 = Angle(122.932, unit=u.degree)
    _ngp_B1950 = FK4Coordinates(192.25, 27.4, unit=(u.degree, u.degree))
    _lon0_B1950 = Angle(123, unit=u.degree)

    def __init__(self, *args, **kwargs):
        super(GalacticCoordinates, self).__init__()

        self._obstime = kwargs.pop('obstime', None)

        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.l = newcoord.l
            self.b = newcoord.b
            self._distance = newcoord._distance
        else:
            super(GalacticCoordinates, self)._initialize_latlon('l', 'b', False, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} l={1:.5f} deg, b={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.l.degrees,
                          self.b.degrees, diststr)

    @property
    def lonangle(self):
        return self.l

    @property
    def latangle(self):
        return self.b


@transformations.coordinate_alias('horizontal')
class HorizontalCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the Horizontal or "az/el" system.

    Parameters
    ----------
    {params}
    equinox : `~astropy.time.Time`, optional
        The equinox for these coordinates.  Defaults to J2000.
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `HorizontalCoordinates` and used
    as this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='az', latnm='el'))

    def __init__(self, *args, **kwargs):
        super(HorizontalCoordinates, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_J2000)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.az = newcoord.az
            self.el = newcoord.el
            self._distance = newcoord._distance
        else:
            super(HorizontalCoordinates, self)._initialize_latlon('az', 'el', False, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} az={1:.5f} deg, el={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.az.degrees,
                          self.el.degrees, diststr)

    @property
    def lonangle(self):
        return self.az

    @property
    def latangle(self):
        return self.el

    @property
    def equinox(self):
        return self._equinox

    @property
    def obstime(self):
        if self._obstime is None:
            return self.equinox
        else:
            return self._obstime


#<--------------------------------transformations------------------------------>
# ICRS to/from FK5
@transformations.static_transform_matrix(ICRSCoordinates, FK5Coordinates)
def icrs_to_fk5():
    """
    B-matrix from USNO circular 179
    """
    from .angles import rotation_matrix

    eta0 = -19.9 / 3600000.
    xi0 = 9.1 / 3600000.
    da0 = -22.9 / 3600000.

    m1 = rotation_matrix(-eta0, 'x')
    m2 = rotation_matrix(xi0, 'y')
    m3 = rotation_matrix(da0, 'z')

    return m1 * m2 * m3


# can't be static because the equinox is needed
@transformations.dynamic_transform_matrix(FK5Coordinates, ICRSCoordinates)
def fk5_to_icrs(fk5c):
    from .earth_orientation import _precess_from_J2000_Capitaine

    pmat = _precess_from_J2000_Capitaine(fk5c.equinox.jyear).T

    # transpose gets equinox -> J2000
    fk5toicrsmat = icrs_to_fk5().T

    return fk5toicrsmat * pmat


# ICRS to/from FK4
# these transformations are very slightly prioritized >1 (lower priority number means
# better path) to prefer the FK5 path over FK4 when possible
# can't be static because the equinox is needed
@transformations.dynamic_transform_matrix(FK4Coordinates, ICRSCoordinates, priority=1.01)
def fk4_to_icrs(fk4c, bypassprec=False):
    # bypassprec is here to make icrs_to_fk4 work more easily - it skips
    # adding the equinox precession part
    from .earth_orientation import _precession_matrix_besselian

    # B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
    B = np.mat([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
                [0.0111814832391717, 0.9999374848933135, -0.0000271625947142],
                [0.0048590037723143, -0.0000271702937440, 0.9999881946023742]])

    if fk4c.obstime.byear != 1950:
        # note this is *julian century*, not besselian
        T = (fk4c.obstime.jyear - 1950) / 100

        # add in correction terms for FK4 rotating system - Murray 89 eqn 29
        B[0, 0] += -2.6455262e-9 * T
        B[0, 1] += -1.1539918689e-6 * T
        B[0, 2] += 2.1111346190e-6 * T
        B[1, 0] += 1.1540628161e-6 * T
        B[1, 1] += -1.29042997e-8 * T
        B[1, 2] += 2.36021478e-8 * T
        B[2, 0] += -2.1112979048e-6 * T
        B[2, 1] += -5.6024448e-9 * T
        B[2, 2] += 1.02587734e-8 * T

    if bypassprec or fk4c.equinox.byear == 1950:
        return B  # if B1950, no precession is needed - B takes us to J2000
    else:
        return B * _precession_matrix_besselian(fk4c.equinox.byear, 1950)


# can't be static because the equinox is needed
@transformations.dynamic_transform_matrix(ICRSCoordinates, FK4Coordinates, priority=1.01)
def icrs_to_fk4(icrs):
    # need inverse instead of transpose because Murray's matrix is *not* a true
    # rotation matrix
    # also note that icrs works here
    return fk4_to_icrs(icrs, bypassprec=True).I


# GalacticCoordinates to/from FK4/FK5
# can't be static because the equinox is needed
@transformations.dynamic_transform_matrix(FK5Coordinates, GalacticCoordinates)
def _fk5_to_gal(fk5coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precess_from_J2000_Capitaine

    # needed mainly to support inverse from galactic
    jequinox = 2000 if fk5coords.equinox is None else fk5coords.equinox.jyear

    mat1 = rotation_matrix(180 - GalacticCoordinates._lon0_J2000.degrees, 'z')
    mat2 = rotation_matrix(90 - GalacticCoordinates._ngp_J2000.dec.degrees, 'y')
    mat3 = rotation_matrix(GalacticCoordinates._ngp_J2000.ra.degrees, 'z')
    # transpose gets equinox -> J2000
    matprec = _precess_from_J2000_Capitaine(jequinox).T
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(GalacticCoordinates, FK5Coordinates)
def _gal_to_fk5(galcoords):
    return _fk5_to_gal(galcoords).T


@transformations.dynamic_transform_matrix(FK4Coordinates, GalacticCoordinates, priority=1.02)
def _fk4_to_gal(fk4coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precession_matrix_besselian

    # needed mainly to support inverse from galactic
    bequinox = 1950 if fk4coords.equinox is None else fk4coords.equinox.byear

    mat1 = rotation_matrix(180 - GalacticCoordinates._lon0_B1950.degrees, 'z')
    mat2 = rotation_matrix(90 - GalacticCoordinates._ngp_B1950.dec.degrees, 'y')
    mat3 = rotation_matrix(GalacticCoordinates._ngp_B1950.ra.degrees, 'z')
    matprec = _precession_matrix_besselian(bequinox, 1950)
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(GalacticCoordinates, FK4Coordinates, priority=1.02)
def _gal_to_fk4(galcoords):
    return _fk4_to_gal(galcoords).T


def _make_transform_graph_docs():
    """
    Generates a string for use with the coordinate package's docstring
    to show the available transforms and coordinate systems
    """
    from inspect import isclass
    from textwrap import dedent

    from .transformations import master_transform_graph

    coosys = [item for item in globals().values()
              if isclass(item) and issubclass(item, SphericalCoordinatesBase)]
    coosys.remove(SphericalCoordinatesBase)
    graphstr = master_transform_graph.to_dot_graph(addnodes=coosys)

    docstr = """
    The diagram below shows all of the coordinate systems built into the
    `~astropy.coordinates` package, their aliases (usable for converting
    other coordinates to them using attribute-style access) and the
    pre-defined transformations between them.  The user is free to
    override any of these transformations by defining new trasnformation
    between these systems, but the pre-defined transformations should be
    sufficient for typical usage.

    The graph also indicates the priority for each transformation as a
    number next to the arrow.  These priorities are used to decide the
    preferred order when two trasnformation paths have the same number
    of steps.  These priorities are defined such that path with a
    *smaller* total priority are favored over larger.
    E.g., the path from `ICRSCoordinates` to `GalacticCoordinates` goes
    through `FK5Coordinates` because the total path length is 2 instead
    of 2.03.


    .. graphviz::

    """

    return dedent(docstr) + '    ' + graphstr.replace('\n', '\n    ')
_transform_graph_docs = _make_transform_graph_docs()
