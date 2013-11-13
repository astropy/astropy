# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the implementations of specific coordinate systems
and the conversions between them.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import inspect
import textwrap

import numpy as np

from ..extern import six
from .angles import Angle
from .coordsystems import SphericalCoordinatesBase
from ..time import Time
from . import transformations
from .. import units as u

__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'AltAz']

# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')


#<--------------Coordinate definitions; transformations are below-------------->
@transformations.coordinate_alias('icrs')
class ICRS(SphericalCoordinatesBase):
    """
    A coordinate in the ICRS.

    If you're looking for "J2000" coordinates, and aren't sure if you
    want to use this or `FK5`, you probably want to use ICRS.
    It's more well-defined as a catalog coordinate and is an inertial
    system.


    Parameters
    ----------
    {params}
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to ICRS and used as this
    coordinate.

    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(ICRS, self).__init__()

        self._obstime = kwargs.pop('obstime', None)

        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(ICRS, self)._initialize_latlon('ra', 'dec', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'RA'
    _repr_lat_name = 'Dec'

    # Default format for to_string
    _default_string_style = 'hmsdms'

    @property
    def ra(self):
        return self._lonangle

    @property
    def dec(self):
        return self._latangle

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
class FK5(SphericalCoordinatesBase):
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
    can be provided, and will be converted to `FK5` and used as this
    coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK5, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_J2000)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')
        if not self._equinox.isscalar:
            raise ValueError('specified equinox is not a scalar')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(FK5, self)._initialize_latlon('ra', 'dec', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'RA'
    _repr_lat_name = 'Dec'

    # Default format for to_string
    _default_string_style = 'hmsdms'

    @property
    def ra(self):
        return self._lonangle

    @property
    def dec(self):
        return self._latangle

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
        newcoord : FK5
            The new coordinate
        """
        from .earth_orientation import precession_matrix_Capitaine

        pmat = precession_matrix_Capitaine(self._equinox, newequinox)

        v = [self.x.value, self.y.value, self.z.value]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance.unit, equinox=newequinox)
        else:
            return self.__class__(x=x, y=y, z=z, equinox=newequinox)


@transformations.coordinate_alias('fk4')
class FK4(SphericalCoordinatesBase):
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
    can be provided, and will be converted to `FK4` and used as this
    coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK4, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_B1950)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')
        if not self._equinox.isscalar:
            raise ValueError('specified equinox is not a scalar')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(FK4, self)._initialize_latlon('ra', 'dec', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'RA'
    _repr_lat_name = 'Dec'

    # Default format for to_string
    _default_string_style = 'hmsdms'

    @property
    def ra(self):
        return self._lonangle

    @property
    def dec(self):
        return self._latangle

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
        newcoord : FK4
            The new coordinate
        """
        from .earth_orientation import _precession_matrix_besselian

        pmat = _precession_matrix_besselian(self._equinox.byear, newequinox.byear)

        v = [self.x.value, self.y.value, self.z.value]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance.unit, equinox=newequinox)
        else:
            return self.__class__(x=x, y=y, z=z, equinox=newequinox)


@transformations.coordinate_alias('fk4_no_e')
class FK4NoETerms(SphericalCoordinatesBase):
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
    can be provided, and will be converted to `FK4NoETerms` and used as this
    coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK4NoETerms, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_B1950)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')
        if not self._equinox.isscalar:
            raise ValueError('specified equinox is not a scalar')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(FK4NoETerms, self)._initialize_latlon('ra', 'dec', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'RA'
    _repr_lat_name = 'Dec'

    # Default format for to_string
    _default_string_style = 'hmsdms'

    @property
    def ra(self):
        return self._lonangle

    @property
    def dec(self):
        return self._latangle

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
        newcoord : FK4NoETerms
            The new coordinate
        """
        from .earth_orientation import _precession_matrix_besselian

        pmat = _precession_matrix_besselian(self._equinox.byear, newequinox.byear)

        v = [self.x.value, self.y.value, self.z.value]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance.unit, equinox=newequinox)
        else:
            return self.__class__(x=x, y=y, z=z, equinox=newequinox)


@transformations.coordinate_alias('galactic')
class Galactic(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates.

    .. note::

        Transformations from Galactic Coordinates to other systems are
        not well-defined because of ambiguities in the definition of
        Galactic Coordinates. See
        `Lie et al. 2011 <http://dx.doi.org/10.1051/0004-6361/201014961>`
        for more details on this.  Here, we use the
        `Reid & Brunthaler 2004 <http://dx.doi.org/10.1086/424960>`
        definition for converting to/from FK5, and assume the IAU
        definition applies for converting to FK4 *without* e-terms.

    Parameters
    ----------
    {params}
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `Galactic` and
    used as this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='l', latnm='b'))

    # North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    # transformations to/from FK4/5

    # These are from Reid & Brunthaler 2004
    _ngp_J2000 = FK5(192.859508, 27.128336, unit=(u.degree, u.degree))
    _lon0_J2000 = Angle(122.932, unit=u.degree)
    # These are from the IAU's definition of galactic coordinates
    _ngp_B1950 = FK4(192.25, 27.4, unit=(u.degree, u.degree))
    _lon0_B1950 = Angle(123, unit=u.degree)

    def __init__(self, *args, **kwargs):
        super(Galactic, self).__init__()

        self._obstime = kwargs.pop('obstime', None)

        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(Galactic, self)._initialize_latlon('l', 'b', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'l'
    _repr_lat_name = 'b'

    # Default format for to_string
    _default_string_style = 'dms'

    @property
    def l(self):
        return self._lonangle

    @property
    def b(self):
        return self._latangle


@transformations.coordinate_alias('altaz')
class AltAz(SphericalCoordinatesBase):
    """
    A coordinate in the altitude/azimuth or "horizontal" system.

    Parameters
    ----------
    {params}
    equinox : `~astropy.time.Time`, optional
        The equinox for these coordinates.  Defaults to J2000.
    obstime : `~astropy.time.Time` or None
        The time of observation for this coordinate.  If None, it will be taken
        to be the same as the `equinox`.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `AltAz` and used
    as this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(lonnm='az', latnm='el'))

    def __init__(self, *args, **kwargs):
        super(AltAz, self).__init__()

        self._equinox = kwargs.pop('equinox', _EQUINOX_J2000)
        self._obstime = kwargs.pop('obstime', None)

        if not isinstance(self._equinox, Time):
            raise TypeError('specified equinox is not a Time object')
        if self._obstime is not None and not isinstance(self._obstime, Time):
            raise TypeError('specified obstime is not None or a Time object')
        if not self._equinox.isscalar:
            raise ValueError('specified equinox is not a scalar')

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self._latangle = newcoord._latangle
            self._lonangle = newcoord._lonangle
            self._distance = newcoord._distance
        else:
            super(AltAz, self)._initialize_latlon('az', 'alt', args, kwargs)

    #strings used for making __repr__ work
    _repr_lon_name = 'az'
    _repr_lat_name = 'alt'

    # Default format for to_string
    _default_string_style = 'dms'

    @property
    def az(self):
        return self._lonangle

    @property
    def alt(self):
        return self._latangle

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
@transformations.static_transform_matrix(ICRS, FK5)
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
@transformations.dynamic_transform_matrix(FK5, ICRS)
def fk5_to_icrs(fk5c):
    from .earth_orientation import _precess_from_J2000_Capitaine

    pmat = _precess_from_J2000_Capitaine(fk5c.equinox.jyear).T

    # transpose gets equinox -> J2000
    fk5toicrsmat = icrs_to_fk5().T

    return fk5toicrsmat * pmat

# FK4-NO-E to/from FK4

# In the present framework, we include two coordinate classes for FK4
# coordinates - one including the E-terms of aberration (FK4), and
# one not including them (FK4NoETerms). In the following functions,
# we describe the transformation between these two.

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


@transformations.transform_function(FK4, FK4NoETerms, priority=1)
def fk4_to_fk4_no_e(fk4c):

    # Extract cartesian vector
    c = fk4c.cartesian
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration
    eterms_a = np.asarray(fk4_e_terms(fk4c.equinox))
    r = r - eterms_a.reshape(3, 1) + np.dot(eterms_a, r) * r

    # Find new distance (for re-normalization)
    d_new = np.sqrt(np.sum(r ** 2))

    # Renormalize
    r = r * d_orig / d_new

    subshape = c.shape[1:]
    x = r[0].reshape(subshape)
    y = r[1].reshape(subshape)
    z = r[2].reshape(subshape)

    newunit = None if fk4c.distance is None else fk4c.distance.unit
    return FK4NoETerms(x=x, y=y, z=z, unit=newunit, equinox=fk4c.equinox)


@transformations.transform_function(FK4NoETerms, FK4, priority=1)
def fk4_no_e_to_fk4(fk4c):

    # Extract cartesian vector
    c = fk4c.cartesian
    r = np.asarray(c.reshape((3, c.size // 3)))

    # Find distance (for re-normalization)
    d_orig = np.sqrt(np.sum(r ** 2))

    # Apply E-terms of aberration
    eterms_a = np.asarray(fk4_e_terms(fk4c.equinox))
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

    newunit = None if fk4c.distance is None else fk4c.distance.unit
    return FK4(x=x, y=y, z=z, unit=newunit, equinox=fk4c.equinox)

# FK5 to/from FK4

# These transformations are very slightly prioritized >1 (lower priority number means
# better path) to prefer the FK5 path over FK4 when possible
#can't be static because the equinox is needed

# B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
B1950_TO_J2000_M = \
    np.mat([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
            [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
            [0.0048590037723143, -0.0000271702937440,  0.9999881946023742]])

FK4_CORR = \
    np.mat([[-0.0026455262, -1.1539918689, +2.1111346190],
            [+1.1540628161, -0.0129042997, +0.0236021478],
            [-2.1112979048, -0.0056024448, +0.0102587734]]) * 1.e-6

# This transformation can't be static because the observation date is needed.
@transformations.dynamic_transform_matrix(FK4NoETerms, FK5, priority=1)
def fk4_no_e_to_fk5(fk4c, skip_precession=False):

    # Add in correction terms for FK4 rotating system - Murray 89 eqn 29
    # Note this is *julian century*, not besselian
    T = (fk4c.obstime.jyear - 1950.) / 100.

    B = B1950_TO_J2000_M + FK4_CORR * T

    # If equinox is not B1950, need to precess first
    if skip_precession or fk4c.equinox == _EQUINOX_B1950:
        return B
    else:
        from .earth_orientation import _precession_matrix_besselian
        return B * _precession_matrix_besselian(fk4c.equinox.byear, 1950)

# This transformation can't be static because the observation date is needed.
@transformations.dynamic_transform_matrix(FK5, FK4NoETerms, priority=1)
def fk5_to_fk4_no_e(fk5c):

    # Get transposed matrix from FK4 -> FK5 assuming equinox B1950 -> J2000
    B = fk4_no_e_to_fk5(fk5c, skip_precession=True).T

    # If equinox is not B1950, need to precess first
    if fk5c.equinox == _EQUINOX_J2000:
        return B
    else:
        from .earth_orientation import precession_matrix_Capitaine
        return B * precession_matrix_Capitaine(fk5c.equinox, _EQUINOX_J2000)

# Galactic to/from FK4/FK5
# can't be static because the equinox is needed
@transformations.dynamic_transform_matrix(FK5, Galactic)
def _fk5_to_gal(fk5coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precess_from_J2000_Capitaine

    # needed mainly to support inverse from galactic
    jequinox = 2000 if fk5coords.equinox is None else fk5coords.equinox.jyear

    mat1 = rotation_matrix(180 - Galactic._lon0_J2000.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_J2000.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_J2000.ra.degree, 'z')
    # transpose gets equinox -> J2000
    matprec = _precess_from_J2000_Capitaine(jequinox).T
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(Galactic, FK5)
def _gal_to_fk5(galcoords):
    return _fk5_to_gal(galcoords).T


@transformations.dynamic_transform_matrix(FK4NoETerms, Galactic, priority=1)
def _fk4_to_gal(fk4coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precession_matrix_besselian

    # needed mainly to support inverse from galactic
    bequinox = 1950 if fk4coords.equinox is None else fk4coords.equinox.byear

    mat1 = rotation_matrix(180 - Galactic._lon0_B1950.degree, 'z')
    mat2 = rotation_matrix(90 - Galactic._ngp_B1950.dec.degree, 'y')
    mat3 = rotation_matrix(Galactic._ngp_B1950.ra.degree, 'z')
    matprec = _precession_matrix_besselian(bequinox, 1950)
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(Galactic, FK4NoETerms, priority=1)
def _gal_to_fk4(galcoords):
    return _fk4_to_gal(galcoords).T


def _make_transform_graph_docs():
    """
    Generates a string for use with the coordinate package's docstring
    to show the available transforms and coordinate systems
    """

    from .transformations import master_transform_graph

    isclass = inspect.isclass
    coosys = [item for item in list(six.itervalues(globals()))
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
    *smaller* total priority are favored over one with larger.


    .. graphviz::

    """

    return textwrap.dedent(docstr) + '    ' + graphstr.replace('\n', '\n    ')
_transform_graph_docs = _make_transform_graph_docs()
