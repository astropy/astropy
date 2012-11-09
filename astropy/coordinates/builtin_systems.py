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


_epoch_j2000 = Time('J2000', scale='utc')


#<--------------Coordinate definitions; transformations are below-------------->
@transformations.coordinate_alias('icrs')
class ICRSCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the ICRS.

    If you're looking for "J2000" coordinates, this is probably what you
    want; ICRS is better defined and is within a few microarcsec of
    J2000. The ICRS is defined in reference to this single epoch.


    Parameters
    ----------
    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to ICRSCoordinates and used as this
    coordinate.

    Attributes
    ----------
    ra : `~astropy.coordinates.angle.RA`
        The right ascension of this coordinate.
    dec : `~astropy.coordinates.angle.Dec`
        The declination of this coordinate.

    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def epoch(self):
        return _epoch_j2000


@transformations.coordinate_alias('fk5')
class FK5Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK5 system.

    Parameters
    ----------
    {params}
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to J2000.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK5Coordinates` and used as this
    coordinate.

    Attributes
    ----------
    ra : `~astropy.coordinates.angle.RA`
        The right ascension of this coordinate.
    dec : `~astropy.coordinates.angle.Dec`
        The declination of this coordinate.
    epoch : `~astropy.time.Time`
        The epoch of this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK5Coordinates, self).__init__()

        self._epoch = kwargs.pop('epoch', _epoch_j2000)

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK5Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def epoch(self):
        return self._epoch

    def precess_to(self, newepoch):
        """
        Precesses the coordinates from their current `epoch` to a new epoch and
        returns the resulting coordinate.

        Parameters
        ----------
        newepoch : `~astropy.time.Time`
            The epoch to precess these coordinates to.

        Returns
        -------
        newcoord : FK5Coordinates
            The new coordinate
        """
        from .earth_orientation import precession_matrix_Capitaine

        pmat = precession_matrix_Capitaine(self._epoch, newepoch)

        v = [self.x, self.y, self.z]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance.unit, epoch=newepoch)
        else:
            return self.__class__(x=x, y=y, z=z, epoch=newepoch)


@transformations.coordinate_alias('fk4')
class FK4Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK4 system.


    Parameters
    ----------
    {params}
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to B1950.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK4Coordinates` and used as this
    coordinate.

    Attributes
    ----------
    ra : `~astropy.coordinates.angle.RA`
        The right ascension of this coordinate.
    dec : `~astropy.coordinates.angle.Dec`
        The declination of this coordinate.
    epoch : `~astropy.time.Time`
        The epoch of this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))

    def __init__(self, *args, **kwargs):
        super(FK4Coordinates, self).__init__()

        self._epoch = kwargs.pop('epoch', Time('B1950', scale='utc'))

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK4Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} RA={1:.5f} deg, Dec={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.ra.degrees,
                          self.dec.degrees, diststr)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec

    @property
    def epoch(self):
        return self._epoch

    def precess_to(self, newepoch):
        """
        Precesses the coordinates from their current `epoch` to a new epoch.

        Parameters
        ----------
        newepoch : `~astropy.time.Time`
            The epoch to precess these coordinates to.

        Returns
        -------
        newcoord : FK4Coordinates
            The new coordinate
        """
        from .earth_orientation import _precession_matrix_besselian

        pmat = _precession_matrix_besselian(self._epoch.byear, newepoch.byear)

        v = [self.x, self.y, self.z]
        x, y, z = np.dot(pmat.A, v)

        if self.distance is not None:
            return self.__class__(x=x, y=y, z=z, unit=self.distance.unit, epoch=newepoch)
        else:
            return self.__class__(x=x, y=y, z=z, epoch=newepoch)


@transformations.coordinate_alias('galactic')
class GalacticCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates.

    Parameters
    ----------
    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `GalacticCoordinates` and
    used as this coordinate.

    Attributes
    ----------
    l : `~astropy.coordinates.angle.Angle`
        The galactic longitude of this coordinate.
    b : `~astropy.coordinates.angle.Angle`
        The galactic latitude of this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='l', latnm='b'))

    #North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
    #transformations to/from FK4/5
    _ngp_J2000 = FK5Coordinates(192.859508, 27.128336, unit=u.degree)
    _long0_J2000 = Angle(122.932, unit=u.degree)
    _ngp_B1950 = FK4Coordinates(192.25, 27.4, unit=u.degree)
    _long0_B1950 = Angle(123, unit=u.degree)

    def __init__(self, *args, **kwargs):
        super(GalacticCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.l = newcoord.l
            self.b = newcoord.b
            self._distance = newcoord._distance
        else:
            super(GalacticCoordinates, self)._initialize_latlong('l', 'b', False, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} l={1:.5f} deg, b={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.l.degrees,
                          self.b.degrees, diststr)
    @property
    def longangle(self):
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
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to J200.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `HorizontalCoordinates` and used
    as this coordinate.

    Attributes
    ----------
    az : `~astropy.coordinates.angle.RA`
        The azimuth of this coordinate.
    el : `~astropy.coordinates.angle.Dec`
        The elevation/altitude of this coordinate.
    epoch : `~astropy.time.Time`
        The epoch of this coordinate.
    """
    __doc__ = __doc__.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='az', latnm='el'))

    def __init__(self, *args, **kwargs):
        super(HorizontalCoordinates, self).__init__()

        self._epoch = kwargs.pop('epoch', _epoch_j2000)

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.az = newcoord.az
            self.el = newcoord.el
            self._distance = newcoord._distance
        else:
            super(HorizontalCoordinates, self)._initialize_latlong('az', 'el', False, args, kwargs)

    def __repr__(self):
        if self.distance is not None:
            diststr = ', Distance={0:.2g} {1!s}'.format(self.distance._value, self.distance._unit)
        else:
            diststr = ''

        msg = "<{0} el={1:.5f} deg, az={2:.5f} deg{3}>"
        return msg.format(self.__class__.__name__, self.el.degrees,
                          self.az.degrees, diststr)

    @property
    def longangle(self):
        return self.az

    @property
    def latangle(self):
        return self.el

    @property
    def epoch(self):
        return self._epoch


#<--------------------------------transformations------------------------------>
#ICRS to/from FK5
@transformations.static_transform_matrix(ICRSCoordinates, FK5Coordinates)
def icrs_to_fk5():
    """
    B-matrix from USNO circular 179
    """
    from .angles import rotation_matrix

    eta0 = -19.9 / 3600000
    xi0 = 9.1 / 3600000
    da0 = -22.9 / 3600000

    m1 = rotation_matrix(-eta0, 'x')
    m2 = rotation_matrix(xi0, 'y')
    m3 = rotation_matrix(da0, 'z')

    return m1 * m2 * m3


#can't be static because the epoch is needed
@transformations.dynamic_transform_matrix(FK5Coordinates, ICRSCoordinates)
def fk5_to_icrs(fk5c):
    from .earth_orientation import _precess_from_J2000_Capitaine

    pmat = _precess_from_J2000_Capitaine(fk5c.epoch.jyear).T

    #transpose gets epoch -> J2000
    fk5toicrsmat = icrs_to_fk5().T

    return fk5toicrsmat * pmat


#ICRS to/from FK4
# these transformations are very slightly prioritized >1 (lower priority number means
# better path) to prefer the FK5 path over FK4 when possible
#can't be static because the epoch is needed
@transformations.dynamic_transform_matrix(FK4Coordinates, ICRSCoordinates, priority=1.01)
def fk4_to_icrs(fk4c):
    from .earth_orientation import _precession_matrix_besselian

    #B1950->J2000 matrix from Murray 1989 A&A 218,325 eqn 28
    B = np.mat([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
                [0.0111814832391717,  0.9999374848933135, -0.0000271625947142],
                [0.0048590037723143, -0.0000271702937440,  0.9999881946023742]])

    if fk4c.epoch is not None and fk4c.epoch.byear != 1950:
        #not this is *julian century*, not besselian
        T = (fk4c.epoch.jyear - 1950) / 100

        #now add in correction terms for FK4 rotating system - Murray 89 eqn 29
        B[0, 0] += -2.6455262e-9 * T
        B[0, 1] += -1.1539918689e-6 * T
        B[0, 2] += 2.1111346190e-6 * T
        B[1, 0] += 1.1540628161e-6 * T
        B[1, 1] += -1.29042997e-8 * T
        B[1, 2] += 2.36021478e-8 * T
        B[2, 0] += -2.1112979048e-6 * T
        B[2, 1] += -5.6024448e-9 * T
        B[2, 2] += 1.02587734e-8 * T

        PB = _precession_matrix_besselian(fk4c.epoch.byear, 1950)

        return B * PB
    else:
        return B


#can't be static because the epoch is needed
@transformations.dynamic_transform_matrix(ICRSCoordinates, FK4Coordinates, priority=1.01)
def icrs_to_fk4(fk4c):
    from .earth_orientation import _precession_matrix_besselian

    # need inverse instead of transpose because Murray's matrix is *not* a true
    # rotation matrix
    return fk4_to_icrs(fk4c).I


#GalacticCoordinates to/from FK4/FK5
#can't be static because the epoch is needed
@transformations.dynamic_transform_matrix(FK5Coordinates, GalacticCoordinates)
def _fk5_to_gal(fk5coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precess_from_J2000_Capitaine

    # needed mainly to support inverse from galactic
    jepoch = 2000 if fk5coords.epoch is None else fk5coords.epoch.jyear

    mat1 = rotation_matrix(180 - GalacticCoordinates._long0_J2000.degrees, 'z')
    mat2 = rotation_matrix(90 - GalacticCoordinates._ngp_J2000.dec.degrees, 'y')
    mat3 = rotation_matrix(GalacticCoordinates._ngp_J2000.ra.degrees, 'z')
    #transpose gets epoch -> J2000
    matprec = _precess_from_J2000_Capitaine(jepoch).T
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(GalacticCoordinates, FK5Coordinates)
def _gal_to_fk5(galcoords):
    return _fk5_to_gal(galcoords).T


@transformations.dynamic_transform_matrix(FK4Coordinates, GalacticCoordinates, priority=1.02)
def _fk4_to_gal(fk4coords):
    from .angles import rotation_matrix
    from .earth_orientation import _precession_matrix_besselian

    # needed mainly to support inverse from galactic
    bepoch = 1950 if fk4coords.epoch is None else fk4coords.epoch.byear

    mat1 = rotation_matrix(180 - GalacticCoordinates._long0_B1950.degrees, 'z')
    mat2 = rotation_matrix(90 - GalacticCoordinates._ngp_B1950.dec.degrees, 'y')
    mat3 = rotation_matrix(GalacticCoordinates._ngp_B1950.ra.degrees, 'z')
    matprec = _precession_matrix_besselian(bepoch, 1950)
    return mat1 * mat2 * mat3 * matprec


@transformations.dynamic_transform_matrix(GalacticCoordinates, FK4Coordinates, priority=1.02)
def _gal_to_fk4(galcoords):
    return _fk4_to_gal(galcoords).T
