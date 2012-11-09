# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains the implementations of specific coordinate systems
and the conversions between them.
"""

from .angles import RA, Dec, Angle
from .coordsystems import SphericalCoordinatesBase
from .import transformations
from .. import units as u

__all__ = ['ICRSCoordinates', 'FK5Coordinates', 'FK4Coordinates',
           'GalacticCoordinates', 'HorizontalCoordinates'
          ]


#<--------------Coordinate definitions; transformations are below-------------->
@transformations.coordinate_alias('icrs')
class ICRSCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the ICRS system.

    If you're looking for "J2000" coordinates, this is probably what you want;
    ICRS is better defined and is within a few microarcsec of J2000.


    Paramaters
    ----------
    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to ICRSCoordinates and used as this
    coordinate.

    """.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        super(ICRSCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(ICRSCoordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


@transformations.coordinate_alias('fk5')
class FK5Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK5 system.

    Paramaters
    ----------
    {params}
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to J2000.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK5Coordinates` and used as this
    coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        from ..time import Time

        super(FK5Coordinates, self).__init__()

        self.epoch = kwargs.pop('epoch', Time('J2000', scale='utc'))

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK5Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


@transformations.coordinate_alias('fk4')
class FK4Coordinates(SphericalCoordinatesBase):
    """
    A coordinate in the FK4 system.


    Paramaters
    ----------
    {params}
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to B1950.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `FK4Coordinates` and used as this
    coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='ra', latnm='dec'))
    def __init__(self, *args, **kwargs):
        from ..time import Time

        super(FK4Coordinates, self).__init__()

        self.epoch = kwargs.pop('epoch', Time('B1950', scale='utc'))

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.ra = newcoord.ra
            self.dec = newcoord.dec
            self._distance = newcoord._distance
        else:
            super(FK4Coordinates, self)._initialize_latlong('ra', 'dec', True, args, kwargs)

    @property
    def longangle(self):
        return self.ra

    @property
    def latangle(self):
        return self.dec


@transformations.coordinate_alias('galactic')
class GalacticCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in Galactic Coordinates.

    Paramaters
    ----------
    {params}

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `GalacticCoordinates` and
    used as this coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='l', latnm='b'))

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

    @property
    def longangle(self):
        return self.l

    @property
    def latangle(self):
        return self.b


@transformations.coordinate_alias('horizontal')
class HorizontalCoordinates(SphericalCoordinatesBase):
    """
    A coordinate in the Horizontal or "alt/az" system.

    Paramaters
    ----------
    {params}
    epoch : `~astropy.time.Time`, optional
        The epoch for these coordinates.  Defaults to J200.

    Alternatively, a single argument that is any kind of spherical coordinate
    can be provided, and will be converted to `HorizontalCoordinates` and used
    as this coordinate.
    """.format(params=SphericalCoordinatesBase._init_docstring_param_templ.format(longnm='az', latnm='alt'))
    def __init__(self, *args, **kwargs):
        from ..time import Time

        super(HorizontalCoordinates, self).__init__()

        self.epoch = kwargs.pop('epoch', Time('J2000', scale='utc'))

        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], SphericalCoordinatesBase):
            newcoord = args[0].transform_to(self.__class__)
            self.az = newcoord.az
            self.alt = newcoord.alt
            self._distance = newcoord._distance
        else:
            super(HorizontalCoordinates, self)._initialize_latlong('az', 'alt', False, args, kwargs)

    @property
    def longangle(self):
        return self.az

    @property
    def latangle(self):
        return self.alt


#<--------------------------------transformations------------------------------>
#ICRS to/from FK5
# @CoordinateSystem.registerTransform(ICRSCoordinates,'self',transtype='smatrix')
# def _fromICRS(icrsc):
#     """
#     B-matrix from USNO circular 179
#     """
#     from ..utils import rotation_matrix

#     eta0 = -19.9/3600000
#     xi0 = 9.1/3600000
#     da0 = -22.9/3600000
#     B = rotation_matrix(-eta0,'x') *\
#         rotation_matrix(xi0,'y') *\
#         rotation_matrix(da0,'z')

#     epoch = icrsc.epoch
#     if icrsc.epoch is None:
#         return B
#     else:
#         return FK5Coordinates._precessionMatrixJ(2000,icrsc.epoch)*B

# @CoordinateSystem.registerTransform('self',ICRSCoordinates,transtype='smatrix')
# def _toICRS(fk5c):
#     return FK5Coordinates._fromICRS(fk5c).T


# #ICRS to/from FK4
# # these transformations are very slightly prioritized >1 (lower priority number means
# # better path) to prefer the FK5 path over FK4 when possible
# @CoordinateSystem.registerTransform('self',FK5Coordinates,transtype='smatrix')
#     def _toFK5(fk4c):
#         from ..obstools import epoch_to_jd,jd_to_epoch


#         #B1950->J2000 matrix from Murray 1989 A&A 218,325
#         B = np.mat([[0.9999256794956877,-0.0111814832204662,-0.0048590038153592],
#                     [0.0111814832391717,0.9999374848933135,-0.0000271625947142],
#                     [0.0048590037723143,-0.0000271702937440,0.9999881946023742]])

#         if fk4c.epoch is not None and fk4c.epoch != 1950:
#             jd = epoch_to_jd(fk4c.epoch,False)
#             jepoch = jd_to_epoch(jd)
#             T = (jepoch - 1950)/100

#             #now add in correction terms for FK4 rotating system
#             B[0,0] += -2.6455262e-9*T
#             B[0,1] += -1.1539918689e-6*T
#             B[0,2] += 2.1111346190e-6*T
#             B[1,0] += 1.1540628161e-6*T
#             B[1,1] += -1.29042997e-8*T
#             B[1,2] += 2.36021478e-8*T
#             B[2,0] += -2.1112979048e-6*T
#             B[2,1] += -5.6024448e-9*T
#             B[2,2] += 1.02587734e-8*T

#             PB = FK4Coordinates._precessionMatrixB(fk4c.epoch,1950)

#             return B*PB
#         else:
#             return B

#     @CoordinateSystem.registerTransform(FK5Coordinates,'self',transtype='smatrix')
#     def _fromFK5(fk5c):
#         #need inverse because Murray's matrix is *not* a true rotation matrix
#         return FK4Coordinates._toFK5(fk5c).I


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


@transformations.dynamic_transform_matrix(FK4Coordinates, GalacticCoordinates)
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


@transformations.dynamic_transform_matrix(GalacticCoordinates, FK4Coordinates)
def _gal_to_fk4(galcoords):
    return _fk4_to_gal(galcoords).T
