# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import numpy as np

from astropy import units as u
from astropy.utils.state import ScienceState
from astropy.utils.decorators import format_doc
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.coordinates.angles import Angle
from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product, matrix_transpose
from astropy.coordinates import representation as r
from astropy.coordinates.baseframe import (BaseCoordinateFrame,
                                           frame_transform_graph,
                                           base_doc)
from astropy.coordinates.attributes import (CoordinateAttribute,
                                            QuantityAttribute,
                                            DifferentialAttribute)
from astropy.coordinates.transformations import AffineTransform
from astropy.coordinates.errors import ConvertError

from .icrs import ICRS

__all__ = ['Galactocentric']


# Measured by minimizing the difference between a plane of coordinates along
#   l=0, b=[-90,90] and the Galactocentric x-z plane
# This is not used directly, but accessed via `get_roll0`.  We define it here to
# prevent having to create new Angle objects every time `get_roll0` is called.
_ROLL0 = Angle(58.5986320306*u.degree)


class galactocentric_frame_defaults(ScienceState):
    """This class controls the global setting of default values for the frame
    attributes in the `~astropy.coordinates.Galactocentric` frame, which may be
    updated in future versions of ``astropy``. Note that when using
    `~astropy.coordinates.Galactocentric`, changing values here will not affect
    any attributes that are set explicitly by passing values in to the
    `~astropy.coordinates.Galactocentric` initializer. Modifying these defaults
    will only affect the frame attribute values when using the frame as, e.g.,
    ``Galactocentric`` or ``Galactocentric()`` with no explicit arguments.

    This class controls the parameter settings by specifying a string name,
    which can be one of:

    - 'pre-v4.0': The current default value, which sets the default frame
      attribute values to their original (pre-astropy-v4.0) values.
    - 'v4.0': The attribute values as updated in Astropy version 4.0.
    - 'latest': An alias of the most recent parameter set (currently: 'v4.0')

    See :ref:`astropy-coordinates-galactocentric-defaults` for more information.

    Examples
    --------
    The default `~astropy.coordinates.Galactocentric` frame parameters can be
    modified globally::

        >>> from astropy.coordinates import galactocentric_frame_defaults
        >>> _ = galactocentric_frame_defaults.set('v4.0')
        >>> Galactocentric() # doctest: +FLOAT_CMP
        <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            (266.4051, -28.936175)>, galcen_distance=8.122 kpc, galcen_v_sun=(12.9, 245.6, 7.78) km / s, z_sun=20.8 pc, roll=0.0 deg)>
        >>> _ = galactocentric_frame_defaults.set('pre-v4.0')
        >>> Galactocentric() # doctest: +FLOAT_CMP
        <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            (266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

    The default parameters can also be updated by using this class as a context
    manager::

        >>> with galactocentric_frame_defaults.set('pre-v4.0'):
        ...     print(Galactocentric()) # doctest: +FLOAT_CMP
        <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            (266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

    Again, changing the default parameter values will not affect frame
    attributes that are explicitly specified::

        >>> import astropy.units as u
        >>> with galactocentric_frame_defaults.set('pre-v4.0'):
        ...     print(Galactocentric(galcen_distance=8.0*u.kpc)) # doctest: +FLOAT_CMP
        <Galactocentric Frame (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            (266.4051, -28.936175)>, galcen_distance=8.0 kpc, galcen_v_sun=(11.1, 232.24, 7.25) km / s, z_sun=27.0 pc, roll=0.0 deg)>

    """

    # the default is to use the original definition of this frame
    # TODO: change this to 'latest' in v4.1?
    _value = 'pre-v4.0'
    _references = None

    @staticmethod
    def get_solar_params_from_string(arg):
        """Return Galactocentric solar parameters given string names for the
        parameter sets.
        """
        # Resolve the meaning of 'latest': The latest parameter set is from v4.0
        # - update this as newer parameter choices are added
        if arg == 'latest':
            arg = 'v4.0'

        params = dict()
        references = dict()

        # Currently, all versions use the same sky position for Sgr A*:
        params['galcen_coord'] = ICRS(ra=266.4051*u.degree,
                                      dec=-28.936175*u.degree)
        references['galcen_coord'] = \
            'https://ui.adsabs.harvard.edu/abs/2004ApJ...616..872R'

        # The roll angle is the same for both frames:
        params['roll'] = 0 * u.deg

        if arg == 'pre-v4.0':
            params['galcen_distance'] = 8.3 * u.kpc
            references['galcen_distance'] = \
                'https://ui.adsabs.harvard.edu/#abs/2009ApJ...692.1075G'

            params['galcen_v_sun'] = r.CartesianDifferential([11.1,
                                                              220+12.24,
                                                              7.25]*u.km/u.s)
            references['galcen_v_sun'] = \
                ['https://ui.adsabs.harvard.edu/#abs/2010MNRAS.403.1829S',
                 'https://ui.adsabs.harvard.edu/#abs/2015ApJS..216...29B']

            params['z_sun'] = 27.0 * u.pc
            references['z_sun'] = \
                'https://ui.adsabs.harvard.edu/#abs/2001ApJ...553..184C'

        elif arg == 'v4.0':
            params['galcen_distance'] = 8.122 * u.kpc
            references['galcen_distance'] = \
                'https://ui.adsabs.harvard.edu/abs/2018A%26A...615L..15G'

            params['galcen_v_sun'] = r.CartesianDifferential([12.9,
                                                              245.6,
                                                              7.78]*u.km/u.s)
            references['galcen_v_sun'] = \
                ['https://ui.adsabs.harvard.edu/abs/2018RNAAS...2..210D',
                 'https://ui.adsabs.harvard.edu/abs/2018A%26A...615L..15G',
                 'https://ui.adsabs.harvard.edu/abs/2004ApJ...616..872R']

            params['z_sun'] = 20.8 * u.pc
            references['z_sun'] = \
                'https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.1417B'

        else:
            raise ValueError('Invalid string input to retrieve solar '
                             'parameters for Galactocentric frame: "{}"'
                             .format(arg))

        return params, references

    @classmethod
    def validate(cls, value):
        if isinstance(value, str):
            params, refs = cls.get_solar_params_from_string(value)
            cls._references = refs
            return params

        elif isinstance(value, dict):
            return value

        elif isinstance(value, Galactocentric):
            # turn the frame instance into a dict of frame attributes
            attrs = dict()
            for k in value.frame_attributes:
                attrs[k] = getattr(value, k)
            cls._references = value.frame_attribute_references()
            return attrs

        else:
            raise ValueError("Invalid input to retrieve solar parameters for "
                             "Galactocentric frame: input must be a string, "
                             "dict, or Galactocentric instance")


doc_components = """
    x : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`x` position component.
    y : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`y` position component.
    z : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`z` position component.

    v_x : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_x` velocity component.
    v_y : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_y` velocity component.
    v_z : `~astropy.units.Quantity`, optional
        Cartesian, Galactocentric :math:`v_z` velocity component.
"""

doc_footer = """
    Other parameters
    ----------------
    galcen_coord : `ICRS`, optional, must be keyword
        The ICRS coordinates of the Galactic center.
    galcen_distance : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the sun to the Galactic center.
    galcen_v_sun : `~astropy.coordinates.representation.CartesianDifferential`, optional, must be keyword
        The velocity of the sun *in the Galactocentric frame* as Cartesian
        velocity components.
    z_sun : `~astropy.units.Quantity`, optional, must be keyword
        The distance from the sun to the Galactic midplane.
    roll : `~astropy.coordinates.Angle`, optional, must be keyword
        The angle to rotate about the final x-axis, relative to the
        orientation for Galactic. For example, if this roll angle is 0,
        the final x-z plane will align with the Galactic coordinates x-z
        plane. Unless you really know what this means, you probably should
        not change this!

    Examples
    --------
    .. testsetup::

        >>> from astropy.coordinates import galactocentric_frame_defaults
        >>> _ = galactocentric_frame_defaults.set('pre-v4.0')

    To transform to the Galactocentric frame with the default
    frame attributes, pass the uninstantiated class name to the
    ``transform_to()`` method of a coordinate frame or
    `~astropy.coordinates.SkyCoord` object::

        >>> import astropy.units as u
        >>> import astropy.coordinates as coord
        >>> c = coord.ICRS(ra=[158.3122, 24.5] * u.degree,
        ...                dec=[-17.3, 81.52] * u.degree,
        ...                distance=[11.5, 24.12] * u.kpc)
        >>> c.transform_to(coord.Galactocentric) # doctest: +FLOAT_CMP
        <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            ( 266.4051, -28.936175)>, galcen_distance=8.3 kpc, galcen_v_sun=( 11.1,  232.24,  7.25) km / s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [( -9.6083819 ,  -9.40062188,  6.52056066),
             (-21.28302307,  18.76334013,  7.84693855)]>

    To specify a custom set of parameters, you have to include extra keyword
    arguments when initializing the Galactocentric frame object::

        >>> c.transform_to(coord.Galactocentric(galcen_distance=8.1*u.kpc)) # doctest: +FLOAT_CMP
        <Galactocentric Coordinate (galcen_coord=<ICRS Coordinate: (ra, dec) in deg
            ( 266.4051, -28.936175)>, galcen_distance=8.1 kpc, galcen_v_sun=( 11.1,  232.24,  7.25) km / s, z_sun=27.0 pc, roll=0.0 deg): (x, y, z) in kpc
            [( -9.40785924,  -9.40062188,  6.52066574),
             (-21.08239383,  18.76334013,  7.84798135)]>

    Similarly, transforming from the Galactocentric frame to another coordinate frame::

        >>> c = coord.Galactocentric(x=[-8.3, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[0.027, 24.12] * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(  86.22349059, 28.83894138,  4.39157788e-05),
             ( 289.66802652, 49.88763881,  8.59640735e+01)]>

    Or, with custom specification of the Galactic center::

        >>> c = coord.Galactocentric(x=[-8.0, 4.5] * u.kpc,
        ...                          y=[0., 81.52] * u.kpc,
        ...                          z=[21.0, 24120.0] * u.pc,
        ...                          z_sun=21 * u.pc, galcen_distance=8. * u.kpc)
        >>> c.transform_to(coord.ICRS) # doctest: +FLOAT_CMP
        <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
            [(  86.2585249 ,  28.85773187,  2.75625475e-05),
             ( 289.77285255,  50.06290457,  8.59216010e+01)]>


    # TODO: this needs to be audited and removed if we change the default
    # defaults (yes) for v4.1. This is a hack to get around the fact that
    # the doctests here actually change the *class* state:

        >>> galactocentric_frame_defaults._value = 'pre-v4.0'

"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class Galactocentric(BaseCoordinateFrame):
    r"""
    A coordinate or frame in the Galactocentric system.

    This frame allows specifying the Sun-Galactic center distance, the height of
    the Sun above the Galactic midplane, and the solar motion relative to the
    Galactic center. However, as there is no modern standard definition of a
    Galactocentric reference frame, it is important to pay attention to the
    default values used in this class if precision is important in your code.
    The default values of the parameters of this frame are taken from the
    original definition of the frame in 2014. As such, the defaults are somewhat
    out of date relative to recent measurements made possible by, e.g., Gaia.
    The defaults can, however, be changed at runtime by setting the parameter
    set name in `~astropy.coordinates.galactocentric_frame_defaults`.

    The current default parameter set is ``"pre-v4.0"``, indicating that the
    parameters were adopted before ``astropy`` version 4.0. A regularly-updated
    parameter set can instead be used by setting
    ``galactocentric_frame_defaults.set ('latest')``, and other parameter set
    names may be added in future versions. To find out the scientific papers
    that the current default parameters are derived from, use
    ``galcen.frame_attribute_references`` (where ``galcen`` is an instance of
    this frame), which will update even if the default parameter set is changed.

    The position of the Sun is assumed to be on the x axis of the final,
    right-handed system. That is, the x axis points from the position of
    the Sun projected to the Galactic midplane to the Galactic center --
    roughly towards :math:`(l,b) = (0^\circ,0^\circ)`. For the default
    transformation (:math:`{\rm roll}=0^\circ`), the y axis points roughly
    towards Galactic longitude :math:`l=90^\circ`, and the z axis points
    roughly towards the North Galactic Pole (:math:`b=90^\circ`).

    For a more detailed look at the math behind this transformation, see
    the document :ref:`coordinates-galactocentric`.

    The frame attributes are listed under **Other Parameters**.
    """

    default_representation = r.CartesianRepresentation
    default_differential = r.CartesianDifferential

    # frame attributes
    galcen_coord = CoordinateAttribute(frame=ICRS)
    galcen_distance = QuantityAttribute(unit=u.kpc)

    galcen_v_sun = DifferentialAttribute(
        allowed_classes=[r.CartesianDifferential])

    z_sun = QuantityAttribute(unit=u.pc)
    roll = QuantityAttribute(unit=u.deg)

    def __init__(self, *args, **kwargs):
        # Set default frame attribute values based on the ScienceState instance
        # for the solar parameters defined above
        default_params = galactocentric_frame_defaults.get()
        self.frame_attribute_references = \
            galactocentric_frame_defaults._references.copy()

        warn = False
        for k in default_params:
            if k in kwargs:
                # If a frame attribute is set by the user, remove its reference
                self.frame_attribute_references.pop(k, None)

            else:
                # If a parameter is read from the defaults, we might want to
                # warn the user that the defaults will change (see below)
                warn = True

            # Keep the frame attribute if it is set by the user, otherwise use
            # the default value
            kwargs[k] = kwargs.get(k, default_params[k])

        # If the frame defaults have not been updated with the ScienceState
        # class, and the user uses any default parameter value, raise a
        # deprecation warning to inform them that the defaults will change in
        # the future:
        if galactocentric_frame_defaults._value == 'pre-v4.0' and warn:
            docs_link = 'http://docs.astropy.org/en/latest/coordinates/galactocentric.html'
            warnings.warn('In v4.1 and later versions, the Galactocentric '
                          'frame will adopt default parameters that may update '
                          'with time. An updated default parameter set is '
                          'already available through the '
                          'astropy.coordinates.galactocentric_frame_defaults '
                          'ScienceState object, as described in but the '
                          'default is currently still set to the pre-v4.0 '
                          'parameter defaults. The safest way to guard against '
                          'changing default parameters in the future is to '
                          'either (1) specify all Galactocentric frame '
                          'attributes explicitly when using the frame, '
                          'or (2) set the galactocentric_frame_defaults '
                          f'parameter set name explicitly. See {docs_link} for more '
                          'information.',
                          AstropyDeprecationWarning)

        super().__init__(*args, **kwargs)

    @classmethod
    def get_roll0(cls):
        """
        The additional roll angle (about the final x axis) necessary to align
        the final z axis to match the Galactic yz-plane.  Setting the ``roll``
        frame attribute to  -this method's return value removes this rotation,
        allowing the use of the `Galactocentric` frame in more general contexts.
        """
        # note that the actual value is defined at the module level.  We make at
        # a property here because this module isn't actually part of the public
        # API, so it's better for it to be accessable from Galactocentric
        return _ROLL0

# ICRS to/from Galactocentric ----------------------->


def get_matrix_vectors(galactocentric_frame, inverse=False):
    """
    Use the ``inverse`` argument to get the inverse transformation, matrix and
    offsets to go from Galactocentric to ICRS.
    """
    # shorthand
    gcf = galactocentric_frame

    # rotation matrix to align x(ICRS) with the vector to the Galactic center
    mat1 = rotation_matrix(-gcf.galcen_coord.dec, 'y')
    mat2 = rotation_matrix(gcf.galcen_coord.ra, 'z')
    # extra roll away from the Galactic x-z plane
    mat0 = rotation_matrix(gcf.get_roll0() - gcf.roll, 'x')

    # construct transformation matrix and use it
    R = matrix_product(mat0, mat1, mat2)

    # Now need to translate by Sun-Galactic center distance around x' and
    # rotate about y' to account for tilt due to Sun's height above the plane
    translation = r.CartesianRepresentation(gcf.galcen_distance * [1., 0., 0.])
    z_d = gcf.z_sun / gcf.galcen_distance
    H = rotation_matrix(-np.arcsin(z_d), 'y')

    # compute total matrices
    A = matrix_product(H, R)

    # Now we re-align the translation vector to account for the Sun's height
    # above the midplane
    offset = -translation.transform(H)

    if inverse:
        # the inverse of a rotation matrix is a transpose, which is much faster
        #   and more stable to compute
        A = matrix_transpose(A)
        offset = (-offset).transform(A)
        offset_v = r.CartesianDifferential.from_cartesian(
            (-gcf.galcen_v_sun).to_cartesian().transform(A))
        offset = offset.with_differentials(offset_v)

    else:
        offset = offset.with_differentials(gcf.galcen_v_sun)

    return A, offset


def _check_coord_repr_diff_types(c):
    if isinstance(c.data, r.UnitSphericalRepresentation):
        raise ConvertError("Transforming to/from a Galactocentric frame "
                           "requires a 3D coordinate, e.g. (angle, angle, "
                           "distance) or (x, y, z).")

    if ('s' in c.data.differentials and
            isinstance(c.data.differentials['s'],
                       (r.UnitSphericalDifferential,
                        r.UnitSphericalCosLatDifferential,
                        r.RadialDifferential))):
        raise ConvertError("Transforming to/from a Galactocentric frame "
                           "requires a 3D velocity, e.g., proper motion "
                           "components and radial velocity.")


@frame_transform_graph.transform(AffineTransform, ICRS, Galactocentric)
def icrs_to_galactocentric(icrs_coord, galactocentric_frame):
    _check_coord_repr_diff_types(icrs_coord)
    return get_matrix_vectors(galactocentric_frame)


@frame_transform_graph.transform(AffineTransform, Galactocentric, ICRS)
def galactocentric_to_icrs(galactocentric_coord, icrs_frame):
    _check_coord_repr_diff_types(galactocentric_coord)
    return get_matrix_vectors(galactocentric_coord, inverse=True)
