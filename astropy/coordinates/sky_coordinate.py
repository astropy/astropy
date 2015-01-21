from __future__ import (absolute_import, division, print_function, unicode_literals)

import re
import collections
import warnings

import numpy as np

from ..utils.compat.misc import override__dir__
from ..extern import six
from ..extern.six.moves import zip
from ..units import Unit, IrreducibleUnit
from .. import units as u
from ..wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from ..utils.exceptions import AstropyDeprecationWarning

from .distances import Distance
from .baseframe import BaseCoordinateFrame, frame_transform_graph, GenericFrame, _get_repr_cls
from .builtin_frames import ICRS
from .representation import (BaseRepresentation, SphericalRepresentation,
                             UnitSphericalRepresentation)

__all__ = ['SkyCoord']

PLUS_MINUS_RE = re.compile(r'(\+|\-)')
J_PREFIXED_RA_DEC_RE = re.compile(
    r"""J                              # J prefix
    ([0-9]{6,7}\.?[0-9]{0,2})          # RA as HHMMSS.ss or DDDMMSS.ss, optional decimal digits
    ([\+\-][0-9]{6}\.?[0-9]{0,2})\s*$  # Dec as DDMMSS.ss, optional decimal digits
    """, re.VERBOSE)


# Define a convenience mapping.  This is used like a module constants
# but is actually dynamically evaluated.
def FRAME_ATTR_NAMES_SET():
    """Set of all possible frame-specific attributes"""
    out = set()
    for frame_cls in frame_transform_graph.frame_set:
        for attr in frame_cls.get_frame_attr_names().keys():
            out.add(attr)
    return out


class SkyCoord(object):
    """High-level object providing a flexible interface for celestial coordinate
    representation, manipulation, and transformation between systems.

    The `SkyCoord` class accepts a wide variety of inputs for initialization. At
    a minimum these must provide one or more celestial coordinate values with
    unambiguous units.  Inputs may be scalars or lists/tuples/arrays, yielding
    scalar or array coordinates (can be checked via ``SkyCoord.isscalar``).
    Typically one also specifies the coordinate frame, though this is not
    required. The general pattern for spherical representations is::

      SkyCoord(COORD, [FRAME], keyword_args ...)
      SkyCoord(LON, LAT, [FRAME], keyword_args ...)
      SkyCoord(LON, LAT, [DISTANCE], frame=FRAME, unit=UNIT, keyword_args ...)
      SkyCoord([FRAME], <lon_attr>=LON, <lat_attr>=LAT, keyword_args ...)

    It is also possible to input coordinate values in other representations
    such as cartesian or cylindrical.  In this case one includes the keyword
    argument ``representation='cartesian'`` (for example) along with data in
    ``x``, ``y``, and ``z``.

    Examples
    --------
    The examples below illustrate common ways of initializing a `SkyCoord`
    object.  For a complete description of the allowed syntax see the
    full coordinates documentation.  First some imports::

      >>> from astropy.coordinates import SkyCoord  # High-level coordinates
      >>> from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
      >>> from astropy.coordinates import Angle, Latitude, Longitude  # Angles
      >>> import astropy.units as u

    The coordinate values and frame specification can now be provided using
    positional and keyword arguments::

      >>> c = SkyCoord(10, 20, unit="deg")  # defaults to ICRS frame
      >>> c = SkyCoord([1, 2, 3], [-30, 45, 8], "icrs", unit="deg")  # 3 coords

      >>> coords = ["1:12:43.2 +1:12:43", "1 12 43.2 +1 12 43"]
      >>> c = SkyCoord(coords, FK4, unit=(u.deg, u.hourangle), obstime="J1992.21")

      >>> c = SkyCoord("1h12m43.2s +1d12m43s", Galactic)  # Units from string
      >>> c = SkyCoord("galactic", l="1h12m43.2s", b="+1d12m43s")

      >>> ra = Longitude([1, 2, 3], unit=u.deg)  # Could also use Angle
      >>> dec = np.array([4.5, 5.2, 6.3]) * u.deg  # Astropy Quantity
      >>> c = SkyCoord(ra, dec, frame='icrs')
      >>> c = SkyCoord(ICRS, ra=ra, dec=dec, obstime='2001-01-02T12:34:56')

      >>> c = FK4(1 * u.deg, 2 * u.deg)  # Uses defaults for obstime, equinox
      >>> c = SkyCoord(c, obstime='J2010.11', equinox='B1965')  # Override defaults

      >>> c = SkyCoord(w=0, u=1, v=2, unit='kpc', frame='galactic', representation='cartesian')

      >>> c = SkyCoord([ICRS(ra=1*u.deg, dec=2*u.deg), ICRS(ra=3*u.deg, dec=4*u.deg)])

    As shown, the frame can be a `~astropy.coordinates.BaseCoordinateFrame`
    class or the corresponding string alias.  The frame classes that are built in
    to astropy are `ICRS`, `FK5`, `FK4`, `FK4NoETerms`, and `Galactic`.
    The string aliases are simply lower-case versions of the class name, and
    allow for creating a `SkyCoord` object and transforming frames without
    explicitly importing the frame classes.

    Parameters
    ----------
    frame : `~astropy.coordinates.BaseCoordinateFrame` class or string, optional
        Type of coordinate frame this `SkyCoord` should represent. Defaults to
        to ICRS if not given or given as None.
    unit : `~astropy.units.Unit`, string, or tuple of :class:`~astropy.units.Unit` or str, optional
        Units for supplied ``LON`` and ``LAT`` values, respectively.  If
        only one unit is supplied then it applies to both ``LON`` and
        ``LAT``.
    obstime : valid `~astropy.time.Time` initializer, optional
        Time of observation
    equinox : valid `~astropy.time.Time` initializer, optional
        Coordinate frame equinox
    representation : str or Representation class
        Specifies the representation, e.g. 'spherical', 'cartesian', or
        'cylindrical'.  This affects the positional args and other keyword args
        which must correspond to the given representation.
    **keyword_args
        Other keyword arguments as applicable for user-defined coordinate frames.
        Common options include:

        ra, dec : valid `~astropy.coordinates.Angle` initializer, optional
            RA and Dec for frames where ``ra`` and ``dec`` are keys in the
            frame's ``representation_component_names``, including `ICRS`,
            `FK5`, `FK4`, and `FK4NoETerms`.
        l, b : valid `~astropy.coordinates.Angle` initializer, optional
            Galactic ``l`` and ``b`` for for frames where ``l`` and ``b`` are
            keys in the frame's ``representation_component_names``, including
            the `Galactic` frame.
        x, y, z : float or `~astropy.units.Quantity`, optional
            Cartesian coordinates values
        w, u, v : float or `~astropy.units.Quantity`, optional
            Cartesian coordinates values for the Galactic frame.
    """


    # Declare that SkyCoord can be used as a Table column by defining the
    # attribute where column attributes will be stored.
    _astropy_column_attrs = None

    def __init__(self, *args, **kwargs):

        # Parse the args and kwargs to assemble a sanitized and validated
        # kwargs dict for initializing attributes for this object and for
        # creating the internal self._sky_coord_frame object
        args = list(args)  # Make it mutable
        kwargs = self._parse_inputs(args, kwargs)

        # Set internal versions of object state attributes
        for attr in FRAME_ATTR_NAMES_SET():
            setattr(self, '_' + attr, kwargs[attr])

        frame = kwargs['frame']
        coord_kwargs = {}
        if 'representation' in kwargs:
            coord_kwargs['representation'] = _get_repr_cls(kwargs['representation'])
        for attr, value in kwargs.items():
            if value is not None and (attr in frame.representation_component_names
                                      or attr in frame.get_frame_attr_names()):
                coord_kwargs[attr] = value

        # Finally make the internal coordinate object.
        self._sky_coord_frame = frame.__class__(**coord_kwargs)

        if not self._sky_coord_frame.has_data:
            raise ValueError('Cannot create a SkyCoord without data')

    @property
    def frame(self):
        return self._sky_coord_frame

    @property
    def representation(self):
        return self.frame.representation

    @representation.setter
    def representation(self, value):
        self.frame.representation = value

    def __len__(self):
        return len(self.frame)

    def __nonzero__(self):  # Py 2.x
        return self.frame.__nonzero__()

    def __bool__(self):  # Py 3.x
        return self.frame.__bool__()

    def __getitem__(self, item):
        self_frame = self._sky_coord_frame
        try:
            # First turn `self` into a mockup of the thing we want - we can copy
            # this to get all the right attributes
            self._sky_coord_frame = self_frame[item]
            return SkyCoord(self)
        finally:
            # now put back the right frame in self
            self._sky_coord_frame = self_frame

    def _parse_inputs(self, args, kwargs):
        """
        Assemble a validated and sanitized keyword args dict for instantiating a
        SkyCoord and coordinate object from the provided `args`, and `kwargs`.
        """
        valid_kwargs = {}

        # Put the SkyCoord attributes like frame, equinox, obstime, location
        # into valid_kwargs dict.  `Frame` could come from args or kwargs, so
        # set valid_kwargs['frame'] accordingly.  The others must be specified
        # by keyword args or else get a None default.  Pop them off of kwargs
        # in the process.
        frame = valid_kwargs['frame'] = _get_frame(args, kwargs)
        if 'representation' in kwargs:
            valid_kwargs['representation'] = _get_repr_cls(kwargs.pop('representation'))

        for attr in FRAME_ATTR_NAMES_SET():
            valid_kwargs[attr] = kwargs.pop(attr, None)

        # Get units
        units = _get_units(args, kwargs)

        # Grab any frame-specific attr names like `ra` or `l` or `distance` from kwargs
        # and migrate to valid_kwargs.
        valid_kwargs.update(_get_representation_attrs(frame, units, kwargs))

        # Error if anything is still left in kwargs
        if kwargs:
            raise ValueError('Unrecognized keyword argument(s) {0}'
                             .format(', '.join("'{0}'".format(key) for key in kwargs)))

        # Finally deal with the unnamed args.  This figures out what the arg[0] is
        # and returns a dict with appropriate key/values for initializing frame class.
        if args:
            if len(args) == 1:
                # One arg which must be a coordinate.  In this case
                # coord_kwargs will contain keys like 'ra', 'dec', 'distance'
                # along with any frame attributes like equinox or obstime which
                # were explicitly specified in the coordinate object (i.e. non-default).
                coord_kwargs = _parse_coordinate_arg(args[0], frame, units, kwargs)

            elif len(args) <= 3:
                frame_attr_names = frame.representation_component_names.keys()
                repr_attr_names = frame.representation_component_names.values()
                coord_kwargs = {}
                for arg, frame_attr_name, repr_attr_name, unit in zip(args, frame_attr_names,
                                                                      repr_attr_names, units):
                    attr_class = frame.representation.attr_classes[repr_attr_name]
                    coord_kwargs[frame_attr_name] = attr_class(arg, unit=unit)

            else:
                raise ValueError('Must supply no more than three positional arguments, got {}'
                                 .format(len(args)))

            # Copy the coord_kwargs into the final valid_kwargs dict.  For each
            # of the coord_kwargs ensure that there is no conflict with a value
            # specified by the user in the original kwargs.
            for attr, coord_value in coord_kwargs.items():
                if (attr in valid_kwargs
                        and valid_kwargs[attr] is not None
                        and np.any(valid_kwargs[attr] != coord_value)):
                    raise ValueError("Coordinate attribute '{0}'={1!r} conflicts with "
                                     "keyword argument '{0}'={2!r}"
                                     .format(attr, coord_value, valid_kwargs[attr]))
                valid_kwargs[attr] = coord_value

        return valid_kwargs

    def transform_to(self, frame):
        """
        Transform this coordinate to a new frame.

        The frame attributes (e.g. equinox or obstime) for the returned object
        depend on the corresponding attributes of SkyCoord object and the
        supplied ``frame``, with the following precedence:

        1. Non-default value in the supplied frame
        2. Non-default value in the SkyCoord instance
        3. Default value in the supplied frame

        Parameters
        ----------
        frame : str or `BaseCoordinateFrame` class / instance or `SkyCoord` instance
            The frame to transform this coordinate into.

        Returns
        -------
        coord : `SkyCoord`
            A new object with this coordinate represented in the `frame` frame.

        Raises
        ------
        ValueError
            If there is no possible transformation route.
        """
        from astropy.coordinates.errors import ConvertError

        frame_kwargs = {}

        # Frame name (string) or frame class?  Coerce into an instance.
        try:
            frame = _get_frame_class(frame)()
        except:
            pass

        if isinstance(frame, SkyCoord):
            frame = frame.frame  # Change to underlying coord frame instance

        if isinstance(frame, BaseCoordinateFrame):
            new_frame_cls = frame.__class__

            # Set the keyword args for making a new frame instance for the
            # transform.  Frame attributes track whether they were explicitly
            # set by user or are just reflecting default values.  Precedence:
            # 1. Non-default value in the supplied frame instance
            # 2. Non-default value in the self instance
            # 3. Default value in the supplied frame instance
            for attr in FRAME_ATTR_NAMES_SET():
                self_val = getattr(self, attr, None)
                frame_val = getattr(frame, attr, None)
                if frame_val is not None and not frame.is_frame_attr_default(attr):
                    frame_kwargs[attr] = frame_val
                elif self_val is not None and not self.is_frame_attr_default(attr):
                    frame_kwargs[attr] = self_val
                elif frame_val is not None:
                    frame_kwargs[attr] = frame_val
        else:
            raise ValueError('Transform `frame` must be a frame name, class, or instance')

        # Get the composite transform to the new frame
        trans = frame_transform_graph.get_transform(self.frame.__class__, new_frame_cls)
        if trans is None:
            raise ConvertError('Cannot transform from {0} to {1}'
                               .format(self.frame.__class__, new_frame_cls))

        # Make a generic frame which will accept all the frame kwargs that
        # are provided and allow for transforming through intermediate frames
        # which may require one or more of those kwargs.
        generic_frame = GenericFrame(frame_kwargs)

        # Do the transformation, returning a coordinate frame of the desired
        # final type (not generic).
        new_coord = trans(self.frame, generic_frame)

        # Finally make the new SkyCoord object from the `new_coord` and
        # remaining frame_kwargs that are not frame_attributes in `new_coord`.
        # We could remove overlaps here, but the init code is set up to accept
        # overlaps as long as the values are identical (which they must be).
        return self.__class__(new_coord, **frame_kwargs)

    def __getattr__(self, attr):
        """
        Overrides getattr to return coordinates that this can be transformed
        to, based on the alias attr in the master transform graph.
        """
        if '_sky_coord_frame' in self.__dict__:
            if self.frame.name == attr:
                return self  # Should this be a deepcopy of self?

            # Anything in the set of all possible frame_attr_names is handled
            # here. If the attr is relevant for the current frame then delegate
            # to self.frame otherwise get it from self._<attr>.
            if attr in FRAME_ATTR_NAMES_SET():
                if attr in self.frame.get_frame_attr_names():
                    return getattr(self.frame, attr)
                else:
                    return getattr(self, '_' + attr)

            # Some attributes might not fall in the above category but still
            # are available through self._sky_coord_frame.
            if not attr.startswith('_') and hasattr(self._sky_coord_frame, attr):
                return getattr(self._sky_coord_frame, attr)

            # Try to interpret as a new frame for transforming.
            frame_cls = frame_transform_graph.lookup_name(attr)
            if frame_cls is not None and self.frame.is_transformable_to(frame_cls):
                return self.transform_to(attr)

        # Fail
        raise AttributeError("'{0}' object has no attribute '{1}'"
                             .format(self.__class__.__name__, attr))

    def __setattr__(self, attr, val):
        # This is to make anything available through __getattr__ immutable
        if '_sky_coord_frame' in self.__dict__:
            if self.frame.name == attr:
                raise AttributeError("'{0}' is immutable".format(attr))

            if (attr in FRAME_ATTR_NAMES_SET() or
                (not attr.startswith('_') and
                 hasattr(self._sky_coord_frame, attr))):
                setattr(self._sky_coord_frame, attr, val)

            frame_cls = frame_transform_graph.lookup_name(attr)
            if frame_cls is not None and self.frame.is_transformable_to(frame_cls):
                raise AttributeError("'{0}' is immutable".format(attr))

        # Otherwise, do the standard Python attribute setting
        super(SkyCoord, self).__setattr__(attr, val)

    @override__dir__
    def __dir__(self):
        """
        Override the builtin `dir` behavior to include:
        - Transforms available by aliases
        - Attribute / methods of the underlying self.frame object
        """

        # determine the aliases that this can be transformed to.
        dir_values = set()
        for name in frame_transform_graph.get_names():
            frame_cls = frame_transform_graph.lookup_name(name)
            if self.frame.is_transformable_to(frame_cls):
                dir_values.add(name)

        # Add public attributes of self.frame
        dir_values.update(set(attr for attr in dir(self.frame) if not attr.startswith('_')))

        # Add all possible frame attributes
        dir_values.update(FRAME_ATTR_NAMES_SET())

        return dir_values

    def __str__(self):
        return str(self.to_string())

    def __repr__(self):
        clsnm = self.__class__.__name__
        coonm = self.frame.__class__.__name__

        crepr = repr(self.frame)
        frameattrs = ''
        if crepr.find('):') != -1:
            frameattrs = ': '+crepr[crepr.index('(')+1:crepr.index('):')]

        s = '<{clsnm} ({coonm}{frameattrs})'.format(**locals())
        return s + crepr[crepr.index(':'):]

    def to_string(self, style='decimal', **kwargs):
        """
        A string representation of the coordinates.

        The default styles definitions are::

          'decimal': 'lat': {'decimal': True, 'unit': "deg"}
                     'lon': {'decimal': True, 'unit': "deg"}
          'dms': 'lat': {'unit': "deg"}
                 'lon': {'unit': "deg"}
          'hmsdms': 'lat': {'alwayssign': True, 'pad': True, 'unit': "deg"}
                    'lon': {'pad': True, 'unit': "hour"}

        See :meth:`~astropy.coordinates.Angle.to_string` for details and
        keyword arguments (the two angles forming the coordinates are are
        both :class:`~astropy.coordinates.Angle` instances). Keyword
        arguments have precedence over the style defaults and are passed
        to :meth:`~astropy.coordinates.Angle.to_string`.

        Parameters
        ----------
        style : {'hmsdms', 'dms', 'decimal'}
            The formatting specification to use. These encode the three most
            common ways to represent coordinates. The default is `decimal`.
        kwargs
            Keyword args passed to :meth:`~astropy.coordinates.Angle.to_string`.
        """

        sph_coord = self.frame.represent_as(SphericalRepresentation)

        styles = {'hmsdms': {'lonargs': {'unit': u.hour, 'pad': True},
                             'latargs': {'unit': u.degree, 'pad': True, 'alwayssign': True}},
                  'dms': {'lonargs': {'unit': u.degree},
                          'latargs': {'unit': u.degree}},
                  'decimal': {'lonargs': {'unit': u.degree, 'decimal': True},
                              'latargs': {'unit': u.degree, 'decimal': True}}
                  }

        lonargs = {}
        latargs = {}

        if style in styles:
            lonargs.update(styles[style]['lonargs'])
            latargs.update(styles[style]['latargs'])
        else:
            raise ValueError('Invalid style.  Valid options are: {0}'.format(",".join(styles)))

        lonargs.update(kwargs)
        latargs.update(kwargs)

        if np.isscalar(sph_coord.lon.value):
            coord_string = (sph_coord.lon.to_string(**lonargs)
                            + " " +
                            sph_coord.lat.to_string(**latargs))
        else:
            coord_string = []
            for lonangle, latangle in zip(sph_coord.lon, sph_coord.lat):
                coord_string += [(lonangle.to_string(**lonargs)
                                 + " " +
                                 latangle.to_string(**latargs))]

        return coord_string

    # High-level convinience methods
    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Angle`
            The on-sky separation between this and the ``other`` coordinate.

        Notes
        -----
        The separation is calculated using the Vincenty formula, which
        is stable at all locations, including poles and antipodes [1]_.

        .. [1] http://en.wikipedia.org/wiki/Great-circle_distance

        """
        from . import Angle
        from .angle_utilities import angular_separation

        if isinstance(other, SkyCoord):
            self_in_other_system = self.transform_to(other.frame)
        elif isinstance(other, BaseCoordinateFrame) and other.has_data:
            # it's a frame
            self_in_other_system = self.transform_to(other)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        lon1 = self_in_other_system.spherical.lon
        lat1 = self_in_other_system.spherical.lat
        lon2 = other.spherical.lon
        lat2 = other.spherical.lat

        # Get the separation as a Quantity, convert to Angle in degrees
        sep = angular_separation(lon1, lat1, lon2, lat2)
        return Angle(sep, unit=u.degree)

    def separation_3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        Parameters
        ----------
        other : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate to get the separation to.

        Returns
        -------
        sep : `~astropy.coordinates.Distance`
            The real-space distance between these two coordinates.

        Raises
        ------
        ValueError
            If this or the other coordinate do not have distances.
        """

        if isinstance(other, SkyCoord):
            self_in_other_system = self.transform_to(other.frame)
        elif isinstance(other, BaseCoordinateFrame) and other.has_data:
            # it's a frame
            self_in_other_system = self.transform_to(other)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        if self.data.__class__ == UnitSphericalRepresentation:
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')
        if other.data.__class__ == UnitSphericalRepresentation:
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        dx = self_in_other_system.cartesian.x - other.cartesian.x
        dy = self_in_other_system.cartesian.y - other.cartesian.y
        dz = self_in_other_system.cartesian.z - other.cartesian.z

        distval = (dx.value ** 2 + dy.value ** 2 + dz.value ** 2) ** 0.5
        return Distance(distval, dx.unit)

    def match_to_catalog_sky(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest on-sky matches of this coordinate in a set of
        catalog coordinates.

        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The base catalog in which to search for matches. Typically this
            will be a coordinate object that is an array (i.e.,
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is
            desired here, as that is correct for matching one set of
            coordinates to another. The next likely use case is ``2``,
            for matching a coordinate catalog against *itself* (``1``
            is inappropriate because each point will find itself as the
            closest match).

        Returns
        -------
        idx : integer array
            Indices into ``catalogcoord`` to get the matched points for
            each of this object's coordinates. Shape matches this
            object.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the closest match for each
            element in this object in ``catalogcoord``. Shape matches
            this object.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the closest match for each element
            in this object in ``catalogcoord``. Shape matches this
            object.

        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_sky
        """
        from .matching import match_coordinates_sky

        if (isinstance(catalogcoord, (SkyCoord, BaseCoordinateFrame))
                and catalogcoord.has_data):
            self_in_catalog_frame = self.transform_to(catalogcoord)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        res = match_coordinates_sky(self_in_catalog_frame, catalogcoord,
                                    nthneighbor=nthneighbor,
                                    storekdtree='_kdtree_sky')
        return res

    def match_to_catalog_3d(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest 3-dimensional matches of this coordinate to a set
        of catalog coordinates.

        This finds the 3-dimensional closest neighbor, which is only different
        from the on-sky distance if ``distance`` is set in this object or the
        ``catalogcoord`` object.

        Parameters
        ----------
        catalogcoord : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The base catalog in which to search for matches. Typically this
            will be a coordinate object that is an array (i.e.,
            ``catalogcoord.isscalar == False``)
        nthneighbor : int, optional
            Which closest neighbor to search for.  Typically ``1`` is
            desired here, as that is correct for matching one set of
            coordinates to another.  The next likely use case is
            ``2``, for matching a coordinate catalog against *itself*
            (``1`` is inappropriate because each point will find
            itself as the closest match).

        Returns
        -------
        idx : integer array
            Indices into ``catalogcoord`` to get the matched points for
            each of this object's coordinates. Shape matches this
            object.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the closest match for each
            element in this object in ``catalogcoord``. Shape matches
            this object.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the closest match for each element
            in this object in ``catalogcoord``. Shape matches this
            object.

        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_3d
        """
        from .matching import match_coordinates_3d

        if (isinstance(catalogcoord, (SkyCoord, BaseCoordinateFrame))
                and catalogcoord.has_data):
            self_in_catalog_frame = self.transform_to(catalogcoord)
        else:
            raise TypeError('Can only get separation to another SkyCoord or a '
                            'coordinate frame with data')

        res = match_coordinates_3d(self_in_catalog_frame, catalogcoord,
                                   nthneighbor=nthneighbor,
                                   storekdtree='_kdtree_3d')

        return res

    def search_around_sky(self, searcharoundcoords, seplimit):
        """
        Searches for all coordinates in this object around a supplied set of
        points within a given on-sky separation.

        Parameters
        ----------
        searcharoundcoords : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate(s) to search around to try to find matching points in
            this `SkyCoord`.
        seplimit : `~astropy.units.Quantity` with angle units
            The on-sky separation to search within.

        Returns
        -------
        idxsearcharound : integer array
            Indices into ``coords1`` that matches to the corresponding element of
            ``idxself``. Shape matches ``idxself``.
        idxself : integer array
            Indices into ``coords2`` that matches to the corresponding element of
            ``idxsearcharound``. Shape matches ``idxsearcharound``.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.

        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        In the current implementation, the return values are always sorted in
        the same order as the ``searcharoundcoords`` (so ``idxsearcharound`` is
        in ascending order).  This is considered an implementation detail,
        though, so it could change in a future release.

        See Also
        --------
        astropy.coordinates.search_around_sky
        """
        from .matching import search_around_sky

        return search_around_sky(searcharoundcoords, self, seplimit,
                                 storekdtree='_kdtree_sky')

    def search_around_3d(self, searcharoundcoords, distlimit):
        """
        Searches for all coordinates in this object around a supplied set of
        points within a given 3D radius.

        Parameters
        ----------
        searcharoundcoords : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate(s) to search around to try to find matching points in
            this `SkyCoord`.
        distlimit : `~astropy.units.Quantity` with distance units
            The physical radius to search within.

        Returns
        -------
        idxsearcharound : integer array
            Indices into ``coords1`` that matches to the corresponding element of
            ``idxself``. Shape matches ``idxself``.
        idxself : integer array
            Indices into ``coords2`` that matches to the corresponding element of
            ``idxsearcharound``. Shape matches ``idxsearcharound``.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.

        Notes
        -----
        This method requires `SciPy <http://www.scipy.org>`_ to be
        installed or it will fail.

        In the current implementation, the return values are always sorted in
        the same order as the ``searcharoundcoords`` (so ``idxsearcharound`` is
        in ascending order).  This is considered an implementation detail,
        though, so it could change in a future release.

        See Also
        --------
        astropy.coordinates.search_around_3d
        """
        from .matching import search_around_3d

        return search_around_3d(searcharoundcoords, self, distlimit,
                                storekdtree='_kdtree_3d')

    def position_angle(self, other):
        """
        Computes the on-sky position angle (East of North) between this
        `SkyCoord` and another.

        Parameters
        ----------
        other : `SkyCoord`
            The other coordinate to compute the position angle to.  It is
            treated as the "head" of the vector of the position angle.

        Returns
        -------
        pa : `~astropy.coordinates.Angle`
            The (positive) position angle of the vector pointing from ``self``
            to ``other``.  If either ``self`` or ``other`` contain arrays, this
            will be an array following the appropriate `numpy` broadcasting
            rules.

        Examples
        --------

        >>> c1 = SkyCoord(0*u.deg, 0*u.deg)
        >>> c2 = SkyCoord(1*u.deg, 0*u.deg)
        >>> c1.position_angle(c2).degree
        90.0
        >>> c3 = SkyCoord(1*u.deg, 1*u.deg)
        >>> c1.position_angle(c3).degree  # doctest: +FLOAT_CMP
        44.995636455344844
        """
        from . import angle_utilities

        if self.frame.name == other.frame.name:
            other_in_self_frame = other
        else:
            other_in_self_frame = other.frame.transform_to(self.frame)

        slat = self.represent_as(UnitSphericalRepresentation).lat
        slon = self.represent_as(UnitSphericalRepresentation).lon
        olat = other_in_self_frame.represent_as(UnitSphericalRepresentation).lat
        olon = other_in_self_frame.represent_as(UnitSphericalRepresentation).lon

        return angle_utilities.position_angle(slon, slat, olon, olat)

    # WCS pixel to/from sky conversions
    def to_pixel(self, wcs, origin=0, mode='all'):
        """
        Convert this coordinate to pixel coordinates using a `~astropy.wcs.WCS`
        object.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to use for convert
        origin : int
            Whether to return 0 or 1-based pixel coordinates.
        mode : 'all' or 'wcs'
            Whether to do the transformation including distortions (``'all'``) or
            only including only the core WCS transformation (``'wcs'``).

        Returns
        -------
        xp, yp : `numpy.ndarray`
            The pixel coordinates

        See Also
        --------
        from_pixel : to do the inverse operation
        astropy.wcs.utils.skycoord_to_pixel : the implementation of this method
        """
        return skycoord_to_pixel(self, wcs=wcs, origin=origin, mode=mode)

    @classmethod
    def from_pixel(cls, xp, yp, wcs, origin=0, mode='all'):
        """
        Create a new `SkyCoord` from pixel coordinates using an
        `~astropy.wcs.WCS` object.

        Parameters
        ----------
        xp, yp : float or `numpy.ndarray`
            The coordinates to convert.
        wcs : `~astropy.wcs.WCS`
            The WCS to use for convert
        origin : int
            Whether to return 0 or 1-based pixel coordinates.
        mode : 'all' or 'wcs'
            Whether to do the transformation including distortions (``'all'``) or
            only including only the core WCS transformation (``'wcs'``).

        Returns
        -------
        coord : an instance of this class
            A new object with sky coordinates corresponding to the input ``xp``
            and ``yp``.

        See Also
        --------
        to_pixel : to do the inverse operation
        astropy.wcs.utils.pixel_to_skycoord : the implementation of this method
        """
        return pixel_to_skycoord(xp, yp, wcs=wcs, origin=origin, mode=mode, cls=cls)

    # Table interactions
    @classmethod
    def guess_from_table(cls, table, **coord_kwargs):
        """
        A convenience method to create and return a new `SkyCoord` from the data
        in an astropy Table.

        This method matches table columns that start with the case-insensitive
        names of the the components of the requested frames, if they are also
        followed by a non-alphanumeric character. It will also match columns
        that *end* with the component name if a non-alphanumeric character is
        *before* it.

        For example, the first rule means columns with names like
        ``'RA[J2000]'`` or ``'ra'`` will be interpreted as ``ra`` attributes for
        `~astropy.coordinates.ICRS` frames, but ``'RAJ2000'`` or ``'radius'``
        are *not*. Similarly, the second rule applied to the
        `~astropy.coordinates.Galactic` frame means that a column named
        ``'gal_l'`` will be used as the the ``l`` component, but ``gall`` or
        ``'fill'`` will not.

        The definition of alphanumeric here is based on Unicode's definition
        of alphanumeric, except without ``_`` (which is normally considered
        alphanumeric).  So for ASCII, this means the non-alphanumeric characters
        are ``<space>_!"#$%&'()*+,-./:;<=>?@[\]^`{|}~``).

        Parameters
        ----------
        table : astropy.Table
            The table to load data from.
        coord_kwargs
            Any additional keyword arguments are passed directly to this class's
            constructor.

        Returns
        -------
        newsc : same as this class
            The new `SkyCoord` (or subclass) object.
        """
        inital_frame = coord_kwargs.get('frame')
        frame = _get_frame([], coord_kwargs)
        coord_kwargs['frame'] = inital_frame

        comp_kwargs = {}
        for comp_name in frame.representation_component_names:
            # this matches things like 'ra[...]'' but *not* 'rad'.
            # note that the "_" must be in there explicitly, because
            # "alphanumeric" usually includes underscores.
            starts_with_comp = comp_name + r'(\W|\b|_)'
            # this part matches stuff like 'center_ra', but *not*
            # 'aura'
            ends_with_comp = r'.*(\W|\b|_)' + comp_name + r'\b'
            #the final regex ORs together the two patterns
            rex = re.compile('(' +starts_with_comp + ')|(' + ends_with_comp + ')',
                             re.IGNORECASE | re.UNICODE)

            for col_name in table.colnames:
                if rex.match(col_name):
                    if comp_name in comp_kwargs:
                        oldname = comp_kwargs[comp_name].name
                        msg = ('Found at least two matches for  component "{0}"'
                               ': "{1}" and "{2}". Cannot continue with this '
                               'ambiguity.')
                        raise ValueError(msg.format(comp_name, oldname, col_name))
                    comp_kwargs[comp_name] = table[col_name]

        for k, v in comp_kwargs.items():
            if k in coord_kwargs:
                raise ValueError('Found column "{0}" in table, but it was '
                                 'already provided as "{1}" keyword to '
                                 'guess_from_table function.'.format(v.name, k))
            else:
                coord_kwargs[k] = v

        return cls(**coord_kwargs)

    # Name resolve
    @classmethod
    def from_name(cls, name, frame='icrs'):
        """
        Given a name, query the CDS name resolver to attempt to retrieve
        coordinate information for that object. The search database, sesame
        url, and  query timeout can be set through configuration items in
        ``astropy.coordinates.name_resolve`` -- see docstring for
        `~astropy.coordinates.get_icrs_coordinates` for more
        information.

        Parameters
        ----------
        name : str
            The name of the object to get coordinates for, e.g. ``'M42'``.
        frame : str or `BaseCoordinateFrame` class or instance
            The frame to transform the object to.

        Returns
        -------
        coord : SkyCoord
            Instance of the SkyCoord class.
        """

        from .name_resolve import get_icrs_coordinates

        icrs_coord = get_icrs_coordinates(name)
        icrs_sky_coord = cls(icrs_coord)
        if frame in ('icrs', icrs_coord.__class__):
            return icrs_sky_coord
        else:
            return icrs_sky_coord.transform_to(frame)


# <----------------Private utility functions below here------------------------->


def _get_frame_class(frame):
    """
    Get a frame class from the input `frame`, which could be a frame name
    string, or frame class.
    """
    import inspect

    if isinstance(frame, six.string_types):
        frame_names = frame_transform_graph.get_names()
        if frame not in frame_names:
            raise ValueError('Coordinate frame {0} not in allowed values {1}'
                             .format(frame, sorted(frame_names)))
        frame_cls = frame_transform_graph.lookup_name(frame)

    elif inspect.isclass(frame) and issubclass(frame, BaseCoordinateFrame):
        frame_cls = frame

    else:
        raise ValueError('Coordinate frame must be a frame name or frame class')

    return frame_cls


def _get_frame(args, kwargs):
    """
    Determine the coordinate frame from input SkyCoord args and kwargs.  This
    modifies args and/or kwargs in-place to remove the item that provided
    `frame`.  It also infers the frame if an input coordinate was provided and
    checks for conflicts.

    This allows for frame to be specified as a string like 'icrs' or a frame
    class like ICRS, but not an instance ICRS() since the latter could have
    non-default representation attributes which would require a three-way merge.
    """
    frame = kwargs.pop('frame', None)

    if frame is None and len(args) > 1:

        # We do not allow frames to be passed as positional arguments if data
        # is passed separately from frame.

        for arg in args:

            if isinstance(arg, (SkyCoord, BaseCoordinateFrame)):
                raise ValueError("{0} instance cannot be passed as a positional "
                                 "argument for the frame, pass it using the "
                                 "frame= keyword instead.".format(arg.__class__.__name__))

    # If the frame is an instance or SkyCoord, we split up the attributes and
    # make it into a class.

    if isinstance(frame, SkyCoord):
        frame = frame.frame

    if isinstance(frame, BaseCoordinateFrame):

        for attr in frame.get_frame_attr_names():
            if attr in kwargs:
                raise ValueError("cannot specify frame attribute '{0}' directly in SkyCoord since a frame instance was passed in".format(attr))
            else:
                kwargs[attr] = getattr(frame, attr)

        frame = frame.__class__

    if frame is not None:
        # Frame was provided as kwarg so validate and coerce into corresponding frame.
        frame_cls = _get_frame_class(frame)
        frame_specified_explicitly = True
    else:
        # Look for the frame in args
        for arg in args:
            try:
                frame_cls = _get_frame_class(arg)
                frame_specified_explicitly = True
            except ValueError:
                pass
            else:
                args.remove(arg)
                warnings.warn("Passing a frame as a positional argument is now "
                              "deprecated, use the frame= keyword argument "
                              "instead.", AstropyDeprecationWarning)
                break
        else:
            # Not in args nor kwargs - default to icrs
            frame_cls = ICRS
            frame_specified_explicitly = False

    # Check that the new frame doesn't conflict with existing coordinate frame
    # if a coordinate is supplied in the args list.  If the frame still had not
    # been set by this point and a coordinate was supplied, then use that frame.
    for arg in args:
        coord_frame_cls = None
        if isinstance(arg, BaseCoordinateFrame):
            coord_frame_cls = arg.__class__
        elif isinstance(arg, SkyCoord):
            coord_frame_cls = arg.frame.__class__

        if coord_frame_cls is not None:
            if not frame_specified_explicitly:
                frame_cls = coord_frame_cls
            elif frame_cls is not coord_frame_cls:
                raise ValueError("Cannot override frame='{0}' of input coordinate with "
                                 "new frame='{1}'.  Instead transform the coordinate."
                                 .format(coord_frame_cls.__name__, frame_cls.__name__))

    if 'representation' in kwargs:
        frame = frame_cls(representation=_get_repr_cls(kwargs['representation']))
    else:
        frame = frame_cls()

    return frame


def _get_units(args, kwargs):
    """
    Get the longitude unit and latitude unit from kwargs.  Possible enhancement
    is to allow input from args as well.
    """
    if 'unit' not in kwargs:
        units = [None, None, None]

    else:
        units = kwargs.pop('unit')

        if isinstance(units, six.string_types):
            units = [x.strip() for x in units.split(',')]
            # Allow for input like unit='deg' or unit='m'
            if len(units) == 1:
                units = [units[0], units[0], units[0]]
        elif isinstance(units, (Unit, IrreducibleUnit)):
            units = [units, units, units]

        try:
            units = [(Unit(x) if x else None) for x in units]
            units.extend(None for x in range(3 - len(units)))
            if len(units) > 3:
                raise ValueError()
        except:
            raise ValueError('Unit keyword must have one to three unit values as '
                             'tuple or comma-separated string')

    return units


def _parse_coordinate_arg(coords, frame, units, init_kwargs):
    """
    Single unnamed arg supplied.  This must be:
    - Coordinate frame with data
    - Representation
    - SkyCoord
    - List or tuple of:
      - String which splits into two values
      - Iterable with two values
      - SkyCoord, frame, or representation objects.

    Returns a dict mapping coordinate attribute names to values (or lists of
    values)
    """
    is_scalar = False  # Differentiate between scalar and list input
    valid_kwargs = {}  # Returned dict of lon, lat, and distance (optional)

    frame_attr_names = frame.representation_component_names.keys()
    repr_attr_names = frame.representation_component_names.values()
    repr_attr_classes = frame.representation.attr_classes.values()
    n_attr_names = len(repr_attr_names)

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, six.string_types):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `frame` it is checked that any coordinate
        # args have the same frame as explicitly supplied, so don't worry here.

        if not coords.has_data:
            raise ValueError('Cannot initialize from a frame without coordinate data')

        data = coords.data.represent_as(frame.representation)

        values = []  # List of values corresponding to representation attrs
        for repr_attr_name in repr_attr_names:
            # If coords did not have an explicit distance then don't include in initializers.
            if (isinstance(coords.data, UnitSphericalRepresentation) and
                    repr_attr_name == 'distance'):
                continue

            # Get the value from `data` in the eventual representation
            values.append(getattr(data, repr_attr_name))

        for attr in FRAME_ATTR_NAMES_SET():
            value = getattr(coords, attr, None)
            use_value = (isinstance(coords, SkyCoord)
                         or attr not in coords._attr_names_with_defaults)
            if use_value and value is not None:
                valid_kwargs[attr] = value

    elif isinstance(coords, BaseRepresentation):
        data = coords.represent_as(frame.representation)
        values = [getattr(data, repr_attr_name) for repr_attr_name in repr_attr_names]

    elif (isinstance(coords, np.ndarray) and coords.dtype.kind in 'if'
          and coords.ndim == 2 and coords.shape[1] <= 3):
        # 2-d array of coordinate values.  Handle specially for efficiency.
        values = coords.transpose()  # Iterates over repr attrs

    elif isinstance(coords, (collections.Sequence, np.ndarray)):
        # Handles list-like input.

        vals = []
        is_ra_dec_representation = ('ra' in frame.representation_component_names and
                                    'dec' in frame.representation_component_names)
        coord_types = (SkyCoord, BaseCoordinateFrame, BaseRepresentation)
        if any(isinstance(coord, coord_types) for coord in coords):
            # this parsing path is used when there are coordinate-like objects
            # in the list - instead of creating lists of values, we create
            # SkyCoords from the list elements and then combine them.
            scs = [SkyCoord(coord, **init_kwargs) for coord in coords]

            # now check that they're all self-consistent in their frame attributes
            # and frame name
            frames_to_check = [sc.frame.name for sc in scs]
            if len(set(frames_to_check)) > 1:
                raise ValueError("List of inputs have different frames: {0}".format(frames_to_check))
            for fattrnm in FRAME_ATTR_NAMES_SET():
                vals = [getattr(sc, fattrnm) for sc in scs]
                for val in vals[1:]:
                    if val != vals[0]:
                        raise ValueError("List of inputs don't give consistent "
                                         "frame attribute {0}: {1}".format(fattrnm, vals))

                valid_kwargs[fattrnm] = getattr(scs[0], fattrnm)

            # Now combine the values, to be used below
            values = []
            for data_attr_name in frame_attr_names:
                data_vals = []
                for sc in scs:
                    data_val = getattr(sc, data_attr_name)
                    data_vals.append(data_val.reshape(1,) if sc.isscalar else data_val)
                concat_vals = np.concatenate(data_vals)
                # Hack because np.concatenate doesn't fully work with Quantity
                if isinstance(concat_vals, u.Quantity):
                    concat_vals._unit = data_val.unit
                values.append(concat_vals)
        else:
            #none of the elements are "frame-like"
            #turn into a list of lists like [[v1_0, v2_0, v3_0], ... [v1_N, v2_N, v3_N]]
            for coord in coords:
                if isinstance(coord, six.string_types):
                    coord1 = coord.split()
                    if len(coord1) == 6:
                        coord = (' '.join(coord1[:3]), ' '.join(coord1[3:]))
                    elif is_ra_dec_representation:
                        coord = _parse_ra_dec(coord)
                    else:
                        coord = coord1
                vals.append(coord)  # Assumes coord is a sequence at this point

            # Do some basic validation of the list elements: all have a length and all
            # lengths the same
            try:
                n_coords = sorted(set(len(x) for x in vals))
            except:
                raise ValueError('One or more elements of input sequence does not have a length')

            if len(n_coords) > 1:
                raise ValueError('Input coordinate values must have same number of elements, found {0}'
                                 .format(n_coords))
            n_coords = n_coords[0]

            # Must have no more coord inputs than representation attributes
            if n_coords > n_attr_names:
                raise ValueError('Input coordinates have {0} values but '
                                 'representation {1} only accepts {2}'
                                 .format(n_coords, frame.representation.get_name(), n_attr_names))

            # Now transpose vals to get [(v1_0 .. v1_N), (v2_0 .. v2_N), (v3_0 .. v3_N)]
            # (ok since we know it is exactly rectangular).  (Note: can't just use zip(*values)
            # because Longitude et al distinguishes list from tuple so [a1, a2, ..] is needed
            # while (a1, a2, ..) doesn't work.
            values = [list(x) for x in zip(*vals)]

            if is_scalar:
                values = [x[0] for x in values]
    else:
        raise ValueError('Cannot parse coordinates from first argument')

    # Finally we have a list of values from which to create the keyword args
    # for the frame initialization.  Validate by running through the appropriate
    # class initializer and supply units (which might be None).
    try:
        for frame_attr_name, repr_attr_class, value, unit in zip(
                frame_attr_names, repr_attr_classes, values, units):
            valid_kwargs[frame_attr_name] = repr_attr_class(value, unit=unit)
    except Exception as err:
        raise ValueError('Cannot parse first argument data "{0}" for attribute '
                         '{1}'.format(value, frame_attr_name), err)
    return valid_kwargs


def _get_representation_attrs(frame, units, kwargs):
    """
    Find instances of the "representation attributes" for specifying data
    for this frame.  Pop them off of kwargs, run through the appropriate class
    constructor (to validate and apply unit), and put into the output
    valid_kwargs.  "Representation attributes" are the frame-specific aliases
    for the underlying data values in the representation, e.g. "ra" for "lon"
    for many equatorial spherical representations, or "w" for "x" in the
    cartesian representation of Galactic.
    """
    frame_attr_names = frame.representation_component_names.keys()
    repr_attr_classes = frame.representation.attr_classes.values()

    valid_kwargs = {}
    for frame_attr_name, repr_attr_class, unit in zip(frame_attr_names, repr_attr_classes, units):
        value = kwargs.pop(frame_attr_name, None)
        if value is not None:
            valid_kwargs[frame_attr_name] = repr_attr_class(value, unit=unit)

    return valid_kwargs


def _parse_ra_dec(coord_str):
    """
    Parse RA and Dec values from a coordinate string. Currently the
    following formats are supported:

     * space separated 6-value format
     * space separated <6-value format, this requires a plus or minus sign
       separation between RA and Dec
     * sign separated format
     * JHHMMSS.ss+DDMMSS.ss format, with up to two optional decimal digits
     * JDDDMMSS.ss+DDMMSS.ss format, with up to two optional decimal digits

    Parameters
    ----------
    coord_str : str
        Coordinate string to parse.

    Returns
    -------
    coord : str or list of str
        Parsed coordinate values.
    """

    if isinstance(coord_str, six.string_types):
        coord1 = coord_str.split()
    else:
        # This exception should never be raised from SkyCoord
        raise TypeError('coord_str must be a single str')

    if len(coord1) == 6:
        coord = (' '.join(coord1[:3]), ' '.join(coord1[3:]))
    elif len(coord1) > 2:
        coord = PLUS_MINUS_RE.split(coord_str)
        coord = (coord[0], ' '.join(coord[1:]))
    elif len(coord1) == 1:
        match_j = J_PREFIXED_RA_DEC_RE.match(coord_str)
        if match_j:
            coord = match_j.groups()
            if len(coord[0].split('.')[0]) == 7:
                coord = ('{0} {1} {2}'.
                         format(coord[0][0:3], coord[0][3:5], coord[0][5:]),
                         '{0} {1} {2}'.
                         format(coord[1][0:3], coord[1][3:5], coord[1][5:]))
            else:
                coord = ('{0} {1} {2}'.
                         format(coord[0][0:2], coord[0][2:4], coord[0][4:]),
                         '{0} {1} {2}'.
                         format(coord[1][0:3], coord[1][3:5], coord[1][5:]))
        else:
            coord = PLUS_MINUS_RE.split(coord_str)
            coord = (coord[0], ' '.join(coord[1:]))
    else:
        coord = coord1

    return coord
