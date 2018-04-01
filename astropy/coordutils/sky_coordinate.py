
import re
import copy
import warnings
import collections

import numpy as np

from .. import _erfa as erfa
from ..utils.compat.misc import override__dir__
from ..units import Unit, IrreducibleUnit
from .. import units as u
from ..constants import c as speed_of_light
from ..wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from ..utils.exceptions import AstropyDeprecationWarning
from ..utils.data_info import MixinInfo
from ..utils import ShapedLikeNDArray
from ..time import Time

from .distances import Distance
from .angles import Angle
from .baseframe import (BaseCoordinateFrame, frame_transform_graph,
                        GenericFrame, _get_repr_cls, _get_diff_cls,
                        _normalize_representation_type)
from .builtin_frames import ICRS, SkyOffsetFrame
from .representation import (BaseRepresentation, SphericalRepresentation,
                             UnitSphericalRepresentation, SphericalDifferential)

__all__ = ['SkyCoord', 'SkyCoordInfo']

PLUS_MINUS_RE = re.compile(r'(\+|\-)')
J_PREFIXED_RA_DEC_RE = re.compile(
    r"""J                              # J prefix
    ([0-9]{6,7}\.?[0-9]{0,2})          # RA as HHMMSS.ss or DDDMMSS.ss, optional decimal digits
    ([\+\-][0-9]{6}\.?[0-9]{0,2})\s*$  # Dec as DDMMSS.ss, optional decimal digits
    """, re.VERBOSE)


class SkyCoordInfo(MixinInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """
    attrs_from_parent = set(['unit'])  # Unit is read-only
    _supports_indexing = False

    @staticmethod
    def default_format(val):
        repr_data = val.info._repr_data
        formats = ['{0.' + compname + '.value:}' for compname
                   in repr_data.components]
        return ','.join(formats).format(repr_data)

    @property
    def unit(self):
        repr_data = self._repr_data
        unit = ','.join(str(getattr(repr_data, comp).unit) or 'None'
                        for comp in repr_data.components)
        return unit

    @property
    def _repr_data(self):
        if self._parent is None:
            return None

        sc = self._parent
        if (issubclass(sc.representation_type, SphericalRepresentation) and
                isinstance(sc.data, UnitSphericalRepresentation)):
            repr_data = sc.represent_as(sc.data.__class__, in_frame_units=True)
        else:
            repr_data = sc.represent_as(sc.representation_type,
                                        in_frame_units=True)
        return repr_data

    def _represent_as_dict(self):
        obj = self._parent
        attrs = (list(obj.representation_component_names) +
                 list(frame_transform_graph.frame_attributes.keys()))

        # Don't output distance if it is all unitless 1.0
        if 'distance' in attrs and np.all(obj.distance == 1.0):
            attrs.remove('distance')

        self._represent_as_dict_attrs = attrs

        out = super()._represent_as_dict()

        out['representation_type'] = obj.representation_type.get_name()
        out['frame'] = obj.frame.name
        # Note that obj.info.unit is a fake composite unit (e.g. 'deg,deg,None'
        # or None,None,m) and is not stored.  The individual attributes have
        # units.

        return out


class SkyCoord(ShapedLikeNDArray):
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
    argument ``representation_type='cartesian'`` (for example) along with data
    in ``x``, ``y``, and ``z``.

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

      >>> c = SkyCoord(w=0, u=1, v=2, unit='kpc', frame='galactic',
      ...              representation_type='cartesian')

      >>> c = SkyCoord([ICRS(ra=1*u.deg, dec=2*u.deg), ICRS(ra=3*u.deg, dec=4*u.deg)])

    Velocity components (proper motions or radial velocities) can also be
    provided in a similar manner::

      >>> c = SkyCoord(ra=1*u.deg, dec=2*u.deg, radial_velocity=10*u.km/u.s)

      >>> c = SkyCoord(ra=1*u.deg, dec=2*u.deg, pm_ra_cosdec=2*u.mas/u.yr, pm_dec=1*u.mas/u.yr)

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
    representation_type : str or Representation class
        Specifies the representation, e.g. 'spherical', 'cartesian', or
        'cylindrical'.  This affects the positional args and other keyword args
        which must correspond to the given representation.
    copy : bool, optional
        If `True` (default), a copy of any coordinate data is made.  This
        argument can only be passed in as a keyword argument.
    **keyword_args
        Other keyword arguments as applicable for user-defined coordinate frames.
        Common options include:

        ra, dec : valid `~astropy.coordinates.Angle` initializer, optional
            RA and Dec for frames where ``ra`` and ``dec`` are keys in the
            frame's ``representation_component_names``, including `ICRS`,
            `FK5`, `FK4`, and `FK4NoETerms`.
        pm_ra_cosdec, pm_dec  : `~astropy.units.Quantity`, optional
            Proper motion components, in angle per time units.
        l, b : valid `~astropy.coordinates.Angle` initializer, optional
            Galactic ``l`` and ``b`` for for frames where ``l`` and ``b`` are
            keys in the frame's ``representation_component_names``, including
            the `Galactic` frame.
        pm_l_cosb, pm_b : `~astropy.units.Quantity`, optional
            Proper motion components in the `Galactic` frame, in angle per time
            units.
        x, y, z : float or `~astropy.units.Quantity`, optional
            Cartesian coordinates values
        u, v, w : float or `~astropy.units.Quantity`, optional
            Cartesian coordinates values for the Galactic frame.
        radial_velocity : `~astropy.units.Quantity`, optional
            The component of the velocity along the line-of-sight (i.e., the
            radial direction), in velocity units.
    """

    # Declare that SkyCoord can be used as a Table column by defining the
    # info property.
    info = SkyCoordInfo()

    def __init__(self, *args, copy=True, **kwargs):

        # Parse the args and kwargs to assemble a sanitized and validated
        # kwargs dict for initializing attributes for this object and for
        # creating the internal self._sky_coord_frame object
        args = list(args)  # Make it mutable
        kwargs = self._parse_inputs(args, kwargs)

        frame = kwargs['frame']
        frame_attr_names = frame.get_frame_attr_names()

        # these are frame attributes set on this SkyCoord but *not* a part of
        # the frame object this SkyCoord contains
        self._extra_frameattr_names = set()

        for attr in kwargs:
            if (attr not in frame_attr_names and
                attr in frame_transform_graph.frame_attributes):
                # Setting it will also validate it.
                setattr(self, attr, kwargs[attr])

        coord_kwargs = {}

        component_names = frame.representation_component_names
        component_names.update(frame.get_representation_component_names('s'))

        # TODO: deprecate representation, remove this in future
        _normalize_representation_type(kwargs)
        if 'representation_type' in kwargs:
            coord_kwargs['representation_type'] = _get_repr_cls(
                kwargs['representation_type'])

        if 'differential_type' in kwargs:
            coord_kwargs['differential_type'] = _get_diff_cls(kwargs['differential_type'])

        for attr, value in kwargs.items():
            if value is not None and (attr in component_names
                                      or attr in frame_attr_names):
                coord_kwargs[attr] = value

        # Finally make the internal coordinate object.
        self._sky_coord_frame = frame.__class__(copy=copy, **coord_kwargs)

        if not self._sky_coord_frame.has_data:
            raise ValueError('Cannot create a SkyCoord without data')

    @property
    def frame(self):
        return self._sky_coord_frame

    @property
    def representation_type(self):
        return self.frame.representation_type

    @representation_type.setter
    def representation_type(self, value):
        self.frame.representation_type = value

    # TODO: deprecate these in future
    @property
    def representation(self):
        return self.frame.representation

    @representation.setter
    def representation(self, value):
        self.frame.representation = value

    @property
    def shape(self):
        return self.frame.shape

    def _apply(self, method, *args, **kwargs):
        """Create a new instance, applying a method to the underlying data.

        In typical usage, the method is any of the shape-changing methods for
        `~numpy.ndarray` (``reshape``, ``swapaxes``, etc.), as well as those
        picking particular elements (``__getitem__``, ``take``, etc.), which
        are all defined in `~astropy.utils.misc.ShapedLikeNDArray`. It will be
        applied to the underlying arrays in the representation (e.g., ``x``,
        ``y``, and ``z`` for `~astropy.coordinates.CartesianRepresentation`),
        as well as to any frame attributes that have a shape, with the results
        used to create a new instance.

        Internally, it is also used to apply functions to the above parts
        (in particular, `~numpy.broadcast_to`).

        Parameters
        ----------
        method : str or callable
            If str, it is the name of a method that is applied to the internal
            ``components``. If callable, the function is applied.
        args : tuple
            Any positional arguments for ``method``.
        kwargs : dict
            Any keyword arguments for ``method``.
        """
        def apply_method(value):
            if isinstance(value, ShapedLikeNDArray):
                return value._apply(method, *args, **kwargs)
            else:
                if callable(method):
                    return method(value, *args, **kwargs)
                else:
                    return getattr(value, method)(*args, **kwargs)

        # create a new but empty instance, and copy over stuff
        new = super().__new__(self.__class__)
        new._sky_coord_frame = self._sky_coord_frame._apply(method,
                                                            *args, **kwargs)
        new._extra_frameattr_names = self._extra_frameattr_names.copy()
        for attr in self._extra_frameattr_names:
            value = getattr(self, attr)
            if getattr(value, 'size', 1) > 1:
                value = apply_method(value)
            elif method == 'copy' or method == 'flatten':
                # flatten should copy also for a single element array, but
                # we cannot use it directly for array scalars, since it
                # always returns a one-dimensional array. So, just copy.
                value = copy.copy(value)
            setattr(new, '_' + attr, value)

        # Copy other 'info' attr only if it has actually been defined.
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if 'info' in self.__dict__:
            new.info = self.info

        return new

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

        # TODO: possibly remove the below.  The representation/differential
        # information should *already* be stored in the frame object, as it is
        # extracted in _get_frame.  So it may be redundent to include it below.
        # TODO: deprecate representation, remove this in future
        _normalize_representation_type(kwargs)
        if 'representation_type' in kwargs:
            valid_kwargs['representation_type'] = _get_repr_cls(
                kwargs.pop('representation_type'))

        if 'differential_type' in kwargs:
            valid_kwargs['differential_type'] = _get_diff_cls(
                kwargs.pop('differential_type'))

        for attr in frame_transform_graph.frame_attributes:
            if attr in kwargs:
                valid_kwargs[attr] = kwargs.pop(attr)

        # Get units
        units = _get_representation_component_units(args, kwargs)

        # Grab any frame-specific attr names like `ra` or `l` or `distance` from
        # kwargs and migrate to valid_kwargs.
        valid_kwargs.update(_get_representation_attrs(frame, units, kwargs))

        # Error if anything is still left in kwargs
        if kwargs:
            # The next few lines add a more user-friendly error message to a
            # common and confusing situation when the user specifies, e.g.,
            # `pm_ra` when they really should be passing `pm_ra_cosdec`. The
            # extra error should only turn on when the positional representation
            # is spherical, and when the component 'pm_<lon>' is passed.
            pm_message = ''
            if frame.representation_type == SphericalRepresentation:
                frame_names = list(frame.get_representation_component_names().keys())
                lon_name = frame_names[0]
                lat_name = frame_names[1]

                if 'pm_{0}'.format(lon_name) in list(kwargs.keys()):
                    pm_message = ('\n\n By default, most frame classes expect '
                                  'the longitudinal proper motion to include '
                                  'the cos(latitude) term, named '
                                  '`pm_{0}_cos{1}`. Did you mean to pass in '
                                  'this component?'
                                  .format(lon_name, lat_name))

            raise ValueError('Unrecognized keyword argument(s) {0}{1}'
                             .format(', '.join("'{0}'".format(key)
                                               for key in kwargs),
                                     pm_message))

        # Finally deal with the unnamed args.  This figures out what the arg[0]
        # is and returns a dict with appropriate key/values for initializing
        # frame class. Note that differentials are *never* valid args, only
        # kwargs.  So they are not accounted for here (unless they're in a frame
        # or SkyCoord object)
        if args:
            if len(args) == 1:
                # One arg which must be a coordinate.  In this case
                # coord_kwargs will contain keys like 'ra', 'dec', 'distance'
                # along with any frame attributes like equinox or obstime which
                # were explicitly specified in the coordinate object (i.e. non-default).
                coord_kwargs = _parse_coordinate_arg(args[0], frame, units, kwargs)

                # Copy other 'info' attr only if it has actually been defined.
                if 'info' in getattr(args[0], '__dict__', ()):
                    self.info = args[0].info

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

    def transform_to(self, frame, merge_attributes=True):
        """Transform this coordinate to a new frame.

        The precise frame transformed to depends on ``merge_attributes``.
        If `False`, the destination frame is used exactly as passed in.
        But this is often not quite what one wants.  E.g., suppose one wants to
        transform an ICRS coordinate that has an obstime attribute to FK4; in
        this case, one likely would want to use this information. Thus, the
        default for ``merge_attributes`` is `True`, in which the precedence is
        as follows: (1) explicitly set (i.e., non-default) values in the
        destination frame; (2) explicitly set values in the source; (3) default
        value in the destination frame.

        Note that in either case, any explicitly set attributes on the source
        `SkyCoord` that are not part of the destination frame's definition are
        kept (stored on the resulting `SkyCoord`), and thus one can round-trip
        (e.g., from FK4 to ICRS to FK4 without loosing obstime).

        Parameters
        ----------
        frame : str, `BaseCoordinateFrame` class or instance, or `SkyCoord` instance
            The frame to transform this coordinate into.  If a `SkyCoord`, the
            underlying frame is extracted, and all other information ignored.
        merge_attributes : bool, optional
            Whether the default attributes in the destination frame are allowed
            to be overridden by explicitly set attributes in the source
            (see note above; default: `True`).

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
        except Exception:
            pass

        if isinstance(frame, SkyCoord):
            frame = frame.frame  # Change to underlying coord frame instance

        if isinstance(frame, BaseCoordinateFrame):
            new_frame_cls = frame.__class__
            # Get frame attributes, allowing defaults to be overridden by
            # explicitly set attributes of the source if ``merge_attributes``.
            for attr in frame_transform_graph.frame_attributes:
                self_val = getattr(self, attr, None)
                frame_val = getattr(frame, attr, None)
                if (frame_val is not None and not
                    (merge_attributes and frame.is_frame_attr_default(attr))):
                    frame_kwargs[attr] = frame_val
                elif (self_val is not None and
                      not self.is_frame_attr_default(attr)):
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
        for attr in (set(new_coord.get_frame_attr_names()) &
                     set(frame_kwargs.keys())):
            frame_kwargs.pop(attr)
        return self.__class__(new_coord, **frame_kwargs)

    def apply_space_motion(self, new_obstime=None, dt=None):
        """
        Compute the position of the source represented by this coordinate object
        to a new time using the velocities stored in this object and assuming
        linear space motion (including relativistic corrections). This is
        sometimes referred to as an "epoch transformation."

        The initial time before the evolution is taken from the ``obstime``
        attribute of this coordinate.  Note that this method currently does not
        support evolving coordinates where the *frame* has an ``obstime`` frame
        attribute, so the ``obstime`` is only used for storing the before and
        after times, not actually as an attribute of the frame. Alternatively,
        if ``dt`` is given, an ``obstime`` need not be provided at all.

        Parameters
        ----------
        new_obstime : `~astropy.time.Time`, optional
            The time at which to evolve the position to. Requires that the
            ``obstime`` attribute be present on this frame.
        dt : `~astropy.units.Quantity`, `~astropy.time.TimeDelta`, optional
            An amount of time to evolve the position of the source. Cannot be
            given at the same time as ``new_obstime``.

        Returns
        -------
        new_coord : `SkyCoord`
            A new coordinate object with the evolved location of this coordinate
            at the new time.  ``obstime`` will be set on this object to the new
            time only if ``self`` also has ``obstime``.
        """

        if (new_obstime is None and dt is None or
                new_obstime is not None and dt is not None):
            raise ValueError("You must specify one of `new_obstime` or `dt`, "
                             "but not both.")

        # Validate that we have velocity info
        if 's' not in self.frame.data.differentials:
            raise ValueError('SkyCoord requires velocity data to evolve the '
                             'position.')

        if 'obstime' in self.frame.frame_attributes:
            raise NotImplementedError("Updating the coordinates in a frame "
                                      "with explicit time dependence is "
                                      "currently not supported. If you would "
                                      "like this functionality, please open an "
                                      "issue on github:\n"
                                      "https://github.com/astropy/astropy")

        if new_obstime is not None and self.obstime is None:
            # If no obstime is already on this object, raise an error if a new
            # obstime is passed: we need to know the time / epoch at which the
            # the position / velocity were measured initially
            raise ValueError('This object has no associated `obstime`. '
                             'apply_space_motion() must receive a time '
                             'difference, `dt`, and not a new obstime.')

        # Compute t1 and t2, the times used in the starpm call, which *only*
        # uses them to compute a delta-time
        t1 = self.obstime
        if dt is None:
            # self.obstime is not None and new_obstime is not None b/c of above
            # checks
            t2 = new_obstime
        else:
            # new_obstime is definitely None b/c of the above checks
            if t1 is None:
                # MAGIC NUMBER: if the current SkyCoord object has no obstime,
                # assume J2000 to do the dt offset. This is not actually used
                # for anything except a delta-t in starpm, so it's OK that it's
                # not necessarily the "real" obstime
                t1 = Time('J2000')
                new_obstime = None  # we don't actually know the inital obstime
                t2 = t1 + dt
            else:
                t2 = t1 + dt
                new_obstime = t2
        # starpm wants tdb time
        t1 = t1.tdb
        t2 = t2.tdb

        # proper motion in RA should not include the cos(dec) term, see the
        # erfa function eraStarpv, comment (4).  So we convert to the regular
        # spherical differentials.
        icrsrep = self.icrs.represent_as(SphericalRepresentation, SphericalDifferential)
        icrsvel = icrsrep.differentials['s']

        try:
            plx = icrsrep.distance.to_value(u.arcsecond, u.parallax())
        except u.UnitConversionError: # No distance: set to 0 by starpm convention
            plx = 0.

        try:
            rv = icrsvel.d_distance.to_value(u.km/u.s)
        except u.UnitConversionError: # No RV
            rv = 0.

        starpm = erfa.starpm(icrsrep.lon.radian, icrsrep.lat.radian,
                             icrsvel.d_lon.to_value(u.radian/u.yr),
                             icrsvel.d_lat.to_value(u.radian/u.yr),
                             plx, rv, t1.jd1, t1.jd2, t2.jd1, t2.jd2)

        icrs2 = ICRS(ra=u.Quantity(starpm[0], u.radian, copy=False),
                     dec=u.Quantity(starpm[1], u.radian, copy=False),
                     pm_ra=u.Quantity(starpm[2], u.radian/u.yr, copy=False),
                     pm_dec=u.Quantity(starpm[3], u.radian/u.yr, copy=False),
                     distance=Distance(parallax=starpm[4] * u.arcsec, copy=False),
                     radial_velocity=u.Quantity(starpm[5], u.km/u.s, copy=False),
                     differential_type=SphericalDifferential)

        # Update the obstime of the returned SkyCoord, and need to carry along
        # the frame attributes
        frattrs = {attrnm: getattr(self, attrnm)
                   for attrnm in self._extra_frameattr_names}
        frattrs['obstime'] = new_obstime
        return self.__class__(icrs2, **frattrs).transform_to(self.frame)

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
            if attr in frame_transform_graph.frame_attributes:
                if attr in self.frame.get_frame_attr_names():
                    return getattr(self.frame, attr)
                else:
                    return getattr(self, '_' + attr, None)

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

            if not attr.startswith('_') and hasattr(self._sky_coord_frame, attr):
                setattr(self._sky_coord_frame, attr, val)
                return

            frame_cls = frame_transform_graph.lookup_name(attr)
            if frame_cls is not None and self.frame.is_transformable_to(frame_cls):
                raise AttributeError("'{0}' is immutable".format(attr))

        if attr in frame_transform_graph.frame_attributes:
            # All possible frame attributes can be set, but only via a private
            # variable.  See __getattr__ above.
            super().__setattr__('_' + attr, val)
            # Validate it
            frame_transform_graph.frame_attributes[attr].__get__(self)
            # And add to set of extra attributes
            self._extra_frameattr_names |= {attr}

        else:
            # Otherwise, do the standard Python attribute setting
            super().__setattr__(attr, val)

    def __delattr__(self, attr):
        # mirror __setattr__ above
        if '_sky_coord_frame' in self.__dict__:
            if self.frame.name == attr:
                raise AttributeError("'{0}' is immutable".format(attr))

            if not attr.startswith('_') and hasattr(self._sky_coord_frame,
                                                    attr):
                delattr(self._sky_coord_frame, attr)
                return

            frame_cls = frame_transform_graph.lookup_name(attr)
            if frame_cls is not None and self.frame.is_transformable_to(frame_cls):
                raise AttributeError("'{0}' is immutable".format(attr))

        if attr in frame_transform_graph.frame_attributes:
            # All possible frame attributes can be deleted, but need to remove
            # the corresponding private variable.  See __getattr__ above.
            super().__delattr__('_' + attr)
            # Also remove it from the set of extra attributes
            self._extra_frameattr_names -= {attr}

        else:
            # Otherwise, do the standard Python attribute setting
            super().__delattr__(attr)

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
        dir_values.update(frame_transform_graph.frame_attributes.keys())

        return dir_values

    def __repr__(self):
        clsnm = self.__class__.__name__
        coonm = self.frame.__class__.__name__
        frameattrs = self.frame._frame_attrs_repr()
        if frameattrs:
            frameattrs = ': ' + frameattrs

        data = self.frame._data_repr()
        if data:
            data = ': ' + data

        return '<{clsnm} ({coonm}{frameattrs}){data}>'.format(**locals())

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
            for lonangle, latangle in zip(sph_coord.lon.ravel(), sph_coord.lat.ravel()):
                coord_string += [(lonangle.to_string(**lonargs)
                                 + " " +
                                 latangle.to_string(**latargs))]
            if len(sph_coord.shape) > 1:
                coord_string = np.array(coord_string).reshape(sph_coord.shape)

        return coord_string

    def is_equivalent_frame(self, other):
        """
        Checks if this object's frame as the same as that of the ``other``
        object.

        To be the same frame, two objects must be the same frame class and have
        the same frame attributes. For two `SkyCoord` objects, *all* of the
        frame attributes have to match, not just those relevant for the object's
        frame.

        Parameters
        ----------
        other : SkyCoord or BaseCoordinateFrame
            The other object to check.

        Returns
        -------
        isequiv : bool
            True if the frames are the same, False if not.

        Raises
        ------
        TypeError
            If ``other`` isn't a `SkyCoord` or a `BaseCoordinateFrame` or subclass.
        """
        if isinstance(other, BaseCoordinateFrame):
            return self.frame.is_equivalent_frame(other)
        elif isinstance(other, SkyCoord):
            if other.frame.name != self.frame.name:
                return False

            for fattrnm in frame_transform_graph.frame_attributes:
                if np.any(getattr(self, fattrnm) != getattr(other, fattrnm)):
                    return False
            return True
        else:
            # not a BaseCoordinateFrame nor a SkyCoord object
            raise TypeError("Tried to do is_equivalent_frame on something that "
                            "isn't frame-like")

    # High-level convenience methods
    def separation(self, other):
        """
        Computes on-sky separation between this coordinate and another.

        .. note::

            If the ``other`` coordinate object is in a different frame, it is
            first transformed to the frame of this object. This can lead to
            unintutive behavior if not accounted for. Particularly of note is
            that ``self.separation(other)`` and ``other.separation(self)`` may
            not give the same answer in this case.

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

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

        .. [1] https://en.wikipedia.org/wiki/Great-circle_distance

        """
        from . import Angle
        from .angle_utilities import angular_separation

        if not self.is_equivalent_frame(other):
            try:
                other = other.transform_to(self, merge_attributes=False)
            except TypeError:
                raise TypeError('Can only get separation to another SkyCoord '
                                'or a coordinate frame with data')

        lon1 = self.spherical.lon
        lat1 = self.spherical.lat
        lon2 = other.spherical.lon
        lat2 = other.spherical.lat

        # Get the separation as a Quantity, convert to Angle in degrees
        sep = angular_separation(lon1, lat1, lon2, lat2)
        return Angle(sep, unit=u.degree)

    def separation_3d(self, other):
        """
        Computes three dimensional separation between this coordinate
        and another.

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

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
        if not self.is_equivalent_frame(other):
            try:
                other = other.transform_to(self, merge_attributes=False)
            except TypeError:
                raise TypeError('Can only get separation to another SkyCoord '
                                'or a coordinate frame with data')

        if issubclass(self.data.__class__, UnitSphericalRepresentation):
            raise ValueError('This object does not have a distance; cannot '
                             'compute 3d separation.')
        if issubclass(other.data.__class__, UnitSphericalRepresentation):
            raise ValueError('The other object does not have a distance; '
                             'cannot compute 3d separation.')

        return Distance((self.cartesian - other.cartesian).norm())

    def spherical_offsets_to(self, tocoord):
        r"""
        Computes angular offsets to go *from* this coordinate *to* another.

        Parameters
        ----------
        tocoord : `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate to offset to.

        Returns
        -------
        lon_offset : `~astropy.coordinates.Angle`
            The angular offset in the longitude direction (i.e., RA for
            equatorial coordinates).
        lat_offset : `~astropy.coordinates.Angle`
            The angular offset in the latitude direction (i.e., Dec for
            equatorial coordinates).

        Raises
        ------
        ValueError
            If the ``tocoord`` is not in the same frame as this one. This is
            different from the behavior of the `separation`/`separation_3d`
            methods because the offset components depend critically on the
            specific choice of frame.

        Notes
        -----
        This uses the sky offset frame machinery, and hence will produce a new
        sky offset frame if one does not already exist for this object's frame
        class.

        See Also
        --------
        separation : for the *total* angular offset (not broken out into components)

        """
        if not self.is_equivalent_frame(tocoord):
            raise ValueError('Tried to use spherical_offsets_to with two non-matching frames!')

        aframe = self.skyoffset_frame()
        acoord = tocoord.transform_to(aframe)

        dlon = acoord.spherical.lon.view(Angle)
        dlat = acoord.spherical.lat.view(Angle)
        return dlon, dlat

    def match_to_catalog_sky(self, catalogcoord, nthneighbor=1):
        """
        Finds the nearest on-sky matches of this coordinate in a set of
        catalog coordinates.

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

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
            object. Unless both this and ``catalogcoord`` have associated
            distances, this quantity assumes that all sources are at a
            distance of 1 (dimensionless).

        Notes
        -----
        This method requires `SciPy <https://www.scipy.org/>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_sky
        SkyCoord.match_to_catalog_3d
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

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

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
        This method requires `SciPy <https://www.scipy.org/>`_ to be
        installed or it will fail.

        See Also
        --------
        astropy.coordinates.match_coordinates_3d
        SkyCoord.match_to_catalog_sky
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

        This is intended for use on `~astropy.coordinates.SkyCoord` objects
        with coordinate arrays, rather than a scalar coordinate.  For a scalar
        coordinate, it is better to use
        `~astropy.coordinates.SkyCoord.separation`.

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

        Parameters
        ----------
        searcharoundcoords : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinates to search around to try to find matching points in
            this `SkyCoord`. This should be an object with array coordinates,
            not a scalar coordinate object.
        seplimit : `~astropy.units.Quantity` with angle units
            The on-sky separation to search within.

        Returns
        -------
        idxsearcharound : integer array
            Indices into ``self`` that matches to the corresponding element of
            ``idxself``. Shape matches ``idxself``.
        idxself : integer array
            Indices into ``searcharoundcoords`` that matches to the
            corresponding element of ``idxsearcharound``. Shape matches
            ``idxsearcharound``.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.

        Notes
        -----
        This method requires `SciPy <https://www.scipy.org/>`_ (>=0.12.0) to be
        installed or it will fail.

        In the current implementation, the return values are always sorted in
        the same order as the ``searcharoundcoords`` (so ``idxsearcharound`` is
        in ascending order).  This is considered an implementation detail,
        though, so it could change in a future release.

        See Also
        --------
        astropy.coordinates.search_around_sky
        SkyCoord.search_around_3d
        """
        from .matching import search_around_sky

        return search_around_sky(searcharoundcoords, self, seplimit,
                                 storekdtree='_kdtree_sky')

    def search_around_3d(self, searcharoundcoords, distlimit):
        """
        Searches for all coordinates in this object around a supplied set of
        points within a given 3D radius.

        This is intended for use on `~astropy.coordinates.SkyCoord` objects
        with coordinate arrays, rather than a scalar coordinate.  For a scalar
        coordinate, it is better to use
        `~astropy.coordinates.SkyCoord.separation_3d`.

        For more on how to use this (and related) functionality, see the
        examples in :doc:`/coordinates/matchsep`.

        Parameters
        ----------
        searcharoundcoords : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinates to search around to try to find matching points in
            this `SkyCoord`. This should be an object with array coordinates,
            not a scalar coordinate object.
        distlimit : `~astropy.units.Quantity` with distance units
            The physical radius to search within.

        Returns
        -------
        idxsearcharound : integer array
            Indices into ``self`` that matches to the corresponding element of
            ``idxself``. Shape matches ``idxself``.
        idxself : integer array
            Indices into ``searcharoundcoords`` that matches to the
            corresponding element of ``idxsearcharound``. Shape matches
            ``idxsearcharound``.
        sep2d : `~astropy.coordinates.Angle`
            The on-sky separation between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.
        dist3d : `~astropy.units.Quantity`
            The 3D distance between the coordinates. Shape matches
            ``idxsearcharound`` and ``idxself``.

        Notes
        -----
        This method requires `SciPy <https://www.scipy.org/>`_ (>=0.12.0) to be
        installed or it will fail.

        In the current implementation, the return values are always sorted in
        the same order as the ``searcharoundcoords`` (so ``idxsearcharound`` is
        in ascending order).  This is considered an implementation detail,
        though, so it could change in a future release.

        See Also
        --------
        astropy.coordinates.search_around_3d
        SkyCoord.search_around_sky
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

        if not self.is_equivalent_frame(other):
            try:
                other = other.transform_to(self, merge_attributes=False)
            except TypeError:
                raise TypeError('Can only get position_angle to another '
                                'SkyCoord or a coordinate frame with data')

        slat = self.represent_as(UnitSphericalRepresentation).lat
        slon = self.represent_as(UnitSphericalRepresentation).lon
        olat = other.represent_as(UnitSphericalRepresentation).lat
        olon = other.represent_as(UnitSphericalRepresentation).lon

        return angle_utilities.position_angle(slon, slat, olon, olat)

    def skyoffset_frame(self, rotation=None):
        """
        Returns the sky offset frame with this `SkyCoord` at the origin.

        Returns
        -------
        astrframe : `~astropy.coordinates.SkyOffsetFrame`
            A sky offset frame of the same type as this `SkyCoord` (e.g., if
            this object has an ICRS coordinate, the resulting frame is
            SkyOffsetICRS, with the origin set to this object)
        rotation : `~astropy.coordinates.Angle` or `~astropy.units.Quantity` with angle units
            The final rotation of the frame about the ``origin``. The sign of
            the rotation is the left-hand rule. That is, an object at a
            particular position angle in the un-rotated system will be sent to
            the positive latitude (z) direction in the final frame.
        """
        return SkyOffsetFrame(origin=self, rotation=rotation)

    def get_constellation(self, short_name=False, constellation_list='iau'):
        """
        Determines the constellation(s) of the coordinates this `SkyCoord`
        contains.

        Parameters
        ----------
        short_name : bool
            If True, the returned names are the IAU-sanctioned abbreviated
            names.  Otherwise, full names for the constellations are used.
        constellation_list : str
            The set of constellations to use.  Currently only ``'iau'`` is
            supported, meaning the 88 "modern" constellations endorsed by the IAU.

        Returns
        -------
        constellation : str or string array
            If this is a scalar coordinate, returns the name of the
            constellation.  If it is an array `SkyCoord`, it returns an array of
            names.

        Notes
        -----
        To determine which constellation a point on the sky is in, this first
        precesses to B1875, and then uses the Delporte boundaries of the 88
        modern constellations, as tabulated by
        `Roman 1987 <http://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/42>`_.

        See Also
        --------
        astropy.coordinates.get_constellation
        """
        from .funcs import get_constellation

        # because of issue #7028, the conversion to a PrecessedGeocentric
        # system fails in some cases.  Work around is to  drop the velocities.
        # they are not needed here since only position infromation is used
        extra_frameattrs = {nm: getattr(self, nm)
                            for nm in self._extra_frameattr_names}
        novel = SkyCoord(self.realize_frame(self.data.without_differentials()),
                         **extra_frameattrs)
        return get_constellation(novel, short_name, constellation_list)

        # the simpler version below can be used when gh-issue #7028 is resolved
        #return get_constellation(self, short_name, constellation_list)

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

    def radial_velocity_correction(self, kind='barycentric', obstime=None,
                                   location=None):
        """
        Compute the correction required to convert a radial velocity at a given
        time and place on the Earth's Surface to a barycentric or heliocentric
        velocity.

        Parameters
        ----------
        kind : str
            The kind of velocity correction.  Must be 'barycentric' or
            'heliocentric'.
        obstime : `~astropy.time.Time` or None, optional
            The time at which to compute the correction.  If `None`, the
            ``obstime`` frame attribute on the `SkyCoord` will be used.
        location : `~astropy.coordinates.EarthLocation` or None, optional
            The observer location at which to compute the correction.  If
            `None`, the  ``location`` frame attribute on the passed-in
            ``obstime`` will be used, and if that is None, the ``location``
            frame attribute on the `SkyCoord` will be used.

        Raises
        ------
        ValueError
            If either ``obstime`` or ``location`` are passed in (not ``None``)
            when the frame attribute is already set on this `SkyCoord`.
        TypeError
            If ``obstime`` or ``location`` aren't provided, either as arguments
            or as frame attributes.

        Returns
        -------
        vcorr : `~astropy.units.Quantity` with velocity units
            The  correction with a positive sign.  I.e., *add* this
            to an observed radial velocity to get the barycentric (or
            heliocentric) velocity. If m/s precision or better is needed,
            see the notes below.

        Notes
        -----
        The barycentric correction is calculated to higher precision than the
        heliocentric correction and includes additional physics (e.g time dilation).
        Use barycentric corrections if m/s precision is required.

        The algorithm here is sufficient to perform corrections at the mm/s level, but
        care is needed in application. Strictly speaking, the barycentric correction is
        multiplicative and should be applied as::

           sc = SkyCoord(1*u.deg, 2*u.deg)
           vcorr = sc.rv_correction(kind='barycentric', obstime=t, location=loc)
           rv = rv + vcorr + rv * vcorr / consts.c

        If your target is nearby and/or has finite proper motion you may need to account
        for terms arising from this. See Wright & Eastmann (2014) for details.

        The default is for this method to use the builtin ephemeris for
        computing the sun and earth location.  Other ephemerides can be chosen
        by setting the `~astropy.coordinates.solar_system_ephemeris` variable,
        either directly or via ``with`` statement.  For example, to use the JPL
        ephemeris, do::

            sc = SkyCoord(1*u.deg, 2*u.deg)
            with coord.solar_system_ephemeris.set('jpl'):
                rv += sc.rv_correction(obstime=t, location=loc)

        """
        # has to be here to prevent circular imports
        from .solar_system import get_body_barycentric_posvel, get_body_barycentric

        # location validation
        timeloc = getattr(obstime, 'location', None)
        if location is None:
            if self.location is not None:
                location = self.location
                if timeloc is not None:
                    raise ValueError('`location` cannot be in both the '
                                     'passed-in `obstime` and this `SkyCoord` '
                                     'because it is ambiguous which is meant '
                                     'for the radial_velocity_correction.')
            elif timeloc is not None:
                location = timeloc
            else:
                raise TypeError('Must provide a `location` to '
                                'radial_velocity_correction, either as a '
                                'SkyCoord frame attribute, as an attribute on '
                                'the passed in `obstime`, or in the method '
                                'call.')

        elif self.location is not None or timeloc is not None:
            raise ValueError('Cannot compute radial velocity correction if '
                             '`location` argument is passed in and there is '
                             'also a  `location` attribute on this SkyCoord or '
                             'the passed-in `obstime`.')

        # obstime validation
        if obstime is None:
            obstime = self.obstime
            if obstime is None:
                raise TypeError('Must provide an `obstime` to '
                                'radial_velocity_correction, either as a '
                                'SkyCoord frame attribute or in the method '
                                'call.')
        elif self.obstime is not None:
            raise ValueError('Cannot compute radial velocity correction if '
                             '`obstime` argument is passed in and it is '
                             'inconsistent with the `obstime` frame '
                             'attribute on the SkyCoord')

        pos_earth, v_earth = get_body_barycentric_posvel('earth', obstime)
        if kind == 'barycentric':
            v_origin_to_earth = v_earth
        elif kind == 'heliocentric':
            v_sun = get_body_barycentric_posvel('sun', obstime)[1]
            v_origin_to_earth = v_earth - v_sun
        else:
            raise ValueError("`kind` argument to radial_velocity_correction must "
                             "be 'barycentric' or 'heliocentric', but got "
                             "'{}'".format(kind))

        gcrs_p, gcrs_v = location.get_gcrs_posvel(obstime)
        # transforming to GCRS is not the correct thing to do here, since we don't want to
        # include aberration (or light deflection)? Instead, only apply parallax if necessary
        if self.data.__class__ is UnitSphericalRepresentation:
            targcart = self.icrs.cartesian
        else:
            # skycoord has distances so apply parallax
            obs_icrs_cart = pos_earth + gcrs_p
            icrs_cart = self.icrs.cartesian
            targcart = icrs_cart - obs_icrs_cart
            targcart /= targcart.norm()

        if kind == 'barycentric':
            beta_obs = (v_origin_to_earth + gcrs_v) / speed_of_light
            gamma_obs = 1 / np.sqrt(1 - beta_obs.norm()**2)
            gr = location.gravitational_redshift(obstime)
            # barycentric redshift according to eq 28 in Wright & Eastmann (2014),
            # neglecting Shapiro delay and effects of the star's own motion
            zb = gamma_obs * (1 + targcart.dot(beta_obs)) / (1 + gr/speed_of_light) - 1
            return zb * speed_of_light
        else:
            # do a simpler correction ignoring time dilation and gravitational redshift
            # this is adequate since Heliocentric corrections shouldn't be used if
            # cm/s precision is required.
            return targcart.dot(v_origin_to_earth + gcrs_v)

    # Table interactions
    @classmethod
    def guess_from_table(cls, table, **coord_kwargs):
        r"""
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
        are ``<space>_!"#$%&'()*+,-./\:;<=>?@[]^`{|}~``).

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
            # the final regex ORs together the two patterns
            rex = re.compile('(' + starts_with_comp + ')|(' + ends_with_comp + ')',
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

    if isinstance(frame, str):
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
        # Copy any extra attributes if they are not explicitly given.
        for attr in frame._extra_frameattr_names:
            kwargs.setdefault(attr, getattr(frame, attr))
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
        # this catches the "single list passed in" case.  For that case we want
        # to allow the first argument to set the class.  That's OK because
        # _parse_coordinate_arg goes and checks that the frames match between
        # the first and all the others
        if (isinstance(arg, (collections.Sequence, np.ndarray)) and
             len(args) == 1 and len(arg) > 0):
            arg = arg[0]

        coord_frame_obj = coord_frame_cls = None
        if isinstance(arg, BaseCoordinateFrame):
            coord_frame_obj = arg
        elif isinstance(arg, SkyCoord):
            coord_frame_obj = arg.frame
        if coord_frame_obj is not None:
            coord_frame_cls = coord_frame_obj.__class__
            frame_diff = coord_frame_obj.get_representation_cls('s')
            if frame_diff is not None:
                # we do this check because otherwise if there's no default
                # differential (i.e. it is None), the code below chokes. but
                # None still gets through if the user *requests* it
                kwargs.setdefault('differential_type', frame_diff)

        if coord_frame_cls is not None:
            if not frame_specified_explicitly:
                frame_cls = coord_frame_cls
            elif frame_cls is not coord_frame_cls:
                raise ValueError("Cannot override frame='{0}' of input coordinate with "
                                 "new frame='{1}'.  Instead transform the coordinate."
                                 .format(coord_frame_cls.__name__, frame_cls.__name__))

    frame_cls_kwargs = {}

    # TODO: deprecate representation, remove this in future
    _normalize_representation_type(kwargs)

    if 'representation_type' in kwargs:
        frame_cls_kwargs['representation_type'] = _get_repr_cls(
            kwargs['representation_type'])

    if 'differential_type' in kwargs:
        frame_cls_kwargs['differential_type'] = _get_diff_cls(
            kwargs['differential_type'])

    return frame_cls(**frame_cls_kwargs)


def _get_representation_component_units(args, kwargs):
    """
    Get the unit from kwargs for the *representation* components (not the
    differentials).
    """
    if 'unit' not in kwargs:
        units = [None, None, None]

    else:
        units = kwargs.pop('unit')

        if isinstance(units, str):
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
        except Exception:
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

    frame_attr_names = list(frame.representation_component_names.keys())
    repr_attr_names = list(frame.representation_component_names.values())
    repr_attr_classes = list(frame.representation.attr_classes.values())
    n_attr_names = len(repr_attr_names)

    # Turn a single string into a list of strings for convenience
    if isinstance(coords, str):
        is_scalar = True
        coords = [coords]

    if isinstance(coords, (SkyCoord, BaseCoordinateFrame)):
        # Note that during parsing of `frame` it is checked that any coordinate
        # args have the same frame as explicitly supplied, so don't worry here.

        if not coords.has_data:
            raise ValueError('Cannot initialize from a frame without coordinate data')

        data = coords.data.represent_as(frame.representation_type)

        values = []  # List of values corresponding to representation attrs
        repr_attr_name_to_drop = []
        for repr_attr_name in repr_attr_names:
            # If coords did not have an explicit distance then don't include in initializers.
            if (isinstance(coords.data, UnitSphericalRepresentation) and
                    repr_attr_name == 'distance'):
                repr_attr_name_to_drop.append(repr_attr_name)
                continue

            # Get the value from `data` in the eventual representation
            values.append(getattr(data, repr_attr_name))

        # drop the ones that were skipped because they were distances
        for nametodrop in repr_attr_name_to_drop:
            nameidx = repr_attr_names.index(nametodrop)
            del repr_attr_names[nameidx]
            del units[nameidx]
            del frame_attr_names[nameidx]
            del repr_attr_classes[nameidx]

        if coords.data.differentials and 's' in coords.data.differentials:
            orig_vel = coords.data.differentials['s']
            vel = coords.data.represent_as(frame.representation, frame.get_representation_cls('s')).differentials['s']
            for frname, reprname in frame.get_representation_component_names('s').items():
                if (reprname == 'd_distance' and not hasattr(orig_vel, reprname) and
                    'unit' in orig_vel.get_name()):
                    continue
                values.append(getattr(vel, reprname))
                units.append(None)
                frame_attr_names.append(frname)
                repr_attr_names.append(reprname)
                repr_attr_classes.append(vel.attr_classes[reprname])

        for attr in frame_transform_graph.frame_attributes:
            value = getattr(coords, attr, None)
            use_value = (isinstance(coords, SkyCoord)
                         or attr not in coords._attr_names_with_defaults)
            if use_value and value is not None:
                valid_kwargs[attr] = value

    elif isinstance(coords, BaseRepresentation):
        if coords.differentials and 's' in coords.differentials:
            diffs = frame.get_representation_cls('s')
            data = coords.represent_as(frame.representation_type, diffs)
            values = [getattr(data, repr_attr_name) for repr_attr_name in repr_attr_names]
            for frname, reprname in frame.get_representation_component_names('s').items():
                values.append(getattr(data.differentials['s'], reprname))
                units.append(None)
                frame_attr_names.append(frname)
                repr_attr_names.append(reprname)
                repr_attr_classes.append(data.differentials['s'].attr_classes[reprname])

        else:
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

            # Check that all frames are equivalent
            for sc in scs[1:]:
                if not sc.is_equivalent_frame(scs[0]):
                    raise ValueError("List of inputs don't have equivalent "
                                     "frames: {0} != {1}".format(sc, scs[0]))

            # Now use the first to determine if they are all UnitSpherical
            allunitsphrepr = isinstance(scs[0].data, UnitSphericalRepresentation)

            # get the frame attributes from the first coord in the list, because
            # from the above we know it matches all the others.  First copy over
            # the attributes that are in the frame itself, then copy over any
            # extras in the SkyCoord
            for fattrnm in scs[0].frame.frame_attributes:
                valid_kwargs[fattrnm] = getattr(scs[0].frame, fattrnm)
            for fattrnm in scs[0]._extra_frameattr_names:
                valid_kwargs[fattrnm] = getattr(scs[0], fattrnm)

            # Now combine the values, to be used below
            values = []
            for data_attr_name, repr_attr_name in zip(frame_attr_names, repr_attr_names):
                if allunitsphrepr and repr_attr_name == 'distance':
                    # if they are *all* UnitSpherical, don't give a distance
                    continue
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
            # none of the elements are "frame-like"
            # turn into a list of lists like [[v1_0, v2_0, v3_0], ... [v1_N, v2_N, v3_N]]
            for coord in coords:
                if isinstance(coord, str):
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
            except Exception:
                raise ValueError('One or more elements of input sequence does not have a length')

            if len(n_coords) > 1:
                raise ValueError('Input coordinate values must have same number of elements, found {0}'
                                 .format(n_coords))
            n_coords = n_coords[0]

            # Must have no more coord inputs than representation attributes
            if n_coords > n_attr_names:
                raise ValueError('Input coordinates have {0} values but '
                                 'representation {1} only accepts {2}'
                                 .format(n_coords,
                                         frame.representation_type.get_name(),
                                         n_attr_names))

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
            valid_kwargs[frame_attr_name] = repr_attr_class(value, unit=unit,
                                                            copy=False)
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

    This also gets any *differential* kwargs, because they go into the same
    frame initializer later on.
    """
    frame_attr_names = frame.representation_component_names.keys()
    repr_attr_classes = frame.representation_type.attr_classes.values()

    valid_kwargs = {}
    for frame_attr_name, repr_attr_class, unit in zip(frame_attr_names, repr_attr_classes, units):
        value = kwargs.pop(frame_attr_name, None)
        if value is not None:
            valid_kwargs[frame_attr_name] = repr_attr_class(value, unit=unit)

    # also check the differentials.  They aren't included in the units keyword,
    # so we only look for the names.

    differential_type = frame.differential_type
    if differential_type is not None:
        for frame_name, repr_name in frame.get_representation_component_names('s').items():
            diff_attr_class = differential_type.attr_classes[repr_name]
            value = kwargs.pop(frame_name, None)
            if value is not None:
                valid_kwargs[frame_name] = diff_attr_class(value)

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

    if isinstance(coord_str, str):
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
