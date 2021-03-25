# Licensed under a 3-clause BSD style license - see LICENSE.rst

from warnings import warn
import collections
import socket
import json
import urllib.request
import urllib.error
import urllib.parse

import numpy as np
import erfa

from astropy import units as u
from astropy import constants as consts
from astropy.units.quantity import QuantityInfoBase
from astropy.utils.exceptions import AstropyUserWarning
from .angles import Angle, Longitude, Latitude
from .representation import CartesianRepresentation, CartesianDifferential
from .errors import UnknownSiteException
from astropy.utils import data


__all__ = ['EarthLocation']

GeodeticLocation = collections.namedtuple('GeodeticLocation', ['lon', 'lat', 'height'])

# Available ellipsoids (defined in erfam.h, with numbers exposed in erfa).
ELLIPSOIDS = ('WGS84', 'GRS80', 'WGS72')

OMEGA_EARTH = u.Quantity(7.292115855306589e-5, 1./u.s)
"""
Rotational velocity of Earth. In UT1 seconds, this would be 2 pi / (24 * 3600),
but we need the value in SI seconds.
See Explanatory Supplement to the Astronomical Almanac, ed. P. Kenneth Seidelmann (1992),
University Science Books.
"""


def _check_ellipsoid(ellipsoid=None, default='WGS84'):
    if ellipsoid is None:
        ellipsoid = default
    if ellipsoid not in ELLIPSOIDS:
        raise ValueError('Ellipsoid {} not among known ones ({})'
                         .format(ellipsoid, ELLIPSOIDS))
    return ellipsoid


def _get_json_result(url, err_str, use_google):

    # need to do this here to prevent a series of complicated circular imports
    from .name_resolve import NameResolveError
    try:
        # Retrieve JSON response from Google maps API
        resp = urllib.request.urlopen(url, timeout=data.conf.remote_timeout)
        resp_data = json.loads(resp.read().decode('utf8'))

    except urllib.error.URLError as e:
        # This catches a timeout error, see:
        #   http://stackoverflow.com/questions/2712524/handling-urllib2s-timeout-python
        if isinstance(e.reason, socket.timeout):
            raise NameResolveError(err_str.format(msg="connection timed out"))
        else:
            raise NameResolveError(err_str.format(msg=e.reason))

    except socket.timeout:
        # There are some cases where urllib2 does not catch socket.timeout
        # especially while receiving response data on an already previously
        # working request
        raise NameResolveError(err_str.format(msg="connection timed out"))

    if use_google:
        results = resp_data.get('results', [])

        if resp_data.get('status', None) != 'OK':
            raise NameResolveError(err_str.format(msg="unknown failure with "
                                                  "Google API"))

    else:  # OpenStreetMap returns a list
        results = resp_data

    if not results:
        raise NameResolveError(err_str.format(msg="no results returned"))

    return results


class EarthLocationInfo(QuantityInfoBase):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """
    _represent_as_dict_attrs = ('x', 'y', 'z', 'ellipsoid')

    def _construct_from_dict(self, map):
        # Need to pop ellipsoid off and update post-instantiation.  This is
        # on the to-fix list in #4261.
        ellipsoid = map.pop('ellipsoid')
        out = self._parent_cls(**map)
        out.ellipsoid = ellipsoid
        return out

    def new_like(self, cols, length, metadata_conflicts='warn', name=None):
        """
        Return a new EarthLocation instance which is consistent with the
        input ``cols`` and has ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        cols : list
            List of input columns
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : EarthLocation (or subclass)
            Empty instance of this class consistent with ``cols``
        """
        # Very similar to QuantityInfo.new_like, but the creation of the
        # map is different enough that this needs its own rouinte.
        # Get merged info attributes shape, dtype, format, description.
        attrs = self.merge_cols_attributes(cols, metadata_conflicts, name,
                                           ('meta', 'format', 'description'))
        # The above raises an error if the dtypes do not match, but returns
        # just the string representation, which is not useful, so remove.
        attrs.pop('dtype')
        # Make empty EarthLocation using the dtype and unit of the last column.
        # Use zeros so we do not get problems for possible conversion to
        # geodetic coordinates.
        shape = (length,) + attrs.pop('shape')
        data = u.Quantity(np.zeros(shape=shape, dtype=cols[0].dtype),
                          unit=cols[0].unit, copy=False)
        # Get arguments needed to reconstruct class
        map = {key: (data[key] if key in 'xyz' else getattr(cols[-1], key))
               for key in self._represent_as_dict_attrs}
        out = self._construct_from_dict(map)
        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out


class EarthLocation(u.Quantity):
    """
    Location on the Earth.

    Initialization is first attempted assuming geocentric (x, y, z) coordinates
    are given; if that fails, another attempt is made assuming geodetic
    coordinates (longitude, latitude, height above a reference ellipsoid).
    When using the geodetic forms, Longitudes are measured increasing to the
    east, so west longitudes are negative. Internally, the coordinates are
    stored as geocentric.

    To ensure a specific type of coordinates is used, use the corresponding
    class methods (`from_geocentric` and `from_geodetic`) or initialize the
    arguments with names (``x``, ``y``, ``z`` for geocentric; ``lon``, ``lat``,
    ``height`` for geodetic).  See the class methods for details.


    Notes
    -----
    This class fits into the coordinates transformation framework in that it
    encodes a position on the `~astropy.coordinates.ITRS` frame.  To get a
    proper `~astropy.coordinates.ITRS` object from this object, use the ``itrs``
    property.
    """

    _ellipsoid = 'WGS84'
    _location_dtype = np.dtype({'names': ['x', 'y', 'z'],
                                'formats': [np.float64]*3})
    _array_dtype = np.dtype((np.float64, (3,)))

    info = EarthLocationInfo()

    def __new__(cls, *args, **kwargs):
        # TODO: needs copy argument and better dealing with inputs.
        if (len(args) == 1 and len(kwargs) == 0 and
                isinstance(args[0], EarthLocation)):
            return args[0].copy()
        try:
            self = cls.from_geocentric(*args, **kwargs)
        except (u.UnitsError, TypeError) as exc_geocentric:
            try:
                self = cls.from_geodetic(*args, **kwargs)
            except Exception as exc_geodetic:
                raise TypeError('Coordinates could not be parsed as either '
                                'geocentric or geodetic, with respective '
                                'exceptions "{}" and "{}"'
                                .format(exc_geocentric, exc_geodetic))
        return self

    @classmethod
    def from_geocentric(cls, x, y, z, unit=None):
        """
        Location on Earth, initialized from geocentric coordinates.

        Parameters
        ----------
        x, y, z : `~astropy.units.Quantity` or array_like
            Cartesian coordinates.  If not quantities, ``unit`` should be given.
        unit : `~astropy.units.UnitBase` object or None
            Physical unit of the coordinate values.  If ``x``, ``y``, and/or
            ``z`` are quantities, they will be converted to this unit.

        Raises
        ------
        astropy.units.UnitsError
            If the units on ``x``, ``y``, and ``z`` do not match or an invalid
            unit is given.
        ValueError
            If the shapes of ``x``, ``y``, and ``z`` do not match.
        TypeError
            If ``x`` is not a `~astropy.units.Quantity` and no unit is given.
        """
        if unit is None:
            try:
                unit = x.unit
            except AttributeError:
                raise TypeError("Geocentric coordinates should be Quantities "
                                "unless an explicit unit is given.")
        else:
            unit = u.Unit(unit)

        if unit.physical_type != 'length':
            raise u.UnitsError("Geocentric coordinates should be in "
                               "units of length.")

        try:
            x = u.Quantity(x, unit, copy=False)
            y = u.Quantity(y, unit, copy=False)
            z = u.Quantity(z, unit, copy=False)
        except u.UnitsError:
            raise u.UnitsError("Geocentric coordinate units should all be "
                               "consistent.")

        x, y, z = np.broadcast_arrays(x, y, z)
        struc = np.empty(x.shape, cls._location_dtype)
        struc['x'], struc['y'], struc['z'] = x, y, z
        return super().__new__(cls, struc, unit, copy=False)

    @classmethod
    def from_geodetic(cls, lon, lat, height=0., ellipsoid=None):
        """
        Location on Earth, initialized from geodetic coordinates.

        Parameters
        ----------
        lon : `~astropy.coordinates.Longitude` or float
            Earth East longitude.  Can be anything that initialises an
            `~astropy.coordinates.Angle` object (if float, in degrees).
        lat : `~astropy.coordinates.Latitude` or float
            Earth latitude.  Can be anything that initialises an
            `~astropy.coordinates.Latitude` object (if float, in degrees).
        height : `~astropy.units.Quantity` or float, optional
            Height above reference ellipsoid (if float, in meters; default: 0).
        ellipsoid : str, optional
            Name of the reference ellipsoid to use (default: 'WGS84').
            Available ellipsoids are:  'WGS84', 'GRS80', 'WGS72'.

        Raises
        ------
        astropy.units.UnitsError
            If the units on ``lon`` and ``lat`` are inconsistent with angular
            ones, or that on ``height`` with a length.
        ValueError
            If ``lon``, ``lat``, and ``height`` do not have the same shape, or
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        For the conversion to geocentric coordinates, the ERFA routine
        ``gd2gc`` is used.  See https://github.com/liberfa/erfa
        """
        ellipsoid = _check_ellipsoid(ellipsoid, default=cls._ellipsoid)
        # We use Angle here since there is no need to wrap the longitude -
        # gd2gc will just take cos/sin anyway.  And wrapping might fail
        # on readonly input.
        lon = Angle(lon, u.degree, copy=False)
        lat = Latitude(lat, u.degree, copy=False)
        # don't convert to m by default, so we can use the height unit below.
        if not isinstance(height, u.Quantity):
            height = u.Quantity(height, u.m, copy=False)
        # get geocentric coordinates. Have to give one-dimensional array.
        xyz = erfa.gd2gc(getattr(erfa, ellipsoid),
                         lon.to_value(u.radian),
                         lat.to_value(u.radian),
                         height.to_value(u.m))
        self = xyz.ravel().view(cls._location_dtype,
                                cls).reshape(xyz.shape[:-1])
        self._unit = u.meter
        self._ellipsoid = ellipsoid
        return self.to(height.unit)

    @classmethod
    def of_site(cls, site_name):
        """
        Return an object of this class for a known observatory/site by name.

        This is intended as a quick convenience function to get basic site
        information, not a fully-featured exhaustive registry of observatories
        and all their properties.

        Additional information about the site is stored in the ``.info.meta``
        dictionary of sites obtained using this method (see the examples below).

        .. note::
            When this function is called, it will attempt to download site
            information from the astropy data server. If you would like a site
            to be added, issue a pull request to the
            `astropy-data repository <https://github.com/astropy/astropy-data>`_ .
            If a site cannot be found in the registry (i.e., an internet
            connection is not available), it will fall back on a built-in list,
            In the future, this bundled list might include a version-controlled
            list of canonical observatories extracted from the online version,
            but it currently only contains the Greenwich Royal Observatory as an
            example case.


        Parameters
        ----------
        site_name : str
            Name of the observatory (case-insensitive).

        Returns
        -------
        site : This class (a `~astropy.coordinates.EarthLocation` or subclass)
            The location of the observatory.

        Examples
        --------

        >>> from astropy.coordinates import EarthLocation
        >>> keck = EarthLocation.of_site('Keck Observatory')  # doctest: +REMOTE_DATA
        >>> keck.geodetic  # doctest: +REMOTE_DATA +FLOAT_CMP
        GeodeticLocation(lon=<Longitude -155.47833333 deg>, lat=<Latitude 19.82833333 deg>, height=<Quantity 4160. m>)
        >>> keck.info  # doctest: +REMOTE_DATA
        name = W. M. Keck Observatory
        dtype = void192
        unit = m
        class = EarthLocation
        n_bad = 0
        >>> keck.info.meta  # doctest: +REMOTE_DATA
        {'source': 'IRAF Observatory Database', 'timezone': 'US/Aleutian'}

        See Also
        --------
        get_site_names : the list of sites that this function can access
        """  # noqa
        registry = cls._get_site_registry()
        try:
            el = registry[site_name]
        except UnknownSiteException as e:
            raise UnknownSiteException(e.site, 'EarthLocation.get_site_names',
                                       close_names=e.close_names)

        if cls is el.__class__:
            return el
        else:
            newel = cls.from_geodetic(*el.to_geodetic())
            newel.info.name = el.info.name
            return newel

    @classmethod
    def of_address(cls, address, get_height=False, google_api_key=None):
        """
        Return an object of this class for a given address by querying either
        the OpenStreetMap Nominatim tool [1]_ (default) or the Google geocoding
        API [2]_, which requires a specified API key.

        This is intended as a quick convenience function to get easy access to
        locations. If you need to specify a precise location, you should use the
        initializer directly and pass in a longitude, latitude, and elevation.

        In the background, this just issues a web query to either of
        the APIs noted above. This is not meant to be abused! Both
        OpenStreetMap and Google use IP-based query limiting and will ban your
        IP if you send more than a few thousand queries per hour [2]_.

        .. warning::
            If the query returns more than one location (e.g., searching on
            ``address='springfield'``), this function will use the **first**
            returned location.

        Parameters
        ----------
        address : str
            The address to get the location for. As per the Google maps API,
            this can be a fully specified street address (e.g., 123 Main St.,
            New York, NY) or a city name (e.g., Danbury, CT), or etc.
        get_height : bool, optional
            This only works when using the Google API! See the ``google_api_key``
            block below. Use the retrieved location to perform a second query to
            the Google maps elevation API to retrieve the height of the input
            address [3]_.
        google_api_key : str, optional
            A Google API key with the Geocoding API and (optionally) the
            elevation API enabled. See [4]_ for more information.


        Returns
        -------
        location : This class (a `~astropy.coordinates.EarthLocation` or subclass)
            The location of the input address.

        References
        ----------
        .. [1] https://nominatim.openstreetmap.org/
        .. [2] https://developers.google.com/maps/documentation/geocoding/start
        .. [3] https://developers.google.com/maps/documentation/elevation/start
        .. [4] https://developers.google.com/maps/documentation/geocoding/get-api-key

        """

        use_google = google_api_key is not None

        # Fail fast if invalid options are passed:
        if not use_google and get_height:
            raise ValueError(
                'Currently, `get_height` only works when using '
                'the Google geocoding API, which requires passing '
                'a Google API key with `google_api_key`. See: '
                'https://developers.google.com/maps/documentation/geocoding/get-api-key '
                'for information on obtaining an API key.')

        if use_google:  # Google
            pars = urllib.parse.urlencode({'address': address,
                                           'key': google_api_key})
            geo_url = ("https://maps.googleapis.com/maps/api/geocode/json?{}"
                       .format(pars))

        else:  # OpenStreetMap
            pars = urllib.parse.urlencode({'q': address,
                                           'format': 'json'})
            geo_url = ("https://nominatim.openstreetmap.org/search?{}"
                       .format(pars))

        # get longitude and latitude location
        err_str = ("Unable to retrieve coordinates for address '{address}'; "
                   "{{msg}}".format(address=address))
        geo_result = _get_json_result(geo_url, err_str=err_str,
                                      use_google=use_google)

        if use_google:
            loc = geo_result[0]['geometry']['location']
            lat = loc['lat']
            lon = loc['lng']

        else:
            loc = geo_result[0]
            lat = float(loc['lat'])  # strings are returned by OpenStreetMap
            lon = float(loc['lon'])

        if get_height:
            pars = {'locations': f'{lat:.8f},{lon:.8f}',
                    'key': google_api_key}
            pars = urllib.parse.urlencode(pars)
            ele_url = ("https://maps.googleapis.com/maps/api/elevation/json?{}"
                       .format(pars))

            err_str = ("Unable to retrieve elevation for address '{address}'; "
                       "{{msg}}".format(address=address))
            ele_result = _get_json_result(ele_url, err_str=err_str,
                                          use_google=use_google)
            height = ele_result[0]['elevation']*u.meter

        else:
            height = 0.

        return cls.from_geodetic(lon=lon*u.deg, lat=lat*u.deg, height=height)

    @classmethod
    def get_site_names(cls):
        """
        Get list of names of observatories for use with
        `~astropy.coordinates.EarthLocation.of_site`.

        .. note::
            When this function is called, it will first attempt to
            download site information from the astropy data server.  If it
            cannot (i.e., an internet connection is not available), it will fall
            back on the list included with astropy (which is a limited and dated
            set of sites).  If you think a site should be added, issue a pull
            request to the
            `astropy-data repository <https://github.com/astropy/astropy-data>`_ .


        Returns
        -------
        names : list of str
            List of valid observatory names

        See Also
        --------
        of_site : Gets the actual location object for one of the sites names
                  this returns.
        """
        return cls._get_site_registry().names

    @classmethod
    def _get_site_registry(cls, force_download=False, force_builtin=False):
        """
        Gets the site registry.  The first time this either downloads or loads
        from the data file packaged with astropy.  Subsequent calls will use the
        cached version unless explicitly overridden.

        Parameters
        ----------
        force_download : bool or str
            If not False, force replacement of the cached registry with a
            downloaded version. If a str, that will be used as the URL to
            download from (if just True, the default URL will be used).
        force_builtin : bool
            If True, load from the data file bundled with astropy and set the
            cache to that.

        returns
        -------
        reg : astropy.coordinates.sites.SiteRegistry
        """
        # need to do this here at the bottom to avoid circular dependencies
        from .sites import get_builtin_sites, get_downloaded_sites

        if force_builtin and force_download:
            raise ValueError('Cannot have both force_builtin and force_download True')

        if force_builtin:
            reg = cls._site_registry = get_builtin_sites()
        else:
            reg = getattr(cls, '_site_registry', None)
            if force_download or not reg:
                try:
                    if isinstance(force_download, str):
                        reg = get_downloaded_sites(force_download)
                    else:
                        reg = get_downloaded_sites()
                except OSError:
                    if force_download:
                        raise
                    msg = ('Could not access the online site list. Falling '
                           'back on the built-in version, which is rather '
                           'limited. If you want to retry the download, do '
                           '{0}._get_site_registry(force_download=True)')
                    warn(AstropyUserWarning(msg.format(cls.__name__)))
                    reg = get_builtin_sites()
                cls._site_registry = reg

        return reg

    @property
    def ellipsoid(self):
        """The default ellipsoid used to convert to geodetic coordinates."""
        return self._ellipsoid

    @ellipsoid.setter
    def ellipsoid(self, ellipsoid):
        self._ellipsoid = _check_ellipsoid(ellipsoid)

    @property
    def geodetic(self):
        """Convert to geodetic coordinates for the default ellipsoid."""
        return self.to_geodetic()

    def to_geodetic(self, ellipsoid=None):
        """Convert to geodetic coordinates.

        Parameters
        ----------
        ellipsoid : str, optional
            Reference ellipsoid to use.  Default is the one the coordinates
            were initialized with.  Available are: 'WGS84', 'GRS80', 'WGS72'

        Returns
        -------
        (lon, lat, height) : tuple
            The tuple contains instances of `~astropy.coordinates.Longitude`,
            `~astropy.coordinates.Latitude`, and `~astropy.units.Quantity`

        Raises
        ------
        ValueError
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        For the conversion to geodetic coordinates, the ERFA routine
        ``gc2gd`` is used.  See https://github.com/liberfa/erfa
        """
        ellipsoid = _check_ellipsoid(ellipsoid, default=self.ellipsoid)
        self_array = self.to(u.meter).view(self._array_dtype, np.ndarray)
        lon, lat, height = erfa.gc2gd(getattr(erfa, ellipsoid), self_array)
        return GeodeticLocation(
            Longitude(lon * u.radian, u.degree,
                      wrap_angle=180.*u.degree, copy=False),
            Latitude(lat * u.radian, u.degree, copy=False),
            u.Quantity(height * u.meter, self.unit, copy=False))

    @property
    def lon(self):
        """Longitude of the location, for the default ellipsoid."""
        return self.geodetic[0]

    @property
    def lat(self):
        """Latitude of the location, for the default ellipsoid."""
        return self.geodetic[1]

    @property
    def height(self):
        """Height of the location, for the default ellipsoid."""
        return self.geodetic[2]

    # mostly for symmetry with geodetic and to_geodetic.
    @property
    def geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return self.to_geocentric()

    def to_geocentric(self):
        """Convert to a tuple with X, Y, and Z as quantities"""
        return (self.x, self.y, self.z)

    def get_itrs(self, obstime=None):
        """
        Generates an `~astropy.coordinates.ITRS` object with the location of
        this object at the requested ``obstime``.

        Parameters
        ----------
        obstime : `~astropy.time.Time` or None
            The ``obstime`` to apply to the new `~astropy.coordinates.ITRS`, or
            if None, the default ``obstime`` will be used.

        Returns
        -------
        itrs : `~astropy.coordinates.ITRS`
            The new object in the ITRS frame
        """
        # Broadcast for a single position at multiple times, but don't attempt
        # to be more general here.
        if obstime and self.size == 1 and obstime.shape:
            self = np.broadcast_to(self, obstime.shape, subok=True)

        # do this here to prevent a series of complicated circular imports
        from .builtin_frames import ITRS
        return ITRS(x=self.x, y=self.y, z=self.z, obstime=obstime)

    itrs = property(get_itrs, doc="""An `~astropy.coordinates.ITRS` object  with
                                     for the location of this object at the
                                     default ``obstime``.""")

    def get_gcrs(self, obstime):
        """GCRS position with velocity at ``obstime`` as a GCRS coordinate.

        Parameters
        ----------
        obstime : `~astropy.time.Time`
            The ``obstime`` to calculate the GCRS position/velocity at.

        Returns
        --------
        gcrs : `~astropy.coordinates.GCRS` instance
            With velocity included.
        """
        # do this here to prevent a series of complicated circular imports
        from .builtin_frames import GCRS

        itrs = self.get_itrs(obstime)
        # Assume the observatory itself is fixed on the ground.
        # We do a direct assignment rather than an update to avoid validation
        # and creation of a new object.
        zeros = np.broadcast_to(0. * u.km / u.s, (3,) + itrs.shape, subok=True)
        itrs.data.differentials['s'] = CartesianDifferential(zeros)
        return itrs.transform_to(GCRS(obstime=obstime))

    def get_gcrs_posvel(self, obstime):
        """
        Calculate the GCRS position and velocity of this object at the
        requested ``obstime``.

        Parameters
        ----------
        obstime : `~astropy.time.Time`
            The ``obstime`` to calculate the GCRS position/velocity at.

        Returns
        --------
        obsgeoloc : `~astropy.coordinates.CartesianRepresentation`
            The GCRS position of the object
        obsgeovel : `~astropy.coordinates.CartesianRepresentation`
            The GCRS velocity of the object
        """
        # GCRS position
        gcrs_data = self.get_gcrs(obstime).data
        obsgeopos = gcrs_data.without_differentials()
        obsgeovel = gcrs_data.differentials['s'].to_cartesian()
        return obsgeopos, obsgeovel

    def gravitational_redshift(self, obstime,
                               bodies=['sun', 'jupiter', 'moon'],
                               masses={}):
        """Return the gravitational redshift at this EarthLocation.

        Calculates the gravitational redshift, of order 3 m/s, due to the
        requested solar system bodies.

        Parameters
        ----------
        obstime : `~astropy.time.Time`
            The ``obstime`` to calculate the redshift at.

        bodies : iterable, optional
            The bodies (other than the Earth) to include in the redshift
            calculation.  List elements should be any body name
            `get_body_barycentric` accepts.  Defaults to Jupiter, the Sun, and
            the Moon.  Earth is always included (because the class represents
            an *Earth* location).

        masses : dict of str to Quantity, optional
            The mass or gravitational parameters (G * mass) to assume for the
            bodies requested in ``bodies``. Can be used to override the
            defaults for the Sun, Jupiter, the Moon, and the Earth, or to
            pass in masses for other bodies.

        Returns
        --------
        redshift :  `~astropy.units.Quantity`
            Gravitational redshift in velocity units at given obstime.
        """
        # needs to be here to avoid circular imports
        from .solar_system import get_body_barycentric

        bodies = list(bodies)
        # Ensure earth is included and last in the list.
        if 'earth' in bodies:
            bodies.remove('earth')
        bodies.append('earth')
        _masses = {'sun': consts.GM_sun,
                   'jupiter': consts.GM_jup,
                   'moon': consts.G * 7.34767309e22*u.kg,
                   'earth': consts.GM_earth}
        _masses.update(masses)
        GMs = []
        M_GM_equivalency = (u.kg, u.Unit(consts.G * u.kg))
        for body in bodies:
            try:
                GMs.append(_masses[body].to(u.m**3/u.s**2, [M_GM_equivalency]))
            except KeyError:
                raise KeyError(f'body "{body}" does not have a mass!')
            except u.UnitsError as exc:
                exc.args += ('"masses" argument values must be masses or '
                             'gravitational parameters',)
                raise

        positions = [get_body_barycentric(name, obstime) for name in bodies]
        # Calculate distances to objects other than earth.
        distances = [(pos - positions[-1]).norm() for pos in positions[:-1]]
        # Append distance from Earth's center for Earth's contribution.
        distances.append(CartesianRepresentation(self.geocentric).norm())
        # Get redshifts due to all objects.
        redshifts = [-GM / consts.c / distance for (GM, distance) in
                     zip(GMs, distances)]
        # Reverse order of summing, to go from small to big, and to get
        # "earth" first, which gives m/s as unit.
        return sum(redshifts[::-1])

    @property
    def x(self):
        """The X component of the geocentric coordinates."""
        return self['x']

    @property
    def y(self):
        """The Y component of the geocentric coordinates."""
        return self['y']

    @property
    def z(self):
        """The Z component of the geocentric coordinates."""
        return self['z']

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if result.dtype is self.dtype:
            return result.view(self.__class__)
        else:
            return result.view(u.Quantity)

    def __array_finalize__(self, obj):
        super().__array_finalize__(obj)
        if hasattr(obj, '_ellipsoid'):
            self._ellipsoid = obj._ellipsoid

    def __len__(self):
        if self.shape == ():
            raise IndexError('0-d EarthLocation arrays cannot be indexed')
        else:
            return super().__len__()

    def _to_value(self, unit, equivalencies=[]):
        """Helper method for to and to_value."""
        # Conversion to another unit in both ``to`` and ``to_value`` goes
        # via this routine. To make the regular quantity routines work, we
        # temporarily turn the structured array into a regular one.
        array_view = self.view(self._array_dtype, np.ndarray)
        if equivalencies == []:
            equivalencies = self._equivalencies
        new_array = self.unit.to(unit, array_view, equivalencies=equivalencies)
        return new_array.view(self.dtype).reshape(self.shape)
