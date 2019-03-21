.. include:: references.txt

.. _astropy-coordinates-remote:

Usage Tips/Suggestions for Methods That Access Remote Resources
***************************************************************

There are currently two methods that rely on getting remote data to work.

The first is the :class:`~astropy.coordinates.SkyCoord` :meth:`~astropy.coordinates.SkyCoord.from_name` method, which uses
`Sesame <http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates
for a particular named object::

    >>> from astropy.coordinates import SkyCoord
    >>> SkyCoord.from_name("PSR J1012+5307")  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        ( 153.1393271,  53.117343)>

The second is the :class:`~astropy.coordinates.EarthLocation` :meth:`~astropy.coordinates.EarthLocation.of_site` method, which
provides a similar quick way to get an
:class:`~astropy.coordinates.EarthLocation` from an observatory name::

    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('Apache Point Observatory')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <EarthLocation (-1463969.3018517173, -5166673.342234327, 3434985.7120456537) m>

The full list of available observatory names can be obtained with
 :meth:`astropy.coordinates.EarthLocation.get_site_names`.

While these methods are convenient, there are several considerations to take
into account:

* Since these methods access online data, the data may evolve over time (for
  example, the accuracy of coordinates might improve, and new observatories
  may be added). Therefore, this means that a script using these and currently
  running may give a different answer in five years. Therefore, users concerned
  with reproducibility should not use these methods in their final scripts,
  but can instead use them to get the values required, and then hard-code them
  into the scripts. For example, we can check the coordinates of the Kitt
  Peak Observatories using::

    >>> loc = EarthLocation.of_site('Kitt Peak')  # doctest: +REMOTE_DATA

  Note that this command requires an internet connection.

  We can then view the actual Cartesian coordinates for the observatory:

    >>> loc  # doctest: +REMOTE_DATA +FLOAT_CMP
    <EarthLocation (-1994502.6043061386, -5037538.54232911, 3358104.9969029757) m>

  This can then be converted into code::

    >>> loc = EarthLocation(-1994502.6043061386, -5037538.54232911, 3358104.9969029757, unit='m')

  This latter line can then be included in a script and will ensure that the
  results stay the same over time.

* The online data may not be accurate enough for your purposes. If maximum
  accuracy is paramount, we recommend that you determine the celestial or
  Earth coordinates yourself and hard-code these, rather than using the
  convenience methods.

* These methods will not function if an internet connection is not available.
  Therefore, if you need to work on a script while offline, follow the
  instructions in the first bullet point above to hard-code the coordinates
  before going offline.
