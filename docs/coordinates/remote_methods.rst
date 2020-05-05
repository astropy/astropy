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

.. testsetup::

    >>> from astropy.coordinates import EarthLocation
    >>> loc = EarthLocation.of_site('Royal Observatory Greenwich', force_builtin=True)

The second is the :class:`~astropy.coordinates.EarthLocation` :meth:`~astropy.coordinates.EarthLocation.of_site` method, which
provides a similar quick way to get an
:class:`~astropy.coordinates.EarthLocation` from an observatory name::

    >>> from astropy.coordinates import EarthLocation
    >>> loc = EarthLocation.of_site('Royal Observatory Greenwich')  # doctest: +SKIP
    >>> loc  # doctest: +FLOAT_CMP
    <EarthLocation (3980608.90246817, -102.47522911, 4966861.27310068) m>

The full list of available observatory names can be obtained with
 :meth:`astropy.coordinates.EarthLocation.get_site_names`.

.. testsetup::

    >>> loc = EarthLocation.of_site('Keck', force_builtin=True)

While these methods are convenient, there are several considerations to take
into account:

* Since these methods access online data, the data may evolve over time (for
  example, the accuracy of coordinates might improve, and new observatories
  may be added). Therefore, this means that a script using these and currently
  running may give a different answer in five years. Therefore, users concerned
  with reproducibility should not use these methods in their final scripts,
  but can instead use them to get the values required, and then hard-code them
  into the scripts. For example, we can check the coordinates of the
  Keck Observatory using::

    >>> loc = EarthLocation.of_site('Keck')  # doctest: +SKIP

  Note that this command requires an internet connection.

  We can then view the actual Cartesian coordinates for the observatory:

    >>> loc  # doctest: +FLOAT_CMP
    <EarthLocation (-5464487.81759887, -2492806.59108569, 2151240.19451846) m>

  This can then be converted into code::

    >>> loc = EarthLocation(-5464487.81759887, -2492806.59108569, 2151240.19451846, unit='m')

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
