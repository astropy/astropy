.. include:: references.txt

.. _astropy-coordinates-remote:

Methods that access remote resources
------------------------------------

There are currently two methods that rely on getting remote data to work.

The first is the :class:`~astropy.coordinates.SkyCoord` :meth:`~astropy.coordinates.SkyCoord.from_name` method, which uses
`Sesame <http://cds.u-strasbg.fr/cgi-bin/Sesame>`_ to retrieve coordinates
for a particular named object::

    >>> SkyCoord.from_name("M42")  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SkyCoord (ICRS): (ra, dec) in deg
        (83.82208, -5.39111)>

The second is the :class:`~astropy.coordinates.EarthLocation` :meth:`~astropy.coordinates.EarthLocation.of_site` method, which
provides a similar quick way to get an
:class:`~astropy.coordinates.EarthLocation` from an observatory name::

    >>> from astropy.coordinates import EarthLocation
    >>> EarthLocation.of_site('Apache Point Observatory')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <EarthLocation (-1463969.3018517173, -5166673.342234327, 3434985.7120456537) m>
    
The full list of available observatory names can be obtained with
the :class:`~astropy.coordinates.EarthLocation` :meth:`~astropy.coordinates.EarthLocation.get_site_names` method.

While these methods are convenient, there are several considerations to take into account:

* Since these methods access online data, the data may evolve over time (for
  example, the accuracy of coordinates might improve, and new observatories
  may be added). However, this means that a script using these and running
  now may give a different answer in five years. Therefore, users concerned
  with reproducibility should not use these methods in their final scripts,
  but can instead use them to get the values required and then hard-code them
  into the scripts. For example, we can check the coordinates of the Kitt
  Peak Observatories using::

    >>> loc = EarthLocation.of_site('Kitt Peak')

  We can then view the actual cartesian coordinates for the observatory:
  
    >>> loc
    <EarthLocation (-1994502.6043061386, -5037538.54232911, 3358104.9969029757) m>
    
  This can then easily be converted to code::
  
    >>> loc = EarthLocation(-1994502.6043061386, -5037538.54232911, 3358104.9969029757, unit='m')
  
  This latter line can then be included in a script and will ensure that the
  results stay the same over time.

* The online data may not be accurate enough for your purposes. If maximum
  accuracy is paramount, we recommend that you determine the celestial or
  Earth coordinates yourself and hard-code these, rather than use the
  convenience methods.
  
* These methods will not function if an internet connection is not available.
  Therefore, if you need to work on a script while offline, follow the
  instructions in the first bullet point above to hard-code the coordinates
  before going offline.
