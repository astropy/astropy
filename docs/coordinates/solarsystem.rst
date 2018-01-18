.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
************************

`astropy.coordinates` can calculate the |SkyCoord| of some of the major solar
system objects. By default, it uses approximate orbital elements calculated
using built-in `ERFA <https://github.com/liberfa/erfa>`_ routines, but it can
also use more precise ones using the JPL ephemerides (which are derived from
dynamical models).  The default JPL ephemerides (DE430) provide predictions
valid roughly for years between 1550 and 2650. The file is 115 MB and will need
to be downloaded the first time you use this functionality, but will be cached
after that.

.. note::
   Using JPL ephemerides requires that the `jplephem
   <https://pypi.python.org/pypi/jplephem>`_ package be installed. This is
   most easily achieved via ``pip install jplephem``, although whatever
   package management system you use might have it as well.

Three functions are provided; :meth:`~astropy.coordinates.get_body`,
:meth:`~astropy.coordinates.get_moon` and
:meth:`~astropy.coordinates.get_body_barycentric`. The first two functions
return |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame, whilst the
latter returns a `~astropy.coordinates.CartesianRepresentation` of the
barycentric position of a body (i.e in the `~astropy.coordinates.ICRS` frame).

Here is an example of using these functions with built-in ephemerides, i.e.,
without the need to download a large ephemerides file::

  >>> from astropy.time import Time
  >>> from astropy.coordinates import solar_system_ephemeris, EarthLocation
  >>> from astropy.coordinates import get_body_barycentric, get_body, get_moon
  >>> t = Time("2014-09-22 23:22")
  >>> loc = EarthLocation.of_site('greenwich') # doctest: +REMOTE_DATA
  >>> with solar_system_ephemeris.set('builtin'):
  ...     jup = get_body('jupiter', t, loc) # doctest: +REMOTE_DATA +IGNORE_OUTPUT
  >>> jup  # doctest: +FLOAT_CMP +REMOTE_DATA
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.68990863, -550931.91188162, 4961151.73733451) m, obsgeovel=(40.15954083, 287.47876693, -0.04597867) m / s): (ra, dec, distance) in (deg, deg, AU)
      (136.91116209, 17.02935409, 5.94386022)>

Above, we used ``solar_system_ephemeris`` as a context, which sets the default
ephemeris while in the ``with`` clause, and resets it at the end.

To get more precise positions, one could use the ``de430`` ephemeris mentioned
above, but between 1950 and 2050 one could also opt for the ``de432s``
ephemeris, which is stored in a smaller, ~10 MB, file (which will be
downloaded and cached when the ephemeris is set):

.. doctest-requires:: jplephem

  >>> solar_system_ephemeris.set('de432s') # doctest: +REMOTE_DATA, +IGNORE_OUTPUT
  <ScienceState solar_system_ephemeris: 'de432s'>
  >>> get_body('jupiter', t, loc) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.69230491, -550931.90674055, 4961151.73597586) m, obsgeovel=(40.15954083, 287.47863521, -0.0459789) m / s): (ra, dec, distance) in (deg, deg, km)
      (136.90234802, 17.03160667, 8.89196021e+08)>
  >>> get_moon(t, loc) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.69230491, -550931.90674055, 4961151.73597586) m, obsgeovel=(40.15954083, 287.47863521, -0.0459789) m / s): (ra, dec, distance) in (deg, deg, km)
      (165.51849203, 2.32863886, 407229.6503193)>
  >>> get_body_barycentric('moon', t) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in km
      (  1.50107535e+08, -866789.11996916, -418963.55218495)>

For one-off calculations with a given ephemeris, one can also pass it directly
to the various functions:

.. doctest-requires:: jplephem

  >>> get_body_barycentric('moon', t, ephemeris='de432s')
  ... # doctest: +REMOTE_DATA, +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in km
      (  1.50107535e+08, -866789.11996916, -418963.55218495)>
  >>> get_body_barycentric('moon', t, ephemeris='builtin')
  ... # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in km
      (  1.50107513e+08, -866838.51786769, -418988.57509287)>

For a list of the bodies for which positions can be calculated, do:

.. note that we skip the next test if jplephem is not installed because if
.. jplephem was not installed, we didn't change the science state higher up

.. doctest-requires:: jplephem

  >>> solar_system_ephemeris.bodies # doctest: +REMOTE_DATA
  ('sun',
   'mercury',
   'venus',
   'earth-moon-barycenter',
   'earth',
   'moon',
   'mars',
   'jupiter',
   'saturn',
   'uranus',
   'neptune',
   'pluto')
  >>> solar_system_ephemeris.set('builtin')
  <ScienceState solar_system_ephemeris: 'builtin'>
  >>> solar_system_ephemeris.bodies
  ('earth',
   'sun',
   'moon',
   'mercury',
   'venus',
   'earth-moon-barycenter',
   'mars',
   'jupiter',
   'saturn',
   'uranus',
   'neptune')

.. note ::
    While the sun is included in the these ephemerides, it is important to
    recognize that `~astropy.coordinates.get_sun` always uses the built-in,
    polynomial model (as this requires no special download). So it is not safe
    to assume that ``get_body(time, 'sun')`` and ``get_sun(time)`` will give
    the same result.
