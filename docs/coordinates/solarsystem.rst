.. include:: references.txt

.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
------------------------

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
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=( 3949481.689878457, -550931.9118838,  4961151.73733447) m, obsgeovel=( 40.1745933,  288.00078051,  0.) m / s): (ra, dec, distance) in (deg, deg, AU)
    ( 136.91116201,  17.02935408,  5.94386022)>

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
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=( 3949481.68990897, -550931.9118838,  4961151.73733447) m, obsgeovel=( 40.1745933,  288.00078051,  0.) m / s): (ra, dec, distance) in (deg, deg, km)
      ( 136.90234781,  17.03160686,   8.89196019e+08)>
  >>> get_moon(t, loc) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=( 3949481.6899252, -550931.91194065,  4961151.73733445) m, obsgeovel=( 40.1745933,  288.00078051,  0.) m / s): (ra, dec, distance) in (deg, deg, km)
          ( 165.51849193,  2.32863887,  407229.65033585)>
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
      (  1.50107516e+08, -866828.92708201, -418980.15909655)>

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
