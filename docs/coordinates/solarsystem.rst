.. _astropy-coordinates-solarsystem:

Solar System Ephemerides
************************

`astropy.coordinates` can calculate the |SkyCoord| of some of the major solar
system objects. By default, it uses approximate orbital elements calculated
using |PyERFA| routines, but it can
also use more precise ones using the JPL ephemerides (which are derived from
dynamical models). The default JPL ephemerides (DE430) provide predictions
valid roughly for the years between 1550 and 2650. The file is 115 MB and will
need to be downloaded the first time you use this functionality, but will be
cached after that. Other JPL ephemerides can be requested by name for specific
use cases you may have (see the examples below).

.. note::
   Using JPL ephemerides requires that the `jplephem
   <https://pypi.org/project/jplephem/>`_ package be installed. This is
   most conveniently achieved via ``python -m pip install jplephem``, although
   whatever package management system you use might have it as well.

Two functions are provided; :func:`~astropy.coordinates.get_body` and
:func:`~astropy.coordinates.get_body_barycentric`.
The first returns |SkyCoord| objects in the `~astropy.coordinates.GCRS` frame,
while the latter returns a `~astropy.coordinates.CartesianRepresentation` of the
barycentric position of a body (i.e., in the `~astropy.coordinates.ICRS` frame).

Examples
========

..
  EXAMPLE START
  Using the Solar System Ephemerides

Here is an example of using these functions with built-in ephemerides (i.e.,
without the need to download a large ephemerides file)::

  >>> from astropy.time import Time
  >>> from astropy.coordinates import solar_system_ephemeris, EarthLocation
  >>> from astropy.coordinates import get_body_barycentric, get_body
  >>> t = Time("2014-09-22 23:22")
  >>> loc = EarthLocation.of_site('greenwich') # doctest: +REMOTE_DATA
  >>> with solar_system_ephemeris.set('builtin'):
  ...     jup = get_body('jupiter', t, loc) # doctest: +REMOTE_DATA +IGNORE_OUTPUT
  >>> jup  # doctest: +FLOAT_CMP +REMOTE_DATA
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.69182405, -550931.91022387, 4961151.73597633) m, obsgeovel=(40.159527, 287.47873161, -0.04597922) m / s): (ra, dec, distance) in (deg, deg, AU)
      (136.91116253, 17.02935396, 5.94386022)>

Above, we used ``solar_system_ephemeris`` as a context, which sets the default
ephemeris while in the ``with`` clause, and resets it at the end.

To get more precise positions than is possible with the built-in ephemeris
(see :ref:`astropy-coordinates-solarsystem-erfa-precision`), you
could use the ``de430`` ephemeris mentioned above, or, if you only care about
times between 1950 and 2050, opt for the ``de432s`` ephemeris, which is stored
in a smaller, ~10 MB, file (which will be downloaded and cached when the
ephemeris is set):

.. doctest-requires:: jplephem

  >>> solar_system_ephemeris.set('de432s') # doctest: +REMOTE_DATA, +IGNORE_OUTPUT
  <ScienceState solar_system_ephemeris: 'de432s'>
  >>> get_body('jupiter', t, loc) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.69182405, -550931.91022387, 4961151.73597633) m, obsgeovel=(40.159527, 287.47873161, -0.04597922) m / s): (ra, dec, distance) in (deg, deg, km)
      (136.90234846, 17.03160654, 8.89196021e+08)>
  >>> get_body('moon', t, loc) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <SkyCoord (GCRS: obstime=2014-09-22 23:22:00.000, obsgeoloc=(3949481.69182405, -550931.91022387, 4961151.73597633) m, obsgeovel=(40.159527, 287.47873161, -0.04597922) m / s): (ra, dec, distance) in (deg, deg, km)
      (165.51854528, 2.32861794, 407229.55638763)>
  >>> get_body_barycentric('moon', t) # doctest: +REMOTE_DATA, +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in km
      (1.50107535e+08, -866789.11996916, -418963.55218495)>

For one-off calculations with a given ephemeris, you can also pass it directly
to the various functions:

.. doctest-requires:: jplephem

  >>> get_body_barycentric('moon', t, ephemeris='de432s')
  ... # doctest: +REMOTE_DATA, +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in km
      (1.50107535e+08, -866789.11996916, -418963.55218495)>
  >>> get_body_barycentric('moon', t, ephemeris='builtin')
  ... # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in AU
      (1.00340683, -0.00579417, -0.00280064)>

..
  EXAMPLE END

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

.. note::
    While the sun is included in these ephemerides, it is important to
    recognize that `~astropy.coordinates.get_sun` always uses the built-in,
    polynomial model (as this requires no special download). So it is not safe
    to assume that ``get_body(time, 'sun')`` and ``get_sun(time)`` will give
    the same result.

.. note::
    When using JPL ephemerides, be aware that answers may change at levels that
    can be surprising if you are not careful about understanding the frame you
    are in.  See for example the case of the DE440s ephemerides, which is
    described in more detail in
    `astropy PR #11608 <https://github.com/astropy/astropy/pull/11608>`_. So
    it is usually best to stay within the same ephemerides for consistency.

.. _astropy-coordinates-solarsystem-erfa-precision:

Precision of the Built-In Ephemeris
===================================

The algorithm for calculating positions and velocities for planets other than
Earth used by |ERFA| is due to J.L. Simon, P. Bretagnon, J. Chapront,
M. Chapront-Touze, G. Francou and J. Laskar (Bureau des Longitudes, Paris,
France).  From comparisons with JPL ephemeris DE102, they quote the maximum
errors over the interval 1800-2050 below. For more details, see the |PyERFA| routine, `erfa.plan94`.
For the Earth, the rms errors in position and velocity are about 4.6 km and
1.4 mm/s, respectively (see `erfa.epv00`).

.. list-table::

  * - Planet
    - L (arcsec)
    - B (arcsec)
    - R (km)
  * - Mercury
    - 4
    - 1
    - 300
  * - Venus
    - 5
    - 1
    - 800
  * - EMB
    - 6
    - 1
    - 1000
  * - Mars
    - 17
    - 1
    - 7700
  * - Jupiter
    - 71
    - 5
    - 76000
  * - Saturn
    - 81
    - 13
    - 267000
  * - Uranus
    - 86
    - 7
    - 712000
  * - Neptune
    - 11
    - 1
    - 253000
