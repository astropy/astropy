.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-coordinates-performance:

Performance Tips
================

If you are using |SkyCoord| for many different coordinates, you will see much
better performance if you create a single |SkyCoord| with arrays of coordinates
as opposed to creating individual |SkyCoord| objects for each individual
coordinate::

    >>> coord = SkyCoord(ra_array, dec_array, unit='deg')  # doctest: +SKIP

Frame attributes can be arrays too, as long as the coordinate data and all of
the frame attributes have shapes that are compatible according to
:doc:`Numpy broadcasting rules <numpy:user/basics.broadcasting>`:


.. testsetup::

    >>> from astropy.coordinates import FK4
    >>> from astropy import units as u

::

    >>> coord = FK4(1 * u.deg, 2 * u.deg, obstime=["J2000", "J2001"])
    >>> coord.shape
    (2,)

In addition, looping over a |SkyCoord| object can be slow. If you need to
transform the coordinates to a different frame, it is much faster to transform a
single |SkyCoord| with arrays of values as opposed to looping over the
|SkyCoord| and transforming them individually.

Finally, for more advanced users, note that you can use broadcasting to
transform |SkyCoord| objects into frames with vector properties.

..
  EXAMPLE START
  Performance Tips for Transforming SkyCoord Objects

To use broadcasting to transform |SkyCoord| objects into frames with vector
properties::

    >>> from astropy.coordinates import SkyCoord, EarthLocation
    >>> from astropy import coordinates as coord
    >>> from astropy.coordinates import golden_spiral_grid
    >>> from astropy.time import Time
    >>> from astropy import units as u
    >>> import numpy as np

    >>> # 1000 locations in a grid on the sky
    >>> coos = SkyCoord(golden_spiral_grid(size=1000))

    >>> # 300 times over the space of 10 hours
    >>> times = Time.now() + np.linspace(-5, 5, 300)*u.hour

    >>> # note the use of broadcasting so that 300 times are broadcast against 1000 positions
    >>> lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
    >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=lapalma)

    >>> # calculate alt-az of each object at each time.
    >>> aa_coos = coos.transform_to(aa_frame)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS

..
  EXAMPLE END


Broadcasting Over Frame Data and Attributes
-------------------------------------------

..
  EXAMPLE START
  Broadcasting Over Frame Data and Attributes

Frames in `astropy.coordinates` support
:doc:`Numpy broadcasting rules <numpy:user/basics.broadcasting>` over both
frame data and frame attributes. This makes it easy and fast to do positional
astronomy calculations and transformations on sweeps of parameters.

Where this really shines is doing fast observability calculations over arrays.
The following example constructs an `~astropy.coordinates.EarthLocation` array
of length :samp:`{L}`, a `~astropy.coordinates.SkyCoord` array of length
:samp:`{M}`, and a `~astropy.time.Time` array of length :samp:`N`. It uses
Numpy broadcasting rules to evaluate a boolean array of shape
:samp:`({L}, {M}, {N})` that is `True` for those observing locations, times,
and sky coordinates, for which the target is above an altitude limit::

    >>> from astropy.coordinates import EarthLocation, AltAz, SkyCoord
    >>> from astropy.coordinates.angles import uniform_spherical_random_surface
    >>> from astropy.time import Time
    >>> from astropy import units as u
    >>> import numpy as np

    >>> L = 25
    >>> M = 100
    >>> N = 50

    >>> # Earth locations of length L
    >>> c = uniform_spherical_random_surface(L)
    >>> locations = EarthLocation.from_geodetic(c.lon, c.lat)

    >>> # Celestial coordinates of length M
    >>> coords = SkyCoord(uniform_spherical_random_surface(M))

    >>> # Observation times of length N
    >>> obstimes = Time('2023-08-04') + np.linspace(0, 24, N) * u.hour

    >>> # AltAz coordinates of shape (L, M, N)
    >>> frame = AltAz(
    ...     location=locations[:, np.newaxis, np.newaxis],
    ...     obstime=obstimes[np.newaxis, np.newaxis, :])
    >>> altaz = coords[np.newaxis, :, np.newaxis].transform_to(frame)  # doctest: +REMOTE_DATA

    >>> min_altitude = 30 * u.deg
    >>> is_above_altitude_limit = (altaz.alt > min_altitude)  # doctest: +REMOTE_DATA
    >>> is_above_altitude_limit.shape  # doctest: +REMOTE_DATA
    (25, 100, 50)

..
  EXAMPLE END


Improving Performance for Arrays of ``obstime``
-----------------------------------------------

The most expensive operations when transforming between observer-dependent coordinate
frames (e.g. ``AltAz``) and sky-fixed frames (e.g. ``ICRS``) are the calculation
of the orientation and position of Earth.

If |SkyCoord| instances are transformed for a large  number of closely spaced ``obstime``,
these calculations can be sped up by factors up to 100, whilst still keeping micro-arcsecond precision,
by utilizing interpolation instead of calculating Earth orientation parameters for each individual point.

..
  EXAMPLE START
  Improving performance for obstime arrays

To use interpolation for the astrometric values in coordinate transformation, use::

   >>> from astropy.coordinates import SkyCoord, EarthLocation, AltAz
   >>> from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
   >>> from astropy.time import Time
   >>> from time import perf_counter
   >>> import numpy as np
   >>> import astropy.units as u


   >>> # array with 10000 obstimes
   >>> obstime = Time('2010-01-01T20:00') + np.linspace(0, 6, 10000) * u.hour
   >>> location = EarthLocation(lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m)
   >>> frame = AltAz(obstime=obstime, location=location)
   >>> crab = SkyCoord(ra='05h34m31.94s', dec='22d00m52.2s')

   >>> # transform with default transformation and print duration
   >>> t0 = perf_counter()
   >>> crab_altaz = crab.transform_to(frame)  # doctest:+IGNORE_WARNINGS +REMOTE_DATA
   >>> print(f'Transformation took {perf_counter() - t0:.2f} s')  # doctest:+IGNORE_OUTPUT
   Transformation took 1.77 s

   >>> # transform with interpolating astrometric values
   >>> t0 = perf_counter()
   >>> with erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)): # doctest:+REMOTE_DATA
   ...     crab_altaz_interpolated = crab.transform_to(frame)  # doctest:+IGNORE_WARNINGS +REMOTE_DATA
   >>> print(f'Transformation took {perf_counter() - t0:.2f} s')  # doctest:+IGNORE_OUTPUT
   Transformation took 0.03 s

   >>> err = crab_altaz.separation(crab_altaz_interpolated)  # doctest:+IGNORE_WARNINGS +REMOTE_DATA
   >>> print(f'Mean error of interpolation: {err.to(u.microarcsecond).mean():.4f}')  # doctest:+ELLIPSIS +REMOTE_DATA
   Mean error of interpolation: 0.0... uarcsec

   >>> # To set erfa_astrom for a whole session, use it without context manager:
   >>> erfa_astrom.set(ErfaAstromInterpolator(300 * u.s))  # doctest:+SKIP

..
  EXAMPLE END


Here, we look into choosing an appropriate ``time_resolution``.
We will transform a single sky coordinate for lots of observation times from
``ICRS`` to ``AltAz`` and evaluate precision and runtime for different values
for ``time_resolution`` compared to the non-interpolating, default approach.

.. plot::
   :include-source:
   :context: reset

    from time import perf_counter

    import numpy as np
    import matplotlib.pyplot as plt

    from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import Time
    import astropy.units as u

    rng = np.random.default_rng(1337)

    n_coords = 10_000
    time_delta = 1 * u.hour
    # n_coords times randomly distributed over time_delta
    t = Time('2020-01-01T20:00:00') + rng.random(n_coords) * time_delta

    location = EarthLocation(
        lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m
    )

    # A celestial object in ICRS
    # crab = SkyCoord.from_name("Crab Nebula")
    crab = SkyCoord(83.6287, 22.0147, unit="deg")

    # target horizontal coordinate frame
    altaz = AltAz(obstime=t, location=location)

    # the reference transform using no interpolation
    t0 = perf_counter()
    no_interp = crab.transform_to(altaz)
    reference = perf_counter() - t0
    print(f'No Interpolation took {reference:.4f} s')

    # now the interpolating approach for different time resolutions
    resolutions = 10.0**np.arange(-1, 5) * u.s
    times = []
    seps = []

    for resolution in resolutions:
        with erfa_astrom.set(ErfaAstromInterpolator(resolution)):
            t0 = perf_counter()
            interp = crab.transform_to(altaz)
            duration = perf_counter() - t0

        print(
            f'Interpolation with {resolution.value: 9.1f} {str(resolution.unit)}'
            f' resolution took {duration:.4f} s'
            f' ({reference / duration:5.1f}x faster) '
        )
        seps.append(no_interp.separation(interp))
        times.append(duration)

    seps = u.Quantity(seps)

    fig, (ax1, ax2) = plt.subplots(
      nrows=2,
      sharex=True,
      gridspec_kw=dict(height_ratios=[2, 1]),
    )

    ax1.plot(
        resolutions.to_value(u.s),
        seps.mean(axis=1).to_value(u.microarcsecond),
        'o', label='mean',
    )

    for p in [25, 50, 75, 95]:
        ax1.plot(
            resolutions.to_value(u.s),
            np.percentile(seps.to_value(u.microarcsecond), p, axis=1),
            'o', label=f'{p}%', color='C1', alpha=p / 100,
        )

    ax1.set_title(f'Transformation of SkyCoord with {n_coords:,} obstimes over {time_delta}')

    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Angular distance to no interpolation / µas')

    ax2.plot(resolutions.to_value(u.s), reference / np.array(times), 's')
    ax2.set_yscale('log')
    ax2.set_ylabel('Speedup')
    ax2.set_xlabel('time resolution / s')

    ax2.yaxis.grid()
    fig.tight_layout()
