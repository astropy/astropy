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
    >>> from astropy.coordinates.angle_utilities import golden_spiral_grid
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
   >>> location = location = EarthLocation(lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m)
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

    # 100_000 times randomly distributed over 12 hours
    t = Time('2020-01-01T20:00:00') + rng.uniform(0, 1, 10_000) * u.hour

    location = location = EarthLocation(
        lon=-17.89 * u.deg, lat=28.76 * u.deg, height=2200 * u.m
    )

    # A celestial object in ICRS
    crab = SkyCoord.from_name("Crab Nebula")

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

    fig = plt.figure()

    ax1, ax2 = fig.subplots(2, 1, gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

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

    ax1.set_title('Transformation of SkyCoord with 100.000 obstimes over 12 hours')

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
