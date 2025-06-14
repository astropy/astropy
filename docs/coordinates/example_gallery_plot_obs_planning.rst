.. _sphx_glr_generated_examples_coordinates_plot_obs-planning.py:

Determining and plotting the altitude/azimuth of a celestial object
===================================================================

..
  EXAMPLE START
  Determining and plotting the altitude/azimuth of a celestial object

This example demonstrates coordinate transformations and the creation of
visibility curves to assist with observing run planning.

In this example, we make a `~astropy.coordinates.SkyCoord` instance for M33.
The altitude-azimuth coordinates are then found using
`astropy.coordinates.EarthLocation` and `astropy.time.Time` objects.

This example is meant to demonstrate the capabilities of the
`astropy.coordinates` package. For more convenient and/or complex observation
planning, consider the `astroplan <https://astroplan.readthedocs.io/>`_
package.

Let's suppose you are planning to visit picturesque Bear Mountain State Park
in New York, USA. You are bringing your telescope with you (of course), and
someone told you M33 is a great target to observe there. You happen to know
you are free at 11:00 PM local time, and you want to know if it will be up.
Astropy can answer that.

.. plot::
   :include-source:

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    >>> from astropy.time import Time
    >>> from astropy.visualization import quantity_support

    :meth:`astropy.coordinates.SkyCoord.from_name` uses Simbad to resolve object
    names and retrieve coordinates.

    Get the coordinates of M33:

    >>> # m33 = SkyCoord.from_name("M33")
    >>> m33 = SkyCoord(23.46206906, 30.66017511, unit="deg")

    Use `astropy.coordinates.EarthLocation` to provide the location of Bear
    Mountain and set the time to 11pm Eastern Daylight Time (EDT) on 2012 July 12:

    >>> bear_mountain = EarthLocation(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)
    >>> utcoffset = -4 * u.hour  # EDT
    >>> time = Time("2012-7-12 23:00:00") - utcoffset

    :meth:`astropy.coordinates.EarthLocation.get_site_names` can be used to get
    locations of major observatories.

    Use `astropy.coordinates` to find the Alt, Az coordinates of M33 at as
    observed from Bear Mountain at 11pm on 2012 July 12:

    >>> m33altaz = m33.transform_to(AltAz(obstime=time, location=bear_mountain))
    >>> print(f"M33's Altitude = {m33altaz.alt:.2}")
    M33's Altitude = 0.13 deg

    This is helpful since it turns out M33 is barely above the horizon at this
    time. It is more informative to find M33's airmass over the course of
    the night.

    Find the Alt, Az coordinates of M33 at 100 times evenly spaced between 10 PM
    and 7 AM EDT:

    >>> midnight = Time("2012-7-13 00:00:00") - utcoffset
    >>> delta_midnight = np.linspace(-2, 10, 100) * u.hour
    >>> frame_July13night = AltAz(obstime=midnight + delta_midnight, location=bear_mountain)
    >>> m33altazs_July13night = m33.transform_to(frame_July13night)

    Convert Alt, Az to airmass with `~astropy.coordinates.AltAz.secz` attribute:

    >>> m33airmasss_July13night = m33altazs_July13night.secz

    Plot the airmass as a function of time:

    >>> with quantity_support():
    ...     fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ...     ax.plot(delta_midnight, m33airmasss_July13night)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_xlim(-2, 10)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_ylim(1, 4)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_xlabel("Hours from EDT Midnight")  # doctest: +IGNORE_OUTPUT
    ...     ax.set_ylabel("Airmass [Sec(z)]")  # doctest: +IGNORE_OUTPUT
    ...     plt.draw()

    Use  :func:`~astropy.coordinates.get_sun` to find the location of the Sun at 1000
    evenly spaced times between noon on July 12 and noon on July 13:

    >>> delta_midnight = np.linspace(-12, 12, 1000) * u.hour
    >>> times_July12_to_13 = midnight + delta_midnight
    >>> frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bear_mountain)
    >>> sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)

    Do the same with :func:`~astropy.coordinates.get_body` to find when the moon is
    up. Be aware that this will need to download a 10 MB file from the internet
    to get a precise location of the moon.

    >>> moon_July12_to_13 = get_body("moon", times_July12_to_13)
    >>> moonaltazs_July12_to_13 = moon_July12_to_13.transform_to(frame_July12_to_13)

    Find the Alt, Az coordinates of M33 at those same times:

    >>> m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)

    Make a figure illustrating nighttime and the altitudes of M33 and
    the Sun over that time:

    >>> with quantity_support():
    ...     fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ...     ax.plot(delta_midnight, sunaltazs_July12_to_13.alt, color="r", label="Sun")  # doctest: +IGNORE_OUTPUT
    ...     ax.plot(
    ...         delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75] * 3, ls="--", label="Moon"
    ...     )  # doctest: +IGNORE_OUTPUT
    ...     mappable = ax.scatter(
    ...         delta_midnight,
    ...         m33altazs_July12_to_13.alt,
    ...         c=m33altazs_July12_to_13.az.value,
    ...         label="M33",
    ...         lw=0,
    ...         s=8,
    ...         cmap="viridis",
    ...     )
    ...     ax.fill_between(
    ...         delta_midnight,
    ...         0 * u.deg,
    ...         90 * u.deg,
    ...         sunaltazs_July12_to_13.alt < (-0 * u.deg),
    ...         color="0.5",
    ...         zorder=0,
    ...     )  # doctest: +IGNORE_OUTPUT
    ...     ax.fill_between(
    ...         delta_midnight,
    ...         0 * u.deg,
    ...         90 * u.deg,
    ...         sunaltazs_July12_to_13.alt < (-18 * u.deg),
    ...         color="k",
    ...         zorder=0,
    ...     )  # doctest: +IGNORE_OUTPUT
    ...     fig.colorbar(mappable).set_label("Azimuth [deg]")  # doctest: +IGNORE_OUTPUT
    ...     ax.legend(loc="upper left")  # doctest: +IGNORE_OUTPUT
    ...     ax.set_xlim(-12 * u.hour, 12 * u.hour)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_xticks((np.arange(13) * 2 - 12) * u.hour)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_ylim(0 * u.deg, 90 * u.deg)  # doctest: +IGNORE_OUTPUT
    ...     ax.set_xlabel("Hours from EDT Midnight")  # doctest: +IGNORE_OUTPUT
    ...     ax.set_ylabel("Altitude [deg]")  # doctest: +IGNORE_OUTPUT
    ...     ax.grid(visible=True)  # doctest: +IGNORE_OUTPUT
    ...     plt.draw()

..
  EXAMPLE END
