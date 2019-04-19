.. include:: references.txt

.. _astropy-coordinates-apply-space-motion:

Accounting for Space Motion
***************************

The |SkyCoord| object supports updating the position of a source given its space
motion and a time at which to evaluate the new position (or a difference between the coordinate's current time and a new one). This is
done using the :meth:`~astropy.coordinates.SkyCoord.apply_space_motion` method.
As an example, first we'll create a |SkyCoord| object with a specified
``obstime``::

    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from astropy.coordinates import SkyCoord
    >>> c = SkyCoord(l=10*u.degree, b=45*u.degree, distance=100*u.pc,
    ...              pm_l_cosb=34*u.mas/u.yr, pm_b=-117*u.mas/u.yr,
    ...              frame='galactic',
    ...              obstime=Time('1988-12-18 05:11:23.5'))

We can now find the position at some other time, taking the space motion into
account. We can either specify the time difference between the observation time
and the desired time::

    >>> c.apply_space_motion(dt=10. * u.year) # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        ( 10.00013356,  44.999675,  99.99999994)
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        ( 33.99980714, -117.00005604,  0.00034117)>
    >>> c.apply_space_motion(dt=-10. * u.year) # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        ( 9.99986643,  45.000325,  100.00000006)
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        ( 34.00019286, -116.99994395, -0.00034117)>

Or, we can specify the new time to evaluate the position at::

    >>> c.apply_space_motion(new_obstime=Time('2017-12-18 01:12:07.3')) # doctest: +FLOAT_CMP
    <SkyCoord (Galactic): (l, b, distance) in (deg, deg, pc)
        ( 10.00038732,  44.99905754,  99.99999985)
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        ( 33.99944073, -117.00016248,  0.00098937)>

If the |SkyCoord| object has no specified radial velocity (RV), the RV is
assumed to be 0. The new position of the source is determined assuming the
source moves in a straight line with constant velocity in an inertial frame.
There are no plans to support more complex evolution (e.g., non-inertial
frames or more complex evolution), as that is out of scope for the ``astropy``
core package (although it may well be in-scope for a variety of affiliated
packages).

Example: Use velocity to compute sky position at different epochs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we'll use *Gaia* `TGAS
<https://www.cosmos.esa.int/web/gaia/dr1>`_ astrometry for a nearby star to
compute the sky position of the source on the date that the 2MASS survey
observed that region of the sky. The TGAS astrometry is provided on the
reference epoch J2015.0, whereas the 2MASS survey occurred in the late 1990's.
For the star of interest, the proper motion is large enough that there are
appreciable differences in the sky position between the two surveys.

After computing the previous position of the source, we'll then cross-match the
source with the 2MASS catalog to compute *Gaia*-2MASS colors for this object
source.

.. note::

    This example requires accessing data from the *Gaia* TGAS and 2MASS
    catalogs. For convenience and speed below, we have created dictionary
    objects that contain the data. We retrieved the data using the Astropy
    affiliated package `astroquery <https://astroquery.readthedocs.io/>`_ using
    the following queries::

        import astropy.coordinates as coord
        import astropy.units as u
        from astroquery.gaia import Gaia
        from astroquery.vizier import Vizier

        job = Gaia.launch_job("SELECT TOP 1 * FROM gaiadr1.tgas_source \
            WHERE parallax_error < 0.3  AND parallax > 5 AND pmra > 100 \
            ORDER BY random_index")
        result_tgas = job.get_results()[0]

        c_tgas = coord.SkyCoord(ra=result_tgas['ra'] * u.deg,
                                dec=result_tgas['dec'] * u.deg)
        v = Vizier(columns=["**"], catalog="II/246/out")
        result_2mass = v.query_region(c, radius=1*u.arcmin)['II/246/out']

The TGAS data from relevant columns for this source (see queries in Note
above)::

    >>> result_tgas = dict(ra=66.44280212823296,
    ...                    dec=-69.99366255906372,
    ...                    parallax=22.764078749733947,
    ...                    pmra=144.91354358297048,
    ...                    pmdec=5.445648092997134,
    ...                    ref_epoch=2015.0,
    ...                    phot_g_mean_mag=7.657174523348196)

The 2MASS data for all sources within 1 arcminute around the above position
(see queries in Note above)::

    >>> result_2mass = dict(RAJ2000=[66.421970000000002, 66.433521999999996,
    ...                              66.420564999999996, 66.485068999999996,
    ...                              66.467928999999998, 66.440815000000001,
    ...                              66.440454000000003],
    ...                     DEJ2000=[-70.003722999999994, -69.990768000000003,
    ...                              -69.992255999999998, -69.994881000000007,
    ...                              -69.994926000000007, -69.993613999999994,
    ...                              -69.990836999999999],
    ...                     Jmag=[16.35, 13.663, 16.171, 16.184, 16.292,
    ...                           6.6420002, 12.275],
    ...                     Hmag=[15.879, 13.955, 15.154, 15.856, 15.642,
    ...                           6.3660002, 12.185],
    ...                     Kmag=[15.581, 14.238, 14.622, 15.398, 15.123,
    ...                           6.2839999, 12.106],
    ...                     Date=['1998-10-24', '1998-10-24', '1998-10-24',
    ...                           '1998-10-24', '1998-10-24', '1998-10-24',
    ...                           '1998-10-24'])

We'll first create a |SkyCoord| object from the information provided in the TGAS
catalog. Note that we set the ``obstime`` of the object to the reference epoch
provided by the TGAS catalog (J2015.0)::

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord, Distance
    >>> from astropy.time import Time
    >>> c = SkyCoord(ra=result_tgas['ra'] * u.deg,
    ...              dec=result_tgas['dec'] * u.deg,
    ...              distance=Distance(parallax=result_tgas['parallax'] * u.mas),
    ...              pm_ra_cosdec=result_tgas['pmra'] * u.mas/u.yr,
    ...              pm_dec=result_tgas['pmdec'] * u.mas/u.yr,
    ...              obstime=Time(result_tgas['ref_epoch'], format='decimalyear'))

We next create a |SkyCoord| object with the sky positions from the 2MASS
catalog, and an `~astropy.time.Time` object for the date of the 2MASS
observations provided in the 2MASS catalog (for the data in this region the
observation date is the same, so we take only the 0th value)::

    >>> catalog_2mass = SkyCoord(ra=result_2mass['RAJ2000'] * u.deg,
    ...                          dec=result_2mass['DEJ2000'] * u.deg)
    >>> epoch_2mass = Time(result_2mass['Date'][0])

We can now use the :meth:`~astropy.coordinates.SkyCoord.apply_space_motion`
method to compute the position of the TGAS source at another epoch. This uses
the proper motion and parallax information to evolve the position of the source
assuming straight-line motion::

    >>> c_2mass_epoch = c.apply_space_motion(epoch_2mass)

Now that we have the coordinates of the TGAS source at the 2MASS epoch, we can
do the cross-match (see also :ref:`astropy-coordinates-separations-matching`)::

    >>> idx, sep, _ = c_2mass_epoch.match_to_catalog_sky(catalog_2mass) # doctest: +SKIP
    >>> sep[0].to_string() # doctest: +FLOAT_CMP +SKIP
    '0d00m00.2818s'
    >>> idx # doctest: +SKIP
    array(5)

The closest source it found is 0.2818 arcseconds away and corresponds to
row index 5 in the 2MASS catalog. We can then, for example, compute *Gaia*-2MASS
colors::

    >>> G = result_tgas['phot_g_mean_mag']
    >>> J = result_2mass['Jmag'][idx] # doctest: +SKIP
    >>> K = result_2mass['Kmag'][idx] # doctest: +SKIP
    >>> G - J, G - K # doctest: +SKIP
    (1.0151743233481962, 1.3731746233481958)
