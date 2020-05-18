.. include:: references.txt

.. |AstrometricFrame| replace:: `~astropy.coordinates.AstrometricFrame`
.. |ICRS| replace:: `~astropy.coordinates.ICRS`
.. |HCRS| replace:: `~astropy.coordinates.HCRS`

.. _astropy-coordinates-astrometricframe:

Astrometric Versions of Coordinate Frames
*****************************************

|AstrometricFrame| creates an astrometric version of a base coordinate frame,
relative to a specified observer location.
The astrometric location of a body is the location of the body as measured
by an observer, which depends on the motion of that body:

* For a body in the solar system, its coordinate normally represents the
  instantaneous location of that body.  However, due to the finite speed of
  light, the light that reaches an observer at a certain observation time
  was emitted at an earlier time.  If the body is moving, then the observer
  will measure the body to be where it was at that earlier time.  This
  effect is also known as "planetary aberration".
* For a cosmic object, its coordinate normally represents its catalog
  location at a reference epoch (e.g., J2000.0).  If the body is moving
  (i.e., has non-zero "proper motion" and/or "radial velocity"), its
  measured location will evolve over time.  In this case, light travel time
  is already included in the catalog location.

.. warning::

    The astrometric calculation assumes that the body is moving in a
    straight line with its specified velocity.  The accuracy of this
    assumption depends on the body's true motion during the time it takes
    for light to travel the observer-body distance (for solar-system
    bodies) or the time since the reference epoch (for cosmic objects).

The astrometric calculation distinguishes between bodies in the solar
system and cosmic objects by the distance from the solar-system barycenter
(SSB).  The coordinates of bodies less than 1 light-year from the SSB are
treated as instantaneous locations.  The coordinates of bodies greater than
1 light-year from the SSB are treated as locations at the reference epoch.

As these are "astrometric" coordinates rather than "apparent" coordinates,
effects such as stellar aberration and gravitational deflection are not
included, at least not explicitly.  Such effects may be implicitly
included if they are accounted for in the base coordinate frame.

Solar-system bodies
===================

Example using the Venus transit
-------------------------------

During the Venus transit of the Sun as seen from Earth.
Instantaneous position/velocity of Venus.
Use a JPL ephemeris for accuracy::

    >>> from astropy.coordinates import (SkyCoord, solar_system_ephemeris,
    ...                                  get_body_barycentric_posvel,
    ...                                  CartesianDifferential)
    >>> import astropy.units as u
    >>> earth = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU,
    ...                  frame='gcrs', obstime='2012-06-06 04:07:29')
    >>> solar_system_ephemeris.set('de432s')  # doctest: +REMOTE_DATA
    <ScienceState solar_system_ephemeris: 'de432s'>
    >>> venus_pos, venus_vel = get_body_barycentric_posvel('venus',
    ...                                                    earth.obstime)
    >>> venus_vel = venus_vel.represent_as(CartesianDifferential)
    >>> venus = SkyCoord(venus_pos.with_differentials(venus_vel),
    ...                  frame='icrs')
    >>> venus = venus.transform_to(earth)
    >>> venus  # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2012-06-06 04:07:29.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, km)
        (74.22533319, 22.7724864, 43190032.91493926)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-7.87423767e+08, -4.03196311e+08, 0.09882208)>

Using |AstrometricFrame|::

    >>> from astropy.coordinates import AstrometricFrame
    >>> ast_gcrs = AstrometricFrame(origin=earth)
    >>> ast_gcrs
    <AstrometricGCRS Frame (origin=<GCRS Coordinate (obstime=2012-06-06 04:07:29.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, AU)
        (0., 0., 0.)>, ref_epoch=J2000.000)>

    >>> venus_ast = venus.transform_to(ast_gcrs)
    >>> venus_ast  # doctest: +REMOTE_DATA
    <SkyCoord (AstrometricGCRS: origin=<GCRS Coordinate (obstime=2012-06-06 04:07:29.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, AU)
        (0., 0., 0.)>, ref_epoch=J2000.000): (ra, dec, distance) in (deg, deg, km)
        (74.23245792, 22.77361214, 43190037.9814361)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-7.87503393e+08, -4.03424187e+08, 0.09816952)>

Transforming to a non-astrometric frame will "undo" the astrometric calculation.
Convert to a real coordinate using :meth:`~astropy.coordinates.AstrometricFrame.as_base`::

    >>> venus_ast_real = SkyCoord(venus_ast.as_base())
    >>> venus_ast_real  # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2012-06-06 04:07:29.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, km)
        (74.23245792, 22.77361214, 43190037.9814361)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-7.87503393e+08, -4.03424187e+08, 0.09816952)>

Compare against :func:`~astropy.coordinates.get_body`::

    >>> from astropy.coordinates import get_body
    >>> venus_gb = get_body('venus', earth.obstime)
    >>> venus_gb  # doctest: +REMOTE_DATA
    <SkyCoord (GCRS: obstime=2012-06-06 04:07:29.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (ra, dec, distance) in (deg, deg, km)
        (74.23245792, 22.77361214, 43190038.0981755)>

Discrepancies are small::

    >>> venus_gb.separation(venus_ast_real).to('mas')  # doctest: +REMOTE_DATA
    <Angle 0.00325284 mas>
    >>> venus_gb.separation_3d(venus_ast_real).to('km')  # doctest: +REMOTE_DATA
    <Distance 0.11674138 km>

Size of errors for planets
--------------------------

For an observer at Earth, these are the astrometric errors due to the linear-motion approximation:

======= ===================== =================== ===================
Planet  Light travel time (s) Angular error (mas) Position error (km)
======= ===================== =================== ===================
Mercury 300--700              <11                 <14
Venus   150--900              <2.0                <4.3
Mars    300--1200             <0.6                <2.2
Jupiter 2000--3200            <0.05               <1.1
Saturn  4200--5300            <0.02               <1.0
Uranus  9500--10500           <0.01               <0.9
Neptune 14000--15000          <0.01               <0.9
======= ===================== =================== ===================

Cosmic objects
==============

Consider a star with a high proper motion.
The reference epoch will be J2010.0, but we supply that later::

    >>> from astropy.coordinates import ICRS
    >>> star_ref = ICRS(ra=10*u.deg, dec=45*u.deg, distance=10*u.pc,
    ...                 pm_ra_cosdec=10*u.arcsec/u.yr, pm_dec=-20*u.arcsec/u.yr,
    ...                 radial_velocity=100*u.km/u.s)
    >>> star_ref
    <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
        (10., 45., 10.)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (10000., -20000., 100.)>

Create the |AstrometricFrame| version of |ICRS|.
For cosmic objects, the observer location does not matter.
Specify ``obstime`` separately because |ICRS| does not have that frame attribute.
Specify the reference epoch::

    >>> ast_icrs = AstrometricFrame(origin=ICRS(0*u.deg, 0*u.deg, 0*u.AU),
    ...                             obstime=earth.obstime, ref_epoch='J2010')
    >>> ast_icrs
    <AstrometricICRS Frame (origin=<ICRS Coordinate: (ra, dec, distance) in (deg, deg, AU)
        (0., 0., 0.)>, obstime=2012-06-06 04:07:29.000, ref_epoch=J2010.000)>

    >>> star_ref.transform_to(ast_icrs)
    <AstrometricICRS Coordinate (origin=<ICRS Coordinate: (ra, dec, distance) in (deg, deg, AU)
        (0., 0., 0.)>, obstime=2012-06-06 04:07:29.000, ref_epoch=J2010.000): (ra, dec, distance) in (deg, deg, pc)
        (10.00953932, 44.98650579, 10.00024406)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (9997.14835121, -20000.18211581, 100.27901197)>

|AstrometricFrame| calls :meth:`~astropy.coordinates.SkyCoord.apply_space_motion` internally, so using that method will give the same answer (see :ref:`astropy-coordinates-apply-space-motion` for usage)::

    >>> star_ref_sc = SkyCoord(star_ref, obstime=ast_icrs.ref_epoch)
    >>> star_ref_sc.apply_space_motion(new_obstime=ast_icrs.obstime)
    <SkyCoord (ICRS): (ra, dec, distance) in (deg, deg, pc)
        (10.00953932, 44.98650579, 10.00024406)
     (pm_ra, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (14134.77417664, -20000.18211581, 100.27901197)>

The astrometric frame can be based on any coordinate frame, not just |ICRS|.
Use the same Earth observer and observation time as earlier, but now with |HCRS| as the base coordinate frame::

    >>> ast_hcrs = AstrometricFrame(origin=earth.hcrs, ref_epoch='J2010')
    >>> ast_hcrs  # doctest: +REMOTE_DATA
    <AstrometricHCRS Frame (origin=<HCRS Coordinate (obstime=2012-06-06 04:07:29.000): (ra, dec, distance) in (deg, deg, AU)
        (254.4656737, -22.66944306, 1.01475669)>, ref_epoch=J2010.000)>

Need to change the finite-difference time step for the |HCRS| to |HCRS| transformation::

    >>> from astropy.coordinates import HCRS, frame_transform_graph
    >>> hcrs_to_hcrs = frame_transform_graph.get_transform(HCRS, HCRS).transforms[0]
    >>> hcrs_to_hcrs.finite_difference_dt = 1*u.yr

Perform the astrometric transformation::

    >>> star_ref.transform_to(ast_hcrs)  # doctest: +REMOTE_DATA
    <AstrometricHCRS Coordinate (origin=<HCRS Coordinate (obstime=2012-06-06 04:07:29.000): (ra, dec, distance) in (deg, deg, AU)
        (254.4656737, -22.66944306, 1.01475669)>, ref_epoch=J2010.000): (ra, dec, distance) in (deg, deg, pc)
        (10.00953937, 44.98650575, 10.00024407)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (9997.28930526, -20000.02695705, 100.27513164)>

Again, check that :meth:`~astropy.coordinates.SkyCoord.apply_space_motion` gives the same answer::

    >>> star_ref_sc = SkyCoord(star_ref, obstime=ast_hcrs.ref_epoch)
    >>> star_ref_sc.apply_space_motion(new_obstime=ast_hcrs.obstime).hcrs  # doctest: +REMOTE_DATA
    <SkyCoord (HCRS: obstime=2012-06-06 04:07:29.000): (ra, dec, distance) in (deg, deg, pc)
        (10.00953937, 44.98650575, 10.00024407)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (9997.28930526, -20000.02695707, 100.27513164)>
