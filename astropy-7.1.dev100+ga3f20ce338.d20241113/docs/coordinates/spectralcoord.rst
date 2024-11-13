.. _astropy-spectralcoord:

Using the SpectralCoord Class
*****************************

.. warning::

    The |SpectralCoord| class is new in Astropy v4.1 and should be considered
    experimental at this time. Note that we do not fully support cases
    where the observer and target are moving relativistically relative to each
    other, so care should be taken in those cases. It is possible that there
    will be API changes in future versions of Astropy based on user feedback. If
    you have specific ideas for how it might be improved, please  let us know on
    the |astropy-dev mailing list| or at http://feedback.astropy.org.

The |SpectralCoord| class provides an interface for representing and
transforming spectral coordinates such as frequencies, wavelengths, and photon
energies, as well as equivalent Doppler velocities. While the plain |Quantity|
class can also represent these kinds of physical quantities, and allow
conversion via dedicated equivalencies (such as :ref:`u.spectral
<astropy-units-spectral-equivalency>` or the :ref:`u.doppler_*
<astropy-units-doppler-equivalencies>` equivalencies), |SpectralCoord| (which is
a sub-class of |Quantity|) aims to make this more straightforward, and can also
be made aware of the observer and target reference frames, allowing for example
transformation from telescope-centric (or topocentric) frames to e.g.
Barycentric or Local Standard of Rest (LSRK and LSRD) velocity frames.

Creating SpectralCoord Objects
==============================

Since the |SpectralCoord| class is a sub-class of |Quantity|, the simplest way
to initialize it is to provide a value (or values) and a unit, or an existing
|Quantity|::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SpectralCoord
    >>> sc1 = SpectralCoord(34.2, unit='GHz')
    >>> sc1
    <SpectralCoord 34.2 GHz>
    >>> sc2 = SpectralCoord([654.2, 654.4, 654.6] * u.nm)
    >>> sc2
    <SpectralCoord [654.2, 654.4, 654.6] nm>

At this point, we are not making any assumptions about the observer frame, or
the target that is being observed. As we will see in subsequent sections, more
information can be provided when initializing |SpectralCoord| objects, but first
we take a look at simple unit conversions with these objects.

Unit conversion
===============

By default, unit conversions between spectral units will work without having to
specify the :ref:`u.spectral <astropy-units-spectral-equivalency>` equivalency::

    >>> sc2.to(u.micron)
    <SpectralCoord [0.6542, 0.6544, 0.6546] micron>
    >>> sc2.to(u.eV)
    <SpectralCoord [1.89520328, 1.89462406, 1.89404519] eV>
    >>> sc2.to(u.THz)
    <SpectralCoord [458.25811373, 458.11805929, 457.97809044] THz>

As is the case with |Quantity| and the :ref:`Doppler equivalencies
<astropy-units-doppler-equivalencies>`, it is also possible to convert these
absolute spectral coordinates into velocities, assuming a particular rest
frequency or wavelength (such as that of a spectral line). For example, to
convert the above values into velocities relative to the Halpha line at 656.65
nm, assuming the optical Doppler convention, you can do::

    >>> sc3 = sc2.to(u.km / u.s,
    ...              doppler_convention='optical',
    ...              doppler_rest=656.65 * u.nm)
    >>> sc3
    <SpectralCoord
       (doppler_rest=656.65 nm
        doppler_convention=optical)
      [-1118.5433977 , -1027.23373258,  -935.92406746] km / s>

The rest value for the Doppler conversion as well as the convention to use are
stored in the resulting ``sc3`` |SpectralCoord| object. You can then convert
back to frequency without having to specify them again::

    >>> sc3.to(u.THz)
    <SpectralCoord
       (doppler_rest=656.65 nm
        doppler_convention=optical)
      [458.25811373, 458.11805929, 457.97809044] THz>

or you can explicitly specify a different convention or rest value to use::

    >>> sc3.to(u.km / u.s, doppler_convention='relativistic')
    <SpectralCoord
       (doppler_rest=656.65 nm
        doppler_convention=relativistic)
      [-1120.63005892, -1028.99362163,  -937.38499411] km / s>

It is also possible to set ``doppler_convention`` and ``doppler_rest`` from the
start, even when creating a |SpectralCoord| in frequency, energy, or
wavelength::

    >>> sc4 = SpectralCoord(343 * u.GHz,
    ...                     doppler_convention='radio',
    ...                     doppler_rest=342.91 * u.GHz)
    >>> sc4.to(u.km / u.s)
    <SpectralCoord
       (doppler_rest=342.91 GHz
        doppler_convention=radio)
      -78.68338987 km / s>


Reference frame transformations
===============================

If you work with any kind of spectral data, you will often need to determine
and/or apply velocity corrections due to different frames of reference, or apply
or remove the effects of redshift. There are two main ways to do this using the
|SpectralCoord| class:

* You can specify or change the velocity offset or redshift
  between the observer and the target without having to specify the
  absolute observer and target, but rather specify a velocity difference.  For example, that you know that there
  is a velocity difference of 15km/s along the line of sight, or that you are
  observing a galaxy at z=3.2. This can be useful for quick analysis but
  will not determine any frame transformations (e.g. from topocentric to
  barycentric) for you.

* You can specify the absolute position of the observer and the target,
  as well as the date of observation, which means that |SpectralCoord| can
  then compute different frame transformations. If information about the
  observer and target are available, this is the recommended approach,
  although it requires you to specify more information when setting up the
  |SpectralCoord|

In the next two sections we will look at each of these in turn.

Specifying radial velocity or redshift manually
-----------------------------------------------

As an example, we will consider an example of a |SpectralCoord| which represents
frequencies which form the x-axis of a (small) spectrum. We happen to know that
the target that was observed appears to be at a redshift of z=0.5, and we will
assume that any frequency shifts due to the Earth's motion are unimportant. In
the reference frame of the telescope, the spectrometer provides 10 values
between 500 and 900nm::

    >>> import numpy as np
    >>> wavs = SpectralCoord(np.linspace(500, 900, 9) * u.nm, redshift=0.5)
    >>> wavs  # doctest: +FLOAT_CMP
    <SpectralCoord
       (observer to target:
          radial_velocity=115304.79153846153 km / s
          redshift=0.5)
      [500., 550., 600., 650., 700., 750., 800., 850., 900.] nm>

We have set redshift=0.5 here so that we can keep track of what frame of reference
our spectral values are in. The ``radial_velocity`` property gives the recession
velocity equivalent to that redshift, and it is indeed large enough that we don't need
to worry about the rotation of the Earth on itself around the Sun (which would be
at most a ~30km/s contribution).

.. note:: In the context of |SpectralCoord|, we use the full relativistic relation
          between redshift and velocity, i.e. :math:`1 + z = \sqrt{(1 + v/c)/(1 - v/c)}`

We now want to shift the wavelengths so that they would be in the rest frame of
the galaxy. We can do this using the
:meth:`~astropy.coordinates.SpectralCoord.to_rest` method::

    >>> wavs_rest = wavs.to_rest()
    >>> wavs_rest
    <SpectralCoord
       (observer to target:
          radial_velocity=0.0 km / s
          redshift=0.0)
      [333.33333333, 366.66666667, 400.        , 433.33333333, 466.66666667,
       500.        , 533.33333333, 566.66666667, 600.        ] nm>

The wavelengths have decreased by 1/3, which is what we expect for z=0.5. Note
that the ``redshift`` and ``radial_velocity`` properties are now zero, since we
are in the reference frame of the target. We can also use the
:meth:`~astropy.coordinates.SpectralCoord.with_radial_velocity_shift` method to more
generically apply redshift and velocity corrections. The simplest way to use
this method is to give a single value that will be applied to the target - if
this value does not have units, it is interpreted as a redshift::

    >>> wavs_orig = wavs_rest.with_radial_velocity_shift(0.5)
    >>> wavs_orig  # doctest: +FLOAT_CMP
    <SpectralCoord
       (observer to target:
          radial_velocity=115304.79153846153 km / s
          redshift=0.5)
      [500., 550., 600., 650., 700., 750., 800., 850., 900.] nm>

This returns an object equivalent to the one we started with, since we've
re-applied a redshift of 0.5. We could also provide a velocity as a |Quantity|::

    >>> wavs_rest.with_radial_velocity_shift(100000 * u.km / u.s)
    <SpectralCoord
       (observer to target:
          radial_velocity=100000.0 km / s
          redshift=0.41458078170200463)
      [471.52692723, 518.67961996, 565.83231268, 612.9850054 , 660.13769813,
       707.29039085, 754.44308357, 801.5957763 , 848.74846902] nm>

which shifts the values to a frame of reference at a redshift of approximately
0.33 (that is, if the spectrum did contain a contribution from an object at
z=0.33, these would be the rest wavelengths for that object.

Specifying an observer and a target explicitly
----------------------------------------------

.. testsetup::

    >>> from astropy.coordinates import EarthLocation
    >>> location = EarthLocation(2225015.30883296, -5440016.41799762, -2481631.27428014, unit='m')

To use the more advanced functionality in |SpectralCoord|, including the ability
to easily transform between different well-defined velocity frames, you will
need to give it information about the location (and optionally velocity) of
the observer and target. This is done by passing either coordinate frame objects
or |SkyCoord| objects. To take a concrete example, let's assume that we are now
observe the source T Tau using the ALMA telescope. To create an observer object
corresponding to this, we can make use of the |EarthLocation| class::

    >>> from astropy.coordinates import EarthLocation
    >>> location = EarthLocation.of_site('ALMA')  # doctest: +SKIP
    >>> location  # doctest: +FLOAT_CMP
    <EarthLocation (2225015.30883296, -5440016.41799762, -2481631.27428014) m>

The three values in meters are geocentric coordinates, i.e. the 3D coordinates
relative to the center of the Earth. See |EarthLocation| for more details about
the different ways of creating these kinds of objects.

Once you have done this, you will need to convert ``location`` to a coordinate
object using the :meth:`~astropy.coordinates.EarthLocation.get_itrs` method,
which takes the observation time (which is important to know for any kind of
velocity frame transformation)::

    >>> from astropy.time import Time
    >>> alma = location.get_itrs(obstime=Time('2019-04-24T02:32:10'))
    >>> alma  # doctest: +FLOAT_CMP
    <ITRS Coordinate (obstime=2019-04-24T02:32:10.000, location=(0., 0., 0.) km): (x, y, z) in m
        (2225015.30883296, -5440016.41799762, -2481631.27428014)>

ITRS here stands for International Terrestrial Reference System which is a 3D
coordinate frame centered on the Earth's center and rotating with the Earth, so
the observatory will be stationary in this frame of reference.

For the target, the simplest way is to use the |SkyCoord| class::

    >>> from astropy.coordinates import SkyCoord
    >>> ttau = SkyCoord('04h21m59.43s +19d32m06.4', frame='icrs',
    ...                 radial_velocity=23.9 * u.km / u.s,
    ...                 distance=144.321 * u.pc)

In this case we specified a radial velocity and a distance for the target (using
the `T Tauri SIMBAD entry
<https://simbad.unistra.fr/simbad/sim-id?Ident=T+Tauri>`_, but it is also
possible to not specify these, which means the target is assumed to be
stationary in the frame in which it is observed, and are assumed to be at large
distance from the Sun (such that any parallax effects would be unimportant if
relevant). The radial velocity is assumed to be in the frame used to define the
target location, so it is relative to the ICRS origin (the Solar System
barycenter) in the above case.

We now define a set of frequencies corresponding to the channels in which fluxes
have been measured (for the purposes of the example here we will assume we have only
11 frequencies)::

    >>> sc_ttau = SpectralCoord(np.linspace(200, 300, 11) * u.GHz,
    ...                         observer=alma, target=ttau)  # doctest: +IGNORE_WARNINGS
    >>> sc_ttau  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer:  <ITRS Coordinate (obstime=2019-04-24T02:32:10.000, location=(0., 0., 0.) km): (x, y, z) in m
                      (2225015.30883296, -5440016.41799762, -2481631.27428014)
                   (v_x, v_y, v_z) in km / s
                      (0., 0., 0.)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=41.03594953774002 km / s
          redshift=0.00013689056329480032)
      [200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300.] GHz>

We can already see above that |SpectralCoord| has computed the difference in
velocity between the observatory and T Tau, which includes the motion of the
observatory around the Earth, the motion of the Earth around the Solar System
barycenter, and the radial velocity of T Tau relative to the Solar System
barycenter. We can get this value directly with::

    >>> sc_ttau.radial_velocity  # doctest: +FLOAT_CMP +REMOTE_DATA
    <Quantity 41.03594948 km / s>

If you work with any kind of spectral data, you will often need to determine
and/or apply velocity corrections due to different frames of reference. For
example if you have observations of the same object on the sky taken at
different dates, it is common to transform these to a common velocity frame of
reference, so that your spectral coordinates are those that would have applied
if the observer had been stationary relative to e.g. the Solar System
Barycenter. You may also want to transform your spectral coordinates so that
they would be in a frame at rest relative to the local standard of rest (LSR),
the center of the Milky Way, the Local Group, or even the Cosmic Microwave
Background (CMB) dipole.

We can transform our frequencies for the observations of T Tau to different
velocity frames using the
:meth:`~astropy.coordinates.SpectralCoord.with_observer_stationary_relative_to`
method. This method can take the name of an existing coordinate/velocity frame,
a :class:`~astropy.coordinates.BaseCoordinateFrame` instance, or any arbitrary
3D position and velocity coordinate object defined either as a
:class:`~astropy.coordinates.BaseCoordinateFrame` or a |SkyCoord| object. Most
commonly-used frames are accessible using strings. For example to transform to a
velocity frame stationary with respect to the center of the Earth (so removing
the effect of the Earth's rotation), we can use the ``'gcrs'`` which stands for
*Geocentric Celestial Reference System* (GCRS)::

    >>> sc_ttau.with_observer_stationary_relative_to('gcrs')  # doctest: +SKIP
    <SpectralCoord
       (observer: <GCRS Coordinate (obstime=2019-04-24T02:32:10.000, obsgeoloc=(0., 0., 0.) m, obsgeovel=(0., 0., 0.) m / s): (x, y, z) in m
                      (-5878853.86171412, -192921.84773269, -2470794.19765021)
                   (v_x, v_y, v_z) in km / s
                      (4.33251262e-09, 8.96175625e-08, -1.49258412e-08)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=40.674086368345165 km / s
          redshift=0.00013568335316072044)
      [200.00024141, 210.00025348, 220.00026555, 230.00027762, 240.00028969,
       250.00030176, 260.00031383, 270.0003259 , 280.00033797, 290.00035004,
       300.00036211] GHz>

As you can see, the frequencies have changed slightly, which is because we have
removed the Doppler shift caused by the Earth's rotation (this can also be seen
in the ``radial_velocity`` property, which has changed by ~0.35 km/s. To use a
velocity reference frame relative to the Solar System barycenter, which is the
origin of the *International Celestial Reference System* (ICRS) system, we can use::

    >>> sc_ttau.with_observer_stationary_relative_to('icrs')  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <ICRS Coordinate: (x, y, z) in m
                      (-1.25867767e+11, -7.48979688e+10, -3.24757657e+10)
                   (v_x, v_y, v_z) in km / s
                      (0., 0., 0.)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=23.9 km / s
          redshift=7.97249967898761e-05)
      [200.0114322 , 210.01200381, 220.01257542, 230.01314703, 240.01371864,
       250.01429025, 260.01486186, 270.01543347, 280.01600508, 290.01657669,
       300.0171483 ] GHz>

Note that in this case the total radial velocity between the observer and the
target matches what we specified when we set up the target, since it was defined
relative to the ICRS origin (the Solar System barycenter). The observer location
is still as before, but the observer velocity is now ~10-20 km/s in x, y, and z,
which is because the observer is now stationary relative to the barycenter so has
a significant velocity relative to the surface of the Earth.

We can also transform the frequencies to the Kinematic Local Standard of Rest
(LSRK) frame of reference, which is a reference frame commonly used in some
branches of astronomy (such as radio astronomy)::

    >>> sc_ttau.with_observer_stationary_relative_to('lsrk')  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <LSRK Coordinate: (x, y, z) in m
                      (-1.25867767e+11, -7.48979688e+10, -3.24757657e+10)
                   (v_x, v_y, v_z) in km / s
                      (0., 0., 0.)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=12.50698856018455 km / s
          redshift=4.171969349386906e-05)
      [200.01903338, 210.01998505, 220.02093672, 230.02188839, 240.02284006,
       250.02379172, 260.02474339, 270.02569506, 280.02664673, 290.0275984 ,
       300.02855007] GHz>


See :ref:`spectralcoord-common-frames` for a list of common velocity frames
available as strings on the |SpectralCoord| class.

Since we can give any arbitrary |SkyCoord| to the
:meth:`~astropy.coordinates.SpectralCoord.with_observer_stationary_relative_to`
method, we can also specify the target itself, to find the frequencies in the
rest frame of the target::

    >>> sc_ttau_targetframe = sc_ttau.with_observer_stationary_relative_to(sc_ttau.target)  # doctest: +REMOTE_DATA
    >>> sc_ttau_targetframe  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <ICRS Coordinate: (x, y, z) in m
                      (-1.25867767e+11, -7.48979688e+10, -3.24757657e+10)
                   (v_x, v_y, v_z) in km / s
                      (9.34149908, 20.49579745, 7.99178839)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=0.0 km / s
          redshift=0.0)
      [200.02737811, 210.02874702, 220.03011592, 230.03148483, 240.03285374,
       250.03422264, 260.03559155, 270.03696045, 280.03832936, 290.03969826,
       300.04106717] GHz>

The ``radial_velocity``, which is the velocity offset between observer and
target, is now zero.

|SpectralCoord| is intended to be versatile and be useful for representing any spectral
values - not just the x-axis of a spectrum, but also for example the
frequencies of spectral features. For example, if we now consider that we found a
spectral feature that appears to have components at the following frequencies
in the frame of reference of the telescope::

    >>> sc_feat = SpectralCoord([115.26, 115.266, 115.267] * u.GHz,
    ...                         observer=alma, target=ttau)  # doctest: +IGNORE_WARNINGS

We can convert these to the rest frame of the target using::

    >>> sc_feat_rest = sc_feat.with_observer_stationary_relative_to(sc_feat.target)  # doctest: +REMOTE_DATA
    >>> sc_feat_rest  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <ICRS Coordinate: (x, y, z) in m
                      (-1.25867767e+11, -7.48979688e+10, -3.24757657e+10)
                   (v_x, v_y, v_z) in km / s
                      (9.34149908, 20.49579745, 7.99178839)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=0.0 km / s
          redshift=0.0)
      [115.27577801, 115.28177883, 115.28277896] GHz>

The frequencies are very close to the rest frequency of the 12CO J=1-0 molecular line transition,
which is 115.2712018 GHz. However, they are not exactly the same, so if the features we see are
indeed from 12CO, then they are Doppler shifted compared to what we consider the rest frame of
T Tau. We can convert these frequencies to velocities assuming the Doppler shift equation
(in this case with the radio convention)::

    >>> sc_feat_rest.to(u.km / u.s, doppler_convention='radio', doppler_rest=115.27120180 * u.GHz)  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <ICRS Coordinate: (x, y, z) in m
                      (-1.25867767e+11, -7.48979688e+10, -3.24757657e+10)
                   (v_x, v_y, v_z) in km / s
                      (9.34149908, 20.49579745, 7.99178839)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=0.0 km / s
          redshift=0.0
        doppler_rest=115.2712018 GHz
        doppler_convention=radio)
      [-11.90160348, -27.50828539, -30.10939905] km / s>

Note that these resulting velocities are different from the ``radial_velocity``
property (which is still zero here) - the latter is the difference in velocity
between observer and target, while the former are how much the spectral values
are Doppler shifted by relative to the rest frequency or wavelength.

So if the features are indeed from 12CO, they have velocities of approximately -11.9, -27.5 and
-30.1 km/s relative to the T tau rest frame.

.. _spectralcoord-common-frames:

Common velocity frames
======================

Any valid astropy coordinate frame can be passed to the
:meth:`~astropy.coordinates.SpectralCoord.with_observer_stationary_relative_to`
method, including string aliases such as ``icrs``. Below we list some of the
frames commonly used to define spectral coordinates in:

The velocity frames available as constants on the |SpectralCoord| class are:

========================== =================================================
Frame name                 Description
========================== =================================================
``'gcrs'``                 Geocentric frame (defined as stationary relative to the GCRS origin)
``'icrs'``                 Barycentric frame (defined as stationary relative to the ICRS origin)
``'hcrs'``                 Heliocentric frame (defined as stationary relative to the HCRS origin)
``'lsrk``                  Kinematic Local Standard of Rest (LSRK),
                           defined as having a velocity of 20 km/s towards
                           18h +30d (B1900) relative to the Solar System
                           Barycenter [1]_.
``'lsrd'``                 Dynamical Local Standard of Rest (LSRD),
                           defined as having a velocity of U=9 km/s,
                           V=12 km/s, and W=7 km/s in Galactic coordinates
                           (equivalent to 16.552945 km/s towards l=53.13
                           and b=25.02) [2]_.
``'lsr'``                  A more recent definition of the Local Standard
                           of rest, with U=11.1 km/s,
                           V=12.24 km/s, and W=7.25 km/s in Galactic coordinates [3]_.
========================== =================================================

Defining custom velocity frames
===============================

As mentioned in the earlier examples on this page, it is possible to pass any
arbitrary :class:`~astropy.coordinates.BaseCoordinateFrame` or |SkyCoord| object
to the :meth:`~astropy.coordinates.SpectralCoord.with_observer_stationary_relative_to` method,
and the observer will be updated to be stationary relative to those coordinates.
As an example, we can define an object that can be used to define a velocity
frame that moves with the local group of galaxies. There is not a unique definition
of this, but for the purposes of this example we use the IAU 1976-recommended
value which states that the Solar System barycenter is moving at 300 km/s towards
l=90 and b=0 in the velocity frame of the local group of galaxies [4]_. Given
this value, we can define the velocity frame using::

    >>> from astropy.coordinates import Galactic
    >>> localgroup_frame = Galactic(u=0 * u.km, v=0 * u.km, w=0 * u.km,
    ...                             U=0 * u.km / u.s, V=-300 * u.km / u.s, W=0 * u.km / u.s,
    ...                             representation_type='cartesian',
    ...                             differential_type='cartesian')

Note that here we specify the velocity as -300, because what we need here is the
velocity of the local group relative to the Solar System barycenter. With this
object, we can then transform a |SpectralCoord| so that the observer is stationary
in that frame of reference::

    >>> sc_ttau.with_observer_stationary_relative_to(localgroup_frame)  # doctest: +FLOAT_CMP +REMOTE_DATA
    <SpectralCoord
       (observer: <Galactic Coordinate: (u, v, w) in m
                      (8.8038652e+10, -5.31344273e+10, 1.09238291e+11)
                   (U, V, W) in km / s
                      (-1.42108547e-14, -300., 2.84217094e-14)>
        target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, pc)
                    (65.497625, 19.53511111, 144.321)
                 (radial_velocity) in km / s
                    (23.9,)>
        observer to target (computed from above):
          radial_velocity=42.33062895275233 km / s
          redshift=0.00014120974955456056)
      [199.99913628, 209.9990931 , 219.99904991, 229.99900673, 239.99896354,
       249.99892036, 259.99887717, 269.99883398, 279.9987908 , 289.99874761,
       299.99870443] GHz>

References
==========

.. [1] Meeks, M. L. 1976, *Methods of experimental physics. Vol._12.
       Astrophysics. Part C: Radio observations*, Section 6.1 by Gordon, M. A.
       `[ADS] <https://ui.adsabs.harvard.edu/abs/1976mep..book.....M>`__.
.. [2] Delhaye, J. 1965, *Galactic Structure*. Edited by Adriaan Blaauw and
       Maarten Schmidt. Published by the University of Chicago Press, p61
       `[ADS] <https://ui.adsabs.harvard.edu/abs/1965gast.book...61D>`__.
.. [3] Schönrich, R., Binney, J., & Dehnen, W. 2010, MNRAS, 403, 1829
       `[ADS] <https://ui.adsabs.harvard.edu/abs/2010MNRAS.403.1829S>`__.
.. [4] *Transactions of the IAU Vol. XVI B Proceedings of the 16th General
      Assembly, Reports of Meetings of Commissions: Comptes Rendus
      Des Séances Des Commissions, Commission 28*.
      `[DOI] <https://doi.org/10.1017/S0251107X00002406>`__

.. The following frames are defined in FITS WCS and may be added here in future:
..
.. ``GALACTOCENTRIC_KLB1986`` Galactocentric frame defined as having a velocity
..                            of 220 km/s towards l=90 and b=0 relative to
..                            the Solar System Barycenter [3]_.
.. ``LOCALGROUP_IAU1976``     Velocity frame representing the motion of the
..                            Local Group of galaxies, and defined as having a velocity
..                            of 300 km/s towards l=90 and b=0 relative to
..                            the Solar System Barycenter [4]_.
.. ``CMBDIPOL_WMAP1``         Velocity frame representing the motion of the
..                            cosmic microwave background (CMB) dipole based on the
..                            1-year WMAP data, and defined as a temperature
..                            difference of 3.346mK (corresponding to approximately
..                            368 km/s) in the direction of l=263.85, b=48.25 [5]_
.. .. [3] Kerr, F. J., & Lynden-Bell, D. 1986, MNRAS, 221, 1023
..       `[ADS] <https://ui.adsabs.harvard.edu/abs/1986MNRAS.221.1023K>`__.
.. .. [5] Bennett, C. L., Halpern, M., Hinshaw, G., et al. 2003, ApJS, 148, 1
..       `[ADS] <https://ui.adsabs.harvard.edu/abs/2003ApJS..148....1B>`__.
