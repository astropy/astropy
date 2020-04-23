.. include:: references.txt

.. _astropy-spectralcoord:

Using the SpectralCoord High-Level Class
****************************************

The |SpectralCoord| class provides an interface for representing and
transforming spectral coordinates such as frequencies, wavelengths, energies,
Doppler velocities, and so on. While the plain |Quantity| class can also
represent these kinds of physical quantities, and allow conversion via dedicated
equivalencies (such as :ref:`u.spectral <astropy-units-doppler-equivalencies>`
or the :ref:`u.doppler_* <astropy-units-doppler-equivalencies>` equivalencies),
|SpectralCoord| (which is a sub-class of |Quantity|) aims to make this more
straightforward, and can also be made aware of the observer and target reference
frames, allowing for example transformation from telescope-centric (or
topocentric) frames to e.g. Barycentric or Local Standard of Rest (LSRK and
LSRD) velocity frames.

Creating SpectralCoord Objects
==============================

Since the |SpectralCoord| class is a sub-class of |Quantity|, the simplest way
to initialize it is to provide a value and a unit, or an existing |Quantity|::

    >>> from astropy import units as u
    >>> from astropy.coordinates import SpectralCoord
    >>> sc1 = SpectralCoord(34.2, unit='GHz')
    >>> sc1
    <SpectralCoord 34.2 GHz,
         radial_velocity=0.0 km / s,
         redshift=0.0,
         doppler_rest=None,
         doppler_convention=None,
         observer=None,
         target=None>
    >>> sc2 = SpectralCoord(654.2 * u.nm)
    >>> sc2
    <SpectralCoord 654.2 nm,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=None,
        doppler_convention=None,
        observer=None,
        target=None>

As we will see in subsequent sections, more information can be provided when
initializing |SpectralCoord| objects, but first we take a look at simple unit
conversions with these objects.

Unit conversion
===============

By default, unit conversions between spectral units will work without having to
specify the :ref:`u.spectral <astropy-units-doppler-equivalencies>` equivalency::

    >>> sc1.to(u.micron)
    <SpectralCoord 8765.86134503 micron,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=None,
        doppler_convention=None,
        observer=None,
        target=None>
    >>> sc1.to(u.keV)
    <SpectralCoord 1.41439835e-07 keV,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=None,
        doppler_convention=None,
        observer=None,
        target=None>
    >>> sc2.to(u.THz)
    <SpectralCoord 458.25811373 THz,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=None,
        doppler_convention=None,
        observer=None,
        target=None>

As is the case with |Quantity| and the :ref:`Doppler equivalencies
<astropy-units-doppler-equivalencies>`, it is also posible to convert these
absolute spectral coordinates into velocities, assuming a particular rest
frequency or wavelength (such as that of a spectral line). For example, to
convert the above values into velocities relative to the Halpha line at
656.65 nm, assuming the optical Doppler convention, you can do::

    >>> sc3 = sc2.to(u.km / u.s,
    ...              doppler_convention='optical',
    ...              doppler_rest=656.65 * u.nm)
    >>> sc3
    <SpectralCoord -1118.5433977 km / s,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=656.65 nm,
        doppler_convention=optical,
        observer=None,
        target=None>

At this point, the rest value for the Doppler conversion as well as the
convention to use are stored in the resulting |SpectralCoord|. You can then
convert back to frequency without having to specify them again::

    >>> sc3.to(u.THz)
    <SpectralCoord 458.25811373 THz,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=656.65 nm,
        doppler_convention=optical,
        observer=None,
        target=None>

or you can explicitly specify a different convention or rest value to use::

    >>> sc3.to(u.km / u.s, doppler_convention='relativistic')  # doctest: +SKIP

.. TODO: the following example does not work for now

It is also possible to set ``doppler_convention`` and ``doppler_rest`` from the
start, even when creating a |SpectralCoord| in frequency, energy, or
wavelength::

    >>> sc4 = SpectralCoord(343 * u.GHz,
    ...                     doppler_convention='radio',
    ...                     doppler_rest=342.91 * u.GHz)
    >>> sc4.to(u.km / u.s)
    <SpectralCoord -78.68338987 km / s,
        radial_velocity=0.0 km / s,
        redshift=0.0,
        doppler_rest=342.91 GHz,
        doppler_convention=radio,
        observer=None,
        target=None>

Specifying an observer and a target
===================================

To use the more advanced funtionality in |SpectralCoord|, you will next need to
give it information about the location (and optionally velocity) of the observer
and target. This is done by passing either coordinate frame objects or
|SkyCoord| objects. To take a concrete example, let's assume that we are
observing the source T Tau using a telescope at Greenwich observatory. To create
an observer object corresponding to this, we can make use of the |EarthLocation|
class::

    >>> from astropy.coordinates import EarthLocation
    >>> location = EarthLocation.of_site('greenwich')  # doctest: +REMOTE_DATA
    >>> location  # doctest: +REMOTE_DATA
    <EarthLocation (3980608.90246817, -102.47522911, 4966861.27310068) m>

See |EarthLocation| for more details about the different ways of creating
these kinds of objects.

Once you have done this, you will need to convert it to a coordinate object
using the :meth:`~astropy.coordinates.EarthLocation.get_itrs` method, which
takes the observation time (which is important to know for any kind of
spctral transformation)::

    >>> from astropy.time import Time
    >>> greenwich = location.get_itrs(obstime=Time('2019-04-24T02:32:10'))  # doctest: +REMOTE_DATA
    >>> greenwich  # doctest: +REMOTE_DATA
    <ITRS Coordinate (obstime=2019-04-24T02:32:10.000): (x, y, z) in m
        (3980608.90246817, -102.47522911, 4966861.27310068)>

For the target, the simplest way is to use the |SkyCoord| class::

    >>> from astropy.coordinates import SkyCoord
    >>> ttau = SkyCoord('04h21m59.43s +19d32m06.4', frame='icrs',
    ...                 radial_velocity=23.9 * u.km / u.s,
    ...                 distance=144.321 * u.pc)  # doctest: +REMOTE_DATA

.. TODO: make it so distance is optional

In this case we specified a radial velocity and a distance for the target (using
the `T Tauri SIMBAD entry
<http://simbad.u-strasbg.fr/simbad/sim-id?Ident=T+Tauri>`_, but it is also
possible to not specify these, which means the target is assumed to be
stationary in the frame in which it is observed, and are assumed to be at large
distance from the Sun (such that any parallax effects would be unimportant if
relevant). The radial velocity is assumed to be in the frame used to define the
target location, so it is relative to the ICRS origin (the Solar System
barycenter) in this case.

Let's now imagine that we detected a spectral line towards T Tau at 654.2 nm. We
can define a |SpectralCoord| to represent this::

    >>> sc_ttau = SpectralCoord(656.8 * u.nm, observer=greenwich, target=ttau)  # doctest: +IGNORE_WARNINGS +REMOTE_DATA
    >>> sc_ttau  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 656.8 nm,
        radial_velocity=40.65452165815913 km / s,
        redshift=0.00013560888732617527,
        doppler_rest=None,
        doppler_convention=None,
        observer=ITRS,
        target=ICRS>

We can already see above that |SpectralCoord| has computed the difference in
velocity between the observatory and T Tau, which includes the motion of the
observatory around the Earth, the motion of the Earth around the Solar System
barycenter, and the radial velocity of T Tau relative to the Solar System
barycenter. We can get this value directly with::

    >>> sc_ttau.radial_velocity  # doctest: +REMOTE_DATA +FLOAT_CMP
    <Quantity 40.65452431 km / s>

Velocity frame transformations
==============================

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

Continuing the example in `Specifying an observer and a target`_, where we had a
wavelength represented as a |SpectralCoord| with the observer set to an
observatory on the Earth, we can transform this to different velocity frames
using the :meth:`~astropy.coordinates.SpectralCoord.in_observer_velocity_frame`
method. This method can take any arbitrary 3-d position and velocity coordinate
object defined either as a :class:`~astropy.coordinates.BaseCoordinateFrame`
or a |SkyCoord| object, or strings referring to existing celestial coordinate
frames. For example to transform to a velocity frame stationary
with respect to the center of the Earth (so removing the effect of the Earth's
rotation), we can use

    >>> sc_ttau.in_observer_velocity_frame('gcrs')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 656.80004286 nm,
        radial_velocity=40.67408630397957 km / s,
        redshift=0.00013567414795998494,
        doppler_rest=None,
        doppler_convention=None,
        observer=ITRS,
        target=ICRS>

To use a velocity reference frame relative to the Solar System barycenter
we can use::

    >>> sc_ttau.in_observer_velocity_frame('icrs')  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 656.76329337 nm,
        radial_velocity=23.89999997498657 km / s,
        redshift=7.972181866892252e-05,
        doppler_rest=None,
        doppler_convention=None,
        observer=ITRS,
        target=ICRS>

Note that in this case the total radial velocity matches what we specified
when we set up the target, since it was defined relative to the ICRS origin.

For other common velocity frames that don't necessarily follow the origin of
celestial coordinate frames, we provide shortcuts as constants on the
|SpectralCoord| object itself::

    >>> sc_ttau.in_observer_velocity_frame(SpectralCoord.LSRK_GORDON1975)  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 656.73833301 nm,
        radial_velocity=12.506991127466303 km / s,
        redshift=4.171883179084613e-05,
        doppler_rest=None,
        doppler_convention=None,
        observer=ITRS,
        target=ICRS>

See :ref:`spectralcoord-common-frames` for a list of frames available as
constants on the |SpectralCoord| class.

Finally, since we can give any arbitrary |SkyCoord| to the
:meth:`~astropy.coordinates.SpectralCoord.in_observer_velocity_frame` method,
we can also specify the target itself, to find a spectral coordinate in the
rest frame of the target:

    >>> sc_ttau_targetframe = sc_ttau.in_observer_velocity_frame(sc_ttau.target)  # doctest: +REMOTE_DATA
    >>> sc_ttau_targetframe  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 656.71093208 nm,
        radial_velocity=-8.203851870936837e-08 km / s,
        redshift=-2.7365104264687125e-13,
        doppler_rest=None,
        doppler_convention=None,
        observer=ITRS,
        target=ICRS>


At this point, the |SpectralCoord| value is the frequency of the spectral line
as it would be measured in a reference frame moving with the target. It seems
the line is probably Halpha at 656.65 nm! But the value is slightly off, which
might indicate that the Halpha emission is from material that is moving
relative to the star. We can estimate the line of sight velocity of this
material by using the Doppler equivalency:

    >>> sc_ha = sc_ttau_targetframe.to(u.km / u.s,
    ...                                doppler_convention='relativistic',
    ...                                doppler_rest=656.65 * u.nm)  # doctest: +REMOTE_DATA
    >>> sc_ha  # doctest: +REMOTE_DATA +FLOAT_CMP
    <SpectralCoord 27.81714706 km / s,
        radial_velocity=-8.203851870936837e-08 km / s,
        redshift=-2.7365104264687125e-13,
        doppler_rest=656.65 nm,
        doppler_convention=relativistic,
        observer=ITRS,
        target=ICRS>

Note that you can convert this to a plain |Quantity| using::

    >>> sc_ha.quantity  # doctest: +REMOTE_DATA +FLOAT_CMP
    <Quantity 27.81714971 km / s>

This tells us that if the emission is from material that is moving at
approximately 28km/s away from us relative to T Tau, so could be e.g. material
accreting onto the star.

.. _spectralcoord-common-frames:

Common velocity frames
======================

The velocity frames available as constants on the |SpectralCoord| class are:

========================== =================================================
Constant Name              Description
========================== =================================================
``GEOCENTRIC``             Defined as stationary relative to the GCRS origin
``BARYCENTRIC``            Defined as stationary relative to the ICRS origin
``HELIOCENTRIC``           Defined as stationary relative to the HCRS origin
``LSRK_GORDON1975``        Kinematic Local Standard of Rest (LSRK),
                           defined as having a velocity of 20 km/s towards
                           18h +30d (B1900) relative to the Solar System
                           Barycenter [1]_.
``LSRD_DELHAYE1965``       Dynamical Local Standard of Rest (LSRD),
                           defined as having a velocity of U=9 km/s,
                           V=12 km/s, and W=7 km/s in Galactic coordinates
                           (equivalent to 16.552945 km/s towards l=53.13
                           and b=25.02 [2]_.
``GALACTOCENTRIC_KLB1986`` Galactocentric frame defined as having a velocity
                           of 220 km/s towards l=90 and b=0 relative to
                           the Solar System Barycenter [3]_.
``LOCALGROUP_IAU1976``     Velocity frame representing the motion of the
                           Local Group of galaxies, and defined as having a velocity
                           of 300 km/s towards l=90 and b=0 relative to
                           the Solar System Barycenter [4]_.
``CMBDIPOL_WMAP1``         Velocity frame representing the motion of the
                           cosmic microwave background (CMB) dipole based on the
                           1-year WMAP data, and defined as a tempreature
                           difference of 3.346mK (corresponding to approximately
                           368 km/s) in the direction of l=263.85, b=48.25 [5]_
========================== =================================================

References
==========

.. [1] Meeks, M. L. 1976, *Methods of experimental physics. Vol._12.
       Astrophysics. Part C: Radio observations*, Section 6.1 by Gordon, M. A.
       `[ADS] <https://ui.adsabs.harvard.edu/abs/1976mep..book.....M>`__.
.. [2] Delhaye, J. 1965, *Galactic Structure*. Edited by Adriaan Blaauw and
       Maarten Schmidt. Published by the University of Chicago Press, p61
       `[ADS] <https://ui.adsabs.harvard.edu/abs/1965gast.book...61D>`__.
.. [3] Kerr, F. J., & Lynden-Bell, D. 1986, MNRAS, 221, 1023
       `[ADS] <https://ui.adsabs.harvard.edu/abs/1986MNRAS.221.1023K>`__.
.. [4] *Transactions of the IAU Vol. XVI B Proceedings of the 16th General
       Assembly, Reports of Meetings of Commissions: Comptes Rendus
       Des Séances Des Commissions, Commission 28*.
       `[DOI] <https://doi.org/10.1017/S0251107X00002406>`__
.. [5] Bennett, C. L., Halpern, M., Hinshaw, G., et al. 2003, ApJS, 148, 1
       `[ADS] <https://ui.adsabs.harvard.edu/abs/2003ApJS..148....1B>`__.
