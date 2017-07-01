.. include:: references.txt

.. _astropy-coordinates-velocities:

Working with velocity data in Astropy coordinates
*************************************************

.. warning::
    Velocities support, new in Astropy v2.0, is an experimental feature and is
    subject to change based on user feedback.  While we do not expect major API
    changes, the possibility exists based on the precedent of earlier changes
    in the ``coordinates`` subpackage based on user feedback from previous
    versions of Astropy.

.. _astropy-coordinate-custom-frame-with-velocities:

Creating frame objects with velocity data
=========================================

The coordinate frame classes now support storing and transforming velocity data
along with positional coordinate data. Similar to the positional data --- that
use the ``Representation`` classes to abstract away the particular
representation and allow re-representing from, e.g., Cartesian to Spherical
--- the velocity information makes use of ``Differential`` classes to do the
same. For more information about the differential classes, see
:ref:`astropy-coordinates-differentials`. Also like the positional data, the
names of the differential (velocity) components depend on the particular
coordinate frame.

The default differential for most frames is the
`~astropy.coordinates.SphericalCosLatDifferential` class, meaning that most
frames expect velocity data in the form of two proper motion components and/or a
radial velocity. When supported, the proper motion components all begin with
``pm_`` and, by default, the longitudinal component should already include the
``cos(latitude)`` term. For example, the proper motion components for the ICRS
frame are (``pm_ra_cosdec``, ``pm_dec``)::

    >>> coord.ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...            pm_ra_cosdec=11*u.mas/u.yr, pm_dec=-27*u.mas/u.yr)
    <ICRS Coordinate: (ra, dec) in deg
        ( 8.67,  53.09)
     (pm_ra_cosdec, pm_dec) in mas / yr
        ( 11., -27.)>

Like the positional data, velocity data must be passed in as
`~astropy.units.Quantity` objects.

TODO: more stuff...

.. _astropy-coordinate-transform-with-velocities:

Transforming frames with velocities
===================================

overview

Affine Transformations
----------------------

details

Finite Difference Transformations
---------------------------------

details, particularly *examples* of the finite difference problem




``SkyCoord`` support for Velocities
===================================

|skycoord| currently does *not* support velocities as of Astropy v2.0.  This is
an intentional choice, allowing the "power-user" community to provide feedback
on the API and functionality in the frame-level classes before it is adopted in
|skycoord| (currently planned for the next Astropy version, v3.0).

Radial Velocity Corrections
===========================

Sperately from the above, Astropy supports computing barycentric or heliocentric
radial velocity corrections.  While in the future this may simply be a
high-level convenience function using the framework described above, the
current implementation is independent to ensure sufficient accuracy (see the
`~astropy.coordinates.SkyCoord.radial_velocity_correction` API docs for
details).

A usage example of this functionality is::

    sc = SkyCoord()
    MORE_NEEDED = True
