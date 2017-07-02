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
--- the velocity data makes use of ``Differential`` classes to do the
same (for more information about the differential classes, see
:ref:`astropy-coordinates-differentials`). Also like the positional data, the
names of the differential (velocity) components depend on the particular
coordinate frame.

Most frames expect velocity data in the form of two proper motion components
and/or a radial velocity because the default differential for most frames is the
`~astropy.coordinates.SphericalCosLatDifferential` class. When supported, the
proper motion components all begin with ``pm_`` and, by default, the
longitudinal component is expected to already include the ``cos(latitude)``
term. For example, the proper motion components for the ICRS frame are
(``pm_ra_cosdec``, ``pm_dec``)::

    >>> from astropy.coordinates import ICRS
    >>> import astropy.units as u
    >>> ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...      pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr)
    <ICRS Coordinate: (ra, dec) in deg
        ( 8.67,  53.09)
     (pm_ra_cosdec, pm_dec) in mas / yr
        ( 4.8, -15.16)>
    >>> ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...      pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr,
    ...      radial_velocity=23.42*u.km/u.s)
    <ICRS Coordinate: (ra, dec) in deg
        ( 8.67,  53.09)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        ( 4.8, -15.16,  23.42)>

Like the positional data, velocity data must be passed in as
`~astropy.units.Quantity` objects.

TODO: how to change the differential class to specify other velocity components

.. _astropy-coordinate-transform-with-velocities:

Transforming frames with velocities
===================================

Transforming coordinate frame instances that contain velocity data is the same
as transforming position-only frame instances::

    >>> from astropy.coordinates import Galactic
    >>> icrs = ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...             pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr)
    >>> icrs.transform_to(Galactic) # doctest: +FLOAT_CMP
    <Galactic Coordinate: (l, b) in deg
        ( 120.38084191, -9.69872044)
     (pm_l_cosb, pm_b) in mas / yr
        ( 3.78957965, -15.44359693)

However, the details of how the velocity components are transformed depends on
the particular path taken through the frame transform graph. If all frames
between the initial frame and the desired frame support transformations with a
`~astropy.coordinates.BaseAffineTransform` subclass (i.e. are matrix
transformations or affine transformations), then the transformations can be
applied explicitly to the velocity data. If this is not the case, the velocity
transformation is computed numerically by finite-differencing the positional
transformation. See the subsections below for more details about these two
methods.

Affine Transformations
----------------------

Frame transformations that involve a rotation and/or an origin shift and/or
a velocity offset are implemented as affine transformations using the
`~astropy.coordinates.BaseAffineTransform` subclassese:
`~astropy.coordinates.StaticMatrixTransform`,
`~astropy.coordinates.DynamicMatrixTransform`, and
`~astropy.coordinates.AffineTransform`.

Matrix-only transformations (e.g., rotations such as
`~astropy.coordinates.ICRS` to `~astropy.coordinates.Galactic`) can be performed
on proper-motion-only data or full-space, 3D velocities::

    >>> icrs = ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...             pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr,
    ...             radial_velocity=23.42*u.km/u.s) # doctest: +FLOAT_CMP
    >>> icrs.transform_to(Galactic) # doctest: +FLOAT_CMP
    <Galactic Coordinate: (l, b) in deg
        ( 120.38084191, -9.69872044)
     (pm_l_cosb, pm_b, radial_velocity) in (mas / yr, mas / yr, km / s)
        ( 3.78957965, -15.44359693,  23.42)>

The same rotation matrix is applied to both the position vector and the velocity
vector. Any transformation that involves a velocity offset requires all 3D
velocity components (which typically requires specifying a distance as well),
for example, `~astropy.coordinates.ICRS` to `~astropy.coordinates.LSR`::

    >>> from astropy.coordinates import LSR
    >>> icrs = ICRS(ra=8.67*u.degree, dec=53.09*u.degree,
    ...             distance=117*u.pc,
    ...             pm_ra_cosdec=4.8*u.mas/u.yr, pm_dec=-15.16*u.mas/u.yr,
    ...             radial_velocity=23.42*u.km/u.s) # doctest: +FLOAT_CMP
    >>> icrs.transform_to(LSR) # doctest: +FLOAT_CMP
    <LSR Coordinate (v_bary=( 11.1,  12.24,  7.25) km / s): (ra, dec, distance) in (deg, deg, pc)
        ( 8.67,  53.09,  117.)
     (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
        (-24.51315607, -2.67935501,  27.07339176)>

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
