.. include:: references.txt

.. _astropy-coordinates-velocities:

Velocities (Radial and Proper Motions) in Coordinates
*****************************************************

.. warning::
    Velocities support, new in Astropy v2.0, is an experimental feature and is
    subject to change based on user feedback.  While we do not expect major API
    changes, the possibility exists based on the precedent of earlier changes
    in the ``coordinates`` subpackage due to user feedback in previous versions
    of Astropy.

Creating frame objects with velocities
======================================

content

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
