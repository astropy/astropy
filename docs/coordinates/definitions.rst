.. _astropy-coordinates-definitions:

Important Definitions
*********************

For reference, below, we define some key terms as they are used in
`~astropy.coordinates`, due to some ambiguities that exist in the
colloquial use of these terms. Chief among these terms is the concept
of a "coordinate system." To some members of the community, "coordinate
system" means the *representation* of a point in space (e.g., "Cartesian
coordinate system" is different from "Spherical polar coordinate
system"). Another use of "coordinate system" is to mean a unique
reference frame with a particular set of reference points (e.g., "the
ICRS coordinate system" or the "J2000 coordinate system"). This second
meaning is further complicated by the fact that such systems use quite
different ways of defining a frame.

Because of the likelihood of confusion between these meanings of
"coordinate system," `~astropy.coordinates` avoids this term wherever
possible, and instead adopts the following terms (loosely inspired by
the IAU2000 resolutions on celestial coordinate systems):

* A "Coordinate Representation" is a particular way of describing a unique
  point in a vector space. (Here, this means three-dimensional space, but future
  extensions might have different dimensionality, particularly if relativistic
  effects are desired.) Examples include Cartesian coordinates, cylindrical
  polar, or latitude/longitude spherical polar coordinates. Note that this term
  applies to the positions, *not* their velocities or other derivatives (which
  are represented as "differential" classes).

* A "Reference System" is a scheme for orienting points in a space and
  describing how they transforms to other systems. Examples include the ICRS,
  equatorial coordinates with mean equinox, or the WGS84 geoid for
  latitude/longitude on the Earth.

* A "Coordinate Frame," "Reference Frame," or just "Frame" is a specific
  realization of a reference system (e.g., the ICRF, or J2000 equatorial
  coordinates). For some systems, there may be only one meaningful frame, while
  others may have many different frames (differentiated by something like a
  different equinox, or a different set of reference points).

* A "Coordinate" is a combination of all of the above that specifies a unique
  point.
