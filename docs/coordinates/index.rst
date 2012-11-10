*******************************************************
Astronomical Coordinate Systems (`astropy.coordinates`)
*******************************************************

Introduction
============

The `~astropy.coordinates` package provides classes for representing celestial
coordinates, as well as tools for converting between standard systems in a
uniform way.

.. note::
	The current `~astropy.coordinates` framework only accepts scalar
	coordinates, i.e. one coordinate per object.  In the next release it will
	be expanded to accept arrays of coordinates.


Getting Started
===============

Coordinate objects are intantiated with a flexible and natural approach::

    >>> from astropy import coordinates as apc
    >>> from astropy import units as u
    >>> apc.Coordinates(ra=10.68458, dec=41.26917, unit=u.degree)
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>
    >>> apc.ICRSCoordinates('00h42m44.3s +41d16m9s')
    <ICRSCoordinates RA=10.68458 deg, Dec=41.26917 deg>

The individual components of a coordinate are `~astropy.coordinates.angles.Angle`
objects, and their values are acceessed using special attributes::

	>>> c = apc.Coordinates(ra=10.68458, dec=41.26917, unit=u.degree)
    >>> c.ra.hours
    0.7123053333333333
    >>> c.dec.radians
    0.7202828960652683
    >>> c.ra.hms
    (0.0, 42, 44.2992000000001)

To convert to some other coordinate system, the easiest method is to use
attribute-style access with short names for the builtin systems, but explicit
transformations via the `transform_to` method are also available::

    >>> c.galactic
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>
    >>> c.transform_to(apc.GalacticCoordinates)
    <GalacticCoordinates l=121.17422 deg, b=-21.57283 deg>


Using `astropy.coordinates`
===========================

Content


Creating Coordinate Objects
---------------------------

Content


Distances
---------

Content

Transforming Between Systems
----------------------------

Content

Custom Coordinate Classes
-------------------------

Content


See Also
========

Some references particularly useful in understanding subtleties of the
coordinate systems implemented here include:

* `Standards Of Funamental Astronomy <http://www.iausofa.org/>`_
	The definitive implementation of IAU-defined algorithms.  The "SOFA Tools
	for Earth Attitude" document is particularly valuable for understanding
	the latest IAU standards in detail.
* `USNO Circular 179 <http://www.usno.navy.mil/USNO/astronomical-applications/publications/circ-179>`_
	A useful guide to the IAU 2000/2003 work surrounding ICRS/IERS/CIRS and
	related problems in precision coordinate system work.
* Meeus, J. "Astronomical Algorithms"
	A valuable text describing details of a wide range of coordinate-related
	problems and concepts.



Reference/API
=============

.. automodapi:: astropy.coordinates
