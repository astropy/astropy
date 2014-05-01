.. _sgr-example:

Example: Defining a new spherical coordinate system
---------------------------------------------------

This document describes in detail how to subclass and define a custom
spherical coordinate system, as mentioned in :doc:`designing` and the
docstring for `~astropy.coordinates.BaseCoordinateFrame`. In this
example, we will define a coordinate system defined by the plane of orbit of
the Sagittarius Dwarf Galaxy (hereafter Sgr; as defined in Majewski et
al. 2003).  The Sgr coordinate system is often referred to in terms of two
angular coordinates, :math:`\Lambda,B`.

We need to define a subclass of `~astropy.coordinates.BaseCoordinateFrame`
that knows the preferred names and units of the coordinate system angles.
In this case, these are Lambda and Beta. Then we have to define the
transformation from this coordinate system to some other built-in system.
Here we will use spherical Galactic coordinates as defined in
`~astropy.coordinates.Galactic`.

The first step is to create a new class, which we'll call
``Sagittarius`` and make it a subclass of
`~astropy.coordinates.BaseCoordinateFrame`::

    import astropy.coordinates as coord
    from astropy.coordinates import frame_transform_graph
    from astropy.coordinates.angles import rotation_matrix

    @frame_transform_graph.add_coord_name
    class Sagittarius(coord.BaseCoordinateFrame):
        """
        A Heliocentric spherical coordinate system defined by the orbit
        of the Sagittarius dwarf galaxy, as described in
            http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
        and further explained in
            http://www.astro.virginia.edu/~srm4n/Sgr/.

        Parameters
        ----------
        representation : `BaseRepresentation` or None
            A representation object or None to have no data (or use the other keywords)
        Lambda : `Angle`, optional, must be keyword
            The longitude for this object (`Beta` must also be given and `representation`
            must be None).
        Beta : `Angle`, optional, must be keyword
            The Declination for this object (`Lambda` must also be given and
            `representation` must be None).
        distance : `Quantity`, optional, must be keyword
            The Distance for this object along the line-of-sight.
            (`representation` must be None).

        """

        preferred_representation = coord.SphericalRepresentation
        preferred_attr_names = OrderedDict([('Lambda', 'lon'), ('Beta', 'lat'),
                                            ('distance', 'distance')])
        preferred_attr_units = {'Lambda': u.degree, 'Beta': u.degree}

Line by line, the first few are simply imports. Next, the decorator method
``@coord.frame_transform_graph.add_coord_name`` tells the global transform
graph to recognize this as a new coordinate name. Next we define the class
as subclass of `~astropy.coordinates.BaseCoordinateFrame`. Then we include
a descriptive docstring. The final three lines are class-level attributes
that specify the preferred representation -- e.g., spherical, cartesian, etc.
-- the names of the individual coordinates, and the preferred units.

Next we have to define the transformation to some other built-in coordinate
system; we will use Galactic coordinates. We can do this by defining functions
that return transformation matrices, or by simply defining function that accept
a coordinate and returns a new coordinate in the new system. We'll start by
constructing the rotation matrix, using the helper function
``rotation_matrix()`` ::

    # Define the Euler angles (from Law & Majewski 2010)
    phi = np.radians(180+3.75)
    theta = np.radians(90-13.46)
    psi = np.radians(180+14.111534)

    # Generate the rotation matrix using the x-convention (see Goldstein)
    D = rotation_matrix(phi, "z", unit=u.radian)
    C = rotation_matrix(theta, "x", unit=u.radian)
    B = rotation_matrix(psi, "z", unit=u.radian)
    sgr_matrix = np.array(B.dot(C).dot(D))

This is done at the module level, since it will be used by both the
transformation from Sgr to Galactic as well as the inverse from Galactic to
Sgr. Now we can define our first transformation function::

    # Galactic to Sgr coordinates
    @frame_transform_graph.transform(coord.FunctionTransform, coord.Galactic, Sagittarius)
    def galactic_to_sgr(gal_coord, sgr_frame):
        """ Compute the transformation from Galactic spherical to
            heliocentric Sgr coordinates.
        """

        l = np.atleast_1d(gal_coord.l.radian)
        b = np.atleast_1d(gal_coord.b.radian)

        X = cos(b)*cos(l)
        Y = cos(b)*sin(l)
        Z = sin(b)

        # Calculate X,Y,Z,distance in the Sgr system
        Xs, Ys, Zs = sgr_matrix.dot(np.array([X, Y, Z]))
        Zs = -Zs

        # Calculate the angular coordinates lambda,beta
        Lambda = np.arctan2(Ys,Xs)*u.radian
        Lambda[Lambda < 0] = Lambda[Lambda < 0] + 2.*np.pi*u.radian
        Beta = np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs))*u.radian

        return Sagittarius(Lambda=Lambda, Beta=Beta,
                           distance=gal_coord.distance)

The decorator ``@frame_transform_graph.transform(coord.FunctionTransform, coord.Galactic, Sagittarius)`` registers this function with the `~astropy.coordinates`
subpackage as a transformation.  Inside the function, we simply follow the
same procedure as detailed by David Law's `transformation code
<http://www.astro.virginia.edu/~srm4n/Sgr/code.html>`_. Note that in this
case, both coordinate systems are heliocentric, so we can simply copy any
distance from the `~astropy.coordinates.Galactic` object.

We then register the inverse transformation by using the transpose of the
rotation matrix (which is faster to compute than the inverse)::

    # Sgr to Galactic coordinates
    @frame_transform_graph.transform(coord.FunctionTransform, Sagittarius, coord.Galactic)
    def sgr_to_galactic(sgr_coord, gal_frame):
        """ Compute the transformation from heliocentric Sgr coordinates to
            spherical Galactic.
        """
        L = np.atleast_1d(sgr_coord.Lambda.radian)
        B = np.atleast_1d(sgr_coord.Beta.radian)

        Xs = cos(B)*cos(L)
        Ys = cos(B)*sin(L)
        Zs = sin(B)
        Zs = -Zs

        X, Y, Z = sgr_matrix.T.dot(np.array([Xs, Ys, Zs]))

        l = np.arctan2(Y,X)*u.radian
        b = np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z))*u.radian

        if l<0:
            l += 2*np.pi*u.radian

        return coord.Galactic(l=l, b=b, distance=sgr_coord.distance)

Now that we've registered these transformations between ``Sagittarius``
and `~astropy.coordinates.Galactic`, we can transform between *any*
coordinate system and ``Sagittarius`` (as long as the other system has a
path to transform to `~astropy.coordinates.Galactic`). For example, to
transform from ICRS coordinates to ``Sagittarius``, we simply::

    >>> import astropy.units as u
    >>> import astropy.coordinates as coord
    >>> icrs = coord.ICRS(280.161732*u.degree, 11.91934*u.degree)
    >>> icrs.transform_to(Sagittarius)  # doctest: +SKIP
    <Sagittarius Coordinate Lambda=346.81827... deg, Beta=-39.28367... deg>

The complete code for the above example is included below for reference.

See Also
========

* Majewski et al. 2003, "A Two Micron All Sky Survey View of the Sagittarius
  Dwarf Galaxy. I. Morphology of the Sagittarius Core and Tidal Arms",
  http://arxiv.org/abs/astro-ph/0304198
* Law & Majewski 2010, "The Sagittarius Dwarf Galaxy: A Model for Evolution in a
  Triaxial Milky Way Halo", http://arxiv.org/abs/1003.1132
* David Law's Sgr info page http://www.astro.virginia.edu/~srm4n/Sgr/

Complete Code for Example
=========================
.. literalinclude:: sgr-example.py
    :linenos:
    :language: python
