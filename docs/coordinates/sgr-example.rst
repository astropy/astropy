.. _sgr-example:

Example: Defining a new spherical coordinate system
---------------------------------------------------

This document describes in detail how to subclass and define a custom spherical
coordinate system, as mentioned in :doc:`designing` and the docstring for
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase`. In this example,
we will define a coordinate system defined by the plane of orbit of the
Sagittarius Dwarf Galaxy (hereafter Sgr; as defined in Majewski et al. 2003).
The Sgr coordinate system is often referred to in terms of two angular
coordinates, :math:`\Lambda,B`.

We'll start by using a template to generate a descriptive docstring for our
subclass. Then we'll add some glue to the default initializer to recognize that
we want to call the coordinates `Lambda` and `Beta`. Then we have to define the
transformation from this coordinate system to some other built-in system, in
this case we will use spherical Galactic Coordinates as defined in
`~astropy.coordinates.builtin_systems.Galactic`.

The first step is to create a new class for our coordinates, let's call it
`SgrCoordinates` and make it a subclass of
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase`::

    import astropy.coordinates as coord

    class SgrCoordinates(coord.SphericalCoordinatesBase):
        """ A spherical coordinate system defined by the orbit of the
            Sagittarius dwarf galaxy, as described in
                http://adsabs.harvard.edu/abs/2003ApJ...599.1082M
            and further explained in
                http://www.astro.virginia.edu/~srm4n/Sgr/.

    """

We can use some convenience variables defined at the superclass level to
generate a descriptive docstring for our subclass::

    __doc__ = __doc__.format(params=coord.SphericalCoordinatesBase. \
            _init_docstring_param_templ.format(lonnm='Lambda', latnm='Beta'))

This may look scary, but let's break down what it is doing. `__doc__` is a class
attribute defined in the superclass,
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase`, and is just a
string with a placeholder for a docstring specific to this subclass. We specify
that extra, descriptive text by using `format()` with the `params` keyword.
`params` expects a string describing the variable names of our current
coordinate system. We can generate this from a template as well by formatting
another class attribute and telling it what we want to name our coordinates (in
this case, Lambda and Beta). The template string is defined as a class attribute
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase._init_docstring_param_templ`,
and expects the name of the longitude-equivalent (Lambda) and
latitude-equivalent (Beta) coordinates::

    coord.SphericalCoordinatesBase._init_docstring_param_templ
                                  .format(lonnm='Lambda', latnm='Beta')

If this looks strange or complicated, just remember that we're just doing this
to help the user (this will provide them with a nice docstring if they look at
the help for this class). **The bottom line:** if you want to define a custom
spherical coordinate system with two angular coordinates, copy and paste the
below line into the class definition and change `'Lambda'` to the name of your
longitude coordinate, and `'Beta'` to the name of your latitude coordinate.::

    __doc__ = __doc__.format(params=coord.SphericalCoordinatesBase. \
            _init_docstring_param_templ.format(lonnm='Lambda', latnm='Beta'))

Next we can define our subclass' initializer, `__init__()`. We start by calling
`__init__()` on the superclass. Then, we add a catch to see if the user passed
in another `~astropy.coordinates.coordsystems.SphericalCoordinatesBase` object
in. If so, we can just transform those coordinates to `SgrCoordinates` and copy
the attributes directly. Otherwise, we have to tell our subclass how to parse
the arguments and keyword arguments that the user passed in. We do this with the
`_initialize_latlon()` classmethod. This convenience function accepts the name
of the longitude coordinate, the name of the latitude coordinate, and the
arguments and keyword arguments passed to this initializer. This function
will handle parsing the args and kwargs to validate and sanitize the input from
the user.::

    def __init__(self, *args, **kwargs):
        super(SgrCoordinates, self).__init__()

        if len(args) == 1 and len(kwargs) == 0 and
            isinstance(args[0], coord.SphericalCoordinatesBase):

            newcoord = args[0].transform_to(self.__class__)
            self._lonangle = newcoord._lonangle
            self._latangle = newcoord._latangle
            self._distance = newcoord._distance
        else:
            super(SgrCoordinates, self).
                _initialize_latlon('Lambda', 'Beta', args, kwargs)

Next we have to tell the class what to use for the longitude and latitude when
doing coordinate transformations with other coordinate systems. We do that by
defining properties of the class -- in this case `Lambda` and `Beta` -- that
reference the internal names for the longitude and latitude::

    @property
    def Lambda(self):
        return self._lonangle

    @property
    def Beta(self):
        return self._latangle

Finally, we customize the default string formatting by specifying the string
names of the longitude and latitude, as well as the default numerical
representation (degree minute second for longitude and degree
minute second for latitude)::

    # strings used for making __repr__ work
    _repr_lon_name = 'Lambda'
    _repr_lat_name = 'Beta'

    # Default format for to_string
    _default_string_style = 'dmsdms'

Now our coordinate system is set up! You can now create `SgrCoordinates` objects
by passing any valid specifiers accepted by
`~astropy.coordinates.coordsystems.SphericalCoordinatesBase`, but you won't be
able to transform to and from any other coordinate systems because we haven't
yet defined how to do that.

We will define the coordinate transformation as a function that accepts a
`~astropy.coordinates.builtin_systems.Galactic` object and returns an
`SgrCoordinates` object. We could alternatively define the transformation by
specifying a transformation matrix (
`~astropy.coordinates.transformations.static_transform_matrix`), but in this
case the transformation is from a left-handed coordinate system to a right-
handed system, so we need more control over the transformation. We'll start by
constructing the rotation matrix, using the helper function
`astropy.coordinates.angles.rotation_matrix` ::

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
transformation from Sgr to Galactic as well as the inverse from Galactic to Sgr.
Now we can define our first transformation function::

    # Galactic to Sgr coordinates
    @transformations.transform_function(coord.Galactic, SgrCoordinates)
    def galactic_to_sgr(galactic_coord):
        """ Compute the transformation from Galactic spherical to
            heliocentric Sgr coordinates.
        """

        l = np.atleast_1d(galactic_coord.l.radian)
        b = np.atleast_1d(galactic_coord.b.radian)

        X = cos(b)*cos(l)
        Y = cos(b)*sin(l)
        Z = sin(b)

        # Calculate X,Y,Z,distance in the Sgr system
        Xs, Ys, Zs = sgr_matrix.dot(np.array([X, Y, Z]))
        Zs = -Zs

        # Calculate the angular coordinates lambda,beta
        Lambda = np.degrees(np.arctan2(Ys,Xs))
        Lambda[Lambda < 0] = Lambda[Lambda < 0] + 360
        Beta = np.degrees(np.arcsin(Zs/np.sqrt(Xs*Xs+Ys*Ys+Zs*Zs)))

        return SgrCoordinates(Lambda, Beta, distance=galactic_coord.distance,
                              unit=(u.degree, u.degree))

The decorator
`@transformations.transform_function(coord.Galactic, SgrCoordinates)`
registers this function with the `Coordinates` subpackage as a transformation.
Inside the function, we simply follow the same procedure as detailed by David
Law's
`transformation code <http://www.astro.virginia.edu/~srm4n/Sgr/code.html>`. Note
that in this case, both coordinate systems are heliocentric, so we can simply
copy any distance from the
`~astropy.coordinates.builtin_systems.Galactic` object.

We then register the inverse transformation by using the Transpose of the
rotation matrix::

    @transformations.transform_function(SgrCoordinates, coord.Galactic)
    def sgr_to_galactic(sgr_coord):
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

        l = np.degrees(np.arctan2(Y,X))
        b = np.degrees(np.arcsin(Z/np.sqrt(X*X+Y*Y+Z*Z)))

        if l<0:
            l += 360

        return coord.Galactic(l, b, distance=sgr_coord.distance,
                              unit=(u.degree, u.degree))

Now that we've registered these transformations between `SgrCoordinates` and
`~astropy.coordinates.builtin_systems.Galactic`, we can transform
between *any* coordinate system and `SgrCoordinates` (as long as the other
system has a path to transform to
`~astropy.coordinates.builtin_systems.Galactic`). For example, to
transform from ICRS coordinates to `SgrCoordinates`, we simply::

    >>> from astropy import units as u
    >>> import astropy.coordinates as coord
    >>> icrs = coord.ICRS(280.161732, 11.91934, unit=(u.degree,u.degree))
    >>> icrs.transform_to(SgrCoordinates)  # doctest: +SKIP
    <SgrCoordinates Lambda=346.81827 deg, Beta=-39.28367 deg>

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
