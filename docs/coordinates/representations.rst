.. include:: references.txt

.. _astropy-coordinates-representations:

Using and Designing Coordinate Representations
**********************************************

Points in a 3D vector space can be represented in different ways, such as
Cartesian, spherical polar, cylindrical, and so on. These underlie the way
coordinate data in `astropy.coordinates` is represented, as described in the
:ref:`astropy-coordinates-overview`. Below, we describe how you can use them on
their own as a way to convert between different representations, including
ones not built-in, and to do simple vector arithmetic.

The built-in representation classes are:

* `~astropy.coordinates.CartesianRepresentation`: Cartesian
  coordinates ``x``, ``y``, and ``z``.
* `~astropy.coordinates.SphericalRepresentation`: spherical
  polar coordinates represented by a longitude (``lon``), a latitude
  (``lat``), and a distance (``distance``). The latitude is a value ranging
  from -90 to 90 degrees.
* `~astropy.coordinates.UnitSphericalRepresentation`:
  spherical polar coordinates on a unit sphere, represented by a longitude
  (``lon``) and latitude (``lat``).
* `~astropy.coordinates.PhysicsSphericalRepresentation`:
  spherical polar coordinates, represented by an inclination (``theta``) and
  azimuthal angle (``phi``), and radius ``r``. The inclination goes from 0 to
  180 degrees, and is related to the latitude in the
  `~astropy.coordinates.SphericalRepresentation` by
  ``theta = 90 deg - lat``.
* `~astropy.coordinates.CylindricalRepresentation`:
  cylindrical polar coordinates, represented by a cylindrical radius
  (``rho``), azimuthal angle (``phi``), and height (``z``).

.. Note::
   For information about using and changing the representation of
   `~astropy.coordinates.SkyCoord` objects, see the
   :ref:`astropy-skycoord-representations` section.

Instantiating and Converting
============================

Representation classes are instantiated with `~astropy.units.Quantity`
objects::

    >>> from astropy import units as u
    >>> from astropy.coordinates.representation import CartesianRepresentation
    >>> car = CartesianRepresentation(3 * u.kpc, 5 * u.kpc, 4 * u.kpc)
    >>> car  # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) in kpc
        (3., 5., 4.)>

Array `~astropy.units.Quantity` objects can also be passed to
representations. They will have the expected shape, which can be changed using
methods with the same names as those for `~numpy.ndarray`, such as ``reshape``,
``ravel``, etc.::

  >>> x = u.Quantity([[1., 0., 0.], [3., 5., 3.]], u.m)
  >>> y = u.Quantity([[0., 2., 0.], [4., 0., -4.]], u.m)
  >>> z = u.Quantity([[0., 0., 3.], [0., 12., -12.]], u.m)
  >>> car_array = CartesianRepresentation(x, y, z)
  >>> car_array  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [[(1.,  0.,   0.), (0.,  2.,   0.), (0.,  0.,   3.)],
       [(3.,  4.,   0.), (5.,  0.,  12.), (3., -4., -12.)]]>
  >>> car_array.shape
  (2, 3)
  >>> car_array.ravel()  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [(1.,  0.,   0.), (0.,  2.,   0.), (0.,  0.,   3.), (3.,  4.,   0.),
       (5.,  0.,  12.), (3., -4., -12.)]>

Representations can be converted to other representations using the
``represent_as`` method::

    >>> from astropy.coordinates.representation import SphericalRepresentation, CylindricalRepresentation
    >>> sph = car.represent_as(SphericalRepresentation)
    >>> sph  # doctest: +FLOAT_CMP
    <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
        (1.03037683, 0.60126422, 7.07106781)>
    >>> cyl = car.represent_as(CylindricalRepresentation)
    >>> cyl  # doctest: +FLOAT_CMP
    <CylindricalRepresentation (rho, phi, z) in (kpc, rad, kpc)
        (5.83095189, 1.03037683, 4.)>

All representations can be converted to each other without loss of
information, with the exception of
`~astropy.coordinates.UnitSphericalRepresentation`. This class
is used to store the longitude and latitude of points but does not contain
any distance to the points, and assumes that they are located on a unit and
dimensionless sphere::

    >>> from astropy.coordinates.representation import UnitSphericalRepresentation
    >>> sph_unit = car.represent_as(UnitSphericalRepresentation)
    >>> sph_unit  # doctest: +FLOAT_CMP
    <UnitSphericalRepresentation (lon, lat) in rad
        (1.03037683, 0.60126422)>

Converting back to Cartesian, the absolute scaling information has been
removed, and the points are still located on a unit sphere::

    >>> sph_unit = car.represent_as(UnitSphericalRepresentation)
    >>> sph_unit.represent_as(CartesianRepresentation)  # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) [dimensionless]
        (0.42426407, 0.70710678, 0.56568542)>


Array Values and NumPy Array Method Analogs
===========================================

Array `~astropy.units.Quantity` objects can also be passed to representations,
and such representations can be sliced, reshaped, etc., using the same
methods as are available to `~numpy.ndarray`::

  >>> import numpy as np
  >>> x = np.linspace(0., 5., 6)
  >>> y = np.linspace(10., 15., 6)
  >>> z = np.linspace(20., 25., 6)
  >>> car_array = CartesianRepresentation(x * u.m, y * u.m, z * u.m)
  >>> car_array  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [(0., 10., 20.), (1., 11., 21.), (2., 12., 22.),
       (3., 13., 23.), (4., 14., 24.), (5., 15., 25.)]>
  >>> car_array[2]  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      (2., 12., 22.)>
  >>> car_array.reshape(3, 2)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [[(0., 10., 20.), (1., 11., 21.)],
       [(2., 12., 22.), (3., 13., 23.)],
       [(4., 14., 24.), (5., 15., 25.)]]>


.. _astropy-coordinates-representations-arithmetic:

Vector Arithmetic
=================

Representations support basic vector arithmetic such as taking the norm,
multiplying with and dividing by quantities, and taking dot and cross products,
as well as adding, subtracting, summing and taking averages of representations,
and multiplying with matrices.

.. Note:: All arithmetic except the matrix multiplication works with
   non-Cartesian representations as well. For taking the norm, multiplication,
   and division, this uses just the non-angular components, while for the other
   operations the representation is converted to Cartesian internally before
   the operation is done, and the result is converted back to the original
   representation. Hence, for optimal speed it may be best to work using
   Cartesian representations.

To see how the operations work, consider the following examples::

  >>> car_array = CartesianRepresentation([[1., 0., 0.], [3., 5.,  3.]] * u.m,
  ...                                     [[0., 2., 0.], [4., 0., -4.]] * u.m,
  ...                                     [[0., 0., 3.], [0.,12.,-12.]] * u.m)
  >>> car_array  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [[(1.,  0.,  0.), (0.,  2.,   0.), (0.,  0.,   3.)],
       [(3.,  4.,  0.), (5.,  0.,  12.), (3., -4., -12.)]]>
  >>> car_array.norm()  # doctest: +FLOAT_CMP
  <Quantity [[ 1.,  2.,  3.],
             [ 5., 13., 13.]] m>
  >>> car_array / car_array.norm()  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) [dimensionless]
      [[(1.        ,  0.        ,  0.        ),
        (0.        ,  1.        ,  0.        ),
        (0.        ,  0.        ,  1.        )],
       [(0.6       ,  0.8       ,  0.        ),
        (0.38461538,  0.        ,  0.92307692),
        (0.23076923, -0.30769231, -0.92307692)]]>
  >>> (car_array[1] - car_array[0]) / (10. * u.s)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m / s
      [(0.2,  0.4,  0. ), (0.5, -0.2,  1.2), (0.3, -0.4, -1.5)]>
  >>> car_array.sum()  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      (12.,  2.,  3.)>
  >>> car_array.mean(axis=0)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [(2. ,  2.,  0. ), (2.5,  1.,  6. ), (1.5, -2., -4.5)]>

  >>> unit_x = UnitSphericalRepresentation(0.*u.deg, 0.*u.deg)
  >>> unit_y = UnitSphericalRepresentation(90.*u.deg, 0.*u.deg)
  >>> unit_z = UnitSphericalRepresentation(0.*u.deg, 90.*u.deg)
  >>> car_array.dot(unit_x)  # doctest: +FLOAT_CMP
  <Quantity [[1., 0., 0.],
             [3., 5., 3.]] m>
  >>> car_array.dot(unit_y)  # doctest: +FLOAT_CMP
  <Quantity [[ 6.12323400e-17,  2.00000000e+00,  0.00000000e+00],
             [ 4.00000000e+00,  3.06161700e-16, -4.00000000e+00]] m>
  >>> car_array.dot(unit_z)  # doctest: +FLOAT_CMP
  <Quantity [[ 6.12323400e-17,  0.00000000e+00,  3.00000000e+00],
             [ 1.83697020e-16,  1.20000000e+01, -1.20000000e+01]] m>
  >>> car_array.cross(unit_x)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [[(0.,  0.,  0.), (0.,   0., -2.), (0.,   3.,  0.)],
       [(0.,  0., -4.), (0.,  12.,  0.), (0., -12.,  4.)]]>

  >>> from astropy.coordinates.matrix_utilities import rotation_matrix
  >>> rotation = rotation_matrix(90 * u.deg, axis='z')
  >>> rotation  # doctest: +FLOAT_CMP
  array([[ 6.12323400e-17,  1.00000000e+00,  0.00000000e+00],
         [-1.00000000e+00,  6.12323400e-17,  0.00000000e+00],
         [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
  >>> car_array.transform(rotation)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in m
      [[( 6.12323400e-17, -1.00000000e+00,   0.),
        ( 2.00000000e+00,  1.22464680e-16,   0.),
        ( 0.00000000e+00,  0.00000000e+00,   3.)],
       [( 4.00000000e+00, -3.00000000e+00,   0.),
        ( 3.06161700e-16, -5.00000000e+00,  12.),
        (-4.00000000e+00, -3.00000000e+00, -12.)]]>

.. _astropy-coordinates-differentials:

Differentials and Derivatives of Representations
================================================

In addition to positions in 3D space, coordinates also deal with proper motions
and radial velocities, which require a way to represent differentials of
coordinates (i.e., finite realizations) of derivatives. To support this, the
representations all have corresponding ``Differential`` classes, which can hold
offsets or derivatives in terms of the components of the representation class.
Adding such an offset to a representation means the offset is taken in the
direction of the corresponding coordinate. (Although for any representation
other than Cartesian, this is only defined relative to a specific location, as
the unit vectors are not invariant.)

To see how this works, consider the following::

  >>> from astropy.coordinates import SphericalRepresentation, SphericalDifferential
  >>> sph_coo = SphericalRepresentation(lon=0.*u.deg, lat=0.*u.deg,
  ...                                   distance=1.*u.kpc)
  >>> sph_derivative = SphericalDifferential(d_lon=1.*u.arcsec/u.yr,
  ...                                        d_lat=0.*u.arcsec/u.yr,
  ...                                        d_distance=0.*u.km/u.s)
  >>> sph_derivative.to_cartesian(base=sph_coo)  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in arcsec kpc / (rad yr)
      (0., 1., 0.)>

Note how the conversion to Cartesian can only be done using a ``base``, since
otherwise the code cannot know what direction an increase in longitude
corresponds to. For ``lon=0``, this is in the ``y`` direction. Now, to get
the coordinates at two later times::

  >>> sph_coo + sph_derivative * [1., 3600*180/np.pi] * u.yr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      [(4.84813681e-06, 0., 1.        ), (7.85398163e-01, 0., 1.41421356)]>

The above shows how addition is not to longitude itself, but in the direction
of increasing longitude: for the large shift, by the equivalent of one radian,
the distance has increased as well (after all, a source will likely not move
along a curve on the sky!). This also means that the order of operations is
important::

  >>> big_offset = SphericalDifferential(1.*u.radian, 0.*u.radian, 0.*u.kpc)
  >>> sph_coo + big_offset + big_offset  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (1.57079633, 0., 2.)>
  >>> sph_coo + (big_offset + big_offset)  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (1.10714872, 0., 2.23606798)>

Often, you may have just a proper motion or a radial velocity, but not both::

  >>> from astropy.coordinates import UnitSphericalDifferential, RadialDifferential
  >>> radvel = RadialDifferential(1000*u.km/u.s)
  >>> sph_coo + radvel * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0., 0., 2.02271217)>
  >>> pm = UnitSphericalDifferential(1.*u.mas/u.yr, 0.*u.mas/u.yr)
  >>> sph_coo + pm * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0.0048481, 0., 1.00001175)>
  >>> pm + radvel  # doctest: +FLOAT_CMP
  <SphericalDifferential (d_lon, d_lat, d_distance) in (mas / yr, mas / yr, km / s)
      (1., 0., 1000.)>
  >>> sph_coo + (pm + radvel) * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0.00239684, 0., 2.02271798)>

Note in the above that the proper motion is defined strictly as a change in
longitude (i.e., it does not include a ``cos(latitude)`` term). There are
special classes where this term is included::

  >>> from astropy.coordinates import UnitSphericalCosLatDifferential
  >>> sph_lat60 = SphericalRepresentation(lon=0.*u.deg, lat=60.*u.deg,
  ...                                     distance=1.*u.kpc)
  >>> pm = UnitSphericalDifferential(1.*u.mas/u.yr, 0.*u.mas/u.yr)
  >>> pm  # doctest: +FLOAT_CMP
  <UnitSphericalDifferential (d_lon, d_lat) in mas / yr
      (1., 0.)>
  >>> pm_coslat = UnitSphericalCosLatDifferential(1.*u.mas/u.yr, 0.*u.mas/u.yr)
  >>> pm_coslat  # doctest: +FLOAT_CMP
  <UnitSphericalCosLatDifferential (d_lon_coslat, d_lat) in mas / yr
      (1., 0.)>
  >>> sph_lat60 + pm * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0.0048481, 1.04719246, 1.00000294)>
  >>> sph_lat60 + pm_coslat * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0.00969597, 1.0471772, 1.00001175)>

Close inspections shows that indeed the changes are as expected. The systems
with and without ``cos(latitude)`` can be converted to each other, provided you
supply the ``base`` (representation)::

  >>> usph_lat60 = sph_lat60.represent_as(UnitSphericalRepresentation)
  >>> pm_coslat2 = pm.represent_as(UnitSphericalCosLatDifferential,
  ...                              base=usph_lat60)
  >>> pm_coslat2  # doctest: +FLOAT_CMP
  <UnitSphericalCosLatDifferential (d_lon_coslat, d_lat) in mas / yr
      (0.5, 0.)>
  >>> sph_lat60 + pm_coslat2 * 1. * u.Myr  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
      (0.0048481, 1.04719246, 1.00000294)>

.. Note:: At present, the differential classes are generally meant to work with
   first derivatives, but they do not check the units of the inputs to enforce
   this. Passing in second derivatives (e.g., acceleration values with
   acceleration units) will succeed, but any transformations that occur through
   re-representation of the differential will not necessarily be correct.

Attaching ``Differential`` Objects to ``Representation`` Objects
================================================================

.. warning::

    The API for this functionality may change in future versions and should be
    viewed as provisional!

``Differential`` objects can be attached to ``Representation`` objects as a way
to encapsulate related information into a single object. ``Differential``
objects can be passed in to the initializer of any of the built-in
``Representation`` classes. For example, to store a single velocity
differential with a position::

  >>> from astropy.coordinates import representation as r
  >>> dif = r.SphericalDifferential(d_lon=1 * u.mas/u.yr,
  ...                               d_lat=2 * u.mas/u.yr,
  ...                               d_distance=3 * u.km/u.s)
  >>> rep = r.SphericalRepresentation(lon=0.*u.deg, lat=0.*u.deg,
  ...                                 distance=1.*u.kpc,
  ...                                 differentials=dif)
  >>> rep  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (deg, deg, kpc)
      (0., 0., 1.)
   (has differentials w.r.t.: 's')>
  >>> rep.differentials  # doctest: +FLOAT_CMP
  {'s': <SphericalDifferential (d_lon, d_lat, d_distance) in (mas / yr, mas / yr, km / s)
       (1., 2., 3.)>}

The ``Differential`` objects are stored as a Python dictionary on the
``Representation`` object with keys equal to the (string) unit with which the
differential derivatives are taken (converted to SI). For example, in this case
the key is ``'s'`` (second) because the ``Differential`` units are velocities, a
time derivative. Passing a single differential to the ``Representation``
initializer will automatically generate the necessary key and store it in the
differentials dictionary, but a dictionary is required to specify multiple
differentials::

  >>> dif2 = r.SphericalDifferential(d_lon=4 * u.mas/u.yr**2,
  ...                                d_lat=5 * u.mas/u.yr**2,
  ...                                d_distance=6 * u.km/u.s**2)
  >>> rep = r.SphericalRepresentation(lon=0.*u.deg, lat=0.*u.deg,
  ...                                 distance=1.*u.kpc,
  ...                                 differentials={'s': dif, 's2': dif2})
  >>> rep.differentials['s']  # doctest: +FLOAT_CMP
  <SphericalDifferential (d_lon, d_lat, d_distance) in (mas / yr, mas / yr, km / s)
      (1., 2., 3.)>
  >>> rep.differentials['s2']  # doctest: +FLOAT_CMP
  <SphericalDifferential (d_lon, d_lat, d_distance) in (mas / yr2, mas / yr2, km / s2)
      (4., 5., 6.)>

``Differential`` objects can also be attached to a ``Representation`` after
creation::

  >>> rep = r.CartesianRepresentation(x=1 * u.kpc, y=2 * u.kpc, z=3 * u.kpc)
  >>> dif = r.CartesianDifferential(*[1, 2, 3] * u.km/u.s)
  >>> rep = rep.with_differentials(dif)
  >>> rep  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in kpc
      (1., 2., 3.)
   (has differentials w.r.t.: 's')>

This works for array data as well, as long as the shape of the
``Differential`` data is the same as that of the ``Representation``::

  >>> xyz = np.arange(12).reshape(3, 4) * u.au
  >>> d_xyz = np.arange(12).reshape(3, 4) * u.km/u.s
  >>> rep = r.CartesianRepresentation(*xyz)
  >>> dif = r.CartesianDifferential(*d_xyz)
  >>> rep = rep.with_differentials(dif)
  >>> rep  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in AU
      [(0., 4.,  8.), (1., 5.,  9.), (2., 6., 10.), (3., 7., 11.)]
   (has differentials w.r.t.: 's')>

As with a ``Representation`` instance without a differential, to convert the
positional data to a new representation, use the ``.represent_as()``::

  >>> rep.represent_as(r.SphericalRepresentation)  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, AU)
      [(1.57079633, 1.10714872,  8.94427191),
       (1.37340077, 1.05532979, 10.34408043),
       (1.24904577, 1.00685369, 11.83215957),
       (1.16590454, 0.96522779, 13.37908816)]>

However, by passing just the desired representation class, only the
``Representation`` has changed, and the differentials are dropped. To
re-represent both the ``Representation`` and any ``Differential`` objects, you
must specify target classes for the ``Differential`` as well::

  >>> rep2 = rep.represent_as(r.SphericalRepresentation, r.SphericalDifferential)
  >>> rep2  # doctest: +FLOAT_CMP
  <SphericalRepresentation (lon, lat, distance) in (rad, rad, AU)
    [(1.57079633, 1.10714872,  8.94427191),
     (1.37340077, 1.05532979, 10.34408043),
     (1.24904577, 1.00685369, 11.83215957),
     (1.16590454, 0.96522779, 13.37908816)]
   (has differentials w.r.t.: 's')>
  >>> rep2.differentials['s']  # doctest: +FLOAT_CMP
  <SphericalDifferential (d_lon, d_lat, d_distance) in (km rad / (AU s), km rad / (AU s), km / s)
      [( 6.12323400e-17, 1.11022302e-16,  8.94427191),
       (-2.77555756e-17, 5.55111512e-17, 10.34408043),
       ( 0.00000000e+00, 0.00000000e+00, 11.83215957),
       ( 5.55111512e-17, 0.00000000e+00, 13.37908816)]>

Shape-changing operations (e.g., reshapes) are propagated to all
``Differential`` objects because they are guaranteed to have the same shape as
their host ``Representation`` object::

  >>> rep.shape
  (4,)
  >>> rep.differentials['s'].shape
  (4,)
  >>> new_rep = rep.reshape(2, 2)
  >>> new_rep.shape
  (2, 2)
  >>> new_rep.differentials['s'].shape
  (2, 2)

This also works for slicing::

  >>> new_rep = rep[:2]
  >>> new_rep.shape
  (2,)
  >>> new_rep.differentials['s'].shape
  (2,)

Operations on representations that return `~astropy.units.Quantity` objects (as
opposed to other ``Representation`` instances) still work, but only operate on
the positional information, for example::

  >>> rep.norm()  # doctest: +FLOAT_CMP
  <Quantity [ 8.94427191, 10.34408043, 11.83215957, 13.37908816] AU>

Operations that involve combining or scaling representations or pairs of
representation objects that contain differentials will currently fail, but
support for some operations may be added in future versions::

  >>> rep + rep
  Traceback (most recent call last):
  ...
  TypeError: Operation 'add' is not supported when differentials are attached to a CartesianRepresentation.

If you have a ``Representation`` with attached ``Differential`` objects, you
can retrieve a copy of the ``Representation`` without the ``Differential``
object and use this ``Differential``-free object for any arithmetic operation::

  >>> 15 * rep.without_differentials()  # doctest: +FLOAT_CMP
  <CartesianRepresentation (x, y, z) in AU
      [( 0.,  60., 120.), (15.,  75., 135.), (30.,  90., 150.),
       (45., 105., 165.)]>

.. _astropy-coordinates-create-repr:

Creating Your Own Representations
=================================

To create your own representation class, your class must inherit from the
`~astropy.coordinates.BaseRepresentation` class. This base has an ``__init__``
method that will put all arguments components through their initializers,
verify they can be broadcast against each other, and store the components on
``self`` as the name prefixed with '_'. Furthermore, through its metaclass it
provides default properties for the components so that they can be accessed
using ``<instance>.<component>``. For the machinery to work, the following
must be defined:

* ``attr_classes`` class attribute (``OrderedDict``):

  Defines through its keys the names of the components (as well as the default
  order), and through its values defines the class of which they should be
  instances (which should be `~astropy.units.Quantity` or a subclass, or
  anything that can initialize it).

* ``from_cartesian`` class method:

  Takes a `~astropy.coordinates.CartesianRepresentation` object and
  returns an instance of your class.

* ``to_cartesian`` method:

  Returns a `~astropy.coordinates.CartesianRepresentation` object.

* ``__init__`` method (optional):

  If you want more than the basic initialization and checks provided by the
  base representation class, or just an explicit signature, you can define your
  own ``__init__``. In general, it is recommended to stay close to the
  signature assumed by the base representation, ``__init__(self, comp1, comp2,
  comp3, copy=True)``, and use ``super`` to call the base representation
  initializer.

Once you do this, you will then automatically be able to call ``represent_as``
to convert other representations to/from your representation class. Your
representation will also be available for use in |skycoord| and all frame
classes.

A representation class may also have a ``_unit_representation`` attribute
(although it is not required). This attribute points to the appropriate
"unit" representation (i.e., a representation that is dimensionless). This is
probably only meaningful for subclasses of
`~astropy.coordinates.SphericalRepresentation`, where it is assumed that it
will be a subclass of `~astropy.coordinates.UnitSphericalRepresentation`.

Finally, if you wish to also use offsets in your coordinate system, two further
methods should be defined (please see
`~astropy.coordinates.SphericalRepresentation` for an example):

* ``unit_vectors`` method:

  Returns a ``dict`` with a
  `~astropy.coordinates.CartesianRepresentation` of unit vectors in the
  direction of each component.

* ``scale_factors`` method:

  Returns a ``dict`` with a `~astropy.units.Quantity` for each component with
  the appropriate physical scale factor for a unit change in that direction.

And furthermore you should define a ``Differential`` class based on
`~astropy.coordinates.BaseDifferential`. This class only needs to define:

* ``base_representation`` attribute:

  A link back to the representation for which this differential holds.


In pseudo-code, this means that a class will look like::

    class MyRepresentation(BaseRepresentation):

        attr_classes = OrderedDict([('comp1', ComponentClass1),
                                     ('comp2', ComponentClass2),
                                     ('comp3', ComponentClass3)])

	# __init__ is optional
        def __init__(self, comp1, comp2, comp3, copy=True):
            super().__init__(comp1, comp2, comp3, copy=copy)
            ...

        @classmethod
        def from_cartesian(self, cartesian):
            ...
            return MyRepresentation(...)

        def to_cartesian(self):
            ...
            return CartesianRepresentation(...)

	# if differential motion is needed
	def unit_vectors(self):
	    ...
	    return {'comp1': CartesianRepresentation(...),
	            'comp2': CartesianRepresentation(...),
		    'comp3': CartesianRepresentation(...)}

        def scale_factors(self):
	    ...
	    return {'comp1': ...,
	            'comp2': ...,
		    'comp3': ...}

    class MyDifferential(BaseDifferential):
        base_representation = MyRepresentation
