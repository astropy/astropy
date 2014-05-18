.. _astropy-coordinates-representations:

Using and Designing Coordinate Representations
----------------------------------------------

As described in :ref:`astropy-coordinates-overview`, the actual coordinate
data in `~astropy.coordinates` inside frames is represented via
"Representation classes". These can be used to store 3-d coordinates in
various representations, such as cartesian, spherical polar, cylindrical, and
so on. The built-in representation classes are:

* `~astropy.coordinates.CartesianRepresentation`: cartesian
  coordinates ``x``, ``y``, and ``z``
* `~astropy.coordinates.SphericalRepresentation`: spherical
  polar coordinates represented by a longitude (``lon``), a latitude
  (``lat``), and a distance (``distance``). The latitude is a value ranging
  from -90 to 90 degrees.
* `~astropy.coordinates.UnitSphericalRepresentation`:
  spherical polar coordinates on a unit sphere, represented by a longitude
  (``lon``) and latitude (``lat``)
* `~astropy.coordinates.PhysicsSphericalRepresentation`:
  spherical polar coordinates, represented by an inclination (``theta``) and
  azimuthal angle (``phi``), and radius ``r``. The inclination goes from 0 to
  180 degrees, and is related to the latitude in the
  `~astropy.coordinates.SphericalRepresentation` by
  ``theta = 90 deg - lat``.
* `~astropy.coordinates.CylindricalRepresentation`:
  cylindrical polar coordinates, represented by a cylindrical radius
  (``rho``), azimuthal angle (``phi``), and height (``z``).

Instantiating and converting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Representation classes should be instantiated with `~astropy.units.Quantity`
objects:

    >>> from astropy import units as u
    >>> from astropy.coordinates.representation import CartesianRepresentation
    >>> car = CartesianRepresentation(3 * u.kpc, 5 * u.kpc, 4 * u.kpc)
    >>> car
    <CartesianRepresentation x=3.0 kpc, y=5.0 kpc, z=4.0 kpc>

Representations can be converted to other representations using the
``represent_as`` method::

    >>> from astropy.coordinates.representation import SphericalRepresentation, CylindricalRepresentation
    >>> sph = car.represent_as(SphericalRepresentation)
    >>> sph
    <SphericalRepresentation lon=1.03037682652 rad, lat=0.601264216679 rad, distance=7.07106781187 kpc>
    >>> cyl = car.represent_as(CylindricalRepresentation)
    >>> cyl
    <CylindricalRepresentation rho=5.83095189485 kpc, phi=1.03037682652 rad, z=4.0 kpc>

All representations can be converted to each other without loss of
information, with the exception of
`~astropy.coordinates.UnitSphericalRepresentation`. This class
is used to store the longitude and latitude of points but does not contain
any distance to the points, and assumes that they are located on a unit and
dimensionless sphere::

    >>> from astropy.coordinates.representation import UnitSphericalRepresentation
    >>> sph_unit = car.represent_as(UnitSphericalRepresentation)
    >>> sph_unit
    <UnitSphericalRepresentation lon=1.03037682652 rad, lat=0.601264216679 rad>

Converting back to cartesian, the absolute scaling information has been
removed, and the points are still located on a unit sphere:

    >>> sph_unit = c.represent_as(UnitSphericalRepresentation)
    >>> sph_unit.represent_as(CartesianRepresentation)
    <CartesianRepresentation x=0.424264068712 , y=0.707106781187 , z=0.565685424949 >

Array values
^^^^^^^^^^^^

Array `~astropy.units.Quantity` objects can also be passed to
representations::

  >>> x = np.random.random(100)
  >>> y = np.random.random(100)
  >>> z = np.random.random(100)
  >>> car_array = CartesianRepresentation(x * u.m, y * u.m, z * u.m)
  >>> car_array
  <CartesianRepresentation (x, y, z) in m
      [(0.7093178096722875, 0.7788407268253749, 0.38429719423939623),
       (0.8434446789489138, 0.4543837868963039, 0.9579319360271565),
       ...
       (0.017984107510863323, 0.8587892718748261, 0.49163383463637955),
       (0.020733247493708196, 0.33551661975140423, 0.279998327062161)]>

Creating your own representations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create your own representation class, your class should:

* Your class should inherit from the ``BaseRepresentation`` class.
* Define ``__init__``
* Define a ``from_cartesian`` class method which should take a
  `~astropy.coordinates.CartesianRepresentation` object and
  return an instance of your class.
* Define a ``to_cartesian`` method which should return a
  `~astropy.coordinates.CartesianRepresentation` object
* Define a ``components`` property that returns a list of the names of the
  coordinate components (such as ``x``, ``lon``, and so on).

In pseudo-code, this means that your class will look like::

    class MyRepresentation(BaseRepresentation):

        def __init__(self, ...):
            ...

        @classmethod
        def from_cartesian(self, cartesian):
            ...
            return MyRepresentation(...)

        def to_cartesian(self):
            ...
            return CartesianRepresentation(...)

        @property
        def components(self):
            return [...]

Once you do this, you will then automatically be able to call
``represent_as`` to convert other representations to/from your representation
class.
