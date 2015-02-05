.. include:: references.txt

.. _astropy-coordinates-representations:

Using and Designing Coordinate Representations
----------------------------------------------

As described in the :ref:`astropy-coordinates-overview`, the actual coordinate
data in `astropy.coordinates` frames is represented via
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

.. Note::
   For information about using and changing the representation of
   `~astropy.coordinates.SkyCoord` objects, see the
   :ref:`astropy-skycoord-representations` section.

Instantiating and converting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Representation classes should be instantiated with `~astropy.units.Quantity`
objects::

    >>> from astropy import units as u
    >>> from astropy.coordinates.representation import CartesianRepresentation
    >>> car = CartesianRepresentation(3 * u.kpc, 5 * u.kpc, 4 * u.kpc)
    >>> car
    <CartesianRepresentation (x, y, z) in kpc
        (3.0, 5.0, 4.0)>

Representations can be converted to other representations using the
``represent_as`` method::

    >>> from astropy.coordinates.representation import SphericalRepresentation, CylindricalRepresentation
    >>> sph = car.represent_as(SphericalRepresentation)
    >>> sph  # doctest: +FLOAT_CMP
    <SphericalRepresentation (lon, lat, distance) in (rad, rad, kpc)
        (1.03037682652, 0.601264216679, 7.07106781187)>
    >>> cyl = car.represent_as(CylindricalRepresentation)
    >>> cyl  # doctest: +FLOAT_CMP
    <CylindricalRepresentation (rho, phi, z) in (kpc, rad, kpc)
        (5.83095189485, 1.03037682652, 4.0)>

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
        (1.03037682652, 0.601264216679)>

Converting back to cartesian, the absolute scaling information has been
removed, and the points are still located on a unit sphere:

    >>> sph_unit = car.represent_as(UnitSphericalRepresentation)
    >>> sph_unit.represent_as(CartesianRepresentation) # doctest: +FLOAT_CMP
    <CartesianRepresentation (x, y, z) in
        (0.424264068712, 0.707106781187, 0.565685424949)>

Array values
^^^^^^^^^^^^

Array `~astropy.units.Quantity` objects can also be passed to
representations::

  >>> import numpy as np
  >>> x = np.random.random(100)
  >>> y = np.random.random(100)
  >>> z = np.random.random(100)
  >>> car_array = CartesianRepresentation(x * u.m, y * u.m, z * u.m)
  >>> car_array  # doctest: +SKIP
  <CartesianRepresentation (x, y, z) in m
      [(0.7093..., 0.7788..., 0.3842...),
       (0.8434..., 0.4543..., 0.9579...),
       ...
       (0.0179..., 0.8587..., 0.4916...),
       (0.0207..., 0.3355..., 0.2799...)]>

.. _astropy-coordinates-create-repr:

Creating your own representations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create your own representation class, your class must inherit from the
``BaseRepresentation`` class.  In addition the following must be defined:

* ``__init__`` method:

  Has a signature like ``__init__(self, comp1, comp2, comp3, copy=True)``
  for inputting the representation component values.

* ``from_cartesian`` class method:

  Takes a `~astropy.coordinates.CartesianRepresentation` object and
  returns an instance of your class.

* ``to_cartesian`` method:

  Returns a `~astropy.coordinates.CartesianRepresentation` object.

* ``components`` property:

  Returns a tuple of the names of the coordinate components (such as ``x``,
  ``lon``, and so on).

* ``attr_classes`` class attribute (``OrderedDict``):

  Defines the initializer class for each component.In most cases this
  class should be derived from `~astropy.units.Quantity`. In particular
  these class initializers must take the value as the first argument and
  accept a ``unit`` keyword which takes a `~astropy.units.Unit`
  initializer or ``None`` to indicate no unit. Also not that the keys of
  this dictionary are treated as the names of the components for this
  representation, with the default ordered given in the order they
  appear as keys.

* ``recommended_units`` dictionary (optional):

  Maps component names to the recommended unit to convert the values of
  that component to.  Can be ``None`` (or missing) to indicate there is
  no preferred unit.  If this dictionary is not defined, no conversion
  of components to particular units will occur.

In pseudo-code, this means that your class will look like::

    class MyRepresentation(BaseRepresentation):

        attr_classes = OrderedDict([('comp1', ComponentClass1),
                                     ('comp2', ComponentClass2),
                                     ('comp3', ComponentClass3)])

        # recommended_units is optional
        recommended_units = {'comp1': u.unit1, 'comp2': u.unit2, 'comp3': u.unit3}

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
            return 'comp1', 'comp2', 'comp3'

Once you do this, you will then automatically be able to call
``represent_as`` to convert other representations to/from your representation
class.  Your representation will also be available for use in |skycoord|
and all frame classes.
