.. _physical_types:

Physical Types
**************

A physical type corresponds to physical quantities with dimensionally
compatible units. For example, the physical type *mass* corresponds to
physical quantities with units that can be converted to kilograms.
Physical types are represented as instances of the |PhysicalType| class.

For a list of physical types, see `astropy.units.physical`.

Accessing Physical Types
========================

.. EXAMPLE START: Accessing Physical Types

Using :func:`~astropy.units.get_physical_type` lets us acquire |PhysicalType|
instances from strings with a name of a physical type, units, |Quantity|
instances, objects that can become quantities (e.g., numbers), and
|PhysicalType| instances.

  >>> import astropy.units as u
  >>> u.get_physical_type('speed')  # from the name of a physical type
  PhysicalType({'speed', 'velocity'})
  >>> u.get_physical_type(u.meter)  # from a unit
  PhysicalType('length')
  >>> u.get_physical_type(1 * u.barn * u.Mpc)  # from a Quantity
  PhysicalType('volume')
  >>> u.get_physical_type(42)  # from a number
  PhysicalType('dimensionless')

The physical type of a unit can be accessed via its
:attr:`~astropy.units.UnitBase.physical_type` attribute::

  >>> u.coulomb.physical_type
  PhysicalType('electrical charge')
  >>> (u.meter ** 2).physical_type
  PhysicalType('area')

.. EXAMPLE END

Using Physical Types
====================

.. EXAMPLE START: Using Physical Types

An equality comparison between a |PhysicalType| and a string will return
`True` if the string is a name of the |PhysicalType|::

  >>> acceleration = u.get_physical_type(u.m / u.s ** 2)
  >>> acceleration == 'acceleration'
  True

Some units may correspond to multiple physical types because compatible
units can be used to quantify different phenomena::

  >>> u.get_physical_type('pressure')
  PhysicalType({'energy density', 'pressure', 'stress'})

We can iterate through the names of a |PhysicalType|::

  >>> for name in u.J.physical_type: print(name)
  energy
  torque
  work

We can test for membership or equality with a string that has the name
of a |PhysicalType|::

  >>> 'energy' == u.J.physical_type
  True
  >>> 'work' in u.J.physical_type
  True

.. EXAMPLE END

Dimensional Analysis
====================

.. EXAMPLE START: Dimensional Analysis With Physical Types

|PhysicalType| instances support multiplication, division,
and exponentiation. Because of this, they can be used for
dimensional analysis::

  >>> length = u.get_physical_type('length')
  >>> time = u.get_physical_type('time')
  >>> length ** 2
  PhysicalType('area')
  >>> 1 / time
  PhysicalType('frequency')

Dimensional analysis can be performed between a |PhysicalType| and a
unit or between a |PhysicalType| and a string with a name of a
|PhysicalType|::

  >>> length ** 2 / u.s
  PhysicalType({'diffusivity', 'kinematic viscosity'})
  >>> length / 'time'
  PhysicalType({'speed', 'velocity'})

.. EXAMPLE END
