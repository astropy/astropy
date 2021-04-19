.. |PhysicalType| replace:: :class:`~astropy.units.PhysicalType`
.. |Quantity| replace:: :class:`~astropy.units.Quantity`

Physical Types
**************

A physical type corresponds to the physical quantities with
dimensionally compatible units.  For example, the physical type *mass*
corresponds to physical quantities that have units of kg. Physical
types are represented as instances of the |PhysicalType| class.

Accessing Physical Types
========================

.. EXAMPLE START: Getting Physical Types

There are two ways to access |PhysicalType| instances. The physical type
of a unit can be accessed via its `~astropy.units.UnitBase.physical_type`
attribute.

  >>> import astropy.units as u
  >>> u.coulomb.physical_type
  PhysicalType('electrical charge')
  >>> (u.meter ** 2).physical_type
  PhysicalType('area')

Using `~astropy.units.get_physical_type` lets us access |PhysicalType|
instances from strings with a name of a physical type, units, |Quantity|
instances, and objects that can become quantities (e.g., numbers).

  >>> u.get_physical_type('speed')  # from the name of a physical type
  PhysicalType('speed')
  >>> u.get_physical_type(u.meter)  # from a unit
  PhysicalType('length')
  >>> u.get_physical_type(1 * u.barn * u.Mpc)  # from a Quantity
  PhysicalType('volume')
  >>> u.get_physical_type(42)  # from a number
  PhysicalType('dimensionless')

.. EXAMPLE END

Dimensional Analysis
====================

.. EXAMPLE START: Dimensional Analysis With Physical Types

`~astropy.units.PhysicalType` instances are useful for dimensional
analysis.

  >>> from astropy import units as u
  >>> length = u.get_physical_type("length")
  >>> time = u.get_physical_type("time")
  >>> length / time
  PhysicalType('speed')

.. EXAMPLE END
