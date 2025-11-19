Unit-Aware Type Annotations
***************************

Python supports static type analysis using the type syntax of `PEP 484
<https://www.python.org/dev/peps/pep-0484/>`_. For a detailed guide on type
hints, function annotations, and other related syntax see the `Real Python Guide
<https://realpython.com/python-type-checking/#type-aliases>`_. Below we describe
how you can be use Quantity type hints and annotations and also include metadata
about the associated units.


We assume the following imports:

::

   >>> import typing as T
   >>> import astropy.units as u
   >>> from astropy.units import Quantity


.. _quantity_type_annotation:

Quantity Type Annotation
========================

A |Quantity| can be used as a type annotation,::

   >>> x: Quantity = 2 * u.km

or as a function annotation.::

   >>> def func(x: Quantity) -> Quantity:
   ...     return x


Preserving Units
^^^^^^^^^^^^^^^^

While the above annotations are useful for annotating the value's type, it
does not inform us of the other most important attribute of a |Quantity|:
the unit.

Unit information may be included by the syntax
``Quantity[unit or "physical_type", shape, numpy.dtype]``.:

   >>> Quantity[u.m]
   typing.Annotated[astropy.units.quantity.Quantity, Unit("m")]
   >>>
   >>> Quantity["length"]
   typing.Annotated[astropy.units.quantity.Quantity, PhysicalType('length')]

See ``typing.Annotated`` for explanation of ``Annotated``

These can also be used on functions

   >>> def func(x: Quantity[u.kpc]) -> Quantity[u.m]:
   ...     return x << u.m


.. _multiple_annotation:

Multiple Annotations
====================

Multiple Quantity and unit-aware |Quantity| annotations are supported using
:class:`~typing.Union` or :class:`~typing.Optional` (including ``|`` operations).

   >>> Quantity[u.m] | None
   typing.Optional[typing.Annotated[astropy.units.quantity.Quantity, Unit("m")]]
   >>>
   >>> Quantity[u.m] | Quantity["time"]
   typing.Union[typing.Annotated[astropy.units.quantity.Quantity, Unit("m")],
                typing.Annotated[astropy.units.quantity.Quantity, PhysicalType('time')]]


.. _typestubs_for_units:

Script for Creating Unit Type Stubs
===================================

Individual unit objects like ``u.deg`` or ``u.m`` are created dynamically
at runtime, which means static type checkers cannot recognize them by default.
To enable proper type checking for these units, type stub (``.pyi``) files
are needed that declare their types.

As an experimental feature, astropy 7.2 includes a script to generate these
stub files automatically:

.. code-block:: bash

   typestubs-for-units

This command will attempt to write the stub files directly into your
installed ``astropy`` package directory. To write to a different location
or see other options, use the ``--help`` flag.

.. warning:: This script is meant to gain information about whether the
    type hints provided this way are sufficient. Once an automatic
    solution has been identified, the script will be removed. In the
    meantime, we welcome feedback.


Type Annotations Module
***********************

.. automodule:: astropy.units.typing
   :members:
   :show-inheritance:
