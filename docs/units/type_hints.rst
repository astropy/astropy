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

Unit Type Stubs
===============

As of ``astropy`` 7.2, type stub files are included for individual unit
objects like ``u.deg`` or ``u.m`` (which are created dynamically at runtime
and thus not easily recognized by static type checkers by default).  These
should appear along with ``astropy``, even if you install from a source
checkout in editable mode.

.. note:: The automatic generation is experimental, and we welcome feedback on
          it. If needed, you can (re)generate the ``*.pyi`` files with the
          ``typestubs-for-units`` script.

Type Annotations Module
***********************

.. automodule:: astropy.units.typing
   :members:
   :show-inheritance:
