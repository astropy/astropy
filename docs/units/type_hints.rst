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
:class:`~typing.Union` or :class:`~typing.Optional`

   >>> T.Union[Quantity[u.m], None]
   typing.Optional[typing.Annotated[astropy.units.quantity.Quantity, Unit("m")]]
   >>>
   >>> T.Union[Quantity[u.m], Quantity["time"]]
   typing.Union[typing.Annotated[astropy.units.quantity.Quantity, Unit("m")],
                typing.Annotated[astropy.units.quantity.Quantity, PhysicalType('time')]]



Type Annotations Module
***********************

.. automodule:: astropy.units.typing
   :members:
   :show-inheritance:
