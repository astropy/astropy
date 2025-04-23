.. currentmodule:: astropy

****************
Astropy Glossary
****************

.. glossary::

   (``n``,)
      A parenthesized number followed by a comma denotes a tuple with one
      element. The trailing comma distinguishes a one-element tuple from a
      parenthesized ``n``.
      This is from NumPy; see https://numpy.org/doc/stable/glossary.html#term-n.

   -like
      ``<Class>-like`` is an instance of the ``Class`` or a valid initializer argument
      for ``Class`` as ``Class(value)``. E.g. :class:`~astropy.units.Quantity`-like
      includes ``"2 * u.km"`` because ``astropy.units.Quantity("2 * u.km")`` works.

   ['physical type']
       The physical type of a quantity can be annotated in square brackets
       following a `~astropy.units.Quantity` (or similar :term:`quantity-like`).

       For example, ``distance : quantity-like ['length']``

   angle-like
      :term:`quantity-like` and a valid initializer for `~astropy.coordinates.Angle`.
      The ``unit`` must be an angular. A string input is interpreted as an angle as
      described in the `~astropy.coordinates.Angle` documentation.

   buffer-like
      Object that implements `Python's buffer protocol
      <https://docs.python.org/3/c-api/buffer.html#bufferobjects>`_.

   coordinate-like
      :class:`~astropy.coordinates.BaseCoordinateFrame` subclass instance, or a
      :class:`~astropy.coordinates.SkyCoord` (or subclass) instance, or a valid
      initializer as described in :ref:`coordinates-initialization-coord`.

   file-like (readable)
      :term:`python:file-like object` object that supports reading with a method ``read``.

      For a formal definition see :class:`~astropy.io.typing.ReadableFileLike`.

   file-like (writeable)
      :term:`python:file-like object` object that supports writing with a method ``write``.

      For a formal definition see :class:`~astropy.io.typing.WriteableFileLike`.

   frame-like
      :class:`~astropy.coordinates.BaseCoordinateFrame` subclass or subclass instance or
      a valid Frame name (string).

   length-like
      :term:`quantity-like` and a valid initializer for
      :class:`~astropy.coordinates.Distance`. The ``unit`` must be a convertible to a
      unit of length.

   number
      Any scalar numeric type. e.g. `float` or `int` or ``numpy.number``.

   quantity-like
      `~astropy.units.Quantity` (or subclass) instance, a number or `array-like
      <https://numpy.org/doc/stable/glossary.html#term-array_like>`_ object, or a string
      which is a valid initializer for `~astropy.units.Quantity`.

      For a formal definition see :obj:`~astropy.units.typing.QuantityLike`.

   table-like
      :class:`~astropy.table.Table` (or subclass) instance or valid initializer for
      :class:`~astropy.table.Table` as described in :ref:`construct_table`. Common types
      include ``dict[list]``, ``list[dict]``, ``list[list]``, and `~numpy.ndarray`
      (structured array).

   time-like
      :class:`~astropy.time.Time` (or subclass) instance or a valid initializer for
      :class:`~astropy.time.Time`, e.g. `str`, array-like[str], `~datetime.datetime`, or
      `~numpy.datetime64`.

   trait type
      In short, a trait type is a class with the following properties:

      - It is a class that can be used as a mixin to add functionality to another class.
      - It should never be instantiated directly.
      - It should not be used as a base class for other classes, but only as a mixin.
      - It can define methods, properties, and attributes -- any of which can be abstract.
      - It can be generic, i.e. it can have type parameters.
      - It can subclass other traits, but should have a linear MRO.

      These are the same set of properties as orthogonal mixin classes, with the added
      emphasis that they can serve as compiled types, if so enabled by a compilation system such as `mypyc <https://mypyc.readthedocs.io/en/latest/>`_.

   unit-like
      :class:`~astropy.units.UnitBase` subclass instance or a valid initializer for
      :class:`~astropy.units.Unit`, e.g., `str` or scalar `~astropy.units.Quantity`.


Optional Packages' Glossary
***************************

.. currentmodule:: matplotlib.pyplot

.. glossary::

   color
      Any valid Matplotlib color.
