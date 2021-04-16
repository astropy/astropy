.. currentmodule:: astropy

****************
Astropy Glossary
****************

.. glossary::

   (`n`,)
      A parenthesized number followed by a comma denotes a tuple with one
      element. The trailing comma distinguishes a one-element tuple from a
      parenthesized ``n``.
      This is from NumPy; see https://numpy.org/doc/stable/glossary.html.

   number
      Any numeric type. eg float or int or any of the ``numpy.number``.

   -like
      Used to indicate on object of that type or that can instantiate the type.
      E.g. :class:`~astropy.units.Quantity`-like includes ``"2 * u.km"``
      because ``astropy.units.Quantity("2 * u.km")`` works.

   unit-like
      Must be an :class:`~astropy.units.UnitBase` (subclass) instance or a
      string or other instance parseable by :class:`~astropy.units.Unit`.

   quantity-like
      Must be an `~astropy.units.Quantity` (or subclass) instance or a string
      parseable by `~astropy.units.Quantity`.
      Note that the interpretation of units in strings depends on the class --
      ``Quantity("180d")`` is 180 **days**, while ``Angle("180d")`` is 180
      **degrees** -- so check the string parses as intended for ``Quantity``.

   angle-like
      :term:`quantity-like`, but interpreted by an angular
      `~astropy.units.SpecificTypeQuantity`, like `~astropy.coordinates.Angle`
      or `~astropy.coordinates.Longitude` or `~astropy.coordinates.Latitude`.
      Note that the interpretation of units in strings depends on the class --
      ``Quantity("180d")`` is 180 days, while ``Angle("180d")`` is 180 degrees
      -- so make sure the string parses as intended for ``Angle``.

   length-like
      :term:`quantity-like`, but interpretable by
      :class:`~astropy.coordinates.Distance`.

   frame-like
      A :class:`~astropy.coordinates.BaseCoordinateFrame` subclass instance or a
      string that can be converted to a Frame by
      :class:`~astropy.coordinates.sky_coordinate_parsers._get_frame_class`.

   coordinate-like
      A Coordinate-type object such as a
      :class:`~astropy.coordinates.BaseCoordinateFrame` subclass instance or a
      :class:`~astropy.coordinates.SkyCoord` (or subclass) instance.

   table-like
      An astropy :class:`~astropy.table.Table` or any object that can
      initialize one. Anything marked as table-like will be processed through
      a :class:`~astropy.table.Table`.

   time-like
      :class:`~astropy.time.Time` or any valid initializer.

   buffer-like
      Anything that implements Python's buffer protocol. See
      https://docs.python.org/3/c-api/buffer.html#bufferobjects

   writable file-like object
      In the context of a :term:`python:file-like object` object, anything
      that supports writing with a method ``write``.

   readable file-like object
      In the context of a :term:`python:file-like object` object, anything
      that supports writing with a method ``read``.


***************************
Optional Packages' Glossary
***************************

.. currentmodule:: matplotlib.pyplot

.. glossary::

   color
      Any valid Matplotlib color.