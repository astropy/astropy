Low-level unit conversion
=========================

Conversion of quantities from one unit to another is handled using the
`Quantity.to <astropy.units.quantity.Quantity.to>` method.  This page
describes some low-level features for handling unit conversion that
are rarely required in user code.

There are two ways of handling conversions between units.

Direct Conversion
-----------------

In this case, given a source and destination unit, the value(s) in the
new units is(are) returned.

  >>> from astropy import units as u
  >>> u.pc.to(u.m, 3.26)
  1.0059308915583043e+17

This converts 3.26 parsecs to meters.

Arrays are permitted as arguments.

  >>> u.h.to(u.s, [1, 2, 5, 10.1])
  array([  3600.,   7200.,  18000.,  36360.])

Incompatible Conversions
------------------------

If you attempt to convert to a incompatible unit, an exception will result:

  >>> cms = u.cm / u.s
  >>> cms.to(u.km)  # doctest: +IGNORE_EXCEPTION_DETAIL
  Traceback (most recent call last):
    ...
  UnitConversionError: 'cm / s' (speed) and 'km' (length) are not convertible

You can check whether a particular conversion is possible using the
`~astropy.units.core.UnitBase.is_equivalent` method::

  >>> u.m.is_equivalent(u.pc)
  True
  >>> u.m.is_equivalent("second")
  False
  >>> (u.m ** 3).is_equivalent(u.l)
  True
