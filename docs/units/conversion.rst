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

Obtaining a Conversion Function
-------------------------------

Finally, one may obtain a function that can be used to convert to the
new unit. Normally this may seem like overkill when all one needs to
do is multiply by a scale factor, but there are cases when the
transformation between units may not be as simple as a single scale
factor, for example when a custom equivalency table is in use.

Conversion to different units involves obtaining a conversion function
and then applying it to the value, or values to be converted.

  >>> cms = u.cm / u.s
  >>> cms_to_kmph = cms.get_converter(u.km / u.hour)
  >>> cms_to_kmph(125.)
  4.5
  >>> cms_to_kmph([1000, 2000])
  array([ 36.,  72.])

Incompatible Conversions
------------------------

If you attempt to convert to a incompatible unit, an exception will result:

  >>> cms.to(u.km)
  Traceback (most recent call last):
    ...
  UnitsError: 'cm / s' (speed) and 'km' (length) are not convertible

You can check whether a particular conversion is possible using the
`~astropy.units.core.UnitBase.is_equivalent` method::

  >>> u.m.is_equivalent(u.pc)
  True
  >>> u.m.is_equivalent("second")
  False
  >>> (u.m ** 3).is_equivalent(u.l)
  True
