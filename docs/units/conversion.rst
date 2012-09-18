Unit Conversion
===============

There are two ways of handling conversions between units.

Direct Conversion
-----------------

In this case, given a source and destination unit, the value(s) in the
new units is(are) returned.

  >>> u.pc.to(u.m, 3.26)
  1.0059317615e+17

This converts 3.26 parsecs to meters.

Arrays are permitted as arguments.

  >>> u.h.to(u.s, [1, 2, 5, 10.1])
  array([  3600.,   7200.,  18000.,  36360.])

Obtaining a Conversion Function
-------------------------------

Finally, one may obtain a function that can be used to convert to the
new unit. Normally this may seem like overkill when all one needs to
do is multiply by a scale factor, but there are cases where it is not
so simple, for example when there are equivalencies in use.

Conversion to different units involves obtaining a conversion function
and then applying it to the value, or values to be converted.

  >>> speed_unit = u.cm / u.s
  >>> speed_converter = speed_unit.get_converter(u.mile / u.hour)
  >>> speed_converter(100.)
  2.2366936292054402
  >>> speed_converter([1000, 2000])
  array([ 22.36936292,  44.73872584])

Incompatible Conversions
------------------------

If you attempt to convert to a incompatible unit, an exception will result:

  >>> speed_unit.to(u.mile)
  ...
  UnitsException: 'cm / (s)' and 'mi' are not convertible

You can check whether a particular conversion is possible using the
`is_equivalent` method::

  >>> u.m.is_equivalent(u.foot)
  True
  >>> u.m.is_equivalent("second")
  False
  >>> (u.m**2).is_equivalent(u.acre)
  True
