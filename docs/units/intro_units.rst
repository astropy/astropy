Introduction to ``astropy.units``
=================================

`astropy.units` is a Python package to handle defining and converting
 between units



Using `astropy.units`
---------------------

Standard units are defined in the module as object instances. 
Of course, any redefinition of any the unit instances will break 
unit behavior.

All units are defined in term of basic unit classes. The following 
classes are currently implemented:

Length
Time
Mass (well, not yet)
Charge (well, not yet)

Units that involve combinations of fundamental units are CompoundUnits.

Examples of use
---------------

  >>> from astropy import units as u
  >>> speed_unit = u.cm / u.s

Conversion to different units involves obtaining a conversion function
and then applying it to the value, or values to be converted.

  >>> speed_unit.convert(u.miles/u.hour)(100.)
  2.2366936292054402

If you attempt to convert to a incompatible unit, an exception will result:

  >>> speed_unit.convert(u.miles)(100.)
  ...
  ValueError: new unit is inconsistent

Users are free to define new units, either fundamental or compound

e.g.:

  >>> fortnight = 14 * u.days + "fortnight"

The addition of a string gives the new unit a name that will show up when
the unit is printed.

Using repr:

  >>> fortnight
  Time(scale=1.209600e+06, name='fortnight')

In this case, one sees it is created using the Time object using a scale
factor of 1209600. Implicit is that Time objects use intrinsic units of 
seconds.

Using print:

  >>> print fortnight
  Units: fortnight

Creating a new fundamental unit is simple

  class Humor(unit.Unit):
  	  """
  	  Standard class for all humor units
  	  
  	  Intrinsic humor unit is a titter
  	  """
  	  def __init__(self, name=", scale=1.0):
  	  	  self.name = name
  	  	  self.scale = scale
  	  	  self.intrinsic = "titter"
  
titter = Humor(name="titter")
chuckle = 5 * titter + "chuckle"
laugh = 4 * chuckle + "laugh"
guffaw = 3 * laugh + "guffaw"
rofl = 4 * guffaw + "rofl"
death_by_laughing = 10 * rofl "death_by_laughing"

  >>> silly_unit = chuckle
  >>> silly_unit.convert(titter)(3)
  15


Special unit cases
------------------

Blah, blah blah
