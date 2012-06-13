"""
light-weight units module. A simple set of classes for tracking units.

preface
=======

This code is adapted from the units code in pynbody 
(see http://code.google.com/p/pynbody/ )

The changes so far have been to remove dependencies on all other pynbody code.
This involves:

1) requiring Python 2.6 or later (this hasn't been tested with Python 3 yet)
2) disabling the dimensional_project method on units (currently commented out)
3) moving the default units defined in the pynbody config file to within the
units module itself

The original pynbody has a GPL license. Eventual use of this code with
astropy will require converting it to the standard astropy license. In 
principle permission has been granted, but that was subject to seeing if
anyone else had contributed code to units besides Andrew Pontzen.

This is the initial snapshot of the units code with dependencies removed.

We do intend to make extensions to both the user interface and capabilities
described below



units
=====

Making units
------------

You can make units in two ways. Either you can create a string, and
instantiate a Unit like this:

>>> units.Unit("Msol kpc**-3")
>>> units.Unit("2.1e12 m_p cm**-2/3")

Or you can do it within python, using the predefined Unit objects

>>> units.Msol * units.kpc**-3
>>> 2.1e12 * units.m_p * units.cm**(-2,3)

In the last example, either a tuple describing a fraction or a
Fraction instance (from the standard python module fractions) is
acceptable.


Getting conversion ratios
-------------------------

To convert one unit to another, use the ``ratio`` member function:

>>> units.Msol.ratio(units.kg)  
1.99e30
>>> (units.Msol / units.kpc**3).ratio(units.m_p/units.cm**3)
4.04e-8

If the units cannot be converted, a UnitsException is raised:

>>> units.Msol.ratio(units.kpc)
UnitsException

Specifying numerical values
---------------------------

Sometimes it's necessary to specify a numerical value in the course
of a conversion. For instance, consider a comoving distance; this
can be specified in pynbody units as follows:

>>> comoving_kpc = units.kpc * units.a

where units.a represents the scalefactor. We can attempt to convert
this to a physical distance as follows

>>> comoving_kpc.ratio(units.kpc)

but this fails, throwing a UnitsException. On the other hand we
can specify a value for the scalefactor when we request the conversion

>>> comoving_kpc.ratio(units.kpc, a=0.5)
0.5

and the conversion completes with the expected result. The units
module also defines units.h for the dimensionless hubble constant,
which can be used similarly. *By default, all conversions happening
within a specific simulation context should pass in values for
a and h as a matter of routine.*

Any IrreducibleUnit (see below) can have a value specified in this way,
but a and h are envisaged to be the most useful applications.

Defining new base units
-----------------------

The units module is fully extensible: you can define and name your own
units which then integrate with all the standard functions.

.. code-block:: python

   litre = units.NamedUnit("litre",0.001*units.m**3)
   gallon = units.NamedUnit("gallon",0.004546*units.m**3)
   gallon.ratio(litre) # 4.546
   (units.pc**3).ratio(litre) # 2.94e52


You can even define completely new dimensions.

.. code-block:: python

    V = units.IrreducibleUnit("V") # define a volt
    C = units.NamedUnit("C", units.J/V) # define a coulomb
    q = units.NamedUnit("q", 1.60217646e-19*C) # elementary charge
    F = units.NamedUnit("F", C/V) # Farad
    epsilon0 = 8.85418e-12 *F/units.m


>>> (q*V).ratio("eV") 
1.000
>>> ((q**2)/(4*math.pi*epsilon0*units.m**2)).ratio("N")
2.31e-28


"""
from .core import *