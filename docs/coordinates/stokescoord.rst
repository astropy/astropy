.. _stokes-coord:

Using the ``StokesCoord`` Class
*******************************

The `.StokesCoord` class provides a minimal wrapper for interacting with Stokes
parameters as handled by WCS.  In WCS (as described in the FITS 4.0
specification), the stokes parameters are allocated an integer number, each
representing a different parameter, for example "I" is given the number 1 and
"Q" the number 2.  The `.StokesCoord` class uses a mapping between these numbers
and the names of the parameters as given by the FITS WCS conventions, and allows
you to refer to the parameters by their string names rather than their numbers.
For example, the default representation of the `.StokesCoord` object uses the
names rather than the numbers.

  >>> import numpy as np
  >>> from astropy.coordinates import StokesCoord
  >>> StokesCoord([1, 2, 3, 4])
  StokesCoord(['I', 'Q', 'U', 'V'])

These "symbols" as the `.StokesCoord` class refers to them as can also be accessed via the `.StokesCoord.symbol` property::

  >>> stokes = StokesCoord([1, 2, 3, 4])
  >>> stokes.symbol
  array(['I', 'Q', 'U', 'V'], dtype='<U2')

And the `.StokesCoord` class can also be instantiated with symbols rather than numbers:

  >>> StokesCoord("I")
  StokesCoord('I')

Numeric values which are not present in the mapping will be represented by the string ``"?"``::

  >>> StokesCoord([1, 10])
  StokesCoord(['I', '?'])

It is possible to add custom number - symbol mappings see :ref:`mapping-stokes-symbols`.


Comparing to Symbols and Numbers
--------------------------------

It is possible to compare the values in a `.StokesCoord` with their parameter names::

  >>> stokes = StokesCoord([1, 2, 3, 4])
  >>> stokes == "I"
  array([ True, False, False, False])

Comparing `.StokesCoord` classes like this compares their underlying numeric
values, by converting the symbol to a number.
To directly compare the symbols you can do::

  >>> stokes = StokesCoord([1, 2, 3, 4, 5, 6])
  >>> stokes.symbol == "?"
  array([False, False, False, False, True, True])

This comparison matches all values where the equivalent symbol isn't known, in
comparison to comparing to values where the underlying numbers could vary but
all are unknown (so they are not equal).

Note that comparing to symbols does not work with ``numpy.equal`` only the ``==`` operator.

.. _mapping-stokes-symbols:

Mapping Symbols to Numeric Values
=================================

Built-in Symbols
----------------

The mapping between the numbers and the parameter names built into astropy is as follows::

  >>> from astropy.coordinates.polarization import FITS_STOKES_VALUE_SYMBOL_MAP
  >>> for number, symbol in FITS_STOKES_VALUE_SYMBOL_MAP.items():
  ...     print(f"{number:-2}: {symbol.symbol} - {symbol.description}")
     1: I - Standard Stokes unpolarized
     2: Q - Standard Stokes linear
     3: U - Standard Stokes linear
     4: V - Standard Stokes circular
    -1: RR - Right-right circular: <RR*>
    -2: LL - Left-left circular: <LL*>
    -3: RL - Right-left cross-circular: Re(<RL*>))
    -4: LR - Left-right cross-circular: Re(<LR*>)=Im(<RL*>)
    -5: XX - X parallel linear: <XX*>
    -6: YY - Y parallel linear: <YY*>
    -7: XY - XY cross linear: Re(<XY*>)
    -8: YX - YX cross linear: Im(<XY*>)


Adding Custom Symbols
---------------------

It is possible to add custom mappings between numbers and parameters.

This can be done with the `.custom_stokes_symbol_mapping` context manager::

  >>> from astropy.coordinates import custom_stokes_symbol_mapping, StokesSymbol
  >>> with custom_stokes_symbol_mapping({10: StokesSymbol("J", "A custom parameter name")}):
  ...     print(StokesCoord(10))
  StokesCoord('J')
