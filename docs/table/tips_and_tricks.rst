Tips and Tricks
***************

Using QTable to store user-defined units in a FITS file
=======================================================

A `~astropy.table.QTable` can be used to store data as :ref:`user-defined
units <astropy:defining_units>` in a FITS file.

.. EXAMPLE START: Combining Units and Quantities

To output data to the current directory::

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy.table import QTable

  >>> kmph = u.def_unit('kmph', u.km / u.h)
  >>> kmph_array = np.ones((5,5)) * kmph
  >>> qt = QTable()
  >>> qt['speeds'] = kmph_array 
  >>> qt.write('kmph_output.fits',overwrite=True,format='fits')

In order to parse and read the file contents, it is necessary to enable 
the new defined unit as shown in :ref:`this example 
<astropy:enabling_units>`::

  >>> u.add_enabled_units(kmph)
  <astropy.units.core._UnitContext object at 0x...>
  >>> QTable.read('kmph_output.fits', format='fits')
  <QTable length=5>
      speeds [5]
      km / h
     float64
    ----------
    1.0 .. 1.0
    1.0 .. 1.0
    1.0 .. 1.0
    1.0 .. 1.0
    1.0 .. 1.0

.. EXAMPLE END