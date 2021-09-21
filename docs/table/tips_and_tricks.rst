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

  >>> lps = u.def_unit('Lps', u.L / u.s)
  >>> lps_array = np.ones((5,5)) * lps
  >>> qt = QTable()
  >>> qt['speeds'] = lps_array 
  >>> qt.write('Lps_output.fits',overwrite=True,format='fits')

In order to parse and read the file contents, it is necessary to enable 
the new defined unit as shown in :ref:`this example 
<astropy:enabling_units>`::

  >>> u.add_enabled_units(lps)
  <astropy.units.core._UnitContext object at 0x...>
  >>> QTable.read('Lps_output.fits', format='fits')
  <QTable length=5>
     speeds [5]
    1000 cm3 / s
      float64
    ------------
      1.0 .. 1.0
      1.0 .. 1.0
      1.0 .. 1.0
      1.0 .. 1.0
      1.0 .. 1.0

.. EXAMPLE END