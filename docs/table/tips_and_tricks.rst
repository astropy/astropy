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

  >>> my_unit = u.def_unit('my_unit') 
  >>> my_unit_array = np.ones((5,5)) * my_unit
  >>> qt = QTable()
  >>> qt['prop_to_measure'] = my_unit_array 
  >>> qt.write('my_unit_output.fits',overwrite=True,format='fits')

In order to parse and read the file contents, it is necessary to enable 
the new defined unit as shown in :ref:`this example 
<astropy:enabling_units>`::

  >>> u.add_enabled_units(my_unit)
  >>> QTable.read('my_unit_output.fits', format='fits')
  <QTable length=5>
  prop_to_measure [5]
        my_unit
        float64
  -------------------
           1.0 .. 1.0
           1.0 .. 1.0
           1.0 .. 1.0
           1.0 .. 1.0
           1.0 .. 1.0

.. EXAMPLE END