Tips and Tricks
***************

Using QTable to store user-defined units in a FITS file
=======================================================

A `~astropy.table.QTable` can be used to store data as :ref:`user-defined
units <astropy:defining_units>` in a FITS file.

.. EXAMPLE START: Combining Units and Quantities

To output data to the current directory::

  >>> titter = u.def_unit('titter') 
  >>> titter_array = np.ones((5,5)) * titter
  >>> qt = QTable()
  >>> qt['Fun'] = titter_array 
  >>> qt.write('output_titter.fits',overwrite=True,format='fits')

In order to parse and read the file contents, it is necessary to enable 
the new defined unit as shown in :ref:`this example 
<astropy:enabling_units>`::

  >>> u.add_enabled_units(titter)
  >>> QTable.read('output_titter.fits', format='fits')
  <QTable length=5>
   Fun [5]
    titter
   float64
  ----------
  1.0 .. 1.0
  1.0 .. 1.0
  1.0 .. 1.0
  1.0 .. 1.0
  1.0 .. 1.0

.. EXAMPLE END