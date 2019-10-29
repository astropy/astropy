.. note that if this is changed from the default approach of using an *include*
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-table-performance:

Performance Tips
================

Constructing |Table| objects row-by-row using
:meth:`~astropy.table.Table.add_row` can be very slow::

    >>> from astropy.table import Table
    >>> t = Table(names=['a', 'b'])
    >>> for i in range(100):
    ...    t.add_row((1, 2))

If you do need to loop in your code to create the rows, a much faster approach
is to construct a list of rows and then create the Table object at the very
end::

  >>> rows = []
  >>> for i in range(100):
  ...    rows.append((1, 2))
  >>> t = Table(rows=rows, names=['a', 'b'])

Writing a |Table| with |MaskedColumn| to ``.ecsv`` using
:meth:`~astropy.table.Table.write` can be very slow::

    >>> from astropy.table import Table
    >>> import numpy as np
    >>> x = np.arange(10000, dtype=float)
    >>> tm = Table([x], masked=True)
    >>> tm.write('tm.ecsv', overwrite=True) # doctest: +SKIP

If you want to write ``.ecsv`` using :meth:`~astropy.table.Table.write`,
then use ``serialize_method='data_mask'``.
It uses the non-masked version of data and it is faster::

    >>> tm.write('tm.ecsv', overwrite=True, serialize_method='data_mask') # doctest: +SKIP

Read FITS with memmap=True
--------------------------

By default :meth:`~astropy.table.Table.read` will read the whole table into memory, which 
can take a lot of memory and can take a lot of time, depending on the table size and 
file format. In some cases, it is possible to only read a subset of the table by choosing 
the option ``memmap=True``.

For FITS binary tables, the data is stored row by row, and it is possible to read only a 
subset of rows, but reading a full column loads the whole table data into memory:

.. doctest-skip::

    >>> import numpy as np
    >>> from astropy.table import Table
    >>> tbl = Table({'a': np.arange(1e7),
    ...              'b': np.arange(1e7, dtype=float),
    ...              'c': np.arange(1e7, dtype=float)})
    >>> tbl.write('test.fits', overwrite=True)
    >>> table = Table.read('test.fits', memmap=True)  # Very fast, doesn't actually load data
    >>> table2 = tbl[:100]  # Fast, will read only first 100 rows
    >>> print(table2)  # Accessing column data triggers the read
     a    b    c  
    ---- ---- ----
    0.0  0.0  0.0
    1.0  1.0  1.0
    2.0  2.0  2.0
    ...  ...  ...
    98.0 98.0 98.0
    99.0 99.0 99.0
    Length = 100 rows
    >>> col = table['my_column']  # Will load all table into memory

At the moment :meth:`~astropy.table.Table.read` does not support ``memmap=True`` for the HDF5 and 
ASCII file formats.
