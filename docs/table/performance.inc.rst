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
