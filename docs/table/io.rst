
Reading and writing Table objects
===================================

Astropy provides a unified interface for reading and writing data
in different formats.  For many common cases this will 
simplify the process of file I/O and reduce the need to master
the separate details of all the I/O packages within Astropy.  For details and 
examples of using this interface see the :ref:`table_io` 
section.

The :class:`~astropy.table.table.Table` class includes two methods,
:meth:`~astropy.table.table.Table.read` and
:meth:`~astropy.table.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
:ref:`built_in_readers_writers`) and new file formats and extensions can be
registered with the :class:`~astropy.table.table.Table` class (see
:ref:`io_registry`).

