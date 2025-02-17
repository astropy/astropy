.. _table_io:
.. _io-unified:

High-level Unified File I/O
***************************

``astropy`` provides a unified interface for reading and writing data in
different formats. For many common cases this will streamline the process of
file I/O and reduce the need to learn the separate details of all of the I/O
packages within ``astropy``. For details on the implementation see
:ref:`io_registry`.


Table of contents

.. toctree::
    :maxdepth: 2

    unified_image
    unified_table

Getting Help on Readers and Writers
-----------------------------------

Each file format is handled by a specific reader or writer, and each of those
functions will have its own set of arguments.

To get help on the available arguments for each format, use the ``help()`` method of the
appropriate ``read()`` or ``write()`` class method, e.g., `astropy.table.Table.read`.
Each of these calls prints a help document which is divided into two sections, the
generic read/write documentation (common to any call) and the format-specific
documentation. For ASCII tables, the format-specific documentation includes the generic
`astropy.io.ascii` package interface and then a description of the particular ASCII
sub-format.

In the examples below we do not show the long output:

.. doctest-skip::

    >>> from astropy.table import Table
    >>> from astropy.nddata import CCDData
    >>> CCDData.read.help('fits')
    >>> Table.read.help('ascii')
    >>> Table.read.help('ascii.latex')
    >>> Table.write.help('hdf5')
    >>> Table.write.help('csv')
