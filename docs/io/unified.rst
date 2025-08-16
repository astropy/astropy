.. _table_io:
.. _io-unified:

High-level Unified File I/O
***************************

``astropy`` uses the :ref:`I/O Registry <io_registry>` to provide a unified interface
for reading and writing data in different formats. For many common cases this will
streamline the process of file I/O and reduce the need to learn the separate details of
all of the low-level I/O packages within ``astropy``.

.. toctree::
    :maxdepth: 2

    unified_image
    unified_table

**Overview**

The fundamental idea for the unified interface is that each data container class such as
`~astropy.table.Table` or `~astropy.nddata.CCDData` has class methods ``read()``
and ``write()`` that can be used to read and write data.

The first positional argument to these methods specifies the input or output. In
general the input can be a file name, a file-like object or a URL. For some formats,
most notably the :ref:`table_io_ascii` formats, the input can also be a string or
list of strings representing the data. The output can be a file name or a file-like
object.

The file format is specified using the ``format`` keyword argument. This is required
unless the format can be uniquely determined from the file name or file content.

**Example**

The example below shows how to read a table in the specialized DAOphot format and write
it back to FITS format. Notice that FITS is a format where the interface recognizes the
format automatically from the file name, so the ``format`` argument is not needed.

.. testsetup::
    >>> import os
    >>> with open('photometry.dat', 'w') as f: # doctest: +IGNORE_OUTPUT
    ...     f.write("#N ID    XCENTER   YCENTER\n")
    ...     f.write("#U ##    pixel     pixel \n")
    ...     f.write("#F %-9d  %-10.3f   %-10.3f\n")
    ...     f.write("#\n")
    ...     f.write("14       138.538   256.405\n")
    ...     f.write("18       18.114    280.170\n")

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')
    >>> t.write('photometry.fits')

.. testcleanup::

    >>> os.remove('photometry.dat')
    >>> os.remove('photometry.fits')

Each file format is handled by a specific reader or writer, and each of those
functions will have its own set of arguments.

**Getting Help**

To get help on the available arguments for each format, use the ``help()`` method of the
appropriate ``read()`` or ``write()`` class method, e.g., `astropy.table.Table.read`.
In the examples below we do not show the long output:

    >>> from astropy.table import Table
    >>> from astropy.nddata import CCDData
    >>> CCDData.read.help('fits')  # doctest: +IGNORE_OUTPUT
    >>> Table.read.help('ascii')  # doctest: +IGNORE_OUTPUT
    >>> Table.read.help('ascii.latex')  # doctest: +IGNORE_OUTPUT
    >>> Table.write.help('hdf5')  # doctest: +IGNORE_OUTPUT
    >>> Table.write.help('csv')  # doctest: +IGNORE_OUTPUT
