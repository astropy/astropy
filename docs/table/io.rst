.. _read_write_tables:

Reading and Writing Table Objects
*********************************

``astropy`` provides a unified interface for reading and writing data in
different formats. For many common cases this will streamline the process of
file I/O and reduce the need to learn the separate details of all of the I/O
packages within ``astropy``. For details and examples of using this interface
see the :ref:`table_io` section.

Getting Started
===============

The :class:`~astropy.table.Table` class includes two methods,
:meth:`~astropy.table.Table.read` and :meth:`~astropy.table.Table.write`, that
make it possible to read from and write to files. A number of formats are
automatically supported (see :ref:`built_in_readers_writers`) and new file
formats and extensions can be registered with the :class:`~astropy.table.Table`
class (see :ref:`io_registry`).

.. EXAMPLE START: Reading and Writing Table Objects

To use this interface, first import the :class:`~astropy.table.Table` class,
then call the :class:`~astropy.table.Table` :meth:`~astropy.table.Table.read`
method with the name of the file and the file format, for instance
``'ascii.daophot'``::

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')  # doctest: +SKIP

It is possible to load tables directly from the Internet using URLs. For
example, download tables from `VizieR catalogs <https://vizier.unistra.fr/>`_
in CDS format (``'ascii.cds'``)::

    >>> t = Table.read("ftp://cdsarc.unistra.fr/pub/cats/VII/253/snrs.dat",
    ...         readme="ftp://cdsarc.unistra.fr/pub/cats/VII/253/ReadMe",
    ...         format="ascii.cds")  # doctest: +SKIP

.. EXAMPLE END


Similarly, for writing, the format can be explicitly specified::

    >>> import astropy.units as u
    >>> t = Table({'object': ['Mrk 421', 'BL Lac'], 'G': [13.2, 13.8] * u.mag})
    >>> t.write('target_table.tex', format='latex')

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('target_table.tex')

.. testcleanup::

    >>> import pathlib
    >>> pathlib.Path.unlink('target_table.tex')

As for the :meth:`~astropy.table.Table.read` method, the format may
be automatically identified in some cases.

Any additional arguments specified will depend on the format. For examples of
this see the section :ref:`built_in_readers_writers`. This section also
provides the full list of choices for the ``format`` argument.

Supported Formats
=================

The :ref:`table_io` has built-in support for the following data file formats:

* :ref:`table_io_ascii`
* :ref:`table_io_hdf5`
* :ref:`table_io_fits`
* :ref:`table_io_votable`
* :ref:`table_io_parquet`


Reading and Writing Column Objects
==================================

.. EXAMPLE START: Reading and Writing Column Objects

Individual table columns do not have their own functions for reading and writing
but it is easy to select just a single column (here "obstime") from a table for writing::

    >>> from astropy.time import Time
    >>> tab = Table({'name': ['AB Aur', 'SU Aur'],
    ...              'obstime': Time(['2013-05-23T14:23:12', '2011-11-11T11:11:11'])})
    >>> tab[['obstime']].write('obstime.fits')

Note the notation ``[['obstime']]`` in the last line - indexing a table with a list of strings
gives us a new table with the columns given by the strings. Since the inner list has only
one element, the resulting table has only one column.

Then, we can read back that single-column table and extract the column from it::

    >>> col = Table.read('obstime.fits').columns[0]
    >>> type(col)
    <class 'astropy.table.column.Column'>

.. testcleanup::

    >>> pathlib.Path.unlink('obstime.fits')

.. EXAMPLE END
