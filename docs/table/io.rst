.. doctest-skip-all

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
    >>> t = Table.read('photometry.dat', format='ascii.daophot')

It is possible to load tables directly from the Internet using URLs. For
example, download tables from `VizieR catalogs <https://vizier.u-strasbg.fr/>`_
in CDS format (``'ascii.cds'``)::

    >>> t = Table.read("ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/snrs.dat",
    ...         readme="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/ReadMe",
    ...         format="ascii.cds")

.. EXAMPLE END

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified::

    >>> t.write(filename, format='latex')

As for the :meth:`~astropy.table.Table.read` method, the format may
be automatically identified in some cases.

Any additional arguments specified will depend on the format. For examples of
this see the section :ref:`built_in_readers_writers`. This section also
provides the full list of choices for the ``format`` argument.

.. _combining_units_and_quantities:

Combining Units and Quantities
==============================

.. EXAMPLE START: Combining Units and Quantities

A |QTable| can be used to store data as :ref:`user-defined units 
<astropy:defining_units>` in a FITS file.

  >>> import numpy as np
  >>> from astropy.table import QTable
  >>> import astropy.units as u
  >>> lps = u.def_unit('Lps', u.L / u.s)
  >>> qt = QTable()
  >>> qt['speeds'] = np.arange(5) * lps
  >>> qt.write('Lps_output.fits', overwrite=True, format='fits')

In order to parse and read the file contents, it is necessary to enable the
new defined unit by calling :func:`~astropy.units.add_enabled_units`, otherwise
the unit will fail to parse::

  >>> qt2 = QTable.read('Lps_output.fits', format='fits') # doctest: +SHOW_WARNINGS
  Traceback (most recent call last):
  astropy.table.meta.YamlParseError
  >>> u.add_enabled_units(lps)
  <astropy.units.core._UnitContext object at 0x...>
  >>> qt2 = QTable.read('Lps_output.fits', format='fits')
  <QTable length=5>
       speeds
    1000 cm3 / s
      float64
    ------------
             0.0
             1.0
             2.0
             3.0
             4.0

.. EXAMPLE END

Supported Formats
=================

The :ref:`table_io` has built-in support for the following data file formats:

* :ref:`table_io_ascii`
* :ref:`table_io_hdf5`
* :ref:`table_io_fits`
* :ref:`table_io_votable`
* :ref:`table_io_parquet`
