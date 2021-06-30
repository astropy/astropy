.. doctest-skip-all

.. _read_write_cosmologies:

Reading and Writing Cosmology Objects
*************************************

``astropy`` provides a unified interface for reading and writing data in
different formats. For many common cases this will streamline the process of
file I/O and reduce the need to learn the separate details of all of the I/O
packages within ``astropy``. For details and examples of using this interface
see the :ref:`table_io` section.

Getting Started
===============

The :class:`~astropy.cosmology.Cosmology` class includes two methods,
:meth:`~astropy.cosmology.Cosmology.read` and
:meth:`~astropy.cosmology.Cosmology.write`, that make it possible to read from
and write to files. Two file formats are automatically supported:
`'JSON' <https://www.json.org/json-en.html>`_ and
`'ECSV' <https://zenodo.org/record/4792325>`_.

.. EXAMPLE START: Reading and Writing Cosmology Objects

To use this interface, first import the :class:`~astropy.cosmology.Cosmology`
class, then call the :class:`~astropy.cosmology.Cosmology`
:meth:`~astropy.cosmology.Cosmology.read` method with the name of the file and
the file format, for instance
``'ascii.ecsv'``::

    >>> from astropy.cosmology import Cosmology
    >>> t = Cosmology.read('COBE.ecsv', format='ascii.ecsv')

.. EXAMPLE END

For both the JSON and ECSV formats, the format can be automatically detected
from the filename extension so providing the ``format`` keyword is not required::

    >>> cosmo = Cosmology.read('lightbird.json')

Similarly, for writing, the ``format`` is optional if the file name extension is
either ``.json`` or ``.ecsv``, but it can be provided for code clarity::

    >>> cosmo.write(filename, format='json', overwrite=True)

Any additional arguments specified will depend on the format. For examples of
this see the section :ref:`built_in_cosmology_readers_writers`. This section
also provides the full list of built-in choices for the ``format`` argument.


Working with Mappings and Tables
================================

Reading and writing :class:`~astropy.cosmology.Cosmology` objects go through
intermediate representations of a dict or `~astropy.table.QTable` instances.
These intermediate representations are made accessible. The Table representation
in particular is useful for e.g. writing a Cosmology to *latex* for a journal
article.

The simplest way to access mapping and Table Cosmology I/O is with convenience
methods on :meth:`~astropy.cosmology.Cosmology.read` and
:meth:`~astropy.cosmology.Cosmology.write`. For reading from a mapping or Table
instance use ``Cosmology.read.from_mapping`` or ``Cosmology.read.from_table``.
For representing a Cosmology instance as a mapping or Table use
``Cosmology.write.to_mapping`` or ``Cosmology.write.to_table``.

.. EXAMPLE START: Planck18 to Table

    >>> from astropy.cosmology import Planck18
    >>> table = Planck18.write.to_table()
    >>> table
    '<QTable length=1>
      name       H0     Om0    Tcmb0    Neff    m_nu [3]    Ob0
      str8    float64 float64 float64 float64   float64   float64
    --------  ------- ------- ------- ------- ----------- -------
    Planck18    67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897'

.. EXAMPLE END

Equivalent functions to the above methods are also available in
:ref:`cosmology.io <built_in_cosmology_readers_writers>`.

- :func:`~astropy.cosmology.io.from_mapping`
- :func:`~astropy.cosmology.io.to_mapping`
- :func:`~astropy.cosmology.io.from_table`
- :func:`~astropy.cosmology.io.to_table`


.. _built_in_cosmology_readers_writers:

Built-In Cosmology Readers/Writers
==================================

The :class:`~astropy.cosmology.Cosmology` class has built-in support for various
input and output formats.

A full list of the supported formats and corresponding classes is shown in the
table below. The ``Suffix`` column indicates the filename suffix indicating a
particular format.

===========================  ======  ============================================================================================
           Format            Suffix                                          Description
===========================  ======  ============================================================================================
                 ascii.ecsv   .ecsv  :class:`~astropy.io.ascii.Ecsv`: Basic table with Enhanced CSV (supporting metadata)
                       json   .json  JSON mapping
===========================  ======  ============================================================================================


Reference/API
=============

.. automodapi:: astropy.cosmology.io
