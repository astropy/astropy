.. doctest-skip-all

.. _table_io:

Unified file read/write interface
***********************************

Astropy provides a unified interface for reading and writing data in different formats.
For many common cases this will simplify the process of file I/O and reduce the need to
master the separate details of all the I/O packages within Astropy.  This functionality is
still in active development and the number of supported formats will be increasing.  For
details on the implementation see :ref:`io_registry`.

Getting started with Table I/O
==============================

The :class:`~astropy.table.Table` class includes two methods,
:meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write`, that make it possible to read from
and write to files. A number of formats are automatically supported (see
`Built-in table readers/writers`_) and new file formats and extensions can be
registered with the :class:`~astropy.table.Table` class (see
:ref:`io_registry`).

To use this interface, first import the :class:`~astropy.table.Table` class, then
simply call the :class:`~astropy.table.Table`
:meth:`~astropy.table.Table.read` method with the name of the file and
the file format, for instance ``'ascii.daophot'``::

    >>> from astropy.table import Table
    >>> t = Table.read('photometry.dat', format='ascii.daophot')

It is possible to load tables directly from the Internet using URLs. For example,
download tables from Vizier catalogues in CDS format (``'ascii.cds'``)::

    >>> t = Table.read("ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/snrs.dat", 
    ...         readme="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/ReadMe", 
    ...         format="ascii.cds")

For certain file formats, the format can be automatically detected, for
example from the filename extension::

    >>> t = Table.read('table.tex')

Similarly, for writing, the format can be explicitly specified::

    >>> t.write(filename, format='latex')

As for the :meth:`~astropy.table.Table.read` method, the format may
be automatically identified in some cases.

The underlying file handler will also automatically detect various
compressed data formats and transparently uncompress them as far as
supported by the Python installation (see
:meth:`~astropy.utils.data.get_readable_fileobj`).

Any additional arguments specified will depend on the format.  For examples of this see the
section `Built-in table readers/writers`_.  This section also provides the full list of
choices for the ``format`` argument.

.. _built_in_readers_writers:

Built-in table readers/writers
==============================

The :class:`~astropy.table.Table` class has built-in support for various input 
and output formats including :ref:`table_io_ascii`, 
-:ref:`table_io_fits`, :ref:`table_io_hdf5`, and :ref:`table_io_votable`.

A full list of the supported formats and corresponding classes
is shown in the table below.
The ``Write`` column indicates those formats that support write functionality, and 
the ``Suffix`` column indicates the filename suffix indicating a particular format.
If the value of ``Suffix`` is ``auto``, the format is auto-detected from the file itself.
Not all formats support auto-detection.

===========================  =====  ======  ============================================================================================
           Format            Write  Suffix                                          Description                                    
===========================  =====  ======  ============================================================================================
                      ascii    Yes          ASCII table in any supported format (uses guessing)                                     
               ascii.aastex    Yes          :class:`~astropy.io.ascii.AASTex`: AASTeX deluxetable used for AAS journals             
                ascii.basic    Yes          :class:`~astropy.io.ascii.Basic`: Basic table with custom delimiters                    
                  ascii.cds     No          :class:`~astropy.io.ascii.Cds`: CDS format table                                                
     ascii.commented_header    Yes          :class:`~astropy.io.ascii.CommentedHeader`: Column names in a commented line            
                  ascii.csv    Yes    .csv  :class:`~astropy.io.ascii.Csv`: Basic table with comma-separated values                         
              ascii.daophot     No          :class:`~astropy.io.ascii.Daophot`: IRAF DAOphot format table                           
                 ascii.ecsv    Yes   .ecsv  :class:`~astropy.io.ascii.Ecsv`: Basic table with Enhanced CSV (supporting metadata)    
          ascii.fixed_width    Yes          :class:`~astropy.io.ascii.FixedWidth`: Fixed width                                              
ascii.fixed_width_no_header    Yes          :class:`~astropy.io.ascii.FixedWidthNoHeader`: Fixed width with no header                       
 ascii.fixed_width_two_line    Yes          :class:`~astropy.io.ascii.FixedWidthTwoLine`: Fixed width with second header line               
                 ascii.html    Yes   .html  :class:`~astropy.io.ascii.HTML`: HTML table                                             
                 ascii.ipac    Yes          :class:`~astropy.io.ascii.Ipac`: IPAC format table                                              
                ascii.latex    Yes    .tex  :class:`~astropy.io.ascii.Latex`: LaTeX table                                           
            ascii.no_header    Yes          :class:`~astropy.io.ascii.NoHeader`: Basic table with no headers                                
                  ascii.rdb    Yes    .rdb  :class:`~astropy.io.ascii.Rdb`: Tab-separated with a type definition header line                
                  ascii.rst    Yes    .rst  :class:`~astropy.io.ascii.RST`: reStructuredText simple format table
           ascii.sextractor     No          :class:`~astropy.io.ascii.SExtractor`: SExtractor format table                          
                  ascii.tab    Yes          :class:`~astropy.io.ascii.Tab`: Basic table with tab-separated values  
                       fits    Yes    auto  :mod:`~astropy.io.fits`: Flexible Image Transport System file
                       hdf5    Yes    auto  HDF5_: Hierarchical Data Format binary file
                    votable    Yes    auto  :mod:`~astropy.io.votable`: Table format used by Virtual Observatory (VO) initiative
===========================  =====  ======  ============================================================================================

.. _table_io_ascii:

ASCII formats
--------------

The :meth:`~astropy.table.Table.read` and
:meth:`~astropy.table.Table.write` methods can be used to read and write formats
supported by `astropy.io.ascii`.

Use ``format='ascii'`` in order to interface to the generic
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions from `astropy.io.ascii`.  When reading a table this means
that all supported ASCII table formats will be tried in order to successfully
parse the input.  For example::

  >>> t = Table.read('astropy/io/ascii/tests/t/latex1.tex', format='ascii')
  >>> print(t)
  cola colb colc
  ---- ---- ----
     a    1    2
     b    3    4

When writing a table with ``format='ascii'`` the output is a basic
character-delimited file with a single header line containing the
column names.

All additional arguments are passed to the `astropy.io.ascii`
:func:`~astropy.io.ascii.read` and :func:`~astropy.io.ascii.write`
functions. Further details are available in the sections on
:ref:`io_ascii_read_parameters` and :ref:`io_ascii_write_parameters`.  For example, to change
column delimiter and the output format for the ``colc`` column use::

  >>> t.write(sys.stdout, format='ascii', delimiter='|', formats={'colc': '%0.2f'})
  cola|colb|colc
  a|1|2.00
  b|3|4.00


.. note::

   When specifying a specific ASCII table format using the unified interface, the format name is
   prefixed with ``ascii`` in order to identify the format as ASCII-based.  Compare the
   table above to the `astropy.io.ascii` list of :ref:`supported formats <supported_formats>` where the prefix is not
   needed. Therefore the following are equivalent::

     >>> dat = ascii.read('file.dat', format='daophot')
     >>> dat = Table.read('file.dat', format='ascii.daophot')

   For compatibility with astropy version 0.2 and earlier, the following format
   values are also allowed in ``Table.read()``: ``daophot``, ``ipac``, ``html``, ``latex``, and ``rdb``.

.. _table_io_fits:

FITS
----

Reading and writing tables in `FITS <https://fits.gsfc.nasa.gov/>`_ format is
supported with ``format='fits'``. In most cases, existing FITS files should be
automatically identified as such based on the header of the file, but if not,
or if writing to disk, then the format should be explicitly specified.

Reading
^^^^^^^^

If a FITS table file contains only a single table, then it can be read in
with::

    >>> from astropy.table import Table
    >>> t = Table.read('data.fits')

If more than one table is present in the file, you can select the HDU
as follows::

    >>> t = Table.read('data.fits', hdu=3)

In this case if the ``hdu`` argument is omitted then the first table found will be
read in and a warning will be emitted::

    >>> t = Table.read('data.fits')
    WARNING: hdu= was not specified but multiple tables are present, reading in first available table (hdu=1) [astropy.io.fits.connect]

Writing
^^^^^^^^

To write a table ``t`` to a new file::

    >>> t.write('new_table.fits')

If the file already exists and you want to overwrite it, then set the
``overwrite`` keyword::

    >>> t.write('existing_table.fits', overwrite=True)

At this time there is no support for appending an HDU to an existing
file or writing multi-HDU files using the Table interface. Instead one
can use the convenience function
:func:`~astropy.io.fits.table_to_hdu` to create a single
binary table HDU and insert or append that to an existing
:class:`~astropy.io.fits.HDUList`.

Keywords
^^^^^^^^^

The FITS keywords associated with an HDU table are represented in the ``meta``
ordered dictionary attribute of a :ref:`Table <astropy-table>`.  After reading
a table one can view the available keywords in a readable format using::

  >>> for key, value in t.meta.items():
  ...     print('{0} = {1}'.format(key, value))

This does not include the "internal" FITS keywords that are required to specify
the FITS table properties (e.g. ``NAXIS``, ``TTYPE1``). ``HISTORY`` and
``COMMENT`` keywords are treated specially and are returned as a list of
values.

Conversely, the following shows examples of setting user keyword values for a
table ``t``::

  >>> t.meta['MY_KEYWD'] = 'my value'
  >>> t.meta['COMMENT'] = ['First comment', 'Second comment', 'etc']
  >>> t.write('my_table.fits', overwrite=True)

The keyword names (e.g. ``MY_KEYWD``) will be automatically capitalized prior
to writing.

At this time, the ``meta`` attribute of the :class:`~astropy.table.Table` class
is simply an ordered dictionary and does not fully represent the structure of a
FITS header (for example, keyword comments are dropped).

.. _fits_astropy_native:

Astropy native objects (mixin columns)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to store not only standard `~astropy.table.Column` objects to a
FITS table HDU, but also the following Astropy native objects
(:ref:`mixin_columns`) within a `~astropy.table.Table` or `~astropy.table.QTable`:

- `astropy.time.Time`
- `astropy.units.Quantity`

Other mixin columns such as `~astropy.coordinates.SkyCoord` or
`~astropy.coordinates.EarthLocation` are not currently supported due to reasons
including extensive metadata and no precise mapping to the FITS standard.

In general a mixin column may contain multiple data components as well as
object attributes beyond the standard Column attributes like ``format`` or
``description``. Abiding by the rules set by the FITS standard requires mapping
of these data components and object attributes to the appropriate FITS table
columns and keywords.  Thus a well defined protocol has been developed to allow
the storage of these mixin columns in FITS while allowing the object to
"round-trip" through the file with no loss of data or attributes.

Quantity
~~~~~~~~

A `~astropy.units.Quantity` mixin column in a `~astropy.table.QTable` is
represented in a FITS table using the ``TUNITn`` FITS column keyword to
incorporate the unit attribute of Quantity.

    >>> from astropy.table import QTable
    >>> import astropy.units as u
    >>> t = QTable([[1, 2] * u.angstrom)])
    >>> t.write('my_table.fits', overwrite=True)
    >>> qt = QTable.read('my_table.fits')
    >>> qt
    <QTable length=2>
      col0
    Angstrom
    float64
    --------
         1.0
         2.0

Time
~~~~

By default, a `~astropy.time.Time` mixin column within a `~astropy.table.Table`
or `~astropy.table.QTable` will be written to FITS in full precision. This will be
done using the FITS time standard by setting the necessary FITS header keywords.

The default behaviour for reading a FITS table into an `~astropy.table.Table`
has historically been to convert all FITS columns to `~astropy.table.Column`
objects, which have closely matching properties. For some columns, however,
closer native astropy representations are possible, and one can indicate these
should be used by passing ``astropy_native=True`` (for backwards compatibility,
this is not done by default). This will convert columns conforming to the
FITS time standard to `~astropy.time.Time` instances, avoiding any loss of
precision. For example::

    >>> from astropy.time import Time
    >>> from astropy.table import Table
    >>> from astropy.coordinates import EarthLocation
    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd',
    ...               location=EarthLocation(-2446354, 4237210, 4077985, unit='m'))
    >>> t.write('my_table.fits', overwrite=True)
    >>> tm = Table.read('my_table.fits', astropy_native=True)
    >>> tm['a']
    <Time object: scale='tt' format='jd' value=[ 2400100.5  2400200.5]>
    >>> tm['a'].location
    <EarthLocation (-2446354.,  4237210.,  4077985.) m>
    >>> tm['a'] == t['a']
    array([ True,  True], dtype=bool)

The same will work with ``QTable``.

In addition to binary table columns, various global time informational FITS
keywords are treated specially with ``astropy_native=True``.  In particular
the keywords ``DATE``, ``DATE-*`` (ISO-8601 datetime strings) and the ``MJD-*``
(MJD date values) will be returned as ``Time`` objects in the Table ``meta``.
For more details regarding the FITS time paper and the implementation,
refer to :ref:`fits_time_column`.

Since not all FITS readers are able to use the FITS time standard, it is also
possible to store `~astropy.time.Time` instances using the `_time_format`.
For this case, none of the special header keywords associated with the
FITS time standard will be set.  When reading this back into Astropy, the
column will be an ordinary Column instead of a `~astropy.time.Time` object.
See the `Details`_ section below for an example.

Details
~~~~~~~

Time as a dimension in astronomical data presents challenges in its
representation in FITS files. The standard has therefore been extended to
describe rigorously the time coordinate in the ``World Coordinate System``
framework. Refer to `FITS WCS paper IV
<http://adsabs.harvard.edu/abs/2015A%26A...574A..36R/>`_ for details.

Allowing ``Time`` columns to be written as time coordinate
columns in FITS tables thus involves storing time values in a way that
ensures retention of precision and mapping the associated metadata to the
relevant FITS keywords.

In accordance with the standard which states that in binary tables one may use
pairs of doubles, the Astropy Time column is written in such a table as a
vector of two doubles ``(TFORM n = ‘2D’) (jd1, jd2)`` where ``JD = jd1 + jd2``.
This reproduces the time values to double-double precision and is the
"lossless" version, exploiting the higher precision provided in binary tables.
Note that ``jd1`` is always a half-integer or integer, while ``abs(jd2) < 1``.
Round-tripping of Astropy written FITS binary tables containing time coordinate
columns has been partially achieved by mapping selected metadata, ``scale`` and
singular ``location`` of `~astropy.time.Time`, to corresponding keywords.  Note
that the arbitrary metadata allowed in `~astropy.table.Table` objects within
the ``meta`` dict is not written and will be lost.

The FITS standard requires an additional translation layer back into
the desired format. In the example stated above, the Time column ``t['a']``
undergoes the translation ``Astropy Time --> FITS --> Astropy Time`` which
corresponds to the format conversion ``mjd --> (jd1, jd2) --> jd``. Thus,
the final conversion from ``(jd1, jd2)`` requires a software implementation which is
fully compliant with the FITS time standard.

Taking this into consideration, the functionality to read/write Time
from/to FITS can be explicitly turned off, by opting to store the time
representation values in the format specified by the ``format`` attribute
of the `~astropy.time.Time` column, instead of the ``(jd1, jd2)`` format, with
no extra metadata in the header. This is the "lossy" version, but can help
portability. For the above example, the FITS column corresponding
to ``t['a']`` will then store ``[100.0 200.0]`` instead of
``[[ 2400100.5, 0. ], [ 2400200.5, 0. ]]``. This is done by using a special
``info.serialize_method`` attribute, as in the following example::

    >>> from astropy.time import Time
    >>> from astropy.table import Table
    >>> from astropy.coordinates import EarthLocation
    >>> t = Table()
    >>> t['a'] = Time([100.0, 200.0], scale='tt', format='mjd')
    >>> t['a'].info.serialize_method['fits'] = 'formatted_value'
    >>> t.write('my_table.fits', overwrite=True)
    >>> tm = Table.read('my_table.fits')
    >>> tm['a']
    <Column name='a' dtype='float64' length=2>
    100.0
    200.0
    >>> tm['a'] == t['a'].value
    array([ True,  True], dtype=bool)

By default, ``serialize_method['fits']`` in a Time column ``info`` is equal to
``'jd1_jd2'``, that is, Time column will be written in full precision.

.. note::

   The Astropy `~astropy.time.Time` object does not precisely map to the FITS
   time standard.

   * FORMAT

     The FITS format considers only three formats, ISO-8601, JD and MJD.
     Astropy Time allows for many other formats like ``unix`` or ``cxcsec``
     for representing the values.

   * LOCATION

     In Astropy Time, location can be an array which is broadcastable to the
     Time values. In the FITS standard, location is a scalar expressed via
     keywords.

   Hence the ``format`` attribute and a vector ``location`` attribute are not
   stored.  After reading from FITS the user must set the ``format`` as desired.

   Reading FITS files with time coordinate columns which are not written by
   Astropy *may* fail.  Astropy supports only a small subset of the rather
   complicated standard.

.. _table_io_hdf5:

HDF5
----

.. _HDF5: http://www.hdfgroup.org/HDF5/
.. _h5py: http://www.h5py.org/

Reading/writing from/to HDF5_ files is
supported with ``format='hdf5'`` (this requires h5py_
to be installed). However, the ``.hdf5``
file extension is automatically recognized when writing files, and HDF5 files
are automatically identified (even with a different extension) when reading
in (using the first few bytes of the file to identify the format), so in most
cases you will not need to explicitly specify ``format='hdf5'``.

Since HDF5 files can contain multiple tables, the full path to the table
should be specified via the ``path=`` argument when reading and writing.
For example, to read a table called ``data`` from an HDF5 file named
``observations.hdf5``, you can do::

    >>> t = Table.read('observations.hdf5', path='data')

To read a table nested in a group in the HDF5 file, you can do::

    >>> t = Table.read('observations.hdf5', path='group/data')

To write a table to a new file, the path should also be specified::

    >>> t.write('new_file.hdf5', path='updated_data')

It is also possible to write a table to an existing file using ``append=True``::

    >>> t.write('observations.hdf5', path='updated_data', append=True)

As with other formats, the ``overwrite=True`` argument is supported for
overwriting existing files. To overwrite only a single table within an HDF5
file that has multiple datasets, use *both* the ``overwrite=True`` and
``append=True`` arguments.

If the metadata of the table cannot be written directly to the HDF5 file 
(e.g. dictionaries), or if you want to preserve the units and description
of tables and columns, use ``serialize_meta=True``::

    >>> t.write('observations.hdf5', path='updated_data', serialize_meta=True)

The way serialized meta are saved in the HDF5 dataset have changed in Astropy 3.0.
Files in the old format are still read correctly. If for some reason the user wants to *write*
in the old format, they will specify the (deprecated) ``compatibility_mode`` keyword

    >>> t.write('observations.hdf5', path='updated_data', serialize_meta=True, compatibility_mode=True)

Finally, when writing to HDF5 files, the ``compression=`` argument can be
used to ensure that the data is compressed on disk::

    >>> t.write('new_file.hdf5', path='updated_data', compression=True)



.. _table_io_jsviewer:

JSViewer
--------

Provides an interactive HTML export of a Table, like the
:class:`~astropy.io.ascii.HTML` writer but using the DataTables_ library, which
allow to visualize interactively an HTML table (with columns sorting, search,
pagination).

To write a table ``t`` to a new file::

    >>> t.write('new_table.html', format='jsviewer')

Several additional parameters can be used:

- *table_id*: the HTML id of the ``<table>`` tag, defaults to ``'table{id}'``
  where ``id`` is the id of the Table object.
- *max_lines*: maximum number of lines.
- *table_class*: HTML classes added to the ``<table>`` tag, can be useful to
  customize the style of the table.
- *jskwargs*: additional arguments passed to :class:`~astropy.table.JSViewer`.
- *css*: CSS style, default to ``astropy.table.jsviewer.DEFAULT_CSS``.
- *htmldict*: additional arguments passed to :class:`~astropy.io.ascii.HTML`.

.. _Datatables: https://www.datatables.net/



.. _table_io_votable:

VO Tables
-----------

Reading/writing from/to `VO table <http://www.ivoa.net/Documents/VOTable/>`_
files is supported with ``format='votable'``. In most cases, existing VO
tables should be automatically identified as such based on the header of the
file, but if not, or if writing to disk, then the format should be explicitly
specified.

If a VO table file contains only a single table, then it can be read in with::

    >>> t = Table.read('aj285677t3_votable.xml')

If more than one table is present in the file, an error will be raised,
unless the table ID is specified via the ``table_id=`` argument::

    >>> t = Table.read('catalog.xml')
    Traceback (most recent call last):
    ...
    ValueError: Multiple tables found: table id should be set via the table_id= argument. The available tables are twomass, spitzer

    >>> t = Table.read('catalog.xml', table_id='twomass')

To write to a new file, the ID of the table should also be specified (unless
``t.meta['ID']`` is defined)::

    >>> t.write('new_catalog.xml', table_id='updated_table', format='votable')

When writing, the ``compression=True`` argument can be used to force
compression of the data on disk, and the ``overwrite=True`` argument can be
used to overwrite an existing file.
