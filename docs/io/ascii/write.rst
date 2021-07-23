.. include:: references.txt

.. _astropy.io.ascii_write:

Writing Tables
==============

:mod:`astropy.io.ascii` is able to write ASCII tables out to a file or file-like
object using the same class structure and basic user interface as for reading
tables.

The |write| function provides a way to write a data table as a
formatted ASCII table.

Examples
--------

..
  EXAMPLE START
  Writing ASCII Tables Using astropy.io.ascii

To write a formatted ASCII table using the |write| function::

  >>> import numpy as np
  >>> from astropy.io import ascii
  >>> from astropy.table import Table
  >>> data = Table()
  >>> data['x'] = np.array([1, 2, 3], dtype=np.int32)
  >>> data['y'] = data['x'] ** 2
  >>> ascii.write(data, 'values.dat', overwrite=True)  # doctest: +SKIP

The ``values.dat`` file will then contain::

  x y
  1 1
  2 4
  3 9

It is also possible and encouraged to use the write functionality from
:mod:`astropy.io.ascii` through a higher level interface in the :ref:`Data
Tables <astropy-table>` package (see :ref:`table_io` for more details). For
example::

  >>> data.write('values.dat', format='ascii', overwrite=True)  # doctest: +SKIP

For a more reproducible ASCII version of your table, we recommend using the
:ref:`ecsv_format`. This stores all the table meta-data (in particular the
column types and units) to a comment section at the beginning while still
maintaining compatibility with most plain CSV readers. It also allows storing
richer data like `~astropy.coordinates.SkyCoord` or multidimensional or
variable-length columns. For our simple example::

  >>> data.write('values.ecsv', overwrite=True)  # doctest: +SKIP

The ``.ecsv`` extension is recognized and implies using ECSV (equivalent to
``format='ascii.ecsv'``). The ``values.ecsv`` file will then contain::

  # %ECSV 1.0
  # ---
  # datatype:
  # - {name: x, datatype: int32}
  # - {name: y, datatype: int32}
  # schema: astropy-2.0
  x y
  1 1
  2 4
  3 9

Most of the input table :ref:`supported_formats` for
reading are also available for writing. This provides a great deal of
flexibility in the format for writing. The example below writes the data as a
LaTeX table, using the option to send the output to ``sys.stdout`` instead of a
file::

  >>> ascii.write(data, format='latex')  # doctest: +SKIP
  \begin{table}
  \begin{tabular}{cc}
  x & y \\
  1 & 1 \\
  2 & 4 \\
  3 & 9 \\
  \end{tabular}
  \end{table}

There is also a faster Cython engine for writing simple formats,
which is enabled by default for these formats (see :ref:`fast_ascii_io`).
To disable this engine, use the parameter ``fast_writer``::

   >>> ascii.write(data, 'values.csv', format='csv', fast_writer=False)  # doctest: +SKIP

..
  EXAMPLE END

.. Note::

   For most supported formats one can write a masked table and then read it back
   without losing information about the masked table entries. This is
   accomplished by using a blank string entry to indicate a masked (missing)
   value. See the :ref:`replace_bad_or_missing_values` section for more
   information.

.. _io_ascii_write_parameters:

Parameters for ``write()``
--------------------------

The |write| function accepts a number of parameters that specify the detailed
output table format. Each of the :ref:`supported_formats` is handled by a
corresponding Writer class that can define different defaults, so the
descriptions below sometimes mention "typical" default values. This refers to
the :class:`~astropy.io.ascii.Basic` writer and other similar Writer classes.

Some output format Writer classes (e.g., :class:`~astropy.io.ascii.Latex` or
:class:`~astropy.io.ascii.AASTex`) accept additional keywords that can
customize the output further. See the documentation of these classes for
details.

**output**: output specifier
  There are two ways to specify the output for the write operation:

  - Name of a file (string)
  - File-like object (from open(), StringIO, etc.)

**table**: input table
  Any value that is supported for initializing a |Table| object (see
  :ref:`construct_table`). This includes a table with a list of columns, a
  dictionary of columns, or from `numpy` arrays (either structured or
  homogeneous).

**format**: output format (default='basic')
  This specifies the format of the ASCII table to be written, such as a basic
  character delimited table, fixed-format table, or a CDS-compatible table,
  etc. The value of this parameter must be one of the :ref:`supported_formats`.

**delimiter**: column delimiter string
  A one-character string used to separate fields which typically defaults to
  the space character. Other common values might be "," or "|" or "\\t".

**comment**: string defining start of a comment line in output table
  For the :class:`~astropy.io.ascii.Basic` Writer this defaults to "# ".
  Which comments are written and how depends on the format chosen.
  The comments are defined as a list of strings in the input table
  ``meta['comments']`` element. Comments in the metadata of the given
  |Table| will normally be written before the header, although
  :class:`~astropy.io.ascii.CommentedHeader` writes table comments after the
  commented header. To disable writing comments, set ``comment=False``.

**formats**: dict of data type converters
  For each key (column name) use the given value to convert the column data to
  a string. If the format value is string-like, then it is used as a Python
  format statement (e.g., '%0.2f' % value). If it is a callable function, then
  that function is called with a single argument containing the column value to
  be converted. Example::

    astropy.io.ascii.write(table, sys.stdout, formats={'XCENTER': '%12.1f',
                                                 'YCENTER': lambda x: round(x, 1)},

**names**: list of output column names
  Define the complete list of output column names to write for the data table,
  overriding the existing column names.

**include_names**: list of names to include in output
  From the list of column names found from the data table or the ``names``
  parameter, select for output only columns within this list. If not supplied
  then include all names.

**exclude_names**: list of names to exclude from output
  Exclude these names from the list of output columns. This is applied *after*
  the ``include_names`` filtering. If not specified then no columns are excluded.

**fill_values**: list of fill value specifiers
  This can be used to fill missing values in the table or replace values with special meaning.

  See the :ref:`replace_bad_or_missing_values` section for more information on
  the syntax. The syntax is almost the same as when reading a table.
  There is a special value ``astropy.io.ascii.masked`` that is used to say
  "output this string for all masked values in a masked table" (the default is
  to use an empty string ``""``)::

      >>> import sys
      >>> from astropy.table import Table, Column, MaskedColumn
      >>> from astropy.io import ascii
      >>> t = Table([(1, 2), (3, 4)], names=('a', 'b'), masked=True)
      >>> t['a'].mask = [True, False]
      >>> ascii.write(t, sys.stdout)
      a b
      "" 3
      2 4
      >>> ascii.write(t, sys.stdout, fill_values=[(ascii.masked, 'N/A')])
      a b
      N/A 3
      2 4

  Note that when writing a table, all values are converted to strings before
  any value is replaced. Because ``fill_values`` only replaces cells that
  are an exact match to the specification, you need to provide the string
  representation (stripped of whitespace) for each value. For example, in
  the following commands ``-99`` is formatted with two digits after the
  comma, so we need to replace ``-99.00`` and not ``-99``::

      >>> t = Table([(-99, 2), (3, 4)], names=('a', 'b'))
      >>> ascii.write(t, sys.stdout, fill_values = [('-99.00', 'no data')],
      ...             formats={'a': '%4.2f'})
      a b
      "no data" 3
      2.00 4

  Similarly, if you replace a value in a column that has a fixed length format
  (e.g., ``'f4.2'``), then the string you want to replace must have the same
  number of characters. In the example above, ``fill_values=[(' nan',' N/A')]``
  would work.

**fill_include_names**: list of column names, which are affected by ``fill_values``
  If not supplied, then ``fill_values`` can affect all columns.

**fill_exclude_names**: list of column names, which are not affected by ``fill_values``
  If not supplied, then ``fill_values`` can affect all columns.

**fast_writer**: whether to use the fast Cython writer
  If this parameter is ``None`` (which it is by default), |write| will attempt
  to use the faster writer (described in :ref:`fast_ascii_io`) if possible.
  Specifying ``fast_writer=False`` disables this behavior.

**Writer** : Writer class (*deprecated* in favor of ``format``)
  This specifies the top-level format of the ASCII table to be written, such as
  a basic character delimited table, fixed-format table, or a CDS-compatible
  table, etc. The value of this parameter must be a Writer class. For basic
  usage this means one of the built-in :ref:`extension_reader_classes`.
  Note that Reader classes and Writer classes are synonymous; in other
  words, Reader classes can also write, but for historical reasons they are
  often called Reader classes.

.. _cds_mrt_format:

CDS/MRT Format
----------

Both `CDS <http://vizier.u-strasbg.fr/doc/catstd.htx>`_ and
`Machine Readable Table (MRT) <https://journals.aas.org/mrt-standards/>`_ formats consist
of a table description and the table data itself. MRT differs slightly from the CDS
format in table description sections. CDS format includes more detailed description
in the form of ``Abstract``, ``Notes``, ``References`` fields etc. and often has it in a
separate file called ``ReadMe`` file. On the other hand, MRT format includes just the
table ``Title``, ``Authors``  Table Caption and ``Notes`` and always has the ``ReadMe``
section together with the data.

The :class:`~astropy.io.ascii.Cds` writer currently supports writing tables to MRT format.

.. note::

    The metadata of the table, apart from column ``unit``, ``name`` and ``description``,
    will not be written in the output file. This also includes CDS/MRT format specific
    ``ReadMe`` fields. They have to be filled in by hand later.

Examples
""""""""

..
  EXAMPLE START
  Writing CDS/MRT Format Tables Using astropy.io.ascii

The following writes a simple ``astropy`` `~astropy.table.Table` to MRT format.

  >>> from astropy.io import ascii
  >>> from astropy.table import Table
  >>> data = ['names e d s i',
  ...         'HD81809 1E-7 22.25608 +2 67',
  ...         'HD103095 -31.6e5 +27.2500 -9E34 -30']
  >>> table = ascii.read(data)
  >>> table.write('simple_table.dat', format='ascii.cds')  # doctest: +SKIP

The file ``simple_table.dat`` will then be as given below. Notice how ``---`` and
``===`` are used to divide the table into different sections::

  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: simple_table.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1- 8  A8     ---    names   Description of names              
  10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e
  16-23  F8.5   ---    d       [22.25/27.25] Description of d    
  25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     
  33-35  I3     ---    i       [-30/67] Description of i         
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  HD81809  1e-07  22.25608   2e+00  67
  HD103095 -3e+06 27.25000  -9e+34 -30

When the table columns do not contain any units, ``---`` is put in the Byte-By-Byte
description for that column. The unit names are tabulated for columns that do have
a ``Unit`` attribute. Also, a ``?`` is prefixed to the column description in the
Byte-By-Byte for ``Masked`` columns or columns that have null values, indicating
them as such. The example below writes the data containing these attributes, using
the option to send the output to ``sys.stdout`` instead of a file::

  >>> from astropy import units as u
  >>> from astropy.table import MaskedColumn
  >>> table.add_column([5.0, 5.0], name='sameF')
  >>> table.add_column([20, 20], name='sameI')
  >>> col_units = [None, u.C, u.kg, u.m/u.s, u.year, None, None]
  >>> table._set_column_attribute('unit', col_units)
  >>> table.add_row(['Sun', '3.25', '0', '5.3e27', '2', '5.0', '20'],
  ...               mask=[False, True, True, False, True, False, False])
  >>> table['e'] = MaskedColumn(table['e'], mask=[False, True, False])
  >>> table['d'] = MaskedColumn(table['d'], mask=[True, True, False])
  >>> table['magnitude'] = [u.Magnitude(25), u.Magnitude(-9), u.Magnitude(1)]

  >>> ascii.write(table, format='ascii.cds')  # doctest: +SKIP
  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: table.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1- 8  A8     ---    names     Description of names                  
  10-14  E5.1   C      e         [0.0/3.25]? Description of e          
  16-19  F4.1   kg     d         ? Description of d                    
  21-28  E8.2   m.s-1  s         [-9e+34/5.3e+27] Description of s     
  30-32  I3     yr     i         [-30/67]? Description of i            
  34-36  F3.1   ---    sameF     [5.0/5.0] Description of sameF        
  38-39  I2     ---    sameI     [20] Description of sameI             
  41-45  E5.1   mag    magnitude [0.0/3981.08] Description of magnitude
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  HD81809  1e-07       2.0e+00  67 5.0 20 1e-10
  HD103095            -9.0e+34 -30 5.0 20 4e+03
  Sun      3e+00  0.0  5.3e+27     5.0 20 4e-01

Columns that are `~astropy.coordinates.SkyCoord` objects or columns with
values that are such objects are recognized as such, and some predefined labels and
description is used for them. Coordinate columns that have `~astropy.coordinates.SphericalRepresentation` are additionally
sub-divided into their coordinate component columns. Representations that have ``ra``
and ``dec`` components are divided into their ``hour``-``min``-``sec``
and ``deg``-``arcmin``-``arcsec`` components respectively. Whereas, columns with
``SkyCoord`` objects in the ``Galactic`` or any of the ``Ecliptic`` frames are divided
into their latitude(``ELAT``/``GLAT``) and longitude components (``ELON``/``GLAT``) only.
The following example illustrates this.

  >>> from astropy.coordinates import SkyCoord
  >>> table = Table()
  >>> table['Name'] = ['ASASSN-15lh']
  >>> table['coord'] = SkyCoord.from_name('ASASSN-15lh')  # doctest: +REMOTE_DATA
  >>> table.write('coord_cols.dat', format='ascii.cds')   # doctest: +SKIP
  >>> table['coord'] = table['coord'].geocentrictrueecliptic  # doctest: +REMOTE_DATA
  >>> table.write('ecliptic_cols.dat', format='ascii.cds')    # doctest: +SKIP

The original table remains accessible as such, while the file is written from a
modified copy of the table. Thus, the contents of ``coords_cols.dat`` will be::

  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: coords_cols.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1-11  A11    ---    Name    Description of Name     
  13-16  F4.1   h      RAh     Right Ascension (hour)  
  18-20  F3.1   min    RAm     Right Ascension (minute)
  22-39  F18.15 s      RAs     Right Ascension (second)
     41  A1     ---    DE-     Sign of Declination     
  42-45  F5.1   deg    DEd     Declination (degree)    
  47-50  F4.1   arcmin DEm     Declination (arcmin)    
  52-67  F16.13 arcsec DEs     Declination (arcsec)    
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 22.0 2.0 15.450000000007265 -61.0 39.0 34.5999960000006

And the file ``ecliptic_cols.dat`` will look like::

  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: ecliptic_cols.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1-11  A11    ---    Name    Description of Name                        
  13-29  F17.13 deg    ELON    Ecliptic Longitude (geocentrictrueecliptic)
  31-48  F18.14 deg    ELAT    Ecliptic Latitude (geocentrictrueecliptic) 
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 306.2242086500961 -45.62178985082456

Also note that no internal conversion or modification takes place for columns with
`~astropy.time.Time` and related values. They should be converted to regular columns
with proper ``unit`` and ``name`` attribute before writing the table. Thus::

  >>> from astropy.table import Column
  >>> from astropy.time import Time, TimeDelta
  >>> from astropy.timeseries import TimeSeries
  >>> ts = TimeSeries(time_start=Time('2019-1-1'), time_delta=2*u.day, n_samples=1)
  >>> table['time'] = Column(ts.time.decimalyear, unit=u.year)
  >>> table['timeDelta'] = Column(TimeDelta(100.0, format='sec').datetime.seconds,
                                  unit=u.s)

  >>> ascii.write(table, format='ascii.cds')  # doctest: +SKIP
  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: ecliptic_cols.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1-11  A11    ---    Name      Description of Name                        
  13-18  F6.1   yr     time      [2019.0/2019.0] Description of time        
  20-22  I3     s      timeDelta [100] Description of timeDelta             
  24-40  F17.13 deg    ELON      Ecliptic Longitude (geocentrictrueecliptic)
  42-59  F18.14 deg    ELAT      Ecliptic Latitude (geocentrictrueecliptic) 
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 2019.0 100 306.2242086500961 -45.62178985082456

..
  EXAMPLE END

.. attention::

    The CDS writer currently supports automatically writing coordinate component
    columns only for tables with a single coordinate column. For tables with more
    than one coordinate columns, only the first found coordinate column will be
    converted to its component columns and the rest of the coordinate columns will
    be converted to string columns. Thus, it should be taken care that the additional
    coordinate columns are dealt with beforehand using ``SkyCoord`` methods.
