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

.. _cds_mrt_format:

Machine-Readable Table Format
-----------------------------

The American Astronomical Society Journals' `Machine-Readable Table (MRT)
<https://journals.aas.org/mrt-standards/>`_ format consists of single file with
the table description header and the table data itself. MRT is similar to the
`CDS <https://vizier.unistra.fr/doc/catstd.htx>`_ format standard, but differs
in the table description sections and the lack of a separate ``ReadMe`` file.
Astropy does not support writing in the CDS format.

The :class:`~astropy.io.ascii.Mrt` writer supports writing tables to MRT format.

.. note::

    The metadata of the table, apart from column ``unit``, ``name`` and
    ``description``, are not written in the output file. Placeholders for
    the title, authors, and table name fields are put into the output file and
    can be edited after writing.

Examples
""""""""

..
  EXAMPLE START
  Writing MRT Format Tables Using astropy.io.ascii

The command ``ascii.write(format='mrt')`` writes an ``astropy`` `~astropy.table.Table`
to the MRT format. Section dividers ``---`` and ``===`` are used to divide the table
into different sections, with the last section always been the actual data.

As the MRT standard requires,
for columns that have a ``unit`` attribute not set to ``None``,
the unit names are tabulated in the Byte-By-Byte
description of the column. When columns do not contain any units, ``---`` is put instead.
A ``?`` is prefixed to the column description in the Byte-By-Byte for ``Masked``
columns or columns that have null values, indicating them as such.

The example below initializes a table with columns that have a ``unit`` attribute and
has masked values.

  >>> from astropy.io import ascii
  >>> from astropy.table import Table, Column, MaskedColumn
  >>> from astropy import units as u
  >>> table = Table()
  >>> table['Name'] = ['ASASSN-15lh', 'ASASSN-14li']
  >>> # MRT Standard requires all quantities in SI units.
  >>> temperature = [0.0334, 0.297] * u.K
  >>> table['Temperature'] = temperature.to(u.keV, equivalencies=u.temperature_energy())
  >>> table['nH'] = Column([0.025, 0.0188], unit=u.Unit(10**22))
  >>> table['Flux'] = ([2.044 * 10**-11] * u.erg * u.cm**-2).to(u.Jy * u.Unit(10**12))
  >>> table['Flux'] = MaskedColumn(table['Flux'], mask=[True, False])
  >>> table['magnitude'] = [u.Magnitude(25), u.Magnitude(-9)]

Note that for columns with `~astropy.time.Time`, `~astropy.time.TimeDelta` and related values,
the writer does not do any internal conversion or modification. These columns should be
converted to regular columns with proper ``unit`` and ``name`` attribute before writing
the table. Thus::

  >>> from astropy.time import Time, TimeDelta
  >>> from astropy.timeseries import TimeSeries
  >>> ts = TimeSeries(time_start=Time('2019-01-01'), time_delta=2*u.day, n_samples=1)
  >>> table['Obs'] = Column(ts.time.decimalyear, description='Time of Observation')
  >>> table['Cadence'] = Column(TimeDelta(100.0, format='sec').datetime.seconds,
  ...                           unit=u.s)

Columns that are `~astropy.coordinates.SkyCoord` objects or columns with
values that are such objects are recognized as such, and some predefined labels and
description is used for them. Coordinate columns that have `~astropy.coordinates.SphericalRepresentation`
are additionally sub-divided into their coordinate component columns. Representations that have
``ra`` and ``dec`` components are divided into their ``hour``-``min``-``sec``
and ``deg``-``arcmin``-``arcsec`` components respectively. Whereas columns with
``SkyCoord`` objects in the ``Galactic`` or any of the ``Ecliptic`` frames are divided
into their latitude(``ELAT``/``GLAT``) and longitude components (``ELON``/``GLAT``) only.
The original table remains accessible as such, while the file is written from a modified
copy of the table. The new coordinate component columns are appended to the end of the table.

It should be noted that the default precision of the latitude, longitude and seconds (of arc)
columns is set at a default number of 12, 10 and 9 digits after the decimal for ``deg``, ``sec``
and ``arcsec`` values, respectively. This default is set to match a machine precision of 1e-15
relative to the original ``SkyCoord`` those columns were extracted from.
As all other columns, the format can be expliclty set by passing the ``formats`` keyword to the
``write`` function or by setting the ``format`` attribute of individual columns (the latter
will only work for columns that are not decomposed).
To customize the number of significant digits, presicions should therefore be specified in the
``formats`` dictionary for the *output* column names, such as
``formats={'RAs': '07.4f', 'DEs': '06.3f'}`` or ``formats={'GLAT': '+10.6f', 'GLON': '9.6f'}``
for milliarcsecond accuracy. Note that the forms with leading zeros for the seconds and
including the sign for latitudes are recommended for better consistency and readability.

The following code illustrates the above.

  >>> from astropy.coordinates import SkyCoord
  >>> table['coord'] = [SkyCoord.from_name('ASASSN-15lh'),
  ...                   SkyCoord.from_name('ASASSN-14li')]  # doctest: +REMOTE_DATA
  >>> table.write('coord_cols.dat', format='ascii.mrt')     # doctest: +SKIP
  >>> table['coord'] = table['coord'].geocentrictrueecliptic  # doctest: +REMOTE_DATA
  >>> table['Temperature'].format = '.5E' # Set default column format.
  >>> table.write('ecliptic_cols.dat', format='ascii.mrt')    # doctest: +SKIP

After execution, the contents of ``coords_cols.dat`` will be::

  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: table.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1-11  A11     ---    Name        Description of Name
  13-23  E11.6   keV    Temperature [0.0/0.01] Description of Temperature
  25-30  F6.4    10+22  nH          [0.01/0.03] Description of nH
  32-36  F5.3   10+12Jy Flux        ? Description of Flux
  38-42  E5.1    mag    magnitude   [0.0/3981.08] Description of magnitude
  44-49  F6.1    ---    Obs         [2019.0/2019.0] Time of Observation
  51-53  I3      s      Cadence     [100] Description of Cadence
  55-56  I2     h      RAh           Right Ascension (hour)
  58-59  I2     min    RAm           Right Ascension (minute)
  61-73  F13.10 s      RAs           Right Ascension (second)
     75  A1     ---    DE-           Sign of Declination
  76-77  I2     deg    DEd           Declination (degree)
  79-80  I2     arcmin DEm           Declination (arcmin)
  82-93  F12.9  arcsec DEs           Declination (arcsec)
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 2.87819e-09 0.0250       1e-10 2019.0 100 22 02 15.4500000000 -61 39 34.599996000
  ASASSN-14li 2.55935e-08 0.0188 2.044 4e+03 2019.0 100 12 48 15.2244072000 +17 46 26.496624000

And the file ``ecliptic_cols.dat`` will look like::

  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: table.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1- 11  A11     ---    Name        Description of Name
  13- 23  E11.6   keV    Temperature [0.0/0.01] Description of Temperature
  25- 30  F6.4    10+22  nH          [0.01/0.03] Description of nH
  32- 36  F5.3   10+12Jy Flux        ? Description of Flux
  38- 42  E5.1    mag    magnitude   [0.0/3981.08] Description of magnitude
  44- 49  F6.1    ---    Obs         [2019.0/2019.0] Time of Observation
  51- 53  I3      s      Cadence     [100] Description of Cadence
  55- 70  F16.12  deg    ELON        Ecliptic Longitude (geocentrictrueecliptic)
  72- 87  F16.12  deg    ELAT        Ecliptic Latitude (geocentrictrueecliptic)
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 2.87819e-09 0.0250       1e-10 2019.0 100 306.224208650096 -45.621789850825
  ASASSN-14li 2.55935e-08 0.0188 2.044 4e+03 2019.0 100 183.754980099243  21.051410763027

Finally, MRT has some specific naming conventions for columns
(`<https://journals.aas.org/mrt-labels/#reflab>`_). For example, if a column contains
the mean error for the data in a column named ``label``, then this column should be named ``e_label``.
These kinds of relative column naming cannot be enforced by the MRT writer
because it does not know what the column data means and thus, the relation between the
columns cannot be figured out. Therefore, it is up to the user to use ``Table.rename_columns``
to appropriately rename any columns before writing the table to MRT format.
The following example shows a similar situation, using the option to send the output to
``sys.stdout`` instead of a file::

  >>> table['error'] = [1e4, 450] * u.Jy  # Error in the Flux values.
  >>> outtab = table.copy()  # So that changes don't affect the original table.
  >>> outtab.rename_column('error', 'e_Flux')
  >>> # re-order so that related columns are placed next to each other.
  >>> outtab = outtab['Name', 'Obs', 'coord', 'Cadence', 'nH', 'magnitude',
  ...                 'Temperature', 'Flux', 'e_Flux']  # doctest: +REMOTE_DATA

  >>> ascii.write(outtab, format='mrt')  # doctest: +SKIP
  Title:
  Authors:
  Table:
  ================================================================================
  Byte-by-byte Description of file: table.dat
  --------------------------------------------------------------------------------
   Bytes Format Units  Label     Explanations
  --------------------------------------------------------------------------------
   1- 11  A11     ---    Name        Description of Name
  13- 18  F6.1    ---    Obs         [2019.0/2019.0] Time of Observation
  20- 22  I3      s      Cadence     [100] Description of Cadence
  24- 29  F6.4    10+22  nH          [0.01/0.03] Description of nH
  31- 35  E5.1    mag    magnitude   [0.0/3981.08] Description of magnitude
  37- 47  E11.6   keV    Temperature [0.0/0.01] Description of Temperature
  49- 53  F5.3   10+12Jy Flux        ? Description of Flux
  55- 61  F7.1    Jy     e_Flux      [450.0/10000.0] Description of e_Flux
  63- 78  F16.12  deg    ELON        Ecliptic Longitude (geocentrictrueecliptic)
  80- 95  F16.12  deg    ELAT        Ecliptic Latitude (geocentrictrueecliptic)
  --------------------------------------------------------------------------------
  Notes:
  --------------------------------------------------------------------------------
  ASASSN-15lh 2019.0 100 0.0250 1e-10 2.87819e-09       10000.0 306.224208650096 -45.621789850825
  ASASSN-14li 2019.0 100 0.0188 4e+03 2.55935e-08 2.044   450.0 183.754980099243  21.051410763027

..
  EXAMPLE END

.. attention::

    The MRT writer currently supports automatic writing of a single coordinate column
    in ``Tables``. For tables with more than one coordinate column of a given kind
    (e.g. equatorial, galactic or ecliptic), only the first found coordinate column
    will be decomposed into its component columns, and the rest of the coordinate
    columns of the same type will be converted to string columns. Thus users should take
    care that the additional coordinate columns are dealt with (e.g. by converting them
    to unique ``float``-valued columns) before using ``SkyCoord`` methods.
