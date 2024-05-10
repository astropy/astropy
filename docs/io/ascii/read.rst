.. include:: references.txt

.. _astropy.io.ascii_read:

Reading Tables
**************

The majority of commonly encountered ASCII tables can be read with the |read|
function::

  >>> from astropy.io import ascii
  >>> data = ascii.read(table)  # doctest: +SKIP

Here ``table`` is the name of a file, a string representation of a table, or a
list of table lines. The return value (``data`` in this case) is a :ref:`Table
<astropy-table>` object.

Help on the ``read()`` function arguments is available interactively as shown in
this example::

.. doctest-skip::

  >>> ascii.read.help()  # Common help for all formats
  >>> ascii.read.help("html")  # Common help plus "html" format-specific args

By default, |read| will try to `guess the table format <#guess-table-format>`_
by trying all of the supported formats.

.. Warning::

   Guessing the file format is often slow for large files because the reader
   tries parsing the file with every allowed format until one succeeds.
   For large files it is recommended to disable guessing with ``guess=False``.

..
  EXAMPLE START
  Reading ASCII Tables Using astropy.io.ascii

For unusually formatted tables where guessing does not work, give additional
hints about the format::

   >>> lines = ['objID                   & osrcid            & xsrcid       ',
   ...          '----------------------- & ----------------- & -------------',
   ...          '              277955213 & S000.7044P00.7513 & XS04861B6_005',
   ...          '              889974380 & S002.9051P14.7003 & XS03957B7_004']
   >>> data = ascii.read(lines, data_start=2, delimiter='&')
   >>> print(data)
     objID         osrcid          xsrcid
   --------- ----------------- -------------
   277955213 S000.7044P00.7513 XS04861B6_005
   889974380 S002.9051P14.7003 XS03957B7_004

Other examples are as follows::

   >>> data = astropy.io.ascii.read('data/nls1_stackinfo.dbout', data_start=2, delimiter='|')  # doctest: +SKIP
   >>> data = astropy.io.ascii.read('data/simple.txt', quotechar="'")  # doctest: +SKIP
   >>> data = astropy.io.ascii.read('data/simple4.txt', format='no_header', delimiter='|')  # doctest: +SKIP
   >>> data = astropy.io.ascii.read('data/tab_and_space.txt', delimiter=r'\s')  # doctest: +SKIP

If the format of a file is known (e.g., it is a fixed-width table or an IPAC
table), then it is more efficient and reliable to provide a value for the
``format`` argument from one of the values in the :ref:`supported_formats`. For
example::

   >>> data = ascii.read(lines, format='fixed_width_two_line', delimiter='&')

See the :ref:`guess_formats` section for additional details on format guessing.

..
  EXAMPLE END

For simpler formats such as CSV, |read| will automatically try reading with the
Cython/C parsing engine, which is significantly faster than the ordinary Python
implementation (described in :ref:`fast_ascii_io`). If the fast engine fails,
|read| will fall back on the Python reader by default. The argument
``fast_reader`` can be specified to control this behavior. For example, to
disable the fast engine::

   >>> data = ascii.read(lines, format='csv', fast_reader=False)

For reading very large tables see the section on :ref:`chunk_reading` or
use `pandas <https://pandas.pydata.org/>`_ (see Note below).

.. Note::

   Reading a table which contains unicode characters is supported with the
   pure Python readers by specifying the ``encoding`` parameter. The fast
   C-readers do not support unicode. For large data files containing unicode,
   we recommend reading the file using `pandas <https://pandas.pydata.org/>`_
   and converting to a :ref:`Table <astropy-table>` via the :ref:`Table -
   Pandas interface <pandas>`.

The |read| function accepts a number of parameters that specify the detailed
table format. Different formats can define different defaults, so the
descriptions below sometimes mention "typical" default values. This refers to
the :class:`~astropy.io.ascii.Basic` format reader and other similar character-separated formats.

.. _io_ascii_read_parameters:

Parameters for ``read()``
=========================

**table** : input table
  There are four ways to specify the table to be read:

  - Path to a file (string)
  - Single string containing all table lines separated by newlines
  - File-like object with a callable read() method
  - List of strings where each list element is a table line

  The first two options are distinguished by the presence of a newline in the
  string. This assumes that valid file names will not normally contain a
  newline, and a valid table input will at least contain two rows.
  Note that a table read in ``no_header`` format can legitimately consist
  of a single row; in this case passing the string as a list with a single
  item will ensure that it is not interpreted as a file name.

**format** : file format (default='basic')
  This specifies the top-level format of the ASCII table; for example,
  if it is a basic character delimited table, fixed format table, or
  a CDS-compatible table, etc. The value of this parameter must
  be one of the :ref:`supported_formats`.

**guess** : try to guess table format (default=None)
  If set to True, then |read| will try to guess the table format by cycling
  through a number of possible table format permutations and attempting to read
  the table in each case. See the `Guess table format`_ section for further details.

**delimiter** : column delimiter string
  A one-character string used to separate fields which typically defaults to
  the space character. Other common values might be "\\s" (whitespace), "," or
  "|" or "\\t" (tab). A value of "\\s" allows any combination of the tab and
  space characters to delimit columns.

**comment** : regular expression defining a comment line in table
  If the ``comment`` regular expression matches the beginning of a table line
  then that line will be discarded from header or data processing.  For the
  ``basic`` format this defaults to "\\s*#" (any whitespace followed by #).

**quotechar** : one-character string to quote fields containing special characters
  This specifies the quote character and will typically be either the single or
  double quote character. This is can be useful for reading text fields with
  spaces in a space-delimited table. The default is typically the double quote.

**header_start** : line index for the header line
  This includes only significant non-comment lines and counting starts at 0. If
  set to None this indicates that there is no header line and the column names
  will be auto-generated. See `Specifying header and data location`_ for more
  details.

**data_start** : line index for the start of data counting
  This includes only significant non-comment lines and counting starts at 0.
  See `Specifying header and data location`_ for more details.

**data_end** : line index for the end of data
  This includes only significant non-comment lines and can be negative to count
  from end. See `Specifying header and data location`_ for more details.

**encoding**: encoding to read the file (``default=None``)
  When `None` use `locale.getpreferredencoding` as an encoding. This matches
  the default behavior of the built-in `open` when no ``mode`` argument is
  provided.

**converters** : ``dict`` specifying output data types
  See the :ref:`io-ascii-read-converters` section for examples. Each key in the
  dictionary is a column name or else a name matching pattern including
  wildcards. The value is one of:

  - Python data type or numpy dtype such as ``int`` or ``np.float32``
  - list of such types which is tried in order until conversion is successful
  - list of converter tuples (this is not common, but see the
    `~astropy.io.ascii.convert_numpy` function for an example).

**names** : list of names corresponding to each data column
  Define the complete list of names for each data column. This will override
  names found in the header (if it exists). If not supplied then
  use names from the header or auto-generated names if there is no header.

**include_names** : list of names to include in output
  From the list of column names found from the header or the ``names``
  parameter, select for output only columns within this list. If not supplied,
  then include all names.

**exclude_names** : list of names to exclude from output
  Exclude these names from the list of output columns. This is applied *after*
  the ``include_names`` filtering. If not specified then no columns are excluded.

**fill_values** : list of fill value specifiers
  Specify input table entries which should be masked in the output table
  because they are bad or missing. See the `Bad or missing values`_ section
  for more information and examples. The default is that any blank table
  values are treated as missing.

**fill_include_names** : list of column names affected by ``fill_values``
  This is a list of column names (found from the header or the ``names``
  parameter) for all columns where values will be filled. `None` (the default) will
  apply ``fill_values`` to all columns.

**fill_exclude_names** : list of column names not affected by ``fill_values``
  This is a list of column names (found from the header or the ``names``
  parameter) for all columns where values will be **not** be filled.
  This parameter takes precedence over ``fill_include_names``.  A value
  of `None` (default) does not exclude any columns.

**outputter_cls** : Outputter class
  This converts the raw data tables value into the
  output object that gets returned by |read|. The default is
  :class:`~astropy.io.ascii.TableOutputter`, which returns a
  :class:`~astropy.table.Table` object (see :ref:`Data Tables <astropy-table>`).

**inputter_cls** : Inputter class
  This is generally not specified.

**data_splitter_cls** : Splitter class to split data columns

**header_splitter_cls** : Splitter class to split header columns

**fast_reader** : whether to use the C engine
  This can be ``True`` or ``False``, and also be a ``dict`` with options.
  (see :ref:`fast_ascii_io`)

Specifying Header and Data Location
===================================

The three parameters ``header_start``, ``data_start``, and ``data_end`` make it
possible to read a table file that has extraneous non-table data included.
This is a case where you need to help out `astropy.io.ascii` and tell it where
to find the header and data.

When a file is processed into a header and data components, any blank lines
(which might have whitespace characters) and commented lines (starting with the
comment character, typically ``#``) are stripped out *before* the header and
data parsing code sees the table content.

Example
-------

..
  EXAMPLE START
  Specifying Header and Data Location for ASCII Tables

To use the parameters ``header_start``, ``data_start``, and ``data_end``
to read a table with non-table data included, take the file below. The column
on the left is not part of the file but instead shows how `astropy.io.ascii` is
viewing each line and the line count index.  ::

  Index    Table content
  ------ ----------------------------------------------------------------
     -  | # This is the start of my data file
     -  |
     0  | Automatically generated by my_script.py at 2012-01-01T12:13:14
     1  | Run parameters: None
     2  | Column header line:
     -  |
     3  | x y z
     -  |
     4  | Data values section:
     -  |
     5  | 1 2 3
     6  | 4 5 6
     -  |
     7  | Run completed at 2012:01-01T12:14:01

In this case you would have ``header_start=3``, ``data_start=5``, and
``data_end=7``. The convention for ``data_end`` follows the normal Python
slicing convention where to select data rows 5 and 6 you would do
``rows[5:7]``. For ``data_end`` you can also supply a negative index to
count backward from the end, so ``data_end=-1`` (like ``rows[5:-1]``) would
work in this case.

..
  EXAMPLE END

.. _replace_bad_or_missing_values:

Bad or Missing Values
=====================

ASCII data tables can contain bad or missing values. A common case is when a
table contains blank entries with no available data.

Examples
--------

..
  EXAMPLE START
  ASCII Tables with Bad or Missing Values

Take this example of a table with blank entries::

  >>> weather_data = """
  ...   day,precip,type
  ...   Mon,1.5,rain
  ...   Tues,,
  ...   Wed,1.1,snow
  ...   """

By default, |read| will interpret blank entries as being bad/missing and output
a masked Table with those entries masked out by setting the corresponding mask
value set to ``True``::

  >>> dat = ascii.read(weather_data)
  >>> print(dat)
  day  precip type
  ---- ------ ----
   Mon    1.5 rain
  Tues     --   --
   Wed    1.1 snow

If you want to replace the masked (missing) values with particular values, set
the masked column ``fill_value`` attribute and then get the "filled" version of
the table. This looks like the following::

  >>> dat['precip'].fill_value = -999
  >>> dat['type'].fill_value = 'N/A'
  >>> print(dat.filled())
  day  precip type
  ---- ------ ----
   Mon    1.5 rain
  Tues -999.0  N/A
   Wed    1.1 snow

ASCII tables may have other indicators of bad or missing data as well. For
example, a table may contain string values that are not a valid representation
of a number (e.g., ``"..."``), or a table may have special values like ``-999``
that are chosen to indicate missing data. The |read| function has a flexible
system to accommodate these cases by marking specified character sequences in
the input data as "missing data" during the conversion process. Whenever
missing data is found the output will be a masked table.

This is done with the ``fill_values`` keyword argument, which can be set to a
single missing-value specification ``<missing_spec>`` or a list of ``<missing_spec>`` tuples::

  fill_values = <missing_spec> | [<missing_spec1>, <missing_spec2>, ...]
  <missing_spec> = (<match_string>, '0', <optional col name 1>, <optional col name 2>, ...)

When reading a table, the second element of a ``<missing_spec>`` should always
be the string ``'0'``, otherwise you may get unexpected behavior [#f1]_. By
default, the ``<missing_spec>`` is applied to all columns unless column name
strings are supplied. An alternate way to limit the columns is via the
``fill_include_names`` and ``fill_exclude_names`` keyword arguments in |read|.

In the example below we read back the weather table after filling the missing
values in with typical placeholders::

  >>> table = ['day   precip  type',
  ...          ' Mon     1.5  rain',
  ...          'Tues  -999.0   N/A',
  ...          ' Wed     1.1  snow']
  >>> t = ascii.read(table, fill_values=[('-999.0', '0', 'precip'), ('N/A', '0', 'type')])
  >>> print(t)
  day  precip type
  ---- ------ ----
   Mon    1.5 rain
  Tues     --   --
   Wed    1.1 snow

.. note::

   The default in |read| is ``fill_values=('','0')``. This marks blank entries as being
   missing for any data type (int, float, or string). If ``fill_values`` is explicitly
   set in the call to |read| then the default behavior of marking blank entries as missing
   no longer applies. For instance setting ``fill_values=None`` will disable this
   auto-masking without setting any other fill values. This can be useful for a string
   column where one of values happens to be ``""``.


.. [#f1] The requirement to put the ``'0'`` there is the legacy of an old
         interface which is maintained for backward compatibility and also to
         match the format of ``fill_value`` for reading with the format of
         ``fill_value`` used for writing tables. On reading, the second
         element of the ``<missing_spec>`` tuple can actually be an arbitrary
         string value which replaces occurrences of the ``<match_string>``
         string in the input stream prior to type conversion. This ends up
         being the value "behind the mask", which should never be directly
         accessed. Only the value ``'0'`` is neutral when attempting to detect
         the column data type and perform type conversion. For instance if you
         used ``'nan'`` for the ``<match_string>`` value then integer columns
         would wind up as float.

..
  EXAMPLE END

Selecting columns for masking
-----------------------------
The |read| function provides the parameters ``fill_include_names`` and ``fill_exclude_names``
to select which columns will be used in the ``fill_values`` masking process described above.

..
  EXAMPLE START
  Using the ``fill_include_names`` and ``fill_exclude_names`` parameters for ASCII tables

The use of these parameters is not common but in some cases can considerably simplify
the code required to read a table. The following gives a simple example to illustrate how
``fill_include_names`` and ``fill_exclude_names`` can be used
in the most basic and typical cases::

  >>> from astropy.io import ascii
  >>> lines = ['a,b,c,d', '1.0,2.0,3.0,4.0', ',,,']
  >>> ascii.read(lines)
  <Table length=2>
     a       b       c       d
  float64 float64 float64 float64
  ------- ------- ------- -------
      1.0     2.0     3.0     4.0
       --      --      --      --

  >>> ascii.read(lines, fill_include_names=['a', 'c'])
  <Table length=2>
     a     b      c     d
  float64 str3 float64 str3
  ------- ---- ------- ----
      1.0  2.0     3.0  4.0
       --           --

  >>> ascii.read(lines, fill_exclude_names=['a', 'c'])
  <Table length=2>
   a      b     c      d
  str3 float64 str3 float64
  ---- ------- ---- -------
   1.0     2.0  3.0     4.0
            --           --

..
  EXAMPLE END

.. _guess_formats:

Guess Table Format
==================

If the ``guess`` parameter in |read| is set to True, then
|read| will try to guess the table format by cycling through a number of
possible table format permutations and attempting to read the table in each
case. The first format which succeeds will be used to read the table. To
succeed, the table must be successfully parsed by the Reader and satisfy the
following column requirements:

 * At least two table columns.
 * No column names are a float or int number.
 * No column names begin or end with space, comma, tab, single quote, double
   quote, or a vertical bar (|).

These requirements reduce the chance for a false positive where a table is
successfully parsed with the wrong format. A common situation is a table
with numeric columns but no header row, and in this case `astropy.io.ascii` will
auto-assign column names because of the restriction on column names that
look like a number.

Guess Order
-----------

The order of guessing is shown by this Python code::

  for format in ("ecsv", "fixed_width_two_line", "rst", "fast_basic", "basic",
                 "fast_rdb", "rdb", "fast_tab", "tab", "cds", "daophot", "sextractor",
                 "ipac", "latex", "aastex"):
      read(format=format)

  for format in ("commented_header", "fast_basic", "basic", "fast_noheader", ""noheader"):
      for delimiter in ("|", ",", " ", "\\s"):
          for quotechar in ('"', "'"):
              read(format=format, delimiter=delimiter, quotechar=quotechar)

Note that the :class:`~astropy.io.ascii.FixedWidth` derived-readers are not
included in the default guess sequence (this causes problems), so to read such
tables you must explicitly specify the format with the ``format`` keyword. Also
notice that formats compatible with the fast reading engine attempt to use the
fast engine before the ordinary reading engine.

If none of the guesses succeed in reading the table (subject to the column
requirements), a final try is made using just the user-supplied parameters but
without checking the column requirements. In this way, a table with only one
column or column names that look like a number can still be successfully read.

The guessing process respects any values of the format, delimiter, and
quotechar parameters as well as options for the fast reader that were
supplied to the read() function. Any guesses that would conflict are
skipped. For example, the call::

 >>> data = ascii.read(table, format="no_header", quotechar="'")

would only try the four delimiter possibilities, skipping all the conflicting
format and quotechar combinations. Similarly, with any setting of
``fast_reader`` that requires use of the fast engine, only the fast
variants in the format list above will be tried.

Disabling
---------

Guessing can be disabled in two ways::

  import astropy.io.ascii
  data = astropy.io.ascii.read(table)               # guessing enabled by default
  data = astropy.io.ascii.read(table, guess=False)  # disable for this call
  astropy.io.ascii.set_guess(False)                 # set default to False globally
  data = astropy.io.ascii.read(table)               # guessing disabled

Debugging
---------

In order to get more insight into the guessing process and possibly debug if
something is not working as expected, use the
`~astropy.io.ascii.get_read_trace()` function. This returns a traceback of the
attempted read formats for the last call to `~astropy.io.ascii.read()`.

Comments and Metadata
=====================

Any comment lines detected during reading are inserted into the output table
via the ``comments`` key in the table's ``.meta`` dictionary.

Example
-------

..
  EXAMPLE START
  Comments and Metadata in ASCII Tables

Comment lines detected during reading are inserted into the output table as
such::

 >>> table='''# TELESCOPE = 30 inch
 ...          # TARGET = PV Ceph
 ...          # BAND = V
 ...          MJD mag
 ...          55555 12.3
 ...          55556 12.4'''
 >>> dat = ascii.read(table)
 >>> print(dat.meta['comments'])
 ['TELESCOPE = 30 inch', 'TARGET = PV Ceph', 'BAND = V']

While :mod:`astropy.io.ascii` will not do any post-processing on comment lines,
custom post-processing can be accomplished by rereading with the metadata line
comments. Here is one example, where comments are of the form "# KEY = VALUE"::

 >>> header = ascii.read(dat.meta['comments'], delimiter='=',
 ...                     format='no_header', names=['key', 'val'])
 >>> print(header)
    key      val
 --------- -------
 TELESCOPE 30 inch
    TARGET PV Ceph
      BAND       V

..
  EXAMPLE END

.. _io-ascii-read-converters:

Converters for Specifying Dtype
===============================

:mod:`astropy.io.ascii` converts the raw string values from the table into
numeric data types by using converter functions such as the Python ``int`` and
``float`` functions or numpy dtype types such as ``np.float64``.

The default converters are::

    default_converters = [int, float, str]

The default converters for each column can be overridden with the
``converters`` keyword::

  >>> import numpy as np
  >>> converters = {'col1': np.uint,
  ...               'col2': np.float32}
  >>> ascii.read('file.dat', converters=converters)  # doctest: +SKIP

In addition to single column names you can use wildcards via `fnmatch` to
select multiple columns. For example, we can set the format for all columns
with a name starting with "col" to an unsigned integer while applying default
converters to all other columns in the table::

  >>> import numpy as np
  >>> converters = {'col*': np.uint}
  >>> ascii.read('file.dat', converters=converters)  # doctest: +SKIP

..
  EXAMPLE START
  Reading True / False values as boolean type

The value in the converters ``dict`` can also be a list of types, in which case
these will be tried in order. This allows for flexible type conversions. For
example, imagine you get read the following table::

  >>> txt = """\
  ...   a   b    c
  ... --- --- -----
  ...   1 3.5  True
  ...   2 4.0 False"""
  >>> t = ascii.read(txt, format='fixed_width_two_line')

By default the ``True`` and ``False`` values will be interpreted as strings.
However, if you want those values to be read as booleans you can do the
following::

  >>> converters = {'*': [int, float, bool, str]}
  >>> t = ascii.read(txt, format='fixed_width_two_line', converters=converters)
  >>> print(t['c'].dtype)
  bool

..
  EXAMPLE END

Advanced usage
--------------

Internally type conversion uses the :func:`~astropy.io.ascii.convert_numpy`
function which returns a two-element tuple ``(converter_func, converter_type)``.
This two-element tuple can be used as the value in a ``converters`` dict.
The type provided to
:func:`~astropy.io.ascii.convert_numpy` must be a valid `NumPy type
<https://numpy.org/doc/stable/user/basics.types.html>`_ such as
``numpy.int``, ``numpy.uint``, ``numpy.int8``, ``numpy.int64``,
``numpy.float``, ``numpy.float64``, or ``numpy.str``.

It is also possible to directly pass an arbitrary conversion function as the
``converter_func`` element of the two-element tuple.

.. _fortran_style_exponents:

Fortran-Style Exponents
=======================

The :ref:`fast converter <fast_conversion_opts>` available with the C
input parser provides an ``exponent_style`` option to define a custom
character instead of the standard ``'e'`` for exponential formats in
the input file, to read, for example, Fortran-style double precision
numbers like ``'1.495978707D+13'``:

  >>> ascii.read('double.dat', format='basic', guess=False,
  ...            fast_reader={'exponent_style': 'D'})  # doctest: +SKIP

The special setting ``'fortran'`` is provided to allow for the
auto-detection of any valid Fortran exponent character (``'E'``,
``'D'``, ``'Q'``), as well as of triple-digit exponents prefixed with no
character at all (e.g., ``'2.1127123261674622-107'``).
All values and exponent characters in the input data are
case-insensitive; any value other than the default ``'E'`` implies the
automatic setting of ``'use_fast_converter': True``.

.. _advanced_customization:

Advanced Customization
======================

Here we provide a few examples that demonstrate how to extend the base
functionality to handle special cases. To go beyond these examples, the
best reference is to read the code for the existing
:ref:`extension_reader_classes`.

Examples
--------

..
  EXAMPLE START
  Advanced Customization to Extend Base Functionality of astropy.io.ascii

For special cases, these examples demonstrate how to extend the base
functionality of `astropy.io.ascii`.

**Define custom readers by class inheritance**

The most useful way to define a new reader class is by inheritance.
This is the way all of the built-in readers are defined, so there are plenty
of examples in the code.

In most cases, you will define one class to handle the header,
one class that handles the data, and a reader class that ties it all together.
Here is an example from the code that defines a reader that is just like
the basic reader, but header and data start in different lines of the file::

  # Note: NoHeader is already included in astropy.io.ascii for convenience.
  class NoHeaderHeader(BasicHeader):
      """Reader for table header without a header

      Set the start of header line number to `None`, which tells the basic
      reader there is no header line.
      """
      start_line = None

  class NoHeaderData(BasicData):
      """Reader for table data without a header

      Data starts at first uncommented line since there is no header line.
      """
      start_line = 0

  class NoHeader(Basic):
      """Read a table with no header line.  Columns are autonamed using
      header.auto_format which defaults to "col%d".  Otherwise this reader
      the same as the :class:`Basic` class from which it is derived.  Example::

        # Table data
        1 2 "hello there"
        3 4 world
      """
      _format_name = 'no_header'
      _description = 'Basic table with no headers'
      header_class = NoHeaderHeader
      data_class = NoHeaderData

In a slightly more involved case, the implementation can also override some of
the methods in the base class::

  # Note: CommentedHeader is already included in astropy.io.ascii for convenience.
  class CommentedHeaderHeader(BasicHeader):
      """Header class for which the column definition line starts with the
      comment character.  See the :class:`CommentedHeader` class  for an example.
      """
      def process_lines(self, lines):
          """Return only lines that start with the comment regexp.  For these
          lines strip out the matching characters."""
          re_comment = re.compile(self.comment)
          for line in lines:
              match = re_comment.match(line)
              if match:
                  yield line[match.end():]

      def write(self, lines):
          lines.append(self.write_comment + self.splitter.join(self.colnames))


  class CommentedHeader(Basic):
      """Read a file where the column names are given in a line that begins with
      the header comment character. ``header_start`` can be used to specify the
      line index of column names, and it can be a negative index (for example -1
      for the last commented line).  The default delimiter is the <space>
      character.::

        # col1 col2 col3
        # Comment line
        1 2 3
        4 5 6
      """
      _format_name = 'commented_header'
      _description = 'Column names in a commented line'

      header_class = CommentedHeaderHeader
      data_class = NoHeaderData


**Define a custom reader functionally**

Instead of defining a new class, it is also possible to obtain an instance
of a reader, and then to modify the properties of this one reader instance
in a function::

   def read_rdb_table(table):
       reader = astropy.io.ascii.Basic()
       reader.header.splitter.delimiter = '\t'
       reader.data.splitter.delimiter = '\t'
       reader.header.splitter.process_line = None
       reader.data.splitter.process_line = None
       reader.data.start_line = 2

       return reader.read(table)


**Create a custom splitter.process_val function**
::

   # The default process_val() normally just strips whitespace.
   # In addition have it replace empty fields with -999.
   def process_val(x):
       """Custom splitter process_val function: Remove whitespace at the beginning
       or end of value and substitute -999 for any blank entries."""
       x = x.strip()
       if x == '':
           x = '-999'
       return x

   # Create an RDB reader and override the splitter.process_val function
   rdb_reader = astropy.io.ascii.get_reader(reader_cls=astropy.io.ascii.Rdb)
   rdb_reader.data.splitter.process_val = process_val

..
  EXAMPLE END

.. _chunk_reading:

Reading Large Tables in Chunks
==============================

The default process for reading ASCII tables is not memory efficient and may
temporarily require much more memory than the size of the file (up to a factor
of 5 to 10). In cases where the temporary memory requirement exceeds available
memory this can cause significant slowdown when disk cache gets used.

In this situation, there is a way to read the table in smaller chunks which are
limited in size. There are two possible ways to do this:

- Read the table in chunks and aggregate the final table along the way. This
  uses only somewhat more memory than the final table requires.
- Use a Python generator function to return a `~astropy.table.Table` object for
  each chunk of the input table. This allows for scanning through arbitrarily
  large tables since it never returns the final aggregate table.

The chunk reading functionality is most useful for very large tables, so this is
available only for the :ref:`fast_ascii_io` readers. The following formats are
supported: ``tab``, ``csv``, ``no_header``, ``rdb``, and ``basic``. The
``commented_header`` format is not directly supported, but as a workaround one
can read using the ``no_header`` format and explicitly supply the column names
using the ``names`` argument.

In order to read a table in chunks you must provide the ``fast_reader`` keyword
argument with a ``dict`` that includes the ``chunk_size`` key with the value
being the approximate size (in bytes) of each chunk of the input table to read.
In addition, if you provide a ``chunk_generator`` key which is set to
``True``, then instead of returning a single table for the whole input it
returns an iterator that provides a table for each chunk of the input.

Examples
--------

..
  EXAMPLE START
  Reading Large Tables in Chunks with astropy.io.ascii

To read an entire table while limiting peak memory usage:
::

  # Read a large CSV table in 100 Mb chunks.

  tbl = ascii.read('large_table.csv', format='csv', guess=False,
                   fast_reader={'chunk_size': 100 * 1000000})

To read the table in chunks with an iterator, we iterate over a CSV table and
select all rows where the ``Vmag`` column is less than 8.0 (e.g., all stars in
table brighter than 8.0 mag). We collect all of these subtables and then stack
them at the end.
::

  from astropy.table import vstack

  # tbls is an iterator over the chunks (no actual reading done yet)
  tbls = ascii.read('large_table.csv', format='csv', guess=False,
                    fast_reader={'chunk_size': 100 * 1000000,
                                 'chunk_generator': True})

  out_tbls = []

  # At this point the file is actually read in chunks.
  for tbl in tbls:
      bright = tbl['Vmag'] < 8.0
      if np.count_nonzero(bright):
          out_tbls.append(tbl[bright])

  out_tbl = vstack(out_tbls)

.. Note:: **Performance**

  Specifying the ``format`` explicitly and using ``guess=False`` is a good idea
  for large tables. This prevents unnecessary guessing in the typical case
  where the format is already known.

  The ``chunk_size`` should generally be set to the largest value that is
  reasonable given available system memory. There is overhead associated
  with processing each chunk, so the fewer chunks the better.

  ..
    EXAMPLE END

.. _io_ascii_how_to_examples:

How to Find and Fix Problems Reading a Table
============================================

The purpose of this section is to provide a few examples how we can
deal with tables that fail to read.

Obtain the Data Table in a Different Format
-------------------------------------------
Sometimes it is easy to obtain the data in a more structured format that
more clearly defines columns and metadata, e.g. a FITS or VO/XML table, or
an ASCII table that uses a different column separator (e.g. comma instead of
white space) or fixed-width columns.
In that case, the fastest solution can be to simply download or export the
data again in a different format.

Find the Problem
----------------
Usually, `astropy.io.ascii.read` tries many different formats until one
succeeds in reading. If it works, that saves you from finding and
setting right options for reading. However, if it fails to find any combination
of format and format options that correctly parses the file, then you will get
a long exception message which shows every format that was tried and ends
with this advice::

  ************************************************************************
  ** ERROR: Unable to guess table format with the guesses listed above. **
  **                                                                    **
  ** To figure out why the table did not read, use guess=False and      **
  ** fast_reader=False, along with any appropriate arguments to read(). **
  ** In particular specify the format and any known attributes like the **
  ** delimiter.                                                         **
  ************************************************************************

To expand on this a bit, you probably know from looking at the file
what format it is in, which must be one of the :ref:`supported_formats`.
For instance maybe it is a basic space-delimited file but has the header
line as a comment like below, which corresponds to the ``commented_header``
format::

  >>> table = """# name id
  ... Jill 1232
  ... Jack Johnson 456"""

In order to find the actual problem with the reading this file, you would do::

  >>> ascii.read(table, format='commented_header', delimiter=' ', guess=False, fast_reader=False)
  Traceback (most recent call last):
    ...
  astropy.io.ascii.core.InconsistentTableError: Number of header columns (2) inconsistent with data columns (3) at data line 1
  Header values: ['name', 'id']
  Data values: ['Jack', 'Johnson', '456']

At this point you can see that the problem is that the 2nd data line has 3 columns while the header says there should be only 2. You might be initially confused by the ``data line 1`` since the problem was in the 3rd line of the file. There are two things happening here. First, ``data line 1`` refers to the count of data lines and does not include any header lines, blank lines, or commented out lines. Second, the count starts from zero, so that ``1`` is the 2nd data line.
See the :ref:`guess_formats` section for additional details on format guessing.


Make the Table Easier to Read
-----------------------------

Sometimes, the parameters for `astropy.io.ascii.read` to specify, for example
``format``, ``delimiter``, ``comment``, ``quote_char``, ``header_start``,
``data_start``, ``data_end``, and ``encoding`` are not enough.
To read just a single table that has a format close to, but not identical
with, any of the :ref:`supported_formats`, the fastest solution may be to open
that one table file in a text editor to modify it until it does conform to a
format that can be read. On the other hand, if we need to
read tables of that specific format again and again, it is better to find a way
to read them with `~astropy.io.ascii` without modifying every file by hand.


Badly formatted header line
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following table will fail to parse (raising an
`~astropy.io.ascii.InconsistentTableError`) because the header line looks as
if there were three columns, while in fact, there are only two::

  Name spectral type
  Vega A0
  Altair A7

Opening this file in a text editor to fix the format is easy::

  Name "spectral type"
  Vega A0
  Altair A7

or::

  Name spectral_type
  Vega A0
  Altair A7

With either of the above changes you can read the file with no problem using default settings.

..
  EXAMPLE START
  Make a table easier to read

To read the table without editing the files, we need to ignore the badly formatted header line and
pass in the names of the column ourselves.
That can be done without any modification of the table file by setting the ``data_start`` parameter::

   >>> table = """
   ... Star spectral type
   ... Vega A0
   ... Altair A7
   ... """
   >>> ascii.read(table, names=["Star", "spectral type"], data_start=1)
   <Table length=2>
    Star  spectral type
    str6       str2
   ------ -------------
     Vega            A0
   Altair            A7

..
  EXAMPLE END

Badly formatted data line
^^^^^^^^^^^^^^^^^^^^^^^^^

Similar principles apply to badly formatted data lines. Here is a
table where the number of columns is not consistent (``alpha Cen``
should be written as ``"alpha Cen"`` to make clear that the two words
"alpha" and "Cen" are part of the same column)::

  Star SpT
  Vega A0
  alpha Cen G2V+K1

When we try to read that with ``guess=False``, astropy throws an
`astropy.io.ascii.InconsistentTableError`::

  >>> from astropy.io import ascii
  >>> table = '''
  ... Star SpT
  ... Vega A0
  ... alpha Cen G2V+K1
  ... '''
  >>> ascii.read(table, guess=False)
  Traceback (most recent call last):
      ...
  astropy.io.ascii.core.InconsistentTableError: Number of header columns (2) inconsistent with data columns in data line 1

This points us to the line with the problem, here line 1 (starting to count
after the header lines and counting the data lines from 0 as usual in Python). In this table with
just two lines the problem is easy to spot, but for longer tables, the line number is
very helpful. We can now fix that line by hand in the file by adding quotes
around ``"alpha Cen"``. Then we can try to read the table again and see
if it works or if there is a another badly formatted data line.

.. _io-ascii-read-gaia-tables:

Reading Gaia Data Tables
========================

.. note::

    The recommended way to access Gaia is via its
    `astroquery.gaia <https://astroquery.readthedocs.io/en/latest/gaia/gaia.html>`_ module.
    However, if you need to access its data file separately via
    ``astropy``, then read on.

Gaia data tables are available in `ECSV
<https://github.com/astropy/astropy-APEs/blob/main/APE6.rst>`_ format including detailed
metadata for the tables and columns (e.g., column descriptions, units, and data types).
For example the DR3 tables are at http://cdn.gea.esac.esa.int/Gaia/gdr3/gaia_source/.

The DR3 data files are not strictly compliant with the ECSV standard because they use
the marker ``null`` to indicate a missing value instead of the required ``""``. In order
to read these files correctly with the full metadata, we need to tell the ECSV reader
to treat ``null`` as the missing value::

  >>> from astropy.table import QTable
  >>> dat = QTable.read(
  ...     "GaiaSource_000000-003111.csv.gz",
  ...     format="ascii.ecsv",
  ...     fill_values=("null", "0")
  ... )  # doctest: +SKIP
