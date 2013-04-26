.. _astropy.io.ascii_read:

.. include:: references.txt

Reading tables
--------------

The majority of commonly encountered ASCII tables can be easily read with the |read|
function::

  >>> from astropy.io import ascii
  >>> data = ascii.read(table)

where ``table`` is the name of a file, a string representation of a table, or a 
list of table lines.  By default |read| will try to `guess the table format <#guess-table-format>`_
by trying all the supported formats.  If this does not work (for unusually
formatted tables) then one needs give `astropy.io.ascii` additional hints about the
format, for example::

   >>> data = astropy.io.ascii.read('t/nls1_stackinfo.dbout', data_start=2, delimiter='|')
   >>> data = astropy.io.ascii.read('t/simple.txt', quotechar="'")
   >>> data = astropy.io.ascii.read('t/simple4.txt', Reader=ascii.NoHeader, delimiter='|')

The |read| function accepts a number of parameters that specify the detailed
table format.  Different Reader classes can define different defaults, so the
descriptions below sometimes mention "typical" default values.  This refers to
the :class:`~astropy.io.ascii.Basic` reader and other similar Reader classes.

Parameters for ``read()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**table** : input table 
  There are four ways to specify the table to be read:

  - Name of a file (string)
  - Single string containing all table lines separated by newlines
  - File-like object with a callable read() method
  - List of strings where each list element is a table line

  The first two options are distinguished by the presence of a newline in the string.  
  This assumes that valid file names will not normally contain a newline.

**Reader** : Reader class (default= :class:`~astropy.io.ascii.Basic`)
  This specifies the top-level format of the ASCII table, for example
  if it is a basic character delimited table, fixed format table, or
  a CDS-compatible table, etc.  The value of this parameter must
  be a Reader class.  For basic usage this means one of the 
  built-in :ref:`extension_reader_classes`.  

**guess**: try to guess table format (default=True)
  If set to True then |read| will try to guess the table format by cycling
  through a number of possible table format permuations and attemping to read
  the table in each case.  See the `Guess table format`_ section for further details.
  
**delimiter** : column delimiter string
  A one-character string used to separate fields which typically defaults to
  the space character.  Other common values might be "\\s" (whitespace), "," or
  "|" or "\\t" (tab).  A value of "\\s" allows any combination of the tab and
  space characters to delimit columns.

**comment** : regular expression defining a comment line in table
  If the ``comment`` regular expression matches the beginning of a table line then that line
  will be discarded from header or data processing.  For the :class:`~astropy.io.ascii.Basic` Reader this
  defaults to "\\s*#" (any whitespace followed by #).  

**quotechar** : one-character string to quote fields containing special characters
  This specifies the quote character and will typically be either the single or double
  quote character.  This is can be useful for reading text fields with spaces in a space-delimited
  table.  The default is typically the double quote.

**header_start** : line index for the header line not counting comment lines
  This specifies in the line index where the header line will be found.  Comment lines are
  not included in this count and the counting starts from 0 (first non-comment line has index=0).
  If set to None this indicates that there is no header line and the column names
  will be auto-generated.  The default is dependent on the Reader.

**data_start**: line index for the start of data not counting comment lines
  This specifies in the line index where the data lines begin where the counting starts
  from 0 and does not include comment lines.  The default is dependent on the Reader.

**data_end**: line index for the end of data (can be negative to count from end)
  If this is not None then it allows for excluding lines at the end that are not
  valid data lines.  A negative value means to count from the end, so -1 would 
  exclude the last line, -2 the last two lines, and so on.

**converters**: dict of data type converters
  See the `Converters`_ section for more information.

**names**: list of names corresponding to each data column
  Define the complete list of names for each data column.  This will override
  names found in the header (if it exists).  If not supplied then
  use names from the header or auto-generated names if there is no header.

**include_names**: list of names to include in output
  From the list of column names found from the header or the ``names``
  parameter, select for output only columns within this list.  If not supplied
  then include all names.
  
**exclude_names**: list of names to exlude from output
  Exclude these names from the list of output columns.  This is applied *after*
  the ``include_names`` filtering.  If not specified then no columns are excluded.

**fill_values**: fill value specifier of lists
  This can be used to fill missing values in the table or replace strings with special meaning.
  See the `Replace bad or missing values`_ section for more information and examples.

**fill_include_names**: list of column names, which are affected by ``fill_values``.
  If not supplied, then ``fill_values`` can affect all columns.

**fill_exclude_names**: list of column names, which are not affected by ``fill_values``.
  If not supplied, then ``fill_values`` can affect all columns.

**Outputter**: Outputter class
  This converts the raw data tables value into the
  output object that gets returned by |read|.  The default is
  :class:`~astropy.io.ascii.core.TableOutputter`, which returns a
  :class:`~astropy.table.Table` object.

**Inputter**: Inputter class
  This is generally not specified.

**data_Splitter**: Splitter class to split data columns

**header_Splitter**: Splitter class to split header columns

.. _replace_bad_or_missing_values:

Replace bad or missing values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ASCII data tables can contain bad or missing values.  A common case is when a table
contains blank entries with no available data, for example::

  day,rain,snow
  Mon,0.0,1.5
  Tues,,       # <-- Weather station down
  Wed,1.1,0.0

By default |read| will interpret blank entries as being bad/missing and output
a masked Table with those entries masked out by setting the corresponding mask
value set to ``True``.

ASCII tables may also have other indicators of bad or missing data.  For example a table
may contain string values that are not a valid representation of a number, e.g. ``"..."``,
or a table may have special values like ``-999`` that are chosen to indicate missing data.
The |read| function has a flexible system to accomodate these cases by replacing string
values in the input data before they are converted.  This is done with the ``fill_values``
argument which replaces ``<old>`` with ``<new>`` before the type conversion is
done::

  fill_values = <fill_spec> | [<fill_spec1>, <fill_spec2>, ...]
  <fill_spec> = (<old>, <new>, <optional col name 1>, <optional col name 2>, ...)

Within the ``<fill_spec>`` tuple the ``<old>`` and ``<new>`` values must be
strings.  These two values are then followed by zero or more column names.  If
column names are included the replacement is limited to those columns listed.
If no columns are specified then the replacement is done in every column,
subject to filtering by ``fill_include_names`` and ``fill_exclude_names`` (see
below).

The ``fill_values`` parameter in |read| takes a single ``<fill_spec>`` or a
list of ``<fill_spec>`` tuples.  If several ``<fill_spec>`` apply to a single
occurence of ``<old>`` then the first one determines the ``<new>`` value.  For
instance the following will replace an empty data value in the ``x`` or ``y``
columns with "1e38" while empty values in any other column will get "-999"::

  >>> ascii.read(table, fill_values=[('', '1e38', 'x', 'y'), ('', '-999')])

The following shows an example where string information needs to be exchanged before the
conversion to float values happens. Here ``no_rain`` and ``no_snow`` is replaced by
``0.0``::

  >>> table = ['day  rain     snow',    # column names
               #---  -------  --------
               'Mon  3.2      no_snow', 
               'Tue  no_rain  1.1', 
               'Wed  0.3      no_snow']
  >>> print ascii.read(table, fill_values=[('no_rain', '0.0'), ('no_snow', '0.0')])
  [('Mon', 3.2, --) ('Tue', --, 1.1) ('Wed', 0.3, --)]

Sometimes these rules apply only to specific columns in the table. Columns can be selected with
``fill_include_names`` or excluded with ``fill_exclude_names``. Also, column names can be
given directly with fill_values::

  >>> asciidata = ['text,no1,no2', 'text1,1,1.',',2,']
  >>> print ascii.read(asciidata, fill_values = ('', 'nan','no1','no2'), delimiter = ',')
  [('text1', 1, 1.0) ('', 2, --)]

Here, the empty value ``''`` in column ``no2`` is replaced by ``nan``, but the ``text``
column remains unaltered.

If any table elements match the fill specification then |read| returns a masked
`~astropy.table.Table` object with the corresponding elements masked out.

.. note::

   The default in |read| is ``fill_values=('','0')``.  This marks blank entries as being
   missing for any data type (int, float, or string).  If ``fill_values`` is explicitly
   set in the call to |read| then the default behavior of marking blank entries as missing
   no longer applies.  For instance setting ``fill_values=None`` will disable this
   auto-masking without setting any other fill values.  This can be useful for a string
   column where one of values happens to be ``""``.

Guess table format
^^^^^^^^^^^^^^^^^^^^^^
If the ``guess`` parameter in |read| is set to True (which is the default) then
|read| will try to guess the table format by cycling through a number of
possible table format permutations and attemping to read the table in each case.
The first format which succeeds and will be used to read the table. To succeed
the table must be successfully parsed by the Reader and satisfy the following
column requirements:

 * At least two table columns
 * No column names are a float or int number
 * No column names begin or end with space, comma, tab, single quote, double quote, or
   a vertical bar (|). 

These requirements reduce the chance for a false positive where a table is
successfully parsed with the wrong format.  A common situation is a table
with numeric columns but no header row, and in this case ``astropy.io.ascii`` will
auto-assign column names because of the restriction on column names that 
look like a number.

The order of guessing is shown by this Python code::
  
  for Reader in (Rdb, Tab, Cds, Daophot, SExtractor, Ipac):
      read(Reader=Reader)
  for Reader in (CommentedHeader, Basic, NoHeader):
      for delimiter in ("|", ",", " ", "\\s"):
          for quotechar in ('"', "'"):
              read(Reader=Reader, delimiter=delimiter, quotechar=quotechar)

Note that the :class:`~astropy.io.ascii.FixedWidth` derived-readers are not included
in the default guess sequence (this causes problems), so to read such tables
one must explicitly specify the reader class with the ``Reader`` keyword.

If none of the guesses succeed in reading the table (subject to the column
requirements) a final try is made using just the user-supplied parameters but
without checking the column requirements.  In this way a table with only one
column or column names that look like a number can still be successfully read.

The guessing process respects any values of the Reader, delimiter, and
quotechar parameters that were supplied to the read() function.  Any guesses
that would conflict are skipped.  For example the call::

 >>> data = astropy.io.ascii.read(table, Reader=NoHeader, quotechar="'")

would only try the four delimiter possibilities, skipping all the conflicting
Reader and quotechar combinations.

Guessing can be disabled in two ways::

  import astropy.io.ascii
  data = astropy.io.ascii.read(table)               # guessing enabled by default
  data = astropy.io.ascii.read(table, guess=False)  # disable for this call
  astropy.io.ascii.set_guess(False)                 # set default to False globally
  data = astropy.io.ascii.read(table)               # guessing disabled
  
Converters
^^^^^^^^^^^^^^

:mod:`astropy.io.ascii` converts the raw string values from the table into
numeric data types by using converter functions such as the Python ``int`` and
``float`` functions.  For example ``int("5.0")`` will fail while float("5.0")
will succeed and return 5.0 as a Python float.  

The default converters are::

    default_converters = [astropy.io.ascii.convert_numpy(numpy.int),
                          astropy.io.ascii.convert_numpy(numpy.float),
                          astropy.io.ascii.convert_numpy(numpy.str)]

These take advantage of the :func:`~astropy.io.ascii.core.convert_numpy`
function which returns a 2-element tuple ``(converter_func, converter_type)``
as described in the previous section.  The type provided to
:func:`~astropy.io.ascii.core.convert_numpy` must be a valid `numpy type
<http://docs.scipy.org/doc/numpy/user/basics.types.html>`_, for example
``numpy.int``, ``numpy.uint``, ``numpy.int8``, ``numpy.int64``,
``numpy.float``, ``numpy.float64``, ``numpy.str``.

The default converters for each column can be overridden with the
``converters`` keyword::

  >>> converters = {'col1': [astropy.io.ascii.convert_numpy(numpy.uint)],
                    'col2': [astropy.io.ascii.convert_numpy(numpy.float32)]}
  >>> ascii.read('file.dat', converters=converters)

Advanced customization
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we provide a few examples that demonstrate how to extend the base
functionality to handle special cases.  To go beyond these simple examples the
best reference is to read the code for the existing
:ref:`extension_reader_classes`.

**Define a custom reader functionally**
::

   def read_rdb_table(table):
       reader = astropy.io.ascii.Basic()
       reader.header.splitter.delimiter = '\t'
       reader.data.splitter.delimiter = '\t'
       reader.header.splitter.process_line = None  
       reader.data.splitter.process_line = None
       reader.data.start_line = 2

       return reader.read(table)

**Define custom readers by class inheritance**
::

   # Note: Tab and Rdb are already included in astropy.io.ascii for convenience.
   class Tab(astropy.io.ascii.Basic):
       def __init__(self):
           astropy.io.ascii.Basic.__init__(self)
           self.header.splitter.delimiter = '\t'
           self.data.splitter.delimiter = '\t'
           # Don't strip line whitespace since that includes tabs
           self.header.splitter.process_line = None  
           self.data.splitter.process_line = None
           # Don't strip data value spaces since that is significant in TSV tables
           self.data.splitter.process_val = None
           self.data.splitter.skipinitialspace = False

   class Rdb(astropy.io.ascii.Tab):
       def __init__(self):
           astropy.io.ascii.Tab.__init__(self)
           self.data.start_line = 2

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
   rdb_reader = astropy.io.ascii.get_reader(Reader=astropy.io.ascii.Rdb)
   rdb_reader.data.splitter.process_val = process_val
