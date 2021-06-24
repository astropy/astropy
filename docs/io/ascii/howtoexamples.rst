.. include:: references.txt

.. _io_ascii_how_to_examples:

How to fix a table that is not read or written correctly
********************************************************

The purpose of this section is to provide a few examples how we can
deal with tables that fail to read or cannot be written with the
default settings.

Obtain the data table in a different format
===========================================
Sometimes we deal with legacy ACSII table and the data only exisits in this
form. However, in other cases, we generate the ASCII table ourselves, for
example, we export it from a spreadsheet or some other application. In that
case, the easist solution is often to export it in a different format such
as fits or as different ASCII form, such as fixed-width instead of
comma-separated.


Making a table easier to read
=============================

Sometimes the Reader/Writer classes are not enough, even with all the
flexibility to customize, the delimiters, column names, format
converters and other properties. To read just a single table that has
a format close to, but not identical with, any of the `Supported formats`_,
the fastest solution is often to open that one table file
in a text editor to modify it until it does conform to a format that
any of the readers above can already read. On the other hand, if we need to
read table of that specific format again and again, it is better to find a
to read them with `~astropy.io.ascii` without modifying every file by hand.


Badly formatted header line
---------------------------
The following table will fail to parse (raising an
`~astropy.io.InconsistentTableError') because the header line looks as
if there were three columns, while in fact, there are only two::

  Name spectral type
  Vega A0
  Altair A7

Opening this file in a text editor to fix the format is easy::

  Name "spectral type"
  Vega A0
  Altair A7

or
::

  Name spectral_type
  Vega A0
  Altair A7

would both work, e.g.

..
  EXAMPLE START
  Make a table easier to read

::

  >>> from astropy.io import ascii
  >>> table = """
  ... Star "spectral type"
  ... Vega A0
  ... Altair A7
  ... """
  >>> ascii.read(table)
  <Table length=2>
   Star  spectral type
   str6       str2
  ------ -------------
    Vega            A0
  Altair            A7

Similarly, the header could be commented out and the names of the table column
can be set manually::

  >>> from astropy.io import ascii
  >>> table = """
  ... #Star spectral type
  ... Vega A0
  ... Altair A7
  ... """
  >>> ascii.read(table, names=["Star", "spectral type"], format='no_header')
  <Table length=2>
   Star  spectral type
   str6       str2
  ------ -------------
    Vega            A0
  Altair            A7

This last experiment shows us that reading works just fine, if we only
find a way to ignore the badly formatted header line. That can be done without
any modification of the table itself with the ``data_start`` parameter::

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
-------------------------

Similar principles apply to badly formatted data lines. Here is a
table where the number of columns is not consistent (`alpha Cen`
should be written as `"alpha Cen"` to make clear that the two words
"alpha" and "Cen" are part of the same column)::

  Star SpT
  Vega A0
  alpha Cen G2V+K1

When we try to read that, astropy throws an
`astropy.io.ascii.InconsistentTableError` and prints out a long list of formats
that it tried on this table and some advice (some output omitted)::

  >>> from astropy.io import ascii
  >>> table = '''
  ... Star SpT
  ... Vega A0
  ... alpha Cen G2V+K1
  ... '''
  >>> ascii.read(table)
  InconsistentTableError: Number of header columns (2) inconsistent with data columns (3) at data line 1
  Header values: ['Star', 'SpT']
  Data values: ['alpha', 'Cen', 'G2V+K1']
  [...]

This points us to the line with the problem, here line 1 (not counting the
header lines and starting to count at 0 as usual in Python). In this table with
just two lines that is easy to see, but for longer tables, the line number is
very helpful. We can now fix that line by hand in the file or just comment it
out (for most formats: putting `#` as the first character in the line) if we
do not really need that line. Then, we can try to read the table again and see
if it works or if there is a another badly formatted data line.

Change a written table
======================

..
  EXAMPLE START
  Write table with end of line string

If a single table needs to be written in a specific format that
astropy does not support, one of the easiest solutions is to write a
data file that is close to the desired format and then edit it by hand
in a text editor. Often, we can add extra columns or place useful
special values and then "search and replace all" in a text
editor. Here, we want to write a table that says "end of line" at the
end of every line. We can try to do that by simply adding another
column to our table that has the content we are looking for::

  >>> from astropy.table import Table
  >>> from astropy.io import ascii
  >>> table = Table({'a':[1, 2], 'b': ["some string", 4]})
  >>> ascii.write(table)
  a b
  1 3
  2 4
  >>> table['end of line'] = 'end of line'
  >>> ascii.write(table)
  a b "end of line"
  1 3 "end of line"
  2 4 "end of line"

We could now open this file in a text editor search and replace `"`
with ``.  Instead, our new column could have simply been a values that
does not occur in the table otherwise (e.g. "-99999") and then we
could search and replace "-9999" with "end of line".

..
  EXAMPLE END


Write our own reader or writer class
====================================

For more complicated cases or when more than a few tables in a
specific format need to be written, `astropy.io.ascii` can be
customized further by writing our own functions or classes. A few
examples for table reading are in the section
`advanced_customization`_.

Fix the program that wrote the table
====================================

If a table cannot be read or written in a desiered format, we can fix
the program used to read or write the table or at least report the
problem for the developers. While that might not solve our immidiate
problem, it will prevent this problem going forward for us and all
future users. If the table was written with astropy, please report it
to the developers: https://www.astropy.org/contribute.html
