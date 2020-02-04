.. include:: references.txt

.. _fixed_width_gallery:

Fixed-Width Gallery
*******************

Fixed-width tables are those where each column has the same width for every row
in the table. This is commonly used to make tables easy to read for humans or
Fortran codes. It also reduces issues with quoting and special characters,
for example::

  Col1   Col2    Col3 Col4
  ---- --------- ---- ----
   1.2   "hello"    1    a
   2.4 's worlds    2    2

There are a number of common variations in the formatting of fixed-width tables
which :mod:`astropy.io.ascii` can read and write. The most significant
difference is whether there is no header line (:class:`~astropy.io.ascii.FixedWidthNoHeader`), one
header line (:class:`~astropy.io.ascii.FixedWidth`), or two header lines
(:class:`~astropy.io.ascii.FixedWidthTwoLine`). Next, there are variations in
the delimiter character, like whether the delimiter appears on either end
("bookends"), or if there is padding around the delimiter.

Details are available in the class API documentation, but the easiest way to
understand all of the options and their interactions is by example.

Reading
=======

..
  EXAMPLE START
  Reading Fixed-Width Tables

FixedWidth
----------

**Nice, typical, fixed-format table:**
::

  >>> from astropy.io import ascii
  >>> table = """
  ... # comment (with blank line above)
  ... |  Col1  |  Col2   |
  ... |  1.2   | "hello" |
  ... |  2.4   |'s worlds|
  ... """
  >>> ascii.read(table, format='fixed_width')
  <Table length=2>
    Col1     Col2
  float64    str9
  ------- ---------
      1.2   "hello"
      2.4 's worlds

**Typical fixed-format table with col names provided:**
::

  >>> table = """
  ... # comment (with blank line above)
  ... |  Col1  |  Col2   |
  ... |  1.2   | "hello" |
  ... |  2.4   |'s worlds|
  ... """
  >>> ascii.read(table, format='fixed_width', names=['name1', 'name2'])
  <Table length=2>
   name1    name2
  float64    str9
  ------- ---------
      1.2   "hello"
      2.4 's worlds

**Weird input table with data values chopped by col extent:**
::

  >>> table = """
  ...   Col1  |  Col2 |
  ...   1.2       "hello"
  ...   2.4   sdf's worlds
  ... """
  >>> ascii.read(table, format='fixed_width')
  <Table length=2>
    Col1    Col2
  float64   str7
  ------- -------
      1.2    "hel
      2.4 df's wo

**Table with double delimiters:**
::

  >>> table = """
  ... || Name ||   Phone ||         TCP||
  ... |  John  | 555-1234 |192.168.1.10X|
  ... |  Mary  | 555-2134 |192.168.1.12X|
  ... |   Bob  | 555-4527 | 192.168.1.9X|
  ... """
  >>> ascii.read(table, format='fixed_width')
  <Table length=3>
  Name  Phone       TCP
  str4   str8      str12
  ---- -------- ------------
  John 555-1234 192.168.1.10
  Mary 555-2134 192.168.1.12
   Bob 555-4527  192.168.1.9

**Table with space delimiter:**
::

  >>> table = """
  ...  Name  --Phone-    ----TCP-----
  ...  John  555-1234    192.168.1.10
  ...  Mary  555-2134    192.168.1.12
  ...   Bob  555-4527     192.168.1.9
  ... """
  >>> ascii.read(table, format='fixed_width', delimiter=' ')
  <Table length=3>
  Name --Phone- ----TCP-----
  str4   str8      str12
  ---- -------- ------------
  John 555-1234 192.168.1.10
  Mary 555-2134 192.168.1.12
   Bob 555-4527  192.168.1.9

**Table with no header row and auto-column naming:**

Use ``header_start`` and ``data_start`` keywords to indicate no header line.
::

  >>> table = """
  ... |  John  | 555-1234 |192.168.1.10|
  ... |  Mary  | 555-2134 |192.168.1.12|
  ... |   Bob  | 555-4527 | 192.168.1.9|
  ... """
  >>> ascii.read(table, format='fixed_width',
  ...            header_start=None, data_start=0)
  <Table length=3>
  col1   col2       col3
  str4   str8      str12
  ---- -------- ------------
  John 555-1234 192.168.1.10
  Mary 555-2134 192.168.1.12
   Bob 555-4527  192.168.1.9

**Table with no header row and with col names provided:**

Second and third rows also have hanging spaces after final "|". Use
header_start and data_start keywords to indicate no header line.
::

  >>> table = ["|  John  | 555-1234 |192.168.1.10|",
  ...          "|  Mary  | 555-2134 |192.168.1.12|  ",
  ...          "|   Bob  | 555-4527 | 192.168.1.9|  "]
  >>> ascii.read(table, format='fixed_width',
  ...            header_start=None, data_start=0,
  ...            names=('Name', 'Phone', 'TCP'))
  <Table length=3>
  Name  Phone       TCP
  str4   str8      str12
  ---- -------- ------------
  John 555-1234 192.168.1.10
  Mary 555-2134 192.168.1.12
   Bob 555-4527  192.168.1.9


FixedWidthNoHeader
------------------

**Table with no header row and auto-column naming. Use the
``fixed_width_no_header`` format for convenience:**
::

  >>> table = """
  ... |  John  | 555-1234 |192.168.1.10|
  ... |  Mary  | 555-2134 |192.168.1.12|
  ... |   Bob  | 555-4527 | 192.168.1.9|
  ... """
  >>> ascii.read(table, format='fixed_width_no_header')
  <Table length=3>
  col1   col2       col3
  str4   str8      str12
  ---- -------- ------------
  John 555-1234 192.168.1.10
  Mary 555-2134 192.168.1.12
   Bob 555-4527  192.168.1.9

**Table with no delimiter with column start and end values specified:**

This uses the col_starts and col_ends keywords. Note that the
col_ends values are inclusive so a position range of zero to five
will select the first six characters.
::

  >>> table = """
  ... #    5   9     17  18      28    <== Column start / end indexes
  ... #    |   |       ||         |    <== Column separation positions
  ...   John   555- 1234 192.168.1.10
  ...   Mary   555- 2134 192.168.1.12
  ...    Bob   555- 4527  192.168.1.9
  ... """
  >>> ascii.read(table, format='fixed_width_no_header',
  ...                 names=('Name', 'Phone', 'TCP'),
  ...                 col_starts=(0, 9, 18),
  ...                 col_ends=(5, 17, 28),
  ...                 )
  <Table length=3>
  Name   Phone      TCP
  str4    str9     str10
  ---- --------- ----------
  John 555- 1234 192.168.1.
  Mary 555- 2134 192.168.1.
   Bob 555- 4527  192.168.1

**Table with no delimiter with only column start or end values specified:**

If only the col_starts keyword is given, it is assumed that each column
ends where the next column starts, and the final column ends at the same
position as the longest line of data.

Conversely, if only the col_ends keyword is given, it is assumed that the first
column starts at position zero and that each successive column starts
immediately after the previous one.

The two examples below read the same table and produce the same result.
::

  >>> table = """
  ... #1       9        19                <== Column start indexes
  ... #|       |         |                <== Column start positions
  ... #<------><--------><------------->  <== Inferred column positions
  ...   John   555- 1234 192.168.1.10
  ...   Mary   555- 2134 192.168.1.123
  ...    Bob   555- 4527  192.168.1.9
  ...    Bill  555-9875  192.255.255.255
  ... """
  >>> ascii.read(table,
  ...                 format='fixed_width_no_header',
  ...                 names=('Name', 'Phone', 'TCP'),
  ...                 col_starts=(1, 9, 19),
  ...                 )
  <Table length=4>
  Name   Phone         TCP
  str4    str9        str15
  ---- --------- ---------------
  John 555- 1234    192.168.1.10
  Mary 555- 2134   192.168.1.123
   Bob 555- 4527     192.168.1.9
  Bill  555-9875 192.255.255.255

  >>> ascii.read(table,
  ...                 format='fixed_width_no_header',
  ...                 names=('Name', 'Phone', 'TCP'),
  ...                 col_ends=(8, 18, 32),
  ...                 )
  <Table length=4>
  Name   Phone        TCP
  str4    str9       str14
  ---- --------- --------------
  John 555- 1234   192.168.1.10
  Mary 555- 2134  192.168.1.123
   Bob 555- 4527    192.168.1.9
  Bill  555-9875 192.255.255.25


FixedWidthTwoLine
-----------------

**Typical fixed-format table with two header lines with some cruft:**
::

  >>> table = """
  ...   Col1    Col2
  ...   ----  ---------
  ...    1.2xx"hello"
  ...   2.4   's worlds
  ... """
  >>> ascii.read(table, format='fixed_width_two_line')
  <Table length=2>
    Col1     Col2
  float64    str9
  ------- ---------
      1.2   "hello"
      2.4 's worlds

..
  EXAMPLE END

..
  EXAMPLE START
  Reading a reStructuredText Table

**reStructuredText table:**
::

  >>> table = """
  ... ======= ===========
  ...   Col1    Col2
  ... ======= ===========
  ...   1.2   "hello"
  ...   2.4   's worlds
  ... ======= ===========
  ... """
  >>> ascii.read(table, format='fixed_width_two_line',
  ...                 header_start=1, position_line=2, data_end=-1)
  <Table length=2>
    Col1     Col2
  float64    str9
  ------- ---------
      1.2   "hello"
      2.4 's worlds

..
  EXAMPLE END

**Text table designed for humans and test having position line before the header line:**
::

  >>> table = """
  ... +------+----------+
  ... | Col1 |   Col2   |
  ... +------|----------+
  ... |  1.2 | "hello"  |
  ... |  2.4 | 's worlds|
  ... +------+----------+
  ... """
  >>> ascii.read(table, format='fixed_width_two_line', delimiter='+',
  ...                 header_start=1, position_line=0, data_start=3, data_end=-1)
  <Table length=2>
    Col1     Col2
  float64    str9
  ------- ---------
      1.2   "hello"
      2.4 's worlds

Writing
=======

..
  EXAMPLE START
  Writing Fixed-Width Tables

FixedWidth
----------

**Define input values ``dat`` for all write examples:**
::

  >>> table = """
  ... | Col1 |  Col2     |  Col3 | Col4 |
  ... | 1.2  | "hello"   |  1    | a    |
  ... | 2.4  | 's worlds |  2    | 2    |
  ... """
  >>> dat = ascii.read(table, format='fixed_width')

**Write a table as a normal fixed-width table:**
::

  >>> ascii.write(dat, format='fixed_width')
  | Col1 |      Col2 | Col3 | Col4 |
  |  1.2 |   "hello" |    1 |    a |
  |  2.4 | 's worlds |    2 |    2 |

**Write a table as a fixed-width table with no padding:**
::

  >>> ascii.write(dat, format='fixed_width', delimiter_pad=None)
  |Col1|     Col2|Col3|Col4|
  | 1.2|  "hello"|   1|   a|
  | 2.4|'s worlds|   2|   2|

**Write a table as a fixed-width table with no bookend:**
::

  >>> ascii.write(dat, format='fixed_width', bookend=False)
  Col1 |      Col2 | Col3 | Col4
   1.2 |   "hello" |    1 |    a
   2.4 | 's worlds |    2 |    2

**Write a table as a fixed-width table with no delimiter:**
::

  >>> ascii.write(dat, format='fixed_width', bookend=False, delimiter=None)
  Col1       Col2  Col3  Col4
   1.2    "hello"     1     a
   2.4  's worlds     2     2

**Write a table as a fixed-width table with no delimiter and formatting:**
::

  >>> ascii.write(dat, format='fixed_width',
  ...                  formats={'Col1': '%-8.3f', 'Col2': '%-15s'})
  |     Col1 |            Col2 | Col3 | Col4 |
  | 1.200    | "hello"         |    1 |    a |
  | 2.400    | 's worlds       |    2 |    2 |

FixedWidthNoHeader
------------------

**Write a table as a normal fixed-width table:**
::

  >>> ascii.write(dat, format='fixed_width_no_header')
  | 1.2 |   "hello" | 1 | a |
  | 2.4 | 's worlds | 2 | 2 |

**Write a table as a fixed-width table with no padding:**
::

  >>> ascii.write(dat, format='fixed_width_no_header', delimiter_pad=None)
  |1.2|  "hello"|1|a|
  |2.4|'s worlds|2|2|

**Write a table as a fixed-width table with no bookend:**
::

  >>> ascii.write(dat, format='fixed_width_no_header', bookend=False)
  1.2 |   "hello" | 1 | a
  2.4 | 's worlds | 2 | 2

**Write a table as a fixed-width table with no delimiter:**
::

  >>> ascii.write(dat, format='fixed_width_no_header', bookend=False,
  ...                  delimiter=None)
  1.2    "hello"  1  a
  2.4  's worlds  2  2

FixedWidthTwoLine
-----------------

**Write a table as a normal fixed-width table:**
::

  >>> ascii.write(dat, format='fixed_width_two_line')
  Col1      Col2 Col3 Col4
  ---- --------- ---- ----
   1.2   "hello"    1    a
   2.4 's worlds    2    2

**Write a table as a fixed width table with space padding and '=' position_char:**
::

  >>> ascii.write(dat, format='fixed_width_two_line',
  ...                  delimiter_pad=' ', position_char='=')
  Col1        Col2   Col3   Col4
  ====   =========   ====   ====
   1.2     "hello"      1      a
   2.4   's worlds      2      2

**Write a table as a fixed-width table with no bookend:**
::

  >>> ascii.write(dat, format='fixed_width_two_line', bookend=True, delimiter='|')
  |Col1|     Col2|Col3|Col4|
  |----|---------|----|----|
  | 1.2|  "hello"|   1|   a|
  | 2.4|'s worlds|   2|   2|

..
  EXAMPLE END
