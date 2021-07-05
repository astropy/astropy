# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

cds.py:
  Classes to read CDS / Vizier table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu), \
         Suyog Garg (suyog7130@gmail.com)
"""


import fnmatch
import itertools
import re
import os
import math
import numpy as np
from contextlib import suppress

from . import core
from . import fixedwidth

from astropy.units import Unit

from io import StringIO
from astropy.table import Table
from astropy.table import Column, MaskedColumn
from string import Template
from textwrap import wrap, fill

MAX_SIZE_README_LINE = 80
MAX_COL_INTLIMIT = 100000


__doctest_skip__ = ['*']

cdsdicts = {'title': 'Title ?',
            'author': '1st author ?',
            'catalogue': '',
            'date': 'Date ?',
            'abstract': 'Abstract ?',
            'authors': 'Authors ?',
            'bibcode': 'ref ?',
            'keywords': '',
            'tableDescription': '',
            'seealso': ''
            }

ByteByByteTemplate = ["Byte-by-byte Description of file: $file",
"--------------------------------------------------------------------------------",
" Bytes Format Units  Label     Explanations",
"--------------------------------------------------------------------------------",
"$bytebybyte",
"--------------------------------------------------------------------------------"]

ReadMeTemplate = os.path.dirname(os.path.realpath(__file__)) + "/src/ReadMe.template"


class CdsSplitter(fixedwidth.FixedWidthSplitter):
    """
    Contains the join function to left align the CDS columns
    when writing to a file.
    """
    def join(self, vals, widths):
        vals = [val + ' ' * (width - len(val)) for val, width in zip(vals, widths)]
        return self.delimiter.join(vals)


class CdsHeader(core.BaseHeader):
    col_type_map = {'e': core.FloatType,
                    'f': core.FloatType,
                    'i': core.IntType,
                    'a': core.StrType}

    'The ReadMe file to construct header from.'
    readme = None

    def get_type_map_key(self, col):
        match = re.match(r'\d*(\S)', col.raw_type.lower())
        if not match:
            raise ValueError('Unrecognized CDS format "{}" for column "{}"'.format(
                col.raw_type, col.name))
        return match.group(1)

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines`` for a CDS
        header.

        Parameters
        ----------
        lines : list
            List of table lines

        """

        # Read header block for the table ``self.data.table_name`` from the read
        # me file ``self.readme``.
        if self.readme and self.data.table_name:
            in_header = False
            readme_inputter = core.BaseInputter()
            f = readme_inputter.get_lines(self.readme)
            # Header info is not in data lines but in a separate file.
            lines = []
            comment_lines = 0
            for line in f:
                line = line.strip()
                if in_header:
                    lines.append(line)
                    if line.startswith(('------', '=======')):
                        comment_lines += 1
                        if comment_lines == 3:
                            break
                else:
                    match = re.match(r'Byte-by-byte Description of file: (?P<name>.+)$',
                                     line, re.IGNORECASE)
                    if match:
                        # Split 'name' in case in contains multiple files
                        names = [s for s in re.split('[, ]+', match.group('name'))
                                 if s]
                        # Iterate on names to find if one matches the tablename
                        # including wildcards.
                        for pattern in names:
                            if fnmatch.fnmatch(self.data.table_name, pattern):
                                in_header = True
                                lines.append(line)
                                break

            else:
                raise core.InconsistentTableError("Can't find table {} in {}".format(
                    self.data.table_name, self.readme))

        found_line = False

        for i_col_def, line in enumerate(lines):
            if re.match(r'Byte-by-byte Description', line, re.IGNORECASE):
                found_line = True
            elif found_line:  # First line after list of file descriptions
                i_col_def -= 1  # Set i_col_def to last description line
                break

        re_col_def = re.compile(r"""\s*
                                    (?P<start> \d+ \s* -)? \s*
                                    (?P<end>   \d+)        \s+
                                    (?P<format> [\w.]+)     \s+
                                    (?P<units> \S+)        \s+
                                    (?P<name>  \S+)
                                    (\s+ (?P<descr> \S.*))?""",
                                re.VERBOSE)

        cols = []
        for line in itertools.islice(lines, i_col_def + 4, None):
            if line.startswith(('------', '=======')):
                break
            match = re_col_def.match(line)
            if match:
                col = core.Column(name=match.group('name'))
                col.start = int(re.sub(r'[-\s]', '',
                                       match.group('start') or match.group('end'))) - 1
                col.end = int(match.group('end'))
                unit = match.group('units')
                if unit == '---':
                    col.unit = None  # "---" is the marker for no unit in CDS table
                else:
                    col.unit = Unit(unit, format='cds', parse_strict='warn')
                col.description = (match.group('descr') or '').strip()
                col.raw_type = match.group('format')
                col.type = self.get_col_type(col)

                match = re.match(
                    r'(?P<limits>[\[\]] \S* [\[\]])?'  # Matches limits specifier (eg [])
                                                       # that may or may not be present
                    r'\?'  # Matches '?' directly
                    r'((?P<equal>=)(?P<nullval> \S*))?'  # Matches to nullval if and only
                                                         # if '=' is present
                    r'(?P<order>[-+]?[=]?)'  # Matches to order specifier:
                                             # ('+', '-', '+=', '-=')
                    r'(\s* (?P<descriptiontext> \S.*))?',  # Matches description text even
                                                           # even if no whitespace is
                                                           # present after '?'
                    col.description, re.VERBOSE)
                if match:
                    col.description = (match.group('descriptiontext') or '').strip()
                    if issubclass(col.type, core.FloatType):
                        fillval = 'nan'
                    else:
                        fillval = '0'

                    if match.group('nullval') == '-':
                        col.null = '---'
                        # CDS tables can use -, --, ---, or ---- to mark missing values
                        # see https://github.com/astropy/astropy/issues/1335
                        for i in [1, 2, 3, 4]:
                            self.data.fill_values.append(('-' * i, fillval, col.name))
                    else:
                        col.null = match.group('nullval')
                        if (col.null is None):
                            col.null = ''
                        self.data.fill_values.append((col.null, fillval, col.name))

                cols.append(col)
            else:  # could be a continuation of the previous col's description
                if cols:
                    cols[-1].description += line.strip()
                else:
                    raise ValueError(f'Line "{line}" not parsable as CDS header')

        self.names = [x.name for x in cols]
        self.cols = cols

    def __splitFloatFormat(self, value):
        """
        Splits a Float string into different parts to find number
        of digits after decimal and check if the value is in Scientific
        notation.

        Parameters
        ----------
        value : the float value

        Returns
        -------
        fmt: (int, int, int, bool, bool)
            List of values describing the Float sting.
            (size, dec, ent, sign, exp)
            size, length of the given string.
            ent, number of digits before decimal point.
            dec, number of digits after decimal point.
            sign, whether or not given value signed.
            exp, is value in Scientific notation?
        """
        regfloat = re.compile(r"""(?P<sign> [+-]*)
                                  (?P<ent> [^eE.]+)
                                  (?P<deciPt> [.]*)
                                  (?P<decimals> [0-9]*)
                                  (?P<exp> [eE]*-*)[0-9]*""",
                             re.VERBOSE)
        mo = regfloat.match(value)

        if mo is None:
            raise Exception(value + " is not a float number")
        return (len(value),
                len(mo.group('ent')),
                len(mo.group('decimals')),
                mo.group('sign') != "",
                mo.group('exp') != "")

    def __set_column_val_limits(self, col):
        """
        Sets the `col.min` and `col.max` column attributes,
        taking into account columns with Null values.
        """
        col.max = max(col)
        col.min = min(col)
        if col.max is np.ma.core.MaskedConstant:
            col.max = None
        if col.min is np.ma.core.MaskedConstant:
            col.min = None

    def columnFloatFormatter(self, col):
        """
        String formatter for a column containing Float values.
        """
        # maxSize: maximum length of string containing the float value.
        # maxEnt: maximum number of digits places before decimal point.
        # maxDec: maximum number of digits places after decimal point.
        # maxPrec: maximum precision of the column values, sum of maxEnt and maxDec.
        maxSize, maxPrec, maxEnt, maxDec = 1, 0, 1, 0
        sign = False
        fformat = 'F'
        #fmt = [0, 0, 0, 0, False]

        # find max sized value in the col
        for rec in col:
            # skip null values
            if rec is None:
                continue

            # find format of the Float string
            fmt = self.__splitFloatFormat(str(rec))

            # if value is in Scientific notation
            if fmt[4] is True:
                # if the previous column value was in normal Float format
                # set maxSize, maxPrec and maxDec to default.
                if fformat == 'F':
                    maxSize, maxPrec, maxDec = 1, 0, 0
                # designate the column to be in Scientific notation
                fformat = 'E'
            else:
                # move to next column value if
                # current value is not in Scientific notation
                # but the column is designated as such because
                # one of the previous values was.
                if fformat == 'E': continue

            if maxSize < fmt[0]: maxSize = fmt[0]
            if maxEnt < fmt[1]: maxEnt = fmt[1]
            if maxDec < fmt[2]: maxDec = fmt[2]
            if fmt[3]: sign = True

            if maxPrec < fmt[1] + fmt[2]: maxPrec = fmt[1] + fmt[2]

        if fformat == 'E':
            col.meta.size = maxSize
            if sign: col.meta.size += 1
            # number of digits after decimal is replaced by the precision
            # for values in Scientific notation, when writing they Format.
            col.fortran_format = fformat + str(col.meta.size) + "." + str(maxPrec)
            col.format = str(col.meta.size) + "." + str(maxDec) + "e"
        else:
            col.meta.size = maxEnt + maxDec + 1
            if sign: col.meta.size += 1
            col.fortran_format = fformat + str(col.meta.size) + "." + str(maxDec)
            col.format = col.fortran_format[1:] + "f"

    def writeByteByByte(self):
        """
        Writes the Byte-By-Byte description of the table.

        See: http://vizier.u-strasbg.fr/doc/catstd-3.1.htx

        Example::

            --------------------------------------------------------------------------------
            Byte-by-byte Description of file: table.dat
            --------------------------------------------------------------------------------
             Bytes Format Units  Label     Explanations
            --------------------------------------------------------------------------------
             1- 8   A8     ---    names    Description of names
            10-14   E5.1   ---    e       [0.0/0.01]? Description of e
            16-18   F3.0   ---    d       ? Description of d
            20-26   E7.1   ---    s       [-9e+34/2.0] Description of s
            28-30   I3     ---    i       [-30/67] Description of i
            32-34   F3.1   ---    sameF   [5.0/5.0] Description of sameF
            36-37   I2     ---    sameI   [20] Description of sameI

            --------------------------------------------------------------------------------
        """
        # get column widths.
        vals_list = []
        col_str_iters = self.data.str_vals()
        for vals in zip(*col_str_iters):
            vals_list.append(vals)

        for i, col in enumerate(self.cols):
            col.width = max([len(vals[i]) for vals in vals_list])
            if self.start_line is not None:
                col.width = max(col.width, len(col.info.name))
        widths = [col.width for col in self.cols]

        startb = 1
        # set default width of Start Byte, End Byte, Format and
        # Label columns in the ByteByByte table.
        sz = [0, 0, 1, 7]
        l = len(str(sum(widths)))
        if l > sz[0]:
            sz[0] = l
            sz[1] = l
        for column in self.cols:
            if len(column.name) > sz[3]:
                sz[3] = len(column.name)

        buff = ""
        maxDescripSize = 16
        nsplit = sum(sz) + 16
        # format string for Start Byte and End Byte
        fmtb = "{0:" + str(sz[0]) + "d}-{1:" + str(sz[1]) + "d} {2:" + str(sz[2]) + "s}"

        bbb = Table(names=['Bytes', 'Format', 'Units', 'Label', 'Explanations'],
                    dtype=[str]*5)

        for i, col in enumerate(self.cols):
            # check if column is MaskedColumn
            col.hasNull = isinstance(col, MaskedColumn)

            # set CDSColumn type, size and format.
            if np.issubdtype(col.dtype, np.int):
                # integer formatter
                self.__set_column_val_limits(col)
                col.meta.size = len(str(col.max))
                l = len(str(col.min))
                if col.meta.size < l: col.meta.size = l
                col.fortran_format = "I" + str(col.meta.size)
                col.format = ">" + col.fortran_format[1:]

            elif np.issubdtype(col.dtype, np.float):
                # float formatter
                self.__set_column_val_limits(col)
                self.columnFloatFormatter(col)

            else:
                # string formatter
                if col.hasNull:
                    mcol = col
                    mcol.fill_value = ""
                    coltmp = Column(mcol.filled(), dtype=str)
                    col.meta.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', coltmp.dtype.str))
                else:
                    col.meta.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', col.dtype.str))
                col.fortran_format = "A" + str(col.meta.size)
                col.format = str(col.meta.size) + "s"

            endb = col.meta.size + startb - 1

            # set column description
            if col.description is not None:
                description = col.description
            else:
                description = "Description of " + col.name

            # set null flag in column description
            if col.hasNull:
                nullflag = "?"
            else:
                nullflag = ""

            # set column unit
            if col.unit is not None:
                col.meta.unit = col.unit.to_string("cds")
            elif col.name.lower().find("magnitude") > -1:
                # ``col.unit`` will still be ``None``, if the unit of column values
                # is ``Magnitude``, because ``astropy.units.Magnitude`` is actually a class.
                # Unlike other units which are instances of ``astropy.units.Unit``,
                # application of the ``Magnitude`` unit calculates the logarithm
                # of the values. Thus, the only way to check for if the column values
                # have ``Magnitude`` unit is to check the column name.
                col.meta.unit = "mag"
            else:
                col.meta.unit = "---"

            # add col limit values to col description
            limVals = ""
            if col.min and col.max:
                if col.fortran_format[0] == 'I':
                    if abs(col.min) < MAX_COL_INTLIMIT and abs(col.max) < MAX_COL_INTLIMIT:
                        if col.min == col.max:
                            limVals = "[{0}]".format(col.min)
                        else:
                            limVals = "[{0}/{1}]".format(col.min, col.max)
                elif col.fortran_format[0] in ('E','F'):
                    limVals = "[{0}/{1}]".format(math.floor(col.min*100)/100.,
                                                 math.ceil(col.max*100)/100.)

            description = "{0}{1} {2}".format(limVals, nullflag, description)

            # find max description length
            if len(description) > maxDescripSize:
                maxDescripSize = len(description)

            # add ByteByByte row to bbb table
            bbb.add_row([fmtb.format(startb, endb, ""),
                         "" if col.fortran_format is None else col.fortran_format,
                         "" if col.meta.unit is None else col.meta.unit,
                         "" if col.name is None else col.name,
                         description])
            startb = endb + 2

        # properly format bbb columns
        bbbLines = StringIO()
        bbb.write(bbbLines, format='ascii.fixed_width_no_header',
                    delimiter=' ', bookend=False, delimiter_pad=None,
                    formats={'Format':'<6s',
                             'Units':'<6s',
                             'Label':'<'+str(sz[3])+'s',
                             'Explanations':''+str(maxDescripSize)+'s'})

        # get formatted bbb lines
        bbbLines = bbbLines.getvalue().splitlines()

        # wrap line if it is too long
        for newline in bbbLines:
            if len(newline) > MAX_SIZE_README_LINE:
                buff += ("\n").join(wrap(newline,
                                            subsequent_indent=" " * nsplit,
                                            width=MAX_SIZE_README_LINE))
                buff += "\n"
            else:
                buff += newline + "\n"

        # last value of ``endb`` is the max column width after formatting.
        self.linewidth = endb

        # add column notes to ByteByByte
        notes = self.cdsdicts.get('notes', None)
        if notes is not None:
            buff += "-" * MAX_SIZE_README_LINE + "\n"
            for line in notes:
                buff += line + "\n"
            buff += "-" * MAX_SIZE_README_LINE + "\n"
        return buff

    def write(self, lines):
        """
        Writes the Header of the CDS table, aka ReadMe, which
        also contains the ByteByByte description of the table.
        """
        # get ByteByByte description and fill the template
        bbbTemplate = Template('\n'.join(ByteByByteTemplate))
        ByteByByte = bbbTemplate.substitute({'file': 'table.dat',
                                    'bytebybyte': self.writeByteByByte()})

        #-- get index of files --#
        # set width of FileName and Lrecl columns
        sz = [14, 0]
        l = len(str(self.linewidth))
        if l > sz[1]:
            sz[1] = l

        # create the File Index table
        fIndexRows = (["ReadMe",
                       MAX_SIZE_README_LINE,
                       ".",
                       "this file"],
                      ['table',
                       self.linewidth,
                       str(len(self.cols[0])),
                       self.cdsdicts['tableDescription']])
        filesIndex = Table(names=['FileName', 'Lrecl', 'Records', 'Explanations'],
                           rows=fIndexRows)

        # get File Index table rows as formatted lines
        fIndexLines = StringIO()
        filesIndex.write(fIndexLines, format='ascii.fixed_width_no_header',
                            delimiter=' ', bookend=False, delimiter_pad=None,
                            formats={'FileName': str(sz[0])+'s',
                                     'Lrecl': ''+str(sz[1])+'d',
                                     'Records': '>8s',
                                     'Explanations': 's'})
        fIndexLines = fIndexLines.getvalue()

        # fill up the full ReadMe
        with open(ReadMeTemplate) as rmf:
            rmTemplate = Template(rmf.read())
        rmTempVals = self.cdsdicts
        rmTempVals.update({'tablesIndex': fIndexLines,
                           'bytebybyte': ByteByByte})
        ReadMe = rmTemplate.substitute(rmTempVals)
        ReadMe += "-" * MAX_SIZE_README_LINE
        lines.append(ReadMe)


class CdsData(fixedwidth.FixedWidthData):
    """CDS table data reader
    """
    splitter_class = CdsSplitter

    def process_lines(self, lines):
        """Skip over CDS header by finding the last section delimiter"""
        # If the header has a ReadMe and data has a filename
        # then no need to skip, as the data lines do not have header
        # info. The ``read`` method adds the table_name to the ``data``
        # attribute.
        if self.header.readme and self.table_name:
            return lines
        i_sections = [i for i, x in enumerate(lines)
                      if x.startswith(('------', '======='))]
        if not i_sections:
            raise core.InconsistentTableError('No CDS section delimiter found')
        return lines[i_sections[-1]+1:]  # noqa

    def write(self, lines):
        self.splitter.delimiter = ' '
        fixedwidth.FixedWidthData.write(self, lines)


class Cds(core.BaseReader):
    """CDS format table.

    See: http://vizier.u-strasbg.fr/doc/catstd.htx

    Example::

      Table: Table name here
      = ==============================================================================
      Catalog reference paper
          Bibliography info here
      ================================================================================
      ADC_Keywords: Keyword ; Another keyword ; etc

      Description:
          Catalog description here.
      ================================================================================
      Byte-by-byte Description of file: datafile3.txt
      --------------------------------------------------------------------------------
         Bytes Format Units  Label  Explanations
      --------------------------------------------------------------------------------
         1-  3 I3     ---    Index  Running identification number
         5-  6 I2     h      RAh    Hour of Right Ascension (J2000)
         8-  9 I2     min    RAm    Minute of Right Ascension (J2000)
        11- 15 F5.2   s      RAs    Second of Right Ascension (J2000)
      --------------------------------------------------------------------------------
      Note (1): A CDS file can contain sections with various metadata.
                Notes can be multiple lines.
      Note (2): Another note.
      --------------------------------------------------------------------------------
        1 03 28 39.09
        2 04 18 24.11

    **About parsing the CDS format**

    The CDS format consists of a table description and the table data.  These
    can be in separate files as a ``ReadMe`` file plus data file(s), or
    combined in a single file.  Different subsections within the description
    are separated by lines of dashes or equal signs ("------" or "======").
    The table which specifies the column information must be preceded by a line
    starting with "Byte-by-byte Description of file:".

    In the case where the table description is combined with the data values,
    the data must be in the last section and must be preceded by a section
    delimiter line (dashes or equal signs only).

    **Basic usage**

    Use the ``ascii.read()`` function as normal, with an optional ``readme``
    parameter indicating the CDS ReadMe file.  If not supplied it is assumed that
    the header information is at the top of the given table.  Examples::

      >>> from astropy.io import ascii
      >>> table = ascii.read("data/cds.dat")
      >>> table = ascii.read("data/vizier/table1.dat", readme="data/vizier/ReadMe")
      >>> table = ascii.read("data/cds/multi/lhs2065.dat", readme="data/cds/multi/ReadMe")
      >>> table = ascii.read("data/cds/glob/lmxbrefs.dat", readme="data/cds/glob/ReadMe")

    The table name and the CDS ReadMe file can be entered as URLs.  This can be used
    to directly load tables from the Internet.  For example, Vizier tables from the
    CDS::

      >>> table = ascii.read("ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/snrs.dat",
      ...             readme="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/253/ReadMe")

    If the header (ReadMe) and data are stored in a single file and there
    is content between the header and the data (for instance Notes), then the
    parsing process may fail.  In this case you can instruct the reader to
    guess the actual start of the data by supplying ``data_start='guess'`` in the
    call to the ``ascii.read()`` function.  You should verify that the output
    data table matches expectation based on the input CDS file.

    **Using a reader object**

    When ``Cds`` reader object is created with a ``readme`` parameter
    passed to it at initialization, then when the ``read`` method is
    executed with a table filename, the header information for the
    specified table is taken from the ``readme`` file.  An
    ``InconsistentTableError`` is raised if the ``readme`` file does not
    have header information for the given table.

      >>> readme = "data/vizier/ReadMe"
      >>> r = ascii.get_reader(ascii.Cds, readme=readme)
      >>> table = r.read("data/vizier/table1.dat")
      >>> # table5.dat has the same ReadMe file
      >>> table = r.read("data/vizier/table5.dat")

    If no ``readme`` parameter is specified, then the header
    information is assumed to be at the top of the given table.

      >>> r = ascii.get_reader(ascii.Cds)
      >>> table = r.read("data/cds.dat")
      >>> #The following gives InconsistentTableError, since no
      >>> #readme file was given and table1.dat does not have a header.
      >>> table = r.read("data/vizier/table1.dat")
      Traceback (most recent call last):
        ...
      InconsistentTableError: No CDS section delimiter found

    Caveats:

    * The Units and Explanations are available in the column ``unit`` and
      ``description`` attributes, respectively.
    * The other metadata defined by this format is not available in the output table.
    """
    _format_name = 'cds'
    _io_registry_format_aliases = ['cds']
    #_io_registry_can_write = False
    _description = 'CDS format table'

    data_class = CdsData
    header_class = CdsHeader

    def __init__(self, readme=None, cdsdict={}):
        super().__init__()
        self.header.readme = readme

        # The cds dict contains the default values for table meta info
        # required to fill up the ReadMe template.
        self.cdsdicts = cdsdicts
        self.cdsdicts.update(cdsdict)
        self.header.cdsdicts = self.cdsdicts
        #self.data.cdsdicts = self.cdsdicts

    def write(self, table=None):
        self.data.header = self.header
        self.header.position_line = None
        self.header.start_line = None
        self.data.start_line = None
        return super().write(table)

    def read(self, table):
        # If the read kwarg `data_start` is 'guess' then the table may have extraneous
        # lines between the end of the header and the beginning of data.
        if self.data.start_line == 'guess':
            # Replicate the first part of BaseReader.read up to the point where
            # the table lines are initially read in.
            with suppress(TypeError):
                # For strings only
                if os.linesep not in table + '':
                    self.data.table_name = os.path.basename(table)

            self.data.header = self.header
            self.header.data = self.data

            # Get a list of the lines (rows) in the table
            lines = self.inputter.get_lines(table)

            # Now try increasing data.start_line by one until the table reads successfully.
            # For efficiency use the in-memory list of lines instead of `table`, which
            # could be a file.
            for data_start in range(len(lines)):
                self.data.start_line = data_start
                with suppress(Exception):
                    table = super().read(lines)
                    return table
        else:
            return super().read(table)
