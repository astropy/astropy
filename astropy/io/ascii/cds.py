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
from contextlib import suppress

from . import core
from . import fixedwidth

from astropy.units import Unit


from astropy.io import ascii
from astropy.table import Table
from .CDSColumn import CDSColumn
import numpy as np
import sys
from string import Template
from textwrap import wrap, fill
import datetime
import logging
import math

MAX_SIZE_README_LINE = 80
MAX_COL_INTLIMIT = 10000000


__doctest_skip__ = ['*']

cdsdicts = {'title': 'Title ?',
            'author': '1st author ?',
            'catalogue': '',
            'date': 'Date ?',
            'abstract': 'Abstract ?',
            'authors': 'Authors ?',
            'bibcode': 'ref ?',
            'keywords': ''
            }


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


class CdsData(fixedwidth.FixedWidthData):   #core.BaseData):
    """CDS table data reader
    """
    splitter_class = fixedwidth.FixedWidthSplitter

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
        #print(lines)
        #print(self.cols)
        self.splitter.delimiter = ' '
        fixedwidth.FixedWidthData.write(self, lines)


class CdsInputter(core.BaseInputter):
    
    def process_lines(self, lines):
        print(lines)
        return lines


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
    inputter_class = CdsInputter

    def __init__(self, readme=None):
        super().__init__()
        self.header.readme = readme
        self.cdsdicts = cdsdicts
        self.header.position_line = None

    def write(self, table=None):
        #tablemaker = CDSTablesMaker()
        #tablemaker.addTable(table, name='astropyTable')
        #CdsTable = tablemaker.returnTable()
        #return core.BaseReader.write(self, table=CdsTable)

        name = 'table'
        description = 'this is a table'
        cdsTable = CDSAstropyTable(table, name, description)
        #cdsTable.makeCDSTable()

        self.data.header = self.header
        #self.header.data = cdsTable.returnLines()
        lines = cdsTable.returnLines()  #--this is done correctly!
        """
        self.header.start_line = None
        self.data.start_line = None
        return core.BaseData.write(self, lines)
        """
        return core.BaseReader.write(self, table=table)

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


class CDSException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class CDSTable:
    """Manage table
    """

    def __init__(self, table):
        """ Constructor
        :param table: astropy table
        """
        self.logger = logging.getLogger('CDSTable')
        self.__bytebybytetemplate = None
        self.__linewidth = None
        self.__cds_columns = [] # array of CDSColumn
        self.notes = None
        self.table = table
        self.nlines = None
        self.nullvalue = None
        self.init_meta_columns()
        #self.columns = self.table.columns

    def get_column(self, name=None):
        """Get CDS meta columns
        :param name: column name (default is None)
        :return: CDSColumn
        """
        if name is None:
            return self.__cds_columns
        for column in self.__cds_columns:
            if column.name == name:
                return column
        return None

    def setByteByByteTemplate(self, name):
        """Set a User Byte-By-Byte template
        """
        self.__bytebybytetemplate = name

    def getByteByByteTemplate(self):
        """Get the Byte-By-Byte template
        """
        return self.__bytebybytetemplate

    def __writeTable(self, fo):
        col_length = len(self.__cds_columns)

        for nrec in range(len(self.table)):
            fo.write(self.__cds_columns[0].value(nrec))
            for i in range(1, col_length):
                fo.write(" " + self.__cds_columns[i].value(nrec))
            fo.write("\n")

    def returnLines(self):
        for col in self.__cds_columns:
            col.parse()
            if self.nullvalue:
                col.set_null_value(self.nullvalue)

        col_length = len(self.__cds_columns)

        lines = []
        for nrec in range(len(self.table)):
            line = []
            line.append(self.__cds_columns[0].value(nrec))
            for i in range(1, col_length):
                line.append(" " + self.__cds_columns[i].value(nrec))
            lines.append(''.join(line))
        return lines

    def makeCDSTable(self, fd=None):
        """Make the standardized table in ASCII aligned format.
        :param fd: file descriptor (by default, the methods creates a new file with CDSTable.name)
        """
        for col in self.__cds_columns:
            col.parse()
            if self.nullvalue:
                col.set_null_value(self.nullvalue)

        if fd is None:
            fd = open(self.name, "w") #, buffering=2048)
        self.__writeTable(fd)
        fd.close()

    def init_meta_columns(self):
        """Initialize list of CDSColumns  (self.__cds_columns)
        """
        if self.table is None:
            raise CDSException("table is empty")

        if isinstance(self.table, str): return;
        if self.__cds_columns: return self.__cds_columns

        self.__cds_columns = []
        for column in self.table.columns:
            col = CDSColumn(self.table[column])
            self.__cds_columns.append(col)

    def getlinewidth(self):
        """Get ASCII table line width
        :return: the line size
        """
        if self.__linewidth != None:
            return self.__linewidth

        self.__linewidth = 0
        for column in self.__cds_columns:
            self.__linewidth += column.size + 1
        self.__linewidth -= 1
        return self.__linewidth

    def gatherSexagesimalColumn(self):# not used
        """gather/detects sexagesimal columns if in different columns
        :return: True if found
        """
        if isinstance(self.table, Table): return

        array = ""
        for col in self.table.columns:
            if self.__column.dtype.name.startswith("i"):
                array +='i'
            elif self.__column.dtype.name.startswith("f"):
                array +='f'
            elif self.__column.dtype.name.startswith("s"):
                array += 's'
            else:
                array += ' '
        n = array.indexof("iifiif")
        if n < 9:
            return False
        return True



class CDSAstropyTable(CDSTable):
    """Manage astropy table
    """

    def __init__(self, table, name=None, description=None):
        """Constructor
        :param table: astropy table
        :param name: table name in output
        :param description: table description
        """
        if not isinstance(table, Table):
            raise CDSException("input is not Astropy Table")

        CDSTable.__init__(self, table)

        self.__filename = "astropytable"
        if name is None:
            self.name = self.__filename
        else:
            self.name = name
        self.description = description
        self.nlines = len(table)


class CDSTablesMaker:
    """Generate standardized tables and ReadMe
    """

    def __init__(self, out=None, debug=False):
        """Constructor
        :param out: the output file (default: stdout)
        :param debug: True/False
        """
        self.__dir = os.path.dirname(os.path.realpath(__file__))
        self.__tables = []
        self.__ReadMetemplate = self.__dir + "/src/ReadMe.template"

        if out != None:
            sys.stdout = open(out, 'w')

        self.title = cdsdicts['title']
        self.author = cdsdicts['author']
        self.catalogue = cdsdicts['catalogue']
        self.date = cdsdicts['date']
        self.abstract = cdsdicts['abstract']
        self.authors = cdsdicts['authors']
        self.bibcode = cdsdicts['bibcode']
        self.keywords = cdsdicts['keywords']
        self.ref = None
        self.__templatevalue = None
        self.logger = logging.getLogger('CDSTablesMaker')

        if debug is True:
            self.logger.basicConfig(level=logging.DEBUG)

    def addTable(self, table, name=None, description=None, nullvalue=None):
        """Add a Table, memorize the meta-data and generate the standardized CDS table
        :param table: table (type accepted: astropy, numpy, filename, CDSTable)
        :param name: the name used in output
        :param description: table description
        :param nullvalue: set a nullvalue (aplyed for all columns)
        :return: CDStable created
        """
        self.logger.debug("add table")
        if isinstance(table, CDSTable):
            cdstable = table
            if name != None:
                cdstable.name = name
            if description != None:
                cdstable.description = description
        elif isinstance(table, Table):
            cdstable = CDSAstropyTable(table, name, description)
        else:
            raise CDSException("type " + type(table) + " is not accepted (only String or astropy.Table)")

        self.__tables.append(cdstable)
        self.logger.debug("append table ok")

        cdstable.nullvalue = nullvalue
        return cdstable

    def returnTable(self):
        """returns the first table"""
        for table in self.__tables:
            return table.table

    def writeCDSTables(self):
        """Write tables in ASCII format
        """
        for table in self.__tables:
            table.makeCDSTable()
            self.logger.debug("make CDS table " + table.name)

    def getTablesInfo(self):
        """get tables information.
        :return: info (dictionary)
        """
        info = []
        for table in self.__tables:
            info.append({'name': table.name, 'lines': len(table.table), 'width': table.getlinewidth()})
        return info

    def getTables(self):
        """Get the CDSTable list
        :return: list of CDSTable
        """
        return self.__tables

    def printTablesIndex(self, outBuffer=False):
        """Print the tables index
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        sz = [14, 0, 8]
        for tab in self.__tables:
            if len(tab.name) > sz[0]:
                sz[0] = len(tab.name)
            l = len(str(tab.getlinewidth()))
            if l > sz[1]:
                sz[1] = l

        fmtt = "{0:" + str(sz[0]) + "s} {1:" + str(sz[1]) + "d} {2:>" + str(sz[2]) + "s}  {3:s}"
        buff = fmtt.format("ReadMe", MAX_SIZE_README_LINE, ".", "this file") + "\n"
        for tab in self.__tables:
            buff += fmtt.format(self.__strFmt(tab.name),
                                tab.getlinewidth(),
                                str(tab.nlines),
                                self.__strFmt(tab.description))
            buff += "\n"

        if outBuffer: return buff
        sys.stdout.write(buff)

    def __strFmt(self, string):
        if string is None:
            return ""
        else:
            return string

    def __strFmtRa(self, column, fmtb, startb):
        ra = column.getSexaRA()
        n = startb
        nvalue = ""
        if column.hasNull: nvalue = "? "

        buff = fmtb.format(n, n + ra.RAh.size - 1, "",
                           ra.RAh.fortran_format, ra.RAh.unit, ra.RAh.name, nvalue + ra.RAh.description)
        n += ra.RAh.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + ra.RAm.size - 1, "",
                            ra.RAm.fortran_format, ra.RAm.unit, ra.RAm.name, nvalue + ra.RAm.description)
        n += ra.RAm.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + ra.RAs.size - 1, "",
                            ra.RAs.fortran_format, ra.RAs.unit, ra.RAs.name, nvalue + ra.RAs.description)
        return buff

    def __strFmtDe(self, column, fmtb, startb):
        de = column.getSexaDE()
        n = startb
        nvalue = ""
        if column.hasNull: nvalue = "? "
        buff = fmtb.format(n, n + de.DEsign.size - 1, "",
                           de.DEsign.fortran_format, de.DEsign.unit, de.DEsign.name, nvalue + de.DEsign.description)
        n += de.DEsign.size
        buff += '\n'
        buff += fmtb.format(n, n + de.DEd.size - 1, "",
                            de.DEd.fortran_format, de.DEd.unit, de.DEd.name, nvalue + de.DEd.description)
        n += de.DEd.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + de.DEm.size - 1, "",
                            de.DEm.fortran_format, de.DEm.unit, de.DEm.name, nvalue + de.DEm.description)
        n += de.DEm.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + de.DEs.size - 1, "",
                            de.DEs.fortran_format, de.DEs.unit, de.DEs.name, nvalue + de.DEs.description)
        return buff

    def __splitLine(self, line, shift=0):
        """Split line 80 char
           :param shift: add left blank
        """
        if shift > MAX_SIZE_README_LINE:
            shift = 0
        return ("\n" + " " * shift).join(wrap(line, width=MAX_SIZE_README_LINE - shift))

    def __add_authors(self, line, shift=0):
        """Split the line containing the authors without separate given and surname
        :param line: authors list in a line
        :param shift: add left blank
        :return: authors formatted string
        """
        # Find all spaces followed by authors's initials
        space_letter_list = re.findall(" (?:[A-Z]\.-?)+", line)

        # Replace founded spaces by !
        for space_letter in space_letter_list:
            line = line.replace(space_letter, "!" + space_letter.strip())

        # Wrap the text by using spaces as breakpoint and then replace ! by spaces so given name and surname
        # are not separate
        new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift * " ").replace("!", " ")

        if self.author:
            if len(new_line)>0:
                return self.author+", "+new_line
            else:
                return self.author
        return new_line

    def __add_keywords(self, line, shift=0):
        """Split the line containing the authors without separate given and surname
        :param line: keywords list in a line
        :param shift: add left blank
        :return: keywods formatted string
        """
        # Replace all spaces that are NOT precede by ; with !
        line = re.sub("(?<!;) ", "!", line)

        # Wrap the text by using spaces as breakpoint and then replace ! by spaces so there is no break line
        # between two ;
        new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift * " ").replace("!", " ")

        return new_line

    def printByteByByte(self, table, outBuffer=False):
        """Print byte-by-byte
        :param table: the CDSTable
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        if isinstance(table, CDSMRTTable):
            buff = table.getByteByByte()
            if table.notes != None:
                buff += "-" * 80 + "\n"
                for line in table.notes:
                    try:
                        buff += line.strip().encode('ascii', 'replace').decode()
                    except Exception as err:
                        self.logger.error("error detected in Notes " + str(err))
                        buff += "?" * len(line) + "\n"
                    buff += "\n"
                buff += "-" * 80 + "\n"

            if outBuffer: return buff
            sys.stdout.write(buff)
            return

        columns = table.get_column()
        startb = 1
        sz = [0, 0, 1, 7]
        l = len(str(table.getlinewidth()))
        if l > sz[0]:
            sz[0] = l
            sz[1] = l
        self.logger.debug("size sz=" + str(sz))
        fmtb = "{0:" + str(sz[0]) + "d}-{1:" + str(sz[1]) + "d} {2:" + str(sz[2]) + "s}"
        for column in columns:
            if len(column.name) > sz[3]:
                sz[3] = len(column.name)
        fmtb += " {3:6s} {4:6s} {5:" + str(sz[3]) + "s} {6:s}"
        buff = ""
        nsplit = sz[0] + sz[1] + sz[2] + sz[3] + 16

        for column in columns:
            endb = column.size + startb - 1
            if column.formatter.fortran_format[0] == 'R':
                buff += self.__strFmtRa(column, fmtb, startb) + "\n"
            elif column.formatter.fortran_format[0] == 'D':
                buff += self.__strFmtDe(column, fmtb, startb) + "\n"
            else:
                description = column.description
                if column.hasNull:
                    nullflag = "?"
                else:
                    nullflag = ""

                borne = ""
                if column.min and column.max:
                    if column.formatter.fortran_format[0] == 'I':
                        if abs(column.min) < MAX_COL_INTLIMIT and abs(column.max) < MAX_COL_INTLIMIT:
                            if column.min == column.max:
                                borne = "[{0}]".format(column.min)
                            else:
                                borne = "[{0}/{1}]".format(column.min, column.max)
                    elif column.formatter.fortran_format[0] in ('E','F'):
                        borne = "[{0}/{1}]".format(math.floor(column.min*100)/100.,
                                                   math.ceil(column.max*100)/100.)

                description = "{0}{1} {2}".format(borne, nullflag, description)
                newline = fmtb.format(startb, endb, "",
                                      self.__strFmt(column.formatter.fortran_format),
                                      self.__strFmt(column.unit),
                                      self.__strFmt(column.name),
                                      description)

                if len(newline) > MAX_SIZE_README_LINE:
                    buff += ("\n").join(wrap(newline,
                                             subsequent_indent=" " * nsplit,
                                             width=MAX_SIZE_README_LINE))
                    buff += "\n"
                else:
                    buff += newline + "\n"
            startb = endb + 2

        if table.notes != None:
            buff += "-" * 80 + "\n"
            for line in table.notes: buff += line + "\n"
            buff += "-" * 80 + "\n"

        if outBuffer: return buff
        sys.stdout.write(buff)

    def __getByteByByteTemplate(self, table):
        templateValue = {'file': table.name,
                         'bytebybyte': self.printByteByByte(table, outBuffer=True)}

        templatename = table.getByteByByteTemplate()
        if templatename is None: templatename = self.__dir + "/src/bytebybyte.template"
        with open(templatename) as filein:
            src = Template(filein.read())
            return src.substitute(templateValue)

    def setReadmeTemplate(self, templatename, templateValue=None):
        """Set a user ReadMe template
        :param templateName: the template name
        :param templateValue: dictionary to fill added variable in templatename
        """
        self.__ReadMetemplate = templatename
        self.__templatevalue = templateValue

    def putRef(self, catname, title=""):
        """Put a reference.
        :param catname: catalgogue name (string)
        :param title: the title (string)
        """
        if self.ref is None: self.ref = []
        if catname is None: raise Exception("catname is required")
        self.ref.append((catname, title))

    def printRef(self, outBuffer):
        """The "See also" section in ReadMe
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        if self.ref is None or len(self.ref) == 0: return

        buf = ""
        # Find the highest string length in first references column
        max_len = len(max([i[0] for i in self.ref], key=len))

        # Write references and align the double dot symbols
        for ref in self.ref:
            buf += self.__splitLine(" {0:<{max}} : {1}".format(ref[0], ref[1], max=max_len)) + "\n"

        if outBuffer is True: return buf
        sys.stdout.write(buf)

    def makeReadMe(self, out=sys.stdout):
        """Print the ReadMe
        :param out: file descriptor (default sys.stdout)
        """
        templateValue = {'catalogue': self.catalogue,
                         'title': self.__splitLine(self.title),
                         'author': self.author,
                         'date': self.date,
                         'abstract': self.__splitLine(self.abstract, shift=2),
                         'authors': self.__add_authors(self.authors, shift=4),
                         'bibcode': "=" + self.bibcode,
                         'keywords': self.__add_keywords(self.keywords, shift=len("keywords: ")),
                         'tablesIndex': self.printTablesIndex(outBuffer=True),
                         'seealso': self.printRef(outBuffer=True),
                         'bytebybyte': '',
                         'today': datetime.datetime.now().strftime("%d-%b-%Y")}

        if self.__templatevalue is not None:
            for key in self.__templatevalue: templateValue[key] = self.__templatevalue[key]

        buff = ""
        for table in self.__tables:
            buff += self.__getByteByByteTemplate(table)
            buff += "\n"
        templateValue['bytebybyte'] = buff

        with open(self.__ReadMetemplate) as filein:
            src = Template(filein.read())
            result = src.substitute(templateValue)
            out.write(result)

    def toMRT(self):
        """Transform tables into MRT (ASCII aligned table with byte-by-byte header)
        """
        templateValue = {'title': self.__splitLine(self.title),
                         'authors': self.__add_authors(self.authors, shift=4)}
        mrt_template = self.__dir + "/MRT.template"

        for table in self.__tables:
            for col in table.get_column():
                col.parse()
            templateValue['bytebybyte'] = self.__getByteByByteTemplate(table)

            with open(table.name, "w") as fd:
                with open(mrt_template) as filein:
                    src = Template(filein.read())
                    result = src.substitute(templateValue)
                    fd.write(result)
                table.makeCDSTable(fd)
