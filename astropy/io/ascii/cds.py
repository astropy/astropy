# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

cds.py:
  Classes to read CDS / Vizier table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import fnmatch
import itertools
import re

from . import core
from . import fixedwidth


__doctest_skip__ = ['*']


class CdsHeader(core.BaseHeader):
    col_type_map = {'e': core.FloatType,
                    'f': core.FloatType,
                    'i': core.IntType,
                    'a': core.StrType}

    def get_type_map_key(self, col):
        match = re.match(r'\d*(\S)', col.raw_type.lower())
        if not match:
            raise ValueError('Unrecognized CDS format "%s" for column "%s"' % (
                col.raw_type, col.name))
        return match.group(1)

    def __init__(self, readme=None):
        """Initialize ReadMe filename.

        :param readme: The ReadMe file to construct header from.
        :type readme: String

        CDS tables have their header information in a separate file
        named "ReadMe". The ``get_cols`` method will read the contents
        of the ReadMe file given by ``self.readme`` and set the various
        properties needed to read the data file. The data file name
        will be the ``table`` passed to the ``read`` method.
        """
        core.BaseHeader.__init__(self)
        self.readme = readme

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines`` for a CDS
        header.

        :param lines: list of table lines
        :returns: list of table Columns
        """
        # Read header block for the table ``self.data.table_name`` from the read
        # me file ``self.readme``.
        if self.readme and self.data.table_name:
            in_header = False
            f = open(self.readme, "r")
            # Header info is not in data lines but in a separate file.
            lines = []
            comment_lines = 0
            for line in f:
                line = line.strip()
                if in_header:
                    lines.append(line)
                    if line.startswith('------') or line.startswith('======='):
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
                raise core.InconsistentTableError("Can't find table {0} in {1}".format(
                    self.data.table_name, self.readme))
            f.close()

        found_line = False

        for i_col_def, line in enumerate(lines):
            if re.match(r'Byte-by-byte Description', line, re.IGNORECASE):
                found_line = True
            elif found_line: # First line after list of file descriptions
                i_col_def -= 1 # Set i_col_def to last description line
                break

        re_col_def = re.compile(r"""\s*
                                    (?P<start> \d+ \s* -)? \s*
                                    (?P<end>   \d+)        \s+
                                    (?P<format> [\w.]+)     \s+
                                    (?P<units> \S+)        \s+
                                    (?P<name>  \S+)        \s+
                                    (?P<descr> \S.+)""",
                                re.VERBOSE)

        cols = []
        for line in itertools.islice(lines, i_col_def+4, None):
            if line.startswith('------') or line.startswith('======='):
                break
            match = re_col_def.match(line)
            if match:
                col = core.Column(name=match.group('name'))
                col.start = int(re.sub(r'[-\s]', '',
                                       match.group('start') or match.group('end'))) - 1
                col.end = int(match.group('end'))
                col.unit = match.group('units')
                if col.unit == '---':
                    col.unit = None  # "---" is the marker for no unit in CDS table
                col.description = match.group('descr').strip()
                col.raw_type = match.group('format')
                col.type = self.get_col_type(col)

                match = re.match(
                    r'\? (?P<equal> =)? (?P<nullval> \S*)', col.description, re.VERBOSE)
                if match:
                    if issubclass(col.type, core.FloatType):
                        fillval = 'nan'
                    else:
                        fillval = '-999'
                    if match.group('nullval') == '':
                        col.null = ''
                    elif match.group('nullval') == '-':
                        col.null = '---'
                    else:
                        col.null = match.group('nullval')
                    self.data.fill_values.append((col.null, fillval, col.name))

                cols.append(col)
            else:  # could be a continuation of the previous col's description
                if cols:
                    cols[-1].description += line.strip()
                else:
                    raise ValueError('Line "%s" not parsable as CDS header' % line)

        self.names = [x.name for x in cols]

        self.cols = cols


class CdsData(core.BaseData):
    """
    CDS table data reader.

    Attributes
    ----------
    data_line : int or str
        If an int, gives the line number at which table data begins.
        It can also be the string 'guess' if the user intends
        for the reader to loop over data lines in order to find the
        location of table data. If data_line is None (as it is by
        default), then the reader will treat input as strict CDS and
        will simply begin inputting data after the appropriate delimiter.
    """
    splitter_class = fixedwidth.FixedWidthSplitter
    data_line = None

    def __init__(self, data_line=None):
        """Initialize CDS data reader.

        :param data_line: A parameter specifying how CdsData should
        look for the beginning of table data.
        :type data_line: int, str, or None
        """
        core.BaseData.__init__(self)
        self.data_line = data_line
            
    def get_str_vals(self):
        """Return valid rows following the last section delimiter as strings"""
        # If the header has a ReadMe and data has a filename
        # then just use data_lines, as the data lines do not have header
        # info. The ``read`` method adds the table_name to the ``data``
        # attribute.
        if self.header.readme and self.table_name:
            return self.splitter(self.data_lines)
        elif isinstance(self.data_line, int):
            # Split beginning at the specified line (start_line = 1, 2, ...)
            return self.splitter(self.data_lines[self.data_line - 1:])
        elif self.data_line is not None and self.data_line != 'guess':
            return core.InconsistentTableError(
                "Invalid value for start_line: '{0}'".format(self.data_line))
            
        i_sections = [i for (i, x) in enumerate(self.data_lines)
                      if x.startswith('------') or x.startswith('=======')]
        if not i_sections:
            raise core.InconsistentTableError('No CDS section delimiter found')

        bottom_lines = self.data_lines[i_sections[-1] + 1:]

        if self.data_line is None:
            return self.splitter(bottom_lines) # Return data after delimiter
        # Begin guessing the row at which table data starts
        type_map = {core.FloatType:float, core.IntType:int, core.StrType:str}
        self.splitter.cols = self.header.cols

        # Need fill values. Usually done later, so set here for splitter cols
        self._set_fill_values(self.splitter.cols)

        # Check lines below final delimiter for first valid row
        for i, vals in enumerate(self.splitter(bottom_lines)):
            if len(vals) != len(self.header.cols): # Incorrect number of columns
                continue
            try:
                for j, val in enumerate(vals):
                    if val not in self.splitter.cols[j].fill_values:
                        col_type = self.header.get_col_type(self.header.cols[j])
                        type_map[col_type](val)
            except ValueError: # One or more values were of incorrect type
                continue
            return self.splitter(bottom_lines[i:]) # First valid line found
            
        return [] # No data found


class Cds(core.BaseReader):
    """Read a CDS format table.  See http://vizier.u-strasbg.fr/doc/catstd.htx.
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
      >>> table = ascii.read("t/cds.dat")
      >>> table = ascii.read("t/vizier/table1.dat", readme="t/vizier/ReadMe")
      >>> table = ascii.read("t/cds/multi/lhs2065.dat", readme="t/cds/multi/ReadMe")
      >>> table = ascii.read("t/cds/glob/lmxbrefs.dat", readme="t/cds/glob/ReadMe")

    **Using a reader object**

    When ``Cds`` reader object is created with a ``readme`` parameter
    passed to it at initialization, then when the ``read`` method is
    executed with a table filename, the header information for the
    specified table is taken from the ``readme`` file.  An
    ``InconsistentTableError`` is raised if the ``readme`` file does not
    have header information for the given table.

      >>> readme = "t/vizier/ReadMe"
      >>> r = ascii.get_reader(ascii.Cds, readme=readme)
      >>> table = r.read("t/vizier/table1.dat")
      >>> # table5.dat has the same ReadMe file
      >>> table = r.read("t/vizier/table5.dat")

    If no ``readme`` parameter is specified, then the header
    information is assumed to be at the top of the given table.

      >>> r = ascii.get_reader(ascii.Cds)
      >>> table = r.read("t/cds.dat")
      >>> #The following gives InconsistentTableError, since no
      >>> #readme file was given and table1.dat does not have a header.
      >>> table = r.read("t/vizier/table1.dat")
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
    _io_registry_can_write = False
    _description = 'CDS format table'

    def __init__(self, readme=None, data_line=None):
        core.BaseReader.__init__(self)
        self.header = CdsHeader(readme)
        self.data = CdsData(data_line)
        self.data.header = self.header

    def write(self, table=None):
        """Not available for the Cds class (raises NotImplementedError)"""
        raise NotImplementedError
