# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

cds.py:
  Classes to read CDS / Vizier table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""


import fnmatch
import itertools
import re
import os
from contextlib import suppress

from . import core
from . import fixedwidth

from astropy.units import Unit


__doctest_skip__ = ['*']


class CdsHeader(core.BaseHeader):
    _subfmt = 'CDS'

    col_type_map = {'e': core.FloatType,
                    'f': core.FloatType,
                    'i': core.IntType,
                    'a': core.StrType}

    'The ReadMe file to construct header from.'
    readme = None

    def get_type_map_key(self, col):
        match = re.match(r'\d*(\S)', col.raw_type.lower())
        if not match:
            raise ValueError('Unrecognized {} format "{}" for column "{}"'.format(
                self._subfmt, col.raw_type, col.name))
        return match.group(1)

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines`` for a CDS/MRT
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
        else:
            raise ValueError('no line with "Byte-by-byte Description" found')

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
                    col.unit = None  # "---" is the marker for no unit in CDS/MRT table
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
                        # CDS/MRT tables can use -, --, ---, or ---- to mark missing values
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


class CdsData(core.BaseData):
    """CDS table data reader
    """
    _subfmt = 'CDS'
    splitter_class = fixedwidth.FixedWidthSplitter

    def process_lines(self, lines):
        """Skip over CDS/MRT header by finding the last section delimiter"""
        # If the header has a ReadMe and data has a filename
        # then no need to skip, as the data lines do not have header
        # info. The ``read`` method adds the table_name to the ``data``
        # attribute.
        if self.header.readme and self.table_name:
            return lines
        i_sections = [i for i, x in enumerate(lines)
                      if x.startswith(('------', '======='))]
        if not i_sections:
            raise core.InconsistentTableError(f'No {self._subfmt} section delimiter found')
        return lines[i_sections[-1]+1:]  # noqa


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
    _io_registry_can_write = False
    _description = 'CDS format table'

    data_class = CdsData
    header_class = CdsHeader

    def __init__(self, readme=None):
        super().__init__()
        self.header.readme = readme

    def write(self, table=None):
        """Not available for the CDS class (raises NotImplementedError)"""
        raise NotImplementedError

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
