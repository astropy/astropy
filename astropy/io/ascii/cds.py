"""Asciitable: an extensible ASCII table reader and writer.

cds.py:
  Classes to read CDS / Vizier table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * Neither the name of the Smithsonian Astrophysical Observatory nor the
##       names of its contributors may be used to endorse or promote products
##       derived from this software without specific prior written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
## ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
## WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
## DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  
## SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import fnmatch
import itertools
import re

import asciitable.core as core
import asciitable.fixedwidth as fixedwidth

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
            f = open(self.readme,"r")
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
                raise core.InconsistentTableError("Cant' find table {0} in {1}".format(
                        self.data.table_name, self.readme))
            f.close()
                       
        for i_col_def, line in enumerate(lines):
            if re.match(r'Byte-by-byte Description', line, re.IGNORECASE):
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
        for i, line in enumerate(itertools.islice(lines, i_col_def+4, None)):
            if line.startswith('------') or line.startswith('======='):
                break
            match = re_col_def.match(line)
            if match:
                col = core.Column(name=match.group('name'), index=i)
                col.start = int(re.sub(r'[-\s]', '', match.group('start') or match.group('end'))) - 1
                col.end = int(match.group('end'))
                col.units = match.group('units')
                col.descr = match.group('descr')
                col.raw_type = match.group('format')
                col.type = self.get_col_type(col)

                match = re.match(r'\? (?P<equal> =)? (?P<nullval> \S*)', col.descr, re.VERBOSE)
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
                    cols[-1].descr += line.strip()
                else:
                    raise ValueError('Line "%s" not parsable as CDS header' % line)

        self.names = [x.name for x in cols]
        names = set(self.names)
        if self.include_names is not None:
            names.intersection_update(self.include_names)
        if self.exclude_names is not None:
            names.difference_update(self.exclude_names)
            
        self.cols = [x for x in cols if x.name in names]
        self.n_data_cols = len(self.cols)

        # Re-index the cols because the FixedWidthSplitter does NOT return the ignored
        # cols (as is the case for typical delimiter-based splitters)
        for i, col in enumerate(self.cols):
            col.index = i
            

class CdsData(core.BaseData):
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
        i_sections = [i for (i, x) in enumerate(lines)
                      if x.startswith('------') or x.startswith('=======')]
        if not i_sections:
            raise core.InconsistentTableError('No CDS section delimiter found')
        return lines[i_sections[-1]+1 : ]


class Cds(core.BaseReader):
    """Read a CDS format table: http://vizier.u-strasbg.fr/doc/catstd.htx.
    Example::

      Table: Spitzer-identified YSOs: Addendum
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
        1 03 28 39.09

    **Basic usage**

    Use the ``asciitable.read()`` function as normal, with an optional ``readme``
    parameter indicating the CDS ReadMe file.  If not supplied it is assumed that
    the header information is at the top of the given table.  Examples::

      >>> import asciitable
      >>> table = asciitable.read("t/cds.dat")
      >>> table = asciitable.read("t/vizier/table1.dat", readme="t/vizier/ReadMe")
      >>> table = asciitable.read("t/cds/multi/lhs2065.dat", readme="t/cds/multi/ReadMe")
      >>> table = asciitable.read("t/cds/glob/lmxbrefs.dat", readme="t/cds/glob/ReadMe")

    **Using a reader object**

    When ``Cds`` reader object is created with a ``readme`` parameter
    passed to it at initialization, then when the ``read`` method is
    executed with a table filename, the header information for the
    specified table is taken from the ``readme`` file.  An
    ``InconsistentTableError`` is raised if the ``readme`` file does not
    have header information for the given table.
    
      >>> readme = "t/vizier/ReadMe"
      >>> r = asciitable.get_reader(asciitable.Cds, readme=readme)
      >>> table = r.read("t/vizier/table1.dat")
      >>> # table5.dat has the same ReadMe file
      >>> table = r.read("t/vizier/table5.dat")

    If no ``readme`` parameter is specified, then the header
    information is assumed to be at the top of the given table.

      >>> r = asciitable.get_reader(asciitable.Cds)
      >>> table = r.read("t/cds.dat")
      >>> #The following gives InconsistentTableError, since no
      >>> #readme file was given and table1.dat does not have a header.
      >>> table = r.read("t/vizier/table1.dat")
      Traceback (most recent call last):
        ...
      InconsistentTableError: No CDS section delimiter found

    Caveats:

    * Format, Units, and Explanations are available in the ``Reader.cols`` attribute.
    * All of the other metadata defined by this format is ignored.

    Code contribution to enhance the parsing to include metadata in a Reader.meta
    attribute would be welcome.

    """
    def __init__(self, readme=None):
        core.BaseReader.__init__(self)
        self.header = CdsHeader(readme)
        self.data = CdsData()

    def write(self, table=None):
        """Not available for the Cds class (raises NotImplementedError)"""
        raise NotImplementedError

CdsReader = Cds

