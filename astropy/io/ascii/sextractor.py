"""Asciitable: an extensible ASCII table reader and writer.

sextractor.py:
  Classes to read SExtractor table format

Built on daophot.py:
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

import re
from . import core
from . import basic

class Daophot(core.BaseReader):
    """Read a SExtractor file. 
       SExtractor is a package for faint-galaxy photometry.
       Bertin & Arnouts 1996, A&A Supp. 317, 393.
       http://www.astromatic.net/software/sextractor

    Example::

      # 1 NUMBER
      # 2 ALPHA_J2000
      # 3 DELTA_J2000
      # 4 FLUX_RADIUS
      # 7 MAG_AUTO  -- note the skipped numbers since flux_radius has 3 columns

      #K MERGERAD   = INDEF                   scaleunit  %-23.7g  
      #K IRAF = NOAO/IRAFV2.10EXPORT version %-23s
      #K USER = davis name %-23s
      #K HOST = tucana computer %-23s
      #
      #N ID    XCENTER   YCENTER   MAG         MERR          MSKY           NITER    \\
      #U ##    pixels    pixels    magnitudes  magnitudes    counts         ##       \\
      #F %-9d  %-10.3f   %-10.3f   %-12.3f     %-14.3f       %-15.7g        %-6d     
      #
      #N         SHARPNESS   CHI         PIER  PERROR                                \\
      #U         ##          ##          ##    perrors                               \\
      #F         %-23.3f     %-12.3f     %-6d  %-13s
      #
      14       138.538   256.405   15.461      0.003         34.85955       4        \\
      -0.032      0.802       0     No_error

    The keywords defined in the #K records are available via the Daophot reader object::

      reader = asciitable.get_reader(Reader=asciitable.Daophot)
      data = reader.read('t/daophot.dat')
      for keyword in reader.keywords:
          print keyword.name, keyword.value, keyword.units, keyword.format
    
    """
    
    def __init__(self):
        core.BaseReader.__init__(self)
        self.header = SextractorHeader()
        self.inputter = core.ContinuationLinesInputter()
        self.data.splitter.delimiter = ' '
        self.data.start_line = 0
        self.data.comment = r'\s*#'
    
    def read(self, table):
        output = core.BaseReader.read(self, table)
        reader = core._get_reader(Reader=basic.NoHeader, comment=r'(?!#K)', 
                                  names = ['temp1','keyword','temp2','value','unit','format'])
        headerkeywords = reader.read(self.comment_lines)

        for line in headerkeywords:
            self.keywords.append(core.Keyword(line['keyword'], line['value'], 
                                              units=line['unit'], format=line['format']))
        self.table = output
        self.cols = self.header.cols

        return self.table

    def write(self, table=None):
        raise NotImplementedError


class SExtractorHeader(core.BaseHeader):
    """Read the header from a file produced by the IRAF DAOphot routine."""
    def __init__(self):
        core.BaseHeader.__init__(self)
        self.comment = r'^\s*#\s*\S\D.*' # Find lines that dont have "# digit"

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines`` for a SExtractor
        header.  The SExtractor header is specialized so that we just copy the entire BaseHeader
        get_cols routine and modify as needed.

        :param lines: list of table lines
        :returns: list of table Columns
        """

        # This assumes that the columns are listed in order, one per line with a 
        # header comment string of the format: "# 1 ID"
        # However, some may be missing and must be inferred from skipped column numbers
        columns = {}
        re_name_def = re.compile(r'^\s*#\s*([0-9]+).*') # E.g. '# 1 ID'
        for line in lines:
            if not line.startswith('#'):
                break                   # End of header lines
            else:
                match = re_name_def.search(line)
                if match:
                    words = match.group(0).strip().strip('#').split()
                    colnumber = int(words[0]) # First string is the column number
                    colname = words[1]   # second string is the column name
                    columns[colnumber] = colname
        # Handle skipped column numbers
        colnumbers = sorted(columns.iterkeys()):
        previous_column = 0
        for n in colnumbers:
            if k != previous_column + 1:
                for c in range(previous_column+1,k):
                    column_name = columns[previous_column]+"_%d" % (c-previous_column)
                    columns[c] = column_name
            previous_column = n

        # Add the columns in order to self.names
        colnumbers = sorted(columns.iterkeys()):
        self.names = []
        for n in colnumbers:
            self.names.append(columns[n])
        
        if not self.names:
            raise core.InconsistentTableError('No column names found in SExtractor header')
        
        self._set_cols_from_names()

