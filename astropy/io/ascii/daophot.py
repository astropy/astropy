# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.
daophot.py:
  Classes to read DAOphot table format
:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np
import itertools as itt
from collections import defaultdict

from . import core
from . import basic
from . import fixedwidth
from ...utils import OrderedDict
from ...extern.six.moves import zip, map
from .misc import first_true_index, first_false_index, groupmore


class DaophotHeader(core.BaseHeader):
    """Read the header from a file produced by the IRAF DAOphot routine."""
    def __init__(self):
        core.BaseHeader.__init__(self)
        self.comment = r'\s*#K'

        #self.re_digit_extract = re.compile( '\d+' )
        self.re_format = re.compile( '%-?(\d+)\.?\d?[sdfg]' )                        #regex for extracting the format strings
        self.re_header_keyword = re.compile(r'[#]K'
                                            r'\s+ (?P<name> \w+)'
                                            r'\s* = (?P<stuff> .+) $',
                                            re.VERBOSE)

        self.aperture_values = ''

    def __str__():
        if hasattr(self, 'lines'):
            return os.linesep.join( self.lines )
        else:
            return super(core.BaseHeader, self).__str__()

    def parse_col_defs(self, grouped_lines_dict):
        '''
        Parse a series of column defintion lines like below.  There may be several
        such blocks in a single file (where continuation characters have already been
        stripped).
        #N ID    XCENTER   YCENTER   MAG         MERR          MSKY           NITER
        #U ##    pixels    pixels    magnitudes  magnitudes    counts         ##
        #F %-9d  %-10.3f   %-10.3f   %-12.3f     %-14.3f       %-15.7g        %-6d
        '''
        line_ids = ('#N', '#U', '#F')
        coldef_dict = defaultdict(list)
        stripper = lambda s: s[2:].strip(' \\')                                          #function to strip identifier lines
        for defblock in zip( *map(grouped_lines_dict.get, line_ids) ):
            for key, line in zip(line_ids, map(stripper, defblock)):
                coldef_dict[key].append( line.split() )

        #save the original columns so we can use it later to reconstruct the original header for writing
        #self._original_columns = sum(self.coldef_dict['#N'], [])

        if self.data.is_multiline:
            #database contains multi-aperture data.
            #Autogen column names, units, formats from last row of column headers
            last_names, last_units, last_formats = list(zip(*map(coldef_dict.get, line_ids)))[-1]
            N_multiline = len(self.data.first_block)
            for i in np.arange(1, N_multiline+1).astype(str):
                #extra column names eg. RAPERT2, SUM2 etc...
                extended_names = list(map( ''.join, zip(last_names, itt.repeat(i)) ))
                if i=='1':      #Enumerate the names starting at 1
                    coldef_dict['#N'][-1] = extended_names
                else:
                    coldef_dict['#N'].append( extended_names )
                    coldef_dict['#U'].append( last_units )
                    coldef_dict['#F'].append( last_formats )

        #Get column widths from column format specifiers
        get_col_width = lambda s: int( self.re_format.search(s).groups()[0] )
        col_widths = [ [get_col_width(f) for f in formats]
                        for formats in coldef_dict['#F'] ]
        # original data format might be shorter than 80 characters and filled with spaces
        row_widths = np.fromiter( map(sum, col_widths), int )
        row_short = Daophot.table_width - row_widths
        #fix last column widths
        for w,r in zip(col_widths, row_short):
            w[-1] += r

        self.col_widths = col_widths
        
        #merge the multi-line header data into single line data
        coldef_dict = dict( (k, sum(v,[])) for (k,v) in coldef_dict.items() )

        return coldef_dict

    def update_meta(self, lines, meta):
        """
        Extract table-level keywords for DAOphot table.  These are indicated by
        a leading '#K ' prefix.
        """
        table_meta = meta['table']

        #self.lines = self.get_header_lines(lines)
        Nlines = len(self.lines)
        if Nlines > 0:
            #group the header lines according to their line identifiers (#K, #N, #U, #F or just # (spacer line))
            get_line_id = lambda s: s.split(None,1)[0]   #function that grabs the line identifier
            #group lines by the line identifier ('#N', '#U', '#F', '#K') and capture line index
            gid, groups = zip( *groupmore(get_line_id, self.lines, range(Nlines)) )
            #groups of lines and their indeces
            grouped_lines, gix = zip( *groups )
            #dict of line groups keyed by line identifiers
            grouped_lines_dict = gld = dict( zip(gid, grouped_lines) )

            #Update the table_meta keywords if necessary
            if '#K' in grouped_lines_dict:
                keywords = OrderedDict( map(self.extract_keyword_line, grouped_lines_dict['#K']) )
                table_meta['keywords'] = keywords

            coldef_dict = self.parse_col_defs( grouped_lines_dict )

            # If there wasn't a #U defined (not sure of DAOphot specification), then
            # replace the empty line with the right number of ## indicators, which matches
            # the DAOphot "no unit" tag.
            #if not len(coldef_dict['#U'])
            #for i, coldef_line in enumerate(coldef_lines):
                #if not coldef_line:
                    #coldef_lines[i] = ['##'] * len(self.names)

            line_ids = ('#N', '#U', '#F')
            for name, unit, fmt in zip(*map( coldef_dict.get, line_ids )):
                meta['cols'][name] = {  'unit'  : unit,
                                        'format': fmt  }

            self.meta = meta
            self.names = coldef_dict['#N']

    def extract_keyword_line(self, line):
        '''extract info from a header keyword line (#K) '''
        m = self.re_header_keyword.match(line)
        if m:
            vals = m.group('stuff').strip().rsplit(None, 2)
            keyword_dict = { 'units': vals[-2],
                            'format': vals[-1],
                            'value' : (vals[0] if len(vals) > 2 else "") }
            return m.group('name'), keyword_dict

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines`` for a DAOphot
        header.  The DAOphot header is specialized so that we just copy the entire BaseHeader
        get_cols routine and modify as needed.

        :param lines: list of table lines
        :returns: list of table Columns
        """

        if not self.names:
            raise core.InconsistentTableError('No column names found in DAOphot header')

        # Create the list of io.ascii column objects
        self._set_cols_from_names()

        # Set unit and format as needed.
        coldefs = self.meta['cols']
        for col in self.cols:
            unit, fmt = map(coldefs[col.name].get, ('unit', 'format'))
            if unit != '##':
                col.unit = unit
            if fmt != '##':
                col.format = fmt

        # Set column start and end positions.
        col_width = sum( self.col_widths, [] )
        ends = np.cumsum( col_width )
        starts = ends - col_width
        for i, col in enumerate(self.cols):
            col.start, col.end = starts[i], ends[i]
            col.span = col.end-col.start
            if hasattr(col, 'format'):
                if any(x in col.format for x in 'fg'):
                    col.type = core.FloatType
                elif 'd' in col.format:
                    col.type = core.IntType
                elif 's' in col.format:
                    col.type = core.StrType

        # INDEF is the missing value marker
        self.data.fill_values.append(('INDEF', '0'))


class DaophotData(core.BaseData):
    splitter_class = fixedwidth.FixedWidthSplitter
    start_line = 0
    comment = r'\s*#'

    def __init__(self):
        core.BaseData.__init__(self)
        self.is_multiline = False

    def get_data_lines(self, lines):

        # Special case for multiline daophot databases. Extract the aperture values from the first multiline data block
        if self.is_multiline:
            #grab the first column of the special block (aperture values) and recreate the aperture description string
            aplist = next( zip(*map(str.split, self.first_block)) )
            self.header.aperture_values = ', '.join( aplist )

        # Set self.data.data_lines to a slice of lines contain the data rows
        core.BaseData.get_data_lines(self, lines)


class DaophotInputter(core.ContinuationLinesInputter):

    continuation_char = '\\'
    multiline_char = '*'
    replace_char = ' '
    re_multiline     = re.compile( r'(#?)[^\\*#]*(\*?)(\\*) ?$' )


    def search_multiline(self, lines, depth=150):
        '''search lines for special continuation character to determine number of continued rows in a
        datablock.  For efficiency, depth gives the upper limit of lines to search.'''

        #The list of apertures given in the #K APERTURES keyword may not be complete!!
        #This happens if the string description of the aperture list is longer than the
        #field width of the #K APERTURES field.
        #In this case we have to figure out how many apertures there are based on the file
        #structure.

        comment, special, cont = zip( *(self.re_multiline.search(l).groups() for l in lines[:depth]) )

        #find first non-comment line
        data_start = first_false_index( comment )
        if data_start is None: #no data in lines[:depth].  This may be because there is no data in the file, or because the header is really huge.  If the latter, increasing the search depth should help
            return None, None, lines[:depth]

        header_lines = lines[:data_start]

        #find first line ending on special row continuation character '*'
        first_special = first_true_index( special[data_start:depth] )   #indexed relative to data_start
        if first_special is None:  #no special lines
            return None, None, header_lines

        #last line ending on special '*', but not on line continue '/'
        last_special = first_false_index( special[data_start+first_special:depth] )  #index relative to first_special
        #if first_special is None: #no end of special lines within search depth!  increase search depth
            #return self.search_multiline( lines, depth=2*depth )

        markers = np.cumsum( [data_start,first_special,last_special] )  #indexing now relative to line[0]
        multiline_block = lines[ markers[1]:markers[-1] ]               #multiline portion of first data block

        return markers, multiline_block, header_lines

    def process_lines(self, lines):

        markers, block, header = self.search_multiline( lines )
        self.data.is_multiline = not markers is None
        self.data.markers = markers
        self.data.first_block = block
        self.data.header.lines = header         #set the header lines returned by the search as a attribute of the header

        if not markers is None:
            lines = lines[ markers[0]: ]

        continuation_char       = self.continuation_char
        multiline_char          = self.multiline_char
        replace_char            = self.replace_char

        parts = []
        outlines = []
        for i, line in enumerate(lines):
            mo = self.re_multiline.search( line )
            if mo:
                comment, special, cont = mo.groups()
                if comment or cont:
                    line = line.replace(continuation_char, replace_char)
                if special:
                    line = line.replace(multiline_char, replace_char)
                if cont and not comment:
                    parts.append( line )
                if not cont:
                    parts.append(line)
                    outlines.append( ''.join(parts) )
                    parts = []
            else:
                raise ValueError( 'multiline re could not match line %i: %s' %(i,line) )

        return outlines


class Daophot(core.BaseReader):
    """Read a DAOphot file.
    Example::
      
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
      14       138.538     INDEF   15.461      0.003         34.85955       4        \\
                  -0.032      0.802       0     No_error
    
    The keywords defined in the #K records are available via the output table
    ``meta`` attribute::
      
      >>> import os
      >>> from astropy.io import ascii
      >>> filename = os.path.join(ascii.__path__[0], 'tests/t/daophot.dat')
      >>> data = ascii.read(filename)
      >>> for name, keyword in data.meta['keywords'].items():
      ...     print(name, keyword['value'], keyword['units'], keyword['format'])
      ...
      MERGERAD INDEF scaleunit %-23.7g
      IRAF NOAO/IRAFV2.10EXPORT version %-23s
      USER  name %-23s
      ...
    
    The unit and formats are available in the output table columns::
      
      >>> for colname in data.colnames:
      ...     col = data[colname]
      ...     print(colname, col.unit, col.format)
      ...
      ID None %-9d
      XCENTER pixels %-10.3f
      YCENTER pixels %-10.3f
      ...
      
    Any column values of INDEF are interpreted as a missing value and will be
    masked out in the resultant table.
    
    In case of multi-aperture daophot files containing repeated entries for the last
    row of fields, extra unique column names will be created by suffixing
    corresponding field names with numbers starting from 2 to N (where N is the
    total number of apertures).
    For example,
    first aperture radius will be RAPERT and corresponding magnitude will be MAG,
    second aperture radius will be RAPERT2 and corresponding magnitude will be MAG2,
    third aperture radius will be RAPERT3 and corresponding magnitude will be MAG3,
    and so on.
    
    """
    _format_name = 'daophot'
    _io_registry_format_aliases = ['daophot']
    _io_registry_can_write = False
    _description = 'IRAF DAOphot format table'

    header_class = DaophotHeader
    data_class = DaophotData
    inputter_class = DaophotInputter

    table_width = 80

    def __init__(self):
        core.BaseReader.__init__(self)
        #The inputter needs to know about the data (see DaophotInputter.process_lines)
        self.inputter.data = self.data



    def write(self, table=None):
        raise NotImplementedError
