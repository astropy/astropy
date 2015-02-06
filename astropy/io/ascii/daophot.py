# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

daophot.py:
  Classes to read DAOphot table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""

from __future__ import absolute_import, division, print_function

import re
import numpy as np
from . import core
from . import basic
from . import fixedwidth
from ...utils import OrderedDict


class DaophotHeader(core.BaseHeader):
    comment = r'\s*#K'
    aperture_values = ''

    def update_meta(self, lines, meta):
        """
        Extract table-level keywords for DAOphot table.  These are indicated by
        a leading '#K ' prefix.
        """
        table_meta = meta['table']
        # Read keywords as a table embedded in the header comments
        comment_lines = [line for line in lines if line.startswith('#')]
        if len(comment_lines) > 0:
            re_header_keyword = re.compile(r'[#]K'
                                           r'\s+ (?P<name> \w+)'
                                           r'\s* = (?P<stuff> .+) $',
                                           re.VERBOSE)

            table_meta['keywords'] = OrderedDict()
            for line in comment_lines:
                m = re_header_keyword.match(line)
                if m:
                    vals = m.group('stuff').strip().rsplit(None, 2)
                    keyword_dict = {'units': vals[-2],
                                    'format': vals[-1]}
                    keyword_dict['value'] = (vals[0] if len(vals) > 2 else "")
                    table_meta['keywords'][m.group('name')] = keyword_dict
                    if m.group('name') == 'APERTURES':
                        self.aperture_values = keyword_dict['value']

    def get_cols(self, lines):
        """Initialize the header Column objects from the table ``lines`` for a DAOphot
        header.  The DAOphot header is specialized so that we just copy the entire BaseHeader
        get_cols routine and modify as needed.

        Parameters
        ----------
        lines : list
            List of table lines

        """

        # Parse a series of column definition lines like below.  There may be several
        # such blocks in a single file (where continuation characters have already been
        # stripped).
        # #N ID    XCENTER   YCENTER   MAG         MERR          MSKY           NITER
        # #U ##    pixels    pixels    magnitudes  magnitudes    counts         ##
        # #F %-9d  %-10.3f   %-10.3f   %-12.3f     %-14.3f       %-15.7g        %-6d
        coldef_lines = ['', '', '']
        starts = ('#N ', '#U ', '#F ')
        col_width = []
        col_len_def = re.compile(r'[0-9]+')
        re_colformat_def = re.compile(r'#F([^#]+)')
        last_coldef_line = ['', '', '']
        for line in lines:
            if not line.startswith('#'):
                break  # End of header lines
            else:
                formatmatch = re_colformat_def.search(line)
                if formatmatch:
                    form = formatmatch.group(1).split()
                    width = ([int(col_len_def.search(s).group()) for s in form])
                    # original data format might be shorter than 80 characters
                    # and filled with spaces
                    width[-1] = 80 - sum(width[:-1])
                    col_width.extend(width)
                    last_width = width
                for i, start in enumerate(starts):
                    if line.startswith(start):
                        line_stripped = line[2:]
                        coldef_lines[i] = coldef_lines[i] + line_stripped
                        last_coldef_line[i] = line_stripped
                        break

        # We need to check whether daophot file has multiple aperture data, in its keywords
        if (',' in self.aperture_values) or (':' in self.aperture_values):
            apertures=[]
            for aper in self.aperture_values.split(','):
                if ':' in aper:
                    # Generate list of apertures from daophot's closed interval range
                    # syntax ap1:apN:apstep
                    ap1, apN, apstep = (float(i) for i in aper.split(':'))
                    apertures.extend(list(np.arange(ap1, apN, apstep)))
                    if (apN-ap1)%apstep == 0:
                        apertures.append(apN)
                else:
                    apertures.append(float(aper))
            # We shall now append the last header multiple times
            for j in range(1, len(apertures)):
                col_width.extend(last_width)
                coldef_lines[0] = coldef_lines[0] + ' ' + ' '.join([name+str(j+1) for name in last_coldef_line[0].split()])
                for i in range(1, len(coldef_lines)):
                    coldef_lines[i] = coldef_lines[i] + last_coldef_line[i]

        # At this point colddef_lines has three lines corresponding to column
        # names, unit, and format.  Get the column names by splitting the
        # first line on whitespace.
        self.names = coldef_lines[0].split()
        if not self.names:
            raise core.InconsistentTableError('No column names found in DAOphot header')

        ends = np.cumsum(col_width)
        starts = ends - col_width

        # If there wasn't a #U defined (not sure of DAOphot specification), then
        # replace the empty line with the right number of ## indicators, which matches
        # the DAOphot "no unit" tag.
        for i, coldef_line in enumerate(coldef_lines):
            if not coldef_line:
                coldef_lines[i] = '## ' * len(self.names)

        # Read the three lines as a basic table.
        reader = core._get_reader(Reader=basic.Basic, comment=None)
        reader.header.comment = None
        coldefs = reader.read(coldef_lines)

        # Create the list of io.ascii column objects
        self._set_cols_from_names()

        # Set unit and format as needed.
        for col in self.cols:
            if coldefs[col.name][0] != '##':
                col.unit = coldefs[col.name][0]
            if coldefs[col.name][1] != '##':
                col.format = coldefs[col.name][1]

        # Set column start and end positions.
        for i, col in enumerate(self.cols):
            col.start = starts[i]
            col.end = ends[i]
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


class DaophotInputter(core.ContinuationLinesInputter):
    no_continue = r'\s*#'


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


    def write(self, table=None):
        raise NotImplementedError
