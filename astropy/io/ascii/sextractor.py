# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" sextractor.py:
  Classes to read SExtractor table format

Built on daophot.py:
:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu)
"""


import re

from . import core


class SExtractorHeader(core.BaseHeader):
    """Read the header from a file produced by SExtractor."""
    comment = r'^\s*#\s*\S\D.*'  # Find lines that don't have "# digit"

    def get_cols(self, lines):
        """
        Initialize the header Column objects from the table ``lines`` for a SExtractor
        header.  The SExtractor header is specialized so that we just copy the entire BaseHeader
        get_cols routine and modify as needed.

        Parameters
        ----------
        lines : list
            List of table lines

        """

        # This assumes that the columns are listed in order, one per line with a
        # header comment string of the format: "# 1 ID short description [unit]"
        # However, some may be missing and must be inferred from skipped column numbers
        columns = {}
        # E.g. '# 1 ID identification number' (no units) or '# 2 MAGERR magnitude of error [mag]'
        # Updated along with issue #4603, for more robust parsing of unit
        re_name_def = re.compile(r"""^\s* \# \s*             # possible whitespace around #
                                 (?P<colnumber> [0-9]+)\s+   # number of the column in table
                                 (?P<colname> [-\w]+)        # name of the column
                                 # column description, match any character until...
                                 (?:\s+(?P<coldescr> \w .+)
                                 # ...until [non-space][space][unit] or [not-right-bracket][end]
                                 (?:(?<!(\]))$|(?=(?:(?<=\S)\s+\[.+\]))))?
                                 (?:\s*\[(?P<colunit>.+)\])?.* # match units in brackets
                                 """, re.VERBOSE)
        dataline = None
        for line in lines:
            if not line.startswith('#'):
                dataline = line  # save for later to infer the actual number of columns
                break                   # End of header lines
            else:
                match = re_name_def.search(line)
                if match:
                    colnumber = int(match.group('colnumber'))
                    colname = match.group('colname')
                    coldescr = match.group('coldescr')
                    colunit = match.group('colunit')  # If no units are given, colunit = None
                    columns[colnumber] = (colname, coldescr, colunit)
        # Handle skipped column numbers
        colnumbers = sorted(columns)
        # Handle the case where the last column is array-like by append a pseudo column
        # If there are more data columns than the largest column number
        # then add a pseudo-column that will be dropped later.  This allows
        # the array column logic below to work in all cases.
        if dataline is not None:
            n_data_cols = len(dataline.split())
        else:
            # handles no data, where we have to rely on the last column number
            n_data_cols = colnumbers[-1]
        # sextractor column number start at 1.
        columns[n_data_cols + 1] = (None, None, None)
        colnumbers.append(n_data_cols + 1)
        if len(columns) > 1:  # only fill in skipped columns when there is genuine column initially
            previous_column = 0
            for n in colnumbers:
                if n != previous_column + 1:
                    for c in range(previous_column + 1, n):
                        column_name = (columns[previous_column][0]
                                       + f"_{c - previous_column}")
                        column_descr = columns[previous_column][1]
                        column_unit = columns[previous_column][2]
                        columns[c] = (column_name, column_descr, column_unit)
                previous_column = n
        # Add the columns in order to self.names
        colnumbers = sorted(columns)[:-1]  # drop the pseudo column
        self.names = []
        for n in colnumbers:
            self.names.append(columns[n][0])

        if not self.names:
            raise core.InconsistentTableError('No column names found in SExtractor header')

        self.cols = []
        for n in colnumbers:
            col = core.Column(name=columns[n][0])
            col.description = columns[n][1]
            col.unit = columns[n][2]
            self.cols.append(col)


class SExtractorData(core.BaseData):
    start_line = 0
    delimiter = ' '
    comment = r'\s*#'


class SExtractor(core.BaseReader):
    """SExtractor format table.

    SExtractor is a package for faint-galaxy photometry (Bertin & Arnouts
    1996, A&A Supp. 317, 393.)

    See: https://sextractor.readthedocs.io/en/latest/

    Example::

      # 1 NUMBER
      # 2 ALPHA_J2000
      # 3 DELTA_J2000
      # 4 FLUX_RADIUS
      # 7 MAG_AUTO [mag]
      # 8 X2_IMAGE Variance along x [pixel**2]
      # 9 X_MAMA Barycenter position along MAMA x axis [m**(-6)]
      # 10 MU_MAX Peak surface brightness above background [mag * arcsec**(-2)]
      1 32.23222 10.1211 0.8 1.2 1.4 18.1 1000.0 0.00304 -3.498
      2 38.12321 -88.1321 2.2 2.4 3.1 17.0 1500.0 0.00908 1.401

    Note the skipped numbers since flux_radius has 3 columns.  The three
    FLUX_RADIUS columns will be named FLUX_RADIUS, FLUX_RADIUS_1, FLUX_RADIUS_2
    Also note that a post-ID description (e.g. "Variance along x") is optional
    and that units may be specified at the end of a line in brackets.

    """
    _format_name = 'sextractor'
    _io_registry_can_write = False
    _description = 'SExtractor format table'

    header_class = SExtractorHeader
    data_class = SExtractorData
    inputter_class = core.ContinuationLinesInputter

    def read(self, table):
        """
        Read input data (file-like object, filename, list of strings, or
        single string) into a Table and return the result.
        """
        out = super().read(table)
        # remove the comments
        if 'comments' in out.meta:
            del out.meta['comments']
        return out

    def write(self, table):
        raise NotImplementedError
