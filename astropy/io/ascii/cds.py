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
import warnings
import numpy as np
from io import StringIO
from contextlib import suppress

from . import core
from . import fixedwidth

from astropy import units as u

from astropy.table import Table
from astropy.table import Column, MaskedColumn
from string import Template
from textwrap import wrap

MAX_SIZE_README_LINE = 80
MAX_COL_INTLIMIT = 100000


__doctest_skip__ = ['*']


BYTE_BY_BYTE_TEMPLATE = [
    "Byte-by-byte Description of file: $file",
    "--------------------------------------------------------------------------------",
    " Bytes Format Units  Label     Explanations",
    "--------------------------------------------------------------------------------",
    "$bytebybyte",
    "--------------------------------------------------------------------------------"]

MRT_TEMPLATE = [
    "Title:",
    "Authors:",
    "Table:",
    "================================================================================",
    "$bytebybyte",
    "Notes:",
    "--------------------------------------------------------------------------------"]


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
                    col.unit = None  # "---" is the marker for no unit in CDS table
                else:
                    col.unit = u.Unit(unit, format='cds', parse_strict='warn')
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

    def _split_float_format(self, value):
        """
        Splits a Float string into different parts to find number
        of digits after decimal and check if the value is in Scientific
        notation.

        Parameters
        ----------
        value : str
            String containing the float value to split.

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
            raise Exception(f'{value} is not a float number')
        return (len(value),
                len(mo.group('ent')),
                len(mo.group('decimals')),
                mo.group('sign') != "",
                mo.group('exp') != "")

    def _set_column_val_limits(self, col):
        """
        Sets the ``col.min`` and ``col.max`` column attributes,
        taking into account columns with Null values.
        """
        col.max = max(col)
        col.min = min(col)
        if col.max is np.ma.core.MaskedConstant:
            col.max = None
        if col.min is np.ma.core.MaskedConstant:
            col.min = None

    def column_float_formatter(self, col):
        """
        String formatter function for a column containing Float values.
        Checks if the values in the given column are in Scientific notation,
        by spliting the value string. It is assumed that the column either has
        float values or Scientific notation.

        A ``col.formatted_width`` attribute is added to the column. It is not added
        if such an attribute is already present, say when the ``formats`` argument
        is passed to the writer. A properly formatted format string is also added as
        the ``col.format`` attribute.

        Parameters
        ----------
        col : A ``Table.Column`` object.
        """
        # maxsize: maximum length of string containing the float value.
        # maxent: maximum number of digits places before decimal point.
        # maxdec: maximum number of digits places after decimal point.
        # maxprec: maximum precision of the column values, sum of maxent and maxdec.
        maxsize, maxprec, maxent, maxdec = 1, 0, 1, 0
        sign = False
        fformat = 'F'

        # Find maximum sized value in the col
        for val in col.str_vals:
            # Skip null values
            if val is None or val == '':
                continue

            # Find format of the Float string
            fmt = self._split_float_format(val)
            # If value is in Scientific notation
            if fmt[4] is True:
                # if the previous column value was in normal Float format
                # set maxsize, maxprec and maxdec to default.
                if fformat == 'F':
                    maxsize, maxprec, maxdec = 1, 0, 0
                # Designate the column to be in Scientific notation.
                fformat = 'E'
            else:
                # Move to next column value if
                # current value is not in Scientific notation
                # but the column is designated as such because
                # one of the previous values was.
                if fformat == 'E':
                    continue

            if maxsize < fmt[0]:
                maxsize = fmt[0]
            if maxent < fmt[1]:
                maxent = fmt[1]
            if maxdec < fmt[2]:
                maxdec = fmt[2]
            if fmt[3]:
                sign = True

            if maxprec < fmt[1] + fmt[2]:
                maxprec = fmt[1] + fmt[2]

        if fformat == 'E':
            if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                col.formatted_width = maxsize
                if sign:
                    col.formatted_width += 1
            # Number of digits after decimal is replaced by the precision
            # for values in Scientific notation, when writing that Format.
            col.fortran_format = fformat + str(col.formatted_width) + "." + str(maxprec)
            col.format = str(col.formatted_width) + "." + str(maxdec) + "e"
        else:
            if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                col.formatted_width = maxent + maxdec + 1
                if sign:
                    col.formatted_width += 1
            col.fortran_format = fformat + str(col.formatted_width) + "." + str(maxdec)
            col.format = col.fortran_format[1:] + "f"

    def write_byte_by_byte(self):
        """
        Writes the Byte-By-Byte description of the table.

        Columns that are `astropy.coordinates.SkyCoord` or `astropy.time.TimeSeries`
        objects or columns with values that are such objects are recognized as such,
        and some predefined labels and description is used for them.
        See the Vizier CDS Standard documentation in the link below for more details
        on these. An example Byte-By-Byte table is shown here.

        See: http://vizier.u-strasbg.fr/doc/catstd-3.1.htx

        Example::

        --------------------------------------------------------------------------------
        Byte-by-byte Description of file: table.dat
        --------------------------------------------------------------------------------
        Bytes Format Units  Label     Explanations
        --------------------------------------------------------------------------------
         1-  8  A8     ---    names   Description of names
        10- 14  E5.1   ---    e       [-3160000.0/0.01] Description of e
        16- 23  F8.5   ---    d       [22.25/27.25] Description of d
        25- 31  E7.1   ---    s       [-9e+34/2.0] Description of s
        33- 35  I3     ---    i       [-30/67] Description of i
        37- 39  F3.1   ---    sameF   [5.0/5.0] Description of sameF
        41- 42  I2     ---    sameI   [20] Description of sameI
        44- 47  F4.1   h      RAh     Right Ascension (hour)
        49- 51  F3.1   min    RAm     Right Ascension (minute)
        53- 70  F18.15 s      RAs     Right Ascension (second)
            72  A1     ---    DE-     Sign of Declination
        73- 76  F5.1   deg    DEd     Declination (degree)
        78- 81  F4.1   arcmin DEm     Declination (arcmin)
        83- 98  F16.13 arcsec DEs     Declination (arcsec)

        --------------------------------------------------------------------------------
        """
        # Get column widths
        vals_list = []
        col_str_iters = self.data.str_vals()
        for vals in zip(*col_str_iters):
            vals_list.append(vals)

        for i, col in enumerate(self.cols):
            col.width = max([len(vals[i]) for vals in vals_list])
            if self.start_line is not None:
                col.width = max(col.width, len(col.info.name))
        widths = [col.width for col in self.cols]

        startb = 1  # Byte count starts at 1.

        # Set default width of the Bytes count column of the Byte-By-Byte table.
        # This ``byte_count_width`` value helps align byte counts with respect
        # to the hyphen using a format string.
        byte_count_width = len(str(sum(widths) + len(self.cols) - 1))

        # Format string for Start Byte and End Byte
        singlebfmt = "{:" + str(byte_count_width) + "d}"
        fmtb = singlebfmt + "-" + singlebfmt
        # Add trailing single whitespaces to Bytes column for better visibility.
        singlebfmt += " "
        fmtb += " "

        # Set default width of Label and Description Byte-By-Byte columns.
        max_label_width, max_descrip_size = 7, 16

        bbb = Table(names=['Bytes', 'Format', 'Units', 'Label', 'Explanations'],
                    dtype=[str] * 5)

        # Iterate over the columns to write Byte-By-Byte rows.
        for i, col in enumerate(self.cols):
            # Check if column is MaskedColumn
            col.has_null = isinstance(col, MaskedColumn)

            # Check if the column format was set by the passed ``formats`` argument.
            if col.format is None:
                if col.info.format is not None:
                    col.format = col.info.format

            if col.format is not None:
                col.formatted_width = max([len(sval) for sval in col.str_vals])

            # Set CDSColumn type, size and format.
            if np.issubdtype(col.dtype, np.integer):
                # Integer formatter
                self._set_column_val_limits(col)
                if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                    col.formatted_width = max(len(str(col.max)), len(str(col.min)))
                col.fortran_format = "I" + str(col.formatted_width)
                col.format = ">" + col.fortran_format[1:]

            elif np.issubdtype(col.dtype, np.dtype(float).type):
                # Float formatter
                self._set_column_val_limits(col)
                self.column_float_formatter(col)

            else:
                # String formatter, ``np.issubdtype(col.dtype, str)`` is ``True``.
                dtype = col.dtype.str
                if col.has_null:
                    mcol = col
                    mcol.fill_value = ""
                    coltmp = Column(mcol.filled(), dtype=str)
                    dtype = coltmp.dtype.str
                if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                    col.formatted_width = int(re.search(r'(\d+)$', dtype).group(1))
                col.fortran_format = "A" + str(col.formatted_width)
                col.format = str(col.formatted_width) + "s"

            endb = col.formatted_width + startb - 1

            # ``mixin`` columns converted to string valued columns will not have a name
            # attribute. In those cases, a ``Unknown`` column label is put, indicating that
            # such columns can be better formatted with some manipulation before calling
            # the CDS/MRT writer.
            if col.name is None:
                col.name = "Unknown"

            # Set column description.
            if col.description is not None:
                description = col.description
            else:
                description = "Description of " + col.name

            # Set null flag in column description
            nullflag = ""
            if col.has_null:
                nullflag = "?"

            # Set column unit
            if col.unit is not None:
                col_unit = col.unit.to_string("cds")
            elif col.name.lower().find("magnitude") > -1:
                # ``col.unit`` can still be ``None``, if the unit of column values
                # is ``Magnitude``, because ``astropy.units.Magnitude`` is actually a class.
                # Unlike other units which are instances of ``astropy.units.Unit``,
                # application of the ``Magnitude`` unit calculates the logarithm
                # of the values. Thus, the only way to check for if the column values
                # have ``Magnitude`` unit is to check the column name.
                col_unit = "mag"
            else:
                col_unit = "---"

            # Add col limit values to col description
            lim_vals = ""
            if (col.min and col.max and
                    not any(x in col.name for x in ['RA', 'DE', 'LON', 'LAT'])):
                # No col limit values for coordinate columns.
                if col.fortran_format[0] == 'I':
                    if abs(col.min) < MAX_COL_INTLIMIT and abs(col.max) < MAX_COL_INTLIMIT:
                        if col.min == col.max:
                            lim_vals = "[{0}]".format(col.min)
                        else:
                            lim_vals = "[{0}/{1}]".format(col.min, col.max)
                elif col.fortran_format[0] in ('E', 'F'):
                    lim_vals = "[{0}/{1}]".format(math.floor(col.min * 100) / 100.,
                                                  math.ceil(col.max * 100) / 100.)

            if lim_vals != '' or nullflag != '':
                description = "{0}{1} {2}".format(lim_vals, nullflag, description)

            # Find the maximum label and description column widths.
            if len(col.name) > max_label_width:
                max_label_width = len(col.name)
            if len(description) > max_descrip_size:
                max_descrip_size = len(description)

            # Add a row for the Sign of Declination in the bbb table
            if col.name == 'DEd' and col[0] < 0.0:
                bbb.add_row([singlebfmt.format(startb),
                             "A1", "---", "DE-",
                             "Sign of Declination"])
                startb += 1

            # Add Byte-By-Byte row to bbb table
            bbb.add_row([singlebfmt.format(startb) if startb == endb
                         else fmtb.format(startb, endb),
                         "" if col.fortran_format is None else col.fortran_format,
                         col_unit,
                         "" if col.name is None else col.name,
                         description])
            startb = endb + 2

        # Properly format bbb columns
        bbblines = StringIO()
        bbb.write(bbblines, format='ascii.fixed_width_no_header',
                  delimiter=' ', bookend=False, delimiter_pad=None,
                  formats={'Format': '<6s',
                           'Units': '<6s',
                           'Label': '<' + str(max_label_width) + 's',
                           'Explanations': '' + str(max_descrip_size) + 's'})

        # Get formatted bbb lines
        bbblines = bbblines.getvalue().splitlines()

        # ``nsplit`` is the number of whitespaces to prefix to long description
        # lines in order to wrap them. It is the sum of the widths of the
        # previous 4 columns plus the number of single spacing between them.
        # The hyphen in the Bytes column is also counted.
        nsplit = byte_count_width * 2 + 1 + 12 + max_label_width + 4

        # Wrap line if it is too long
        buff = ""
        for newline in bbblines:
            if len(newline) > MAX_SIZE_README_LINE:
                buff += ("\n").join(wrap(newline,
                                         subsequent_indent=" " * nsplit,
                                         width=MAX_SIZE_README_LINE))
                buff += "\n"
            else:
                buff += newline + "\n"

        # Last value of ``endb`` is the sum of column widths after formatting.
        self.linewidth = endb

        # Remove the last extra newline character from Byte-By-Byte.
        buff = buff[:-1]
        return buff

    def write(self, lines):
        """
        Writes the Header of the CDS table, aka ReadMe, which
        also contains the Byte-By-Byte description of the table.
        """
        from astropy.coordinates import SkyCoord

        # list to store indices of columns that are modified.
        to_pop = []

        # For columns that are instances of ``SkyCoord`` and other ``mixin`` columns
        # or whose values are objects of these classes.
        for i, col in enumerate(self.cols):
            # If col is a ``Column`` object but its values are ``SkyCoord`` objects,
            # convert the whole column to ``SkyCoord`` object, which helps in applying
            # SkyCoord methods directly.
            if not isinstance(col, SkyCoord) and isinstance(col[0], SkyCoord):
                try:
                    col = SkyCoord(col)
                except (ValueError, TypeError):
                    # If only the first value of the column is a ``SkyCoord`` object,
                    # the column cannot be converted to a ``SkyCoord`` object.
                    # These columns are converted to ``Column`` object and then converted
                    # to string valued column.
                    if not isinstance(col, Column):
                        col = Column(col)
                    col = Column([str(val) for val in col])
                    self.cols[i] = col
                    continue

            # Replace single ``SkyCoord`` column by its coordinate components.
            if isinstance(col, SkyCoord):
                # If coordinates are given in RA/DEC, divide each them into hour/deg,
                # minute/arcminute, second/arcsecond columns.
                if 'ra' in col.representation_component_names.keys():
                    ra_col, dec_col = col.ra.hms, col.dec.dms
                    coords = [ra_col.h, ra_col.m, ra_col.s,
                              dec_col.d, dec_col.m, dec_col.s]
                    names = ['RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs']
                    coord_units = [u.h, u.min, u.second,
                                   u.deg, u.arcmin, u.arcsec]
                    coord_descrip = ['Right Ascension (hour)', 'Right Ascension (minute)',
                                     'Right Ascension (second)', 'Declination (degree)',
                                     'Declination (arcmin)', 'Declination (arcsec)']
                    for coord, name, coord_unit, descrip in zip(
                            coords, names, coord_units, coord_descrip):
                        # Have Sign of Declination only in the DEd column.
                        if name in ['DEm', 'DEs']:
                            coord_col = Column(list(np.abs(coord)), name=name,
                                               unit=coord_unit, description=descrip)
                        else:
                            coord_col = Column(list(coord), name=name, unit=coord_unit,
                                               description=descrip)
                        # Set default number of digits after decimal point for the
                        # second values.
                        if name in ['RAs', 'DEs']:
                            coord_col.format = '.12f'
                        self.cols.append(coord_col)

                # For all other coordinate types, simply divide into two columns
                # for latitude and longitude resp. with the unit used been as it is.
                else:
                    # Galactic coordinates.
                    if col.name == 'galactic':
                        lon_col = Column(col.l, name='GLON',
                                         description='Galactic Longitude',
                                         unit=col.representation_component_units['l'],
                                         format='.12f')
                        lat_col = Column(col.b, name='GLAT',
                                         description='Galactic Latitude',
                                         unit=col.representation_component_units['b'],
                                         format='.12f')
                        self.cols.append(lon_col)
                        self.cols.append(lat_col)

                    # Ecliptic coordinates, can be any of various available.
                    elif 'ecliptic' in col.name:
                        lon_col = Column(col.lon, name='ELON',
                                         description='Ecliptic Longitude (' + col.name + ')',
                                         unit=col.representation_component_units['lon'],
                                         format='.12f')
                        lat_col = Column(col.lat, name='ELAT',
                                         description='Ecliptic Latitude (' + col.name + ')',
                                         unit=col.representation_component_units['lat'],
                                         format='.12f')
                        self.cols.append(lon_col)
                        self.cols.append(lat_col)

                    # Convert all other ``SkyCoord`` columns that are not in the above three
                    # representations to string valued columns.
                    else:
                        self.cols.append(Column(col.to_string()))

                to_pop.append(i)   # Delete original ``SkyCoord`` column.

            # Convert all other ``mixin`` columns to ``Column`` objects.
            # Parsing these may still lead to errors!
            elif not isinstance(col, Column):
                col = Column(col)
                # If column values are ``object`` types, convert them to string.
                if np.issubdtype(col.dtype, np.dtype(object).type):
                    col = Column([str(val) for val in col])
                self.cols[i] = col

        # Delete original ``SkyCoord`` column, if there were any.
        for i in to_pop:
            self.cols.pop(i)

        # Check for any left over extra coordinate columns.
        if any(x in self.colnames for x in ['RAh', 'DEd', 'ELON', 'GLAT']):
            # If there were any ``SkyCoord`` columns after the first one, then they would
            # have been skipped the division into their component columns. This is done in
            # order to not replace the data in the component columns already obtained.
            # Explicit renaming of the extra coordinate component columns by appending some
            # suffix to their name, so as to distinguish them, is not implemented.
            # Such extra ``SkyCoord`` columns are converted to string valued columns,
            # together with issuance of a warning.
            for i, col in enumerate(self.cols):
                if isinstance(col, SkyCoord):
                    self.cols[i] = Column(col.to_string())
                    message = 'Table already has coordinate system in CDS/MRT-syle columns.' \
                              + f' So column {i} is being skipped with designation' \
                              + ' of an `Unknown` string valued column.'
                    warnings.warn(message, UserWarning)

        # Get Byte-By-Byte description and fill the template
        bbb_template = Template('\n'.join(BYTE_BY_BYTE_TEMPLATE))
        byte_by_byte = bbb_template.substitute({'file': 'table.dat',
                                                'bytebybyte': self.write_byte_by_byte()})

        # Fill up the full ReadMe
        rm_template = Template('\n'.join(MRT_TEMPLATE))
        readme_filled = rm_template.substitute({'bytebybyte': byte_by_byte})
        lines.append(readme_filled)


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

    **Writing**

    Use ``ascii.write(format='cds')`` to  write tables to Machine Readable Table (MRT)
    format, which is currently supported. MRT format differs slightly from the CDS
    format in table description sections of the ReadMe. It also has both the data and
    the ReadMe in a single file.

    Note that the metadata of the table, apart from units, column names and description,
    will not be written. These have to be filled in by hand later.

    See also: :ref:`cds_mrt_format`

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
    _description = 'CDS format table'

    data_class = CdsData
    header_class = CdsHeader

    def __init__(self, readme=None):
        super().__init__()
        self.header.readme = readme

        # ``CdsData`` class inherits from ``FixedWidthData`` which has the
        # default ``start_line`` at 1. For CDS format writing start line
        # should be at 0.
        self.data.start_line = None

    def write(self, table=None):
        # Construct for writing empty table is not yet done.
        if len(table) == 0:
            raise NotImplementedError
        self.data.header = self.header
        self.header.position_line = None
        self.header.start_line = None
        # Create a copy of the ``table``, so that it the copy gets modified and
        # written to the file, while the original table remains as it is.
        table = table.copy()
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
