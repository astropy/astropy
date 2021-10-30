# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Classes to read AAS MRT table format

Ref: https://journals.aas.org/mrt-standards

:Copyright: Smithsonian Astrophysical Observatory (2021)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu), \
         Suyog Garg (suyog7130@gmail.com)
"""

import re
import math
import warnings
import numpy as np
from io import StringIO

from . import core
from . import fixedwidth, cds

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


class MrtSplitter(fixedwidth.FixedWidthSplitter):
    """
    Contains the join function to left align the MRT columns
    when writing to a file.
    """
    def join(self, vals, widths):
        vals = [val + ' ' * (width - len(val)) for val, width in zip(vals, widths)]
        return self.delimiter.join(vals)


class MrtHeader(cds.CdsHeader):
    _subfmt = 'MRT'

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
            lead = ''
            if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                col.formatted_width = maxent + maxdec + 1
                if sign:
                    col.formatted_width += 1
            elif col.format.startswith('0'):
                # Keep leading zero, if already set in format - primarily for `seconds` columns
                # in coordinates; may need extra case if this is to be also supported with `sign`.
                lead = '0'
            col.fortran_format = fformat + str(col.formatted_width) + "." + str(maxdec)
            col.format = lead + col.fortran_format[1:] + "f"

    def write_byte_by_byte(self):
        """
        Writes the Byte-By-Byte description of the table.

        Columns that are `astropy.coordinates.SkyCoord` or `astropy.time.TimeSeries`
        objects or columns with values that are such objects are recognized as such,
        and some predefined labels and description is used for them.
        See the Vizier MRT Standard documentation in the link below for more details
        on these. An example Byte-By-Byte table is shown here.

        See: http://vizier.u-strasbg.fr/doc/catstd-3.1.htx

        Example::

        --------------------------------------------------------------------------------
        Byte-by-byte Description of file: table.dat
        --------------------------------------------------------------------------------
        Bytes Format Units  Label     Explanations
        --------------------------------------------------------------------------------
         1- 8  A8     ---    names   Description of names
        10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e
        16-23  F8.5   ---    d       [22.25/27.25] Description of d
        25-31  E7.1   ---    s       [-9e+34/2.0] Description of s
        33-35  I3     ---    i       [-30/67] Description of i
        37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF
        41-42  I2     ---    sameI   [20] Description of sameI
        44-45  I2     h      RAh     Right Ascension (hour)
        47-48  I2     min    RAm     Right Ascension (minute)
        50-67  F18.15 s      RAs     Right Ascension (second)
           69  A1     ---    DE-     Sign of Declination
        70-71  I2     deg    DEd     Declination (degree)
        73-74  I2     arcmin DEm     Declination (arcmin)
        76-91  F16.13 arcsec DEs     Declination (arcsec)

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

            if col.format is not None:
                col.formatted_width = max([len(sval) for sval in col.str_vals])

            # Set MRTColumn type, size and format.
            if np.issubdtype(col.dtype, np.integer):
                # Integer formatter
                self._set_column_val_limits(col)
                if getattr(col, 'formatted_width', None) is None:  # If ``formats`` not passed.
                    col.formatted_width = max(len(str(col.max)), len(str(col.min)))
                col.fortran_format = "I" + str(col.formatted_width)
                if col.format is None:
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
            # the MRT writer.
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
                    not any(x in col.name for x in ['RA', 'DE', 'LON', 'LAT', 'PLN', 'PLT'])):
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
            if col.name == 'DEd':
                bbb.add_row([singlebfmt.format(startb),
                             "A1", "---", "DE-",
                             "Sign of Declination"])
                col.fortran_format = 'I2'
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
        Writes the Header of the MRT table, aka ReadMe, which
        also contains the Byte-By-Byte description of the table.
        """
        from astropy.coordinates import SkyCoord

        # Recognised ``SkyCoord.name`` forms with their default column names (helio* require SunPy).
        coord_systems = {'galactic': ('GLAT', 'GLON', 'b', 'l'),
                         'ecliptic': ('ELAT', 'ELON', 'lat', 'lon'),      # 'geocentric*ecliptic'
                         'heliographic': ('HLAT', 'HLON', 'lat', 'lon'),  # '_carrington|stonyhurst'
                         'helioprojective': ('HPLT', 'HPLN', 'Ty', 'Tx')}
        eqtnames = ['RAh', 'RAm', 'RAs', 'DEd', 'DEm', 'DEs']

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

            # Replace single ``SkyCoord`` column by its coordinate components if no coordinate
            # columns of the correspoding type exist yet.
            if isinstance(col, SkyCoord):
                # If coordinates are given in RA/DEC, divide each them into hour/deg,
                # minute/arcminute, second/arcsecond columns.
                if ('ra' in col.representation_component_names.keys() and
                        len(set(eqtnames) - set(self.colnames)) == 6):
                    ra_c, dec_c = col.ra.hms, col.dec.dms
                    coords = [ra_c.h.round().astype('i1'), ra_c.m.round().astype('i1'), ra_c.s,
                              dec_c.d.round().astype('i1'), dec_c.m.round().astype('i1'), dec_c.s]
                    coord_units = [u.h, u.min, u.second,
                                   u.deg, u.arcmin, u.arcsec]
                    coord_descrip = ['Right Ascension (hour)', 'Right Ascension (minute)',
                                     'Right Ascension (second)', 'Declination (degree)',
                                     'Declination (arcmin)', 'Declination (arcsec)']
                    for coord, name, coord_unit, descrip in zip(
                            coords, eqtnames, coord_units, coord_descrip):
                        # Have Sign of Declination only in the DEd column.
                        if name in ['DEm', 'DEs']:
                            coord_col = Column(list(np.abs(coord)), name=name,
                                               unit=coord_unit, description=descrip)
                        else:
                            coord_col = Column(list(coord), name=name, unit=coord_unit,
                                               description=descrip)
                        # Set default number of digits after decimal point for the
                        # second values, and deg-min to (signed) 2-digit zero-padded integer.
                        if name == 'RAs':
                            coord_col.format = '013.10f'
                        elif name == 'DEs':
                            coord_col.format = '012.9f'
                        elif name == 'RAh':
                            coord_col.format = '2d'
                        elif name == 'DEd':
                            coord_col.format = '+03d'
                        elif name.startswith(('RA', 'DE')):
                            coord_col.format = '02d'
                        self.cols.append(coord_col)
                    to_pop.append(i)   # Delete original ``SkyCoord`` column.

                # For all other coordinate types, simply divide into two columns
                # for latitude and longitude resp. with the unit used been as it is.

                else:
                    frminfo = ''
                    for frame, latlon in coord_systems.items():
                        if frame in col.name and len(set(latlon[:2]) - set(self.colnames)) == 2:
                            if frame != col.name:
                                frminfo = f' ({col.name})'
                            lon_col = Column(getattr(col, latlon[3]), name=latlon[1],
                                             description=f'{frame.capitalize()} Longitude{frminfo}',
                                             unit=col.representation_component_units[latlon[3]],
                                             format='.12f')
                            lat_col = Column(getattr(col, latlon[2]), name=latlon[0],
                                             description=f'{frame.capitalize()} Latitude{frminfo}',
                                             unit=col.representation_component_units[latlon[2]],
                                             format='+.12f')
                            self.cols.append(lon_col)
                            self.cols.append(lat_col)
                            to_pop.append(i)   # Delete original ``SkyCoord`` column.

                # Convert all other ``SkyCoord`` columns that are not in the above three
                # representations to string valued columns. Those could either be types not
                # supported yet (e.g. 'helioprojective'), or already present and converted.
                # If there were any extra ``SkyCoord`` columns of one kind after the first one,
                # then their decomposition into their component columns has been skipped.
                # This is done in order to not create duplicate component columns.
                # Explicit renaming of the extra coordinate component columns by appending some
                # suffix to their name, so as to distinguish them, is not yet implemented.
                if i not in to_pop:
                    warnings.warn(f"Coordinate system of type '{col.name}' already stored in table "
                                  f"as CDS/MRT-syle columns or of unrecognized type. So column {i} "
                                  f"is being skipped with designation of a string valued column "
                                  f"`{self.colnames[i]}`.", UserWarning)
                    self.cols.append(Column(col.to_string(), name=self.colnames[i]))
                    to_pop.append(i)   # Delete original ``SkyCoord`` column.

            # Convert all other ``mixin`` columns to ``Column`` objects.
            # Parsing these may still lead to errors!
            elif not isinstance(col, Column):
                col = Column(col)
                # If column values are ``object`` types, convert them to string.
                if np.issubdtype(col.dtype, np.dtype(object).type):
                    col = Column([str(val) for val in col])
                self.cols[i] = col

        # Delete original ``SkyCoord`` columns, if there were any.
        for i in to_pop[::-1]:
            self.cols.pop(i)

        # Check for any left over extra coordinate columns.
        if any(x in self.colnames for x in ['RAh', 'DEd', 'ELON', 'GLAT']):
            # At this point any extra ``SkyCoord`` columns should have been converted to string
            # valued columns, together with issuance of a warning, by the coordinate parser above.
            # This test is just left here as a safeguard.
            for i, col in enumerate(self.cols):
                if isinstance(col, SkyCoord):
                    self.cols[i] = Column(col.to_string(), name=self.colnames[i])
                    message = ('Table already has coordinate system in CDS/MRT-syle columns. '
                               f'So column {i} should have been replaced already with '
                               f'a string valued column `{self.colnames[i]}`.')
                    raise core.InconsistentTableError(message)

        # Get Byte-By-Byte description and fill the template
        bbb_template = Template('\n'.join(BYTE_BY_BYTE_TEMPLATE))
        byte_by_byte = bbb_template.substitute({'file': 'table.dat',
                                                'bytebybyte': self.write_byte_by_byte()})

        # Fill up the full ReadMe
        rm_template = Template('\n'.join(MRT_TEMPLATE))
        readme_filled = rm_template.substitute({'bytebybyte': byte_by_byte})
        lines.append(readme_filled)


class MrtData(cds.CdsData):
    """MRT table data reader
    """
    _subfmt = 'MRT'
    splitter_class = MrtSplitter

    def write(self, lines):
        self.splitter.delimiter = ' '
        fixedwidth.FixedWidthData.write(self, lines)


class Mrt(core.BaseReader):
    """AAS MRT (Machine-Readable Table) format table.

    **Reading**
    ::

      >>> from astropy.io import ascii
      >>> table = ascii.read('data.mrt', format='mrt')

    **Writing**

    Use ``ascii.write(table, 'data.mrt', format='mrt')`` to  write tables to
    Machine Readable Table (MRT) format.

    Note that the metadata of the table, apart from units, column names and
    description, will not be written. These have to be filled in by hand later.

    See also: :ref:`cds_mrt_format`.

    Caveats:

    * The Units and Explanations are available in the column ``unit`` and
      ``description`` attributes, respectively.
    * The other metadata defined by this format is not available in the output table.
    """
    _format_name = 'mrt'
    _io_registry_format_aliases = ['mrt']
    _io_registry_can_write = True
    _description = 'MRT format table'

    data_class = MrtData
    header_class = MrtHeader

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
