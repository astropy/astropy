# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer.

mrt.py:
  Classes to read and write AAS Machine Readable Table format

:Copyright: Smithsonian Astrophysical Observatory (2011)
:Author: Tom Aldcroft (aldcroft@head.cfa.harvard.edu), \
         Suyog Garg (suyog7130@gmail.com)
"""


import warnings
import numpy as np

from . import core
from .cds import CdsHeader, CdsData

from astropy import units as u

from astropy.table import Column
from string import Template

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


class MrtHeader(CdsHeader):

    def write(self, lines):
        """
        Writes the Header of the MRT table which
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


class MrtData(CdsData):
    """MRT table data reader
    """


class Mrt(core.BaseReader):
    """MRT format table.

    See: FIXME

    Example::

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

      !!!FIXME!!!

    **Writing**

    Use ``ascii.write(format='mrt')`` to  write tables to Machine Readable Table (MRT)
    format, which is currently supported. MRT format differs slightly from the CDS
    format in table description sections of the ReadMe. It also has both the data and
    the ReadMe in a single file.

    Note that the metadata of the table, apart from units, column names and description,
    will not be written. These have to be filled in by hand later.

    See also: :ref:`cds_mrt_format`

    Caveats:

    * The Units and Explanations are available in the column ``unit`` and
      ``description`` attributes, respectively.
    * The other metadata defined by this format is not available in the output table.
    """
    _format_name = 'mrt'
    _io_registry_format_aliases = ['mrt']
    _description = 'CDS format table'

    data_class = MrtData
    header_class = MrtHeader

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
