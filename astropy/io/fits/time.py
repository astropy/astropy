# Licensed under a 3-clause BSD style license - see LICENSE.rst

from collections import OrderedDict

import numpy as np

from . import Header, Card

from ...coordinates import EarthLocation
from ...table import Table, QTable, Column
from ...time import Time
from ...time.formats import FITS_DEPRECATED_SCALES
from ...units import Quantity


class FitsTime(object):
    """
    A class to maintain the metadata of time columns stored in a FITS binary table
    (header and data) and to convert those time columns into `~astropy.time.Time` columns.

    Parameters
    ----------
    nfields : int
        Number of fields(columns) in the binary table.
    """

    def __init__(self, nfields):
        """
        Set default time metadata.
        """
        # Default time scale UTC
        self.scale = 'UTC'
        # If none of the reference values are provided, MJD = 0.0 must be assumed.
        self.mjd = None
        # Default time unit
        self.unit = 's'
        self.time_columns = OrderedDict()
        for i in range(1, nfields + 1):
            self.time_columns[i] = None

    def set_global_time(self, key, value, comment):
        """
        Set the global time reference frame attributes.

        Parameters
        ----------
        key : string
            FITS global time reference frame keyword.
        value : int, str, list
            value associated with specified keyword.
        comment : str
            comment associated with specified keyword.
        """
        if key == 'TIMESYS':
            self.scale = value
        elif key == 'MJDREF':
            self.mjd = value
        elif key == 'JDREF':
            self.jd = value
        elif key == 'DATEREF':
            self.date = value
        elif key == 'TREFPOS':
            self.pos = value
        elif key == 'TREFDIR':
            self.dir = value
        elif key == 'TIMEUNIT':
            self.unit = value
        elif key == 'TIMEOFFS':
            self.offs = value
        elif key == 'OBSGEO-X':
            self.loc_x = value
        elif key == 'OBSGEO-Y':
            self.loc_y = value
        elif key == 'OBSGEO-Z':
            self.loc_z = value

    def set_column_override(self, key, value, comment):
        """
        Set the time column specific override attributes.

        Parameters
        ----------
        key : string
            FITS time column specific keyword.
        value : int, str, list
            value associated with specified keyword.
        comment : str
            comment associated with specified keyword.
        """
        idx = int(key[-1])
        if self.time_columns[idx] is None:
            self.time_columns[idx] = OrderedDict()
        if key[:-1] == 'TCTYP':
            self.time_columns[idx]['scale'] = value
        elif key[:-1] == 'TRPOS':
            self.time_columns[idx]['pos'] = value
        elif key[:-1] == 'TCUNI':
            self.time_columns[idx]['unit'] = value
        elif key[:-1] == 'TCAPF':
            self.time_columns[idx]['APy_format'] = value

    def convert_to_time(self, table):
        """
        Convert time columns to Astropy Time columns.

        Parameters
        ----------
        table : astropy.table.Table
            The table whose time columns are to be converted.
        """
        for idx, time_col in self.time_columns.items():
            if time_col is not None:
                time_colname = table.colnames[idx - 1]
                if time_col['APy_format'] is not None:
                    table[time_colname] = Time(table[time_colname][:,0], table[time_colname][:,1], format='jd', scale=time_col['scale'].lower())
                    table[time_colname].format = time_col['APy_format'].lower()
                    table[time_colname].location = EarthLocation(Quantity(self.loc_x, unit='m'), Quantity(self.loc_y, unit='m'), Quantity(self.loc_z, unit='m'))
                else:
                    # Still have to complete this to read FITS files not written by astropy
                    pass


def replace_time_table(table):
    """
    Replace Time columns in a Table with non-mixin columns containing
    each element as a vector of two doubles (jd1, jd2) and return a FITS header
    with appropriate time coordinate keywords.
    jd = jd1 + jd2 represents time in the Julian Date format with high-precision.

    Parameters
    ----------
    table : Table
        The table whose Time columns are to be replaced.

    Returns
    -------
    table : Table
        The table with replaced Time columns
    hdr : Header
        `astropy.io.fits` Header object containing Cards associated with the FITS time coordinate
    """
    # Global time coordinate frame keywords
    hdr = Header([Card(keyword='TIMESYS', value='UTC', comment='Default time scale'),
                  Card(keyword='JDREF', value =0.0, comment='Time columns are jd = jd1 + jd2'),
                  Card(keyword='TREFPOS', value='TOPOCENTER', comment='Time reference position')])

    time_cols = table.columns.isinstance(Time)
    for col in time_cols:
            table[col.info.name] = Column(np.empty(col.shape + (2,)))
            table[col.info.name][...,0] = col.jd1
            table[col.info.name][...,1] = col.jd2

            # Get column position(index)
            n = table.get_pos(col.info.name)

            # Time column override keywords
            hdr.append(Card(keyword='TCTYP%d' %n, value=col.scale.upper()))
            hdr.append(Card(keyword='TCUNI%d' %n, value='d'))
            
            # Time column reference positions 
            if col.scale in ('tai', 'tt', 'ut1', 'utc'):
                hdr.append(Card(keyword='TRPOS%d' %n, value='TOPOCENTER'))
                hdr.append(Card(keyword='OBSGEO-X', value=col.location.x.value))
                hdr.append(Card(keyword='OBSGEO-Y', value=col.location.y.value))
                hdr.append(Card(keyword='OBSGEO-Z', value=col.location.z.value))
            elif col.scale == 'tcg':
                hdr.append(Card(keyword='TRPOS%d' %n, value='GEOCENTER'))
            else:
                hdr.append(Card(keyword='TRPOS%d' %n, value='BARYCENTER'))

            # Astropy specific keyword for storing Time format
            hdr.append(Card(keyword='TCAPF%d' %n, value=col.format.upper()))

    return table, hdr
