# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings
from collections import defaultdict, OrderedDict

import numpy as np

from . import Header, Card

from ...coordinates import EarthLocation
from ...table import Column
from ...time import Time
from ...utils.exceptions import AstropyUserWarning


class FITS_time(object):
    """A class to read FITS binary table time columns as `~astropy.time.Time`

    This class reads the metadata associated with time coordinates, as 
    stored in a FITS binary table header, converts time columns into 
    `~astropy.time.Time` columns and reads global reference times as 
    `~astropy.time.Time` instances.
    """

    # Global time reference coordinate keywords
    TIME_KEYWORDS = {'TIMESYS' : 'scale',
                     'MJDREF' : 'mjdref',
                     'JDREF' : 'jdref',
                     'DATEREF' : 'dateref',
                     'TREFPOS' : 'pos',
                     'TREFDIR' : 'dir',
                     'TIMEUNIT' : 'unit',
                     'TIMEOFFS' : 'offs',
                     'OBSGEO-X' : 'loc_x',
                     'OBSGEO-Y' : 'loc_y',
                     'OBSGEO-Z' : 'loc_z',
                     'DATE' : 'date',
                     'DATE-OBS' : 'date_obs',
                     'DATE-END' : 'date_end'}

    # Column-specific time override keywords
    COLUMN_TIME_KEYWORDS = {'TCTYP' : 'scale',
                            'TRPOS' : 'pos',
                            'TCUNI' : 'unit'}

    # Set AstroPy Time global information
    GLOBAL_TIME_INFO = {'TIMESYS' : ('UTC','Default time scale'),
                        'JDREF' : (0.0,'Time columns are jd = jd1 + jd2'),
                        'TREFPOS' : ('TOPOCENTER','Time reference position')}

    # Compatibility of Time Scales and Reference Positions
    TIME_SCALE_REF = {'tai': 'TOPOCENTER',
                      'tt' : 'TOPOCENTER',
                      'ut1' : 'TOPOCENTER',
                      'utc' : 'TOPOCENTER',
                      'tcg' : 'GEOCENTER',
                      'tcb' : 'BARYCENTER',
                      'tdb' : 'BARYCENTER'}

    def __init__(self):
        # Set defaults for global time scale, reference, etc.
        self.global_info = {'scale' : 'UTC',
                            'mjdref' : None,
                            'unit' : 's'}
        # Set default dictionary for time columns
        self.time_columns = defaultdict(OrderedDict)

    @classmethod
    def is_time_column_keyword(cls, keyword):
        """
        Check if the FITS header keyword is a time column-specific keyword.

        Parameters
        ----------
        keyword : str
            FITS keyword.
        """
        try:
            key, idx = re.match(r'([A-Z]+)([0-9]+)', keyword).groups()
            if key in cls.COLUMN_TIME_KEYWORDS:
                return True
            return False
        except:
            return False

    def set_global_time(self, keyword, value, comment=None):
        """
        Set the global time reference frame attributes.

        Parameters
        ----------
        key : str
            FITS global time reference frame keyword.
        value : int, str, list
            value associated with specified keyword.
        comment : str
            comment associated with specified keyword.
        """
        if keyword in self.TIME_KEYWORDS:
            self.global_info[self.TIME_KEYWORDS[keyword]] = value
        else:
            raise ValueError('Illegal global time keyword: {} = {}'.format(keyword, value))

    def set_column_override(self, keyword, value, comment=None):
        """
        Set the time column specific override attributes.

        Parameters
        ----------
        key : str
            FITS time column specific keyword.
        value : int, str, list
            value associated with specified keyword.
        comment : str
            comment associated with specified keyword.
        """
        if self.is_time_column_keyword(keyword):
            key, idx = re.match(r'([A-Z]+)([0-9]+)', keyword).groups()
            self.time_columns[int(idx)][self.COLUMN_TIME_KEYWORDS[key]] = value
        else:
            raise ValueError('Illegal column-specific time keyword: {} = {}'.format(keyword, value))

    def convert_time_columns(self, table):
        """
        Convert time columns to Astropy Time columns.

        Parameters
        ----------
        table : astropy.table.Table
            The table whose time columns are to be converted.
        """
        for idx, time_col in self.time_columns.items():
            time_colname = table.colnames[idx - 1]
            table[time_colname] = Time(table[time_colname][...,0], table[time_colname][...,1],
                                       format='jd', scale=time_col['scale'].lower())
            # Add location if attribute `pos` is set
            try:
                if time_col['pos'] == 'TOPOCENTER':
                    table[time_colname].location = EarthLocation(self.global_info['loc_x'],
                                                                 self.global_info['loc_y'],
                                                                 self.global_info['loc_z'], 
                                                                 unit='m')
            except:
                pass

    def read_time(self, hdr, table):
        """
        Read time coordinates in a FITS Binary Table as Time.

        Parameters
        ----------
        hdr : `~astropy.io.fits.header.Header`
            FITS Header
        table : astropy.table.Table
            The table whose time columns are to be read as Time
        Returns
        -------
        hdr : `~astropy.io.fits.header.Header`
            Modified FITS Header (time metadata removed)
        """
        # Make a "copy" (not just a view) of the input header, since it
        # may get modified.  the data is still a "view" (for now)
        hcopy = hdr.copy(strip=True)

        for key, value, comment in hdr.cards:
            if (key.upper() in self.TIME_KEYWORDS):

                self.set_global_time(key, value, comment)
                hcopy.remove(key)

            elif (self.is_time_column_keyword(key.upper())):

                self.set_column_override(key, value, comment)
                hcopy.remove(key)

        self.convert_time_columns(table)

        return hcopy

    @classmethod
    def replace_time_table(cls, table):
        """
        Replace Time columns in a Table with non-mixin columns containing
        each element as a vector of two doubles (jd1, jd2) and return a FITS 
        header with appropriate time coordinate keywords.
        jd = jd1 + jd2 represents time in the Julian Date format with 
        high-precision.

        Parameters
        ----------
        table : astropy.table.Table
            The table whose Time columns are to be replaced.

        Returns
        -------
        table : astropy.table.Table
            The table with replaced Time columns
        hdr : `~astropy.io.fits.header.Header`
            Header containing Cards associated with the FITS time coordinate
        """
        # Shallow copy of the input table
        newtable = table.copy(copy_data=False)

        # Global time coordinate frame keywords
        hdr = Header([Card(keyword=key, value=val[0], comment=val[1]) 
                      for key, val in cls.GLOBAL_TIME_INFO.items()])

        time_cols = table.columns.isinstance(Time)

        # Geocentric Position
        pos_geo = None

        for col in time_cols:
            jd12 = np.array([col.jd1, col.jd2])
            jd12 = np.rollaxis(jd12, 0, jd12.ndim)
            newtable.replace_column(col.info.name, Column(jd12, unit='d'))

            # Get column position(index)
            n = table.colnames.index(col.info.name) + 1

            # Time column override keywords
            hdr.append(Card(keyword='TCTYP{}'.format(n), value=col.scale.upper()))

            # Time column reference positions
            if col.location is None:
                if pos_geo is not None:
                    warnings.warn(
                        'Time Column "{}" has no specified location, but global Time Position '
                        'is present, which will be the default for this column in '
                        'FITS specification.'.format(col.info.name),
                        AstropyUserWarning)
            else:
                hdr.append(Card(keyword='TRPOS{}'.format(n), value='TOPOCENTER'))
                # Compatibility of Time Scales and Reference Positions
                if cls.TIME_SCALE_REF[col.scale] != 'TOPOCENTER':
                    warnings.warn(
                        'Earth Location "TOPOCENTER" for Time Column "{}" is incompatabile '
                        'with scale "{}".'.format(col.info.name, col.scale.upper()),
                        AstropyUserWarning)
                if col.location.size > 1:
                    raise ValueError('Vectorized Location of Time Column "{}" cannot be written, '
                                     'as it is not yet supported'.format(col.info.name))
                if pos_geo is None:
                    hdr.extend([Card(keyword='OBSGEO-'+ dim.upper(), value=getattr(col.location, dim).value)
                                for dim in ('x', 'y', 'z')])
                    pos_geo = col.location
                elif pos_geo != col.location:
                    raise ValueError('Multiple Time Columns with different geocentric observatory '
                                     'locations ({}, {}, {}) , ({}, {}, {}) are encountered. '
                                     'These are not yet supported.'.format(pos_geo.x, pos_geo.y,
                                     pos_geo.z, col.location.x, col.location.y, col.location.z))

        return newtable, hdr
