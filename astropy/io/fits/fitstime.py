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
                     'MJDREF' : 'mjd',
                     'JDREF' : 'jd',
                     'DATEREF' : 'date',
                     'TREFPOS' : 'pos',
                     'TREFDIR' : 'dir',
                     'TIMEUNIT' : 'unit',
                     'TIMEOFFS' : 'offs',
                     'OBSGEO-X' : 'loc_x',
                     'OBSGEO-Y' : 'loc_y',
                     'OBSGEO-Z' : 'loc_z'}

    # Column-specific time override keywords
    COLUMN_TIME_KEYWORDS = {'TCTYP' : 'scale',
                            'TRPOS' : 'pos',
                            'TCUNI' : 'unit'}

    # AstroPy-specific keywords
    ASTROPY_TIME_KEYWORDS = {'FORMAT' : 'format'}

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
                      'tcb' : 'BARYCENTER'}

    astropy_namespace = 'HIERARCH ASTROPY TIME'

    def __init__(self):
        # Set defaults for global time scale, reference, etc.
        self.global_info = {'scale' : 'UTC',
                            'mjd' : None,
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
        if keyword[5:].isdigit() and keyword[:5] in cls.COLUMN_TIME_KEYWORDS:
            return True
        return False

    def set_global_time(self, key, value, comment=None):
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
        if key in self.TIME_KEYWORDS:
            self.global_info[self.TIME_KEYWORDS[key]] = value
        else:
            raise ValueError('Illegal global time keyword.')

    def set_column_override(self, key, value, comment=None):
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
        if self.is_time_column_keyword(key):
            idx = int(key[5:])
            self.time_columns[idx][self.COLUMN_TIME_KEYWORDS[key[:5]]] = value
        else:
            raise ValueError('Illegal column-specific time keyword.')

    def set_astropy_time(self, key, value, comment=None):
        """
        Set the astropy namespace attributes. 

        Parameters
        ----------
        key : str
            FITS astropy time specific keyword.
        value : int, str, list
            value associated with specified keyword.
        comment : str
            comment associated with specified keyword.
        """
        if key.startswith('ASTROPY TIME'):
            sub_keys = re.split('(\d+)', key)
            if len(sub_keys) > 1 and sub_keys[-1] == '':
                idx = sub_keys[-2]
                self.time_columns[int(idx)][self.ASTROPY_TIME_KEYWORDS[key[13:-len(idx)]]] = value
        else:
            raise ValueError('Illegal astropy time keyword.')

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
            if time_col['format'] is not None:
                table[time_colname] = Time(table[time_colname][:,0], table[time_colname][:,1],
                                           format='jd', scale=time_col['scale'].lower())
                table[time_colname].format = time_col['format'].lower()
                try:
                    if time_col['pos'] == 'TOPOCENTER':
                        table[time_colname].location = EarthLocation(self.global_info['loc_x'],
                                                                     self.global_info['loc_y'],
                                                                     self.global_info['loc_z'], 
                                                                     unit='m')
                except:
                    pass
            else:
                # Still have to complete this to read FITS files not written by astropy
                pass

    def read_time(self, hdr, table):
        """
        Set the time coordinate state of a FITS Binary Table.

        Parameters
        ----------
        hdr : `~astropy.io.fits.header.Header`
            FITS Header
        table : astropy.table.Table
            The table whose time columns are to be read as Time
        """
        for key, value, comment in hdr.cards:
            if (key.upper() in TIME_KEYWORDS):

                self.set_global_time(key, value, comment)
                hdr.remove(key)

            elif (is_time_column_keyword(key.upper())):

                self.set_column_override(key, value, comment)
                hdr.remove(key)

        self.convert_time_columns(table)

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
            hdr.append(Card(keyword='TCTYP%d' %n, value=col.scale.upper()))

            # Astropy specific keyword for storing Time format
            hdr.append(Card(keyword='{0} FORMAT%d'.format(cls.astropy_namespace) %n, value=col.format.upper()))

            # Time column reference positions
            if col.location is None:
                if pos_geo is not None:
                    warnings.warn(
                        'Time Column "{0}" has no specified location, but global Time Position '
                        'is present, which will be the default for this column in '
                        'FITS specification.'.format(col.info.name),
                        AstropyUserWarning)
            else:
                # Compatibility of Time Scales and Reference Positions
                pos_ref = cls.TIME_SCALE_REF[col.scale]
                hdr.append(Card(keyword='TRPOS%d' %n, value=pos_ref))
                if pos_ref == 'TOPOCENTER':
                    if col.location.size > 1:
                        raise ValueError('Vectorized Location of Time Column "{0}" cannot be written, '
                                         'as it is not yet supported'.format(col.info.name))
                    if pos_geo is None:
                        hdr.extend([Card(keyword='OBSGEO-'+ dim.upper(), value=getattr(col.location, dim).value) for dim in ('x', 'y', 'z')])
                        pos_geo = col.location
                    elif pos_geo != col.location:
                        raise ValueError('Multiple Time Columns with different geocentric observatory '
                                         'locations ({0}, {1}, {2}) , ({3}, {4}, {5}) are encountered. '
                                         'These are not yet supported.'.format(pos_geo.x, pos_geo.y,
                                         pos_geo.z, col.location.x, col.location.y, col.location.z))
                else:
                    warnings.warn(
                        'Location cannot be written for "{0}" due to incompatability of '
                        '{1} and TOPOCENTER.'.format(col.info.name, col.scale.upper()),
                        AstropyUserWarning)

        return newtable, hdr
