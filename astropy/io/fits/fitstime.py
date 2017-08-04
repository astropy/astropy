# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings
from collections import defaultdict, OrderedDict

import numpy as np

from . import Header, Card

from ...coordinates import EarthLocation
from ...table import Column
from ...time import Time
from ...time.core import BARYCENTRIC_SCALES
from ...time.formats import FITS_DEPRECATED_SCALES
from ...utils.exceptions import AstropyUserWarning

# The following is based on the FITS WCS Paper IV, "Representations of time
# coordinates in FITS".
# http://adsabs.harvard.edu/abs/2015A%26A...574A..36R


# FITS WCS standard specified "4-3" form for non-linear coordinate types
TCTYP_RE_TYPE = re.compile(r'(?P<type>[A-Z]+)[-]+')
TCTYP_RE_ALGO = re.compile(r'(?P<algo>[A-Z]+)\s*')


# FITS Time standard specified time units
FITS_TIME_UNIT = ['s', 'd', 'a', 'cy', 'min', 'h', 'yr', 'ta', 'Ba']


# Global time reference coordinate keywords
TIME_KEYWORDS = ('TIMESYS', 'MJDREF', 'JDREF', 'DATEREF',
                 'TREFPOS', 'TREFDIR', 'TIMEUNIT', 'TIMEOFFS',
                 'OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z', 'DATE',
                 'DATE-OBS', 'DATE-AVG', 'DATE-BEG', 'DATE-END',
                 'MJD-OBS', 'MJD-AVG', 'MJD-BEG', 'MJD-END')


# Column-specific time override keywords
COLUMN_TIME_KEYWORDS = ('TCTYP', 'TCUNI', 'TRPOS')


# Column-specific keywords regex
COLUMN_TIME_KEYWORD_REGEXP = '({0})[0-9]+'.format(
    '|'.join(COLUMN_TIME_KEYWORDS))


def is_time_column_keyword(keyword):
    """
    Check if the FITS header keyword is a time column-specific keyword.

    Parameters
    ----------
    keyword : str
        FITS keyword.
    """
    return re.match(COLUMN_TIME_KEYWORD_REGEXP, keyword) is not None


# Set astropy time global information
GLOBAL_TIME_INFO = {'TIMESYS': ('UTC', 'Default time scale'),
                    'JDREF': (0.0, 'Time columns are jd = jd1 + jd2'),
                    'TREFPOS': ('TOPOCENTER', 'Time reference position')}


def _verify_global_info(global_info):
    """
    Given the global time reference frame information, verify that
    each global time coordinate attribute will be given a valid value.

    Parameters
    ----------
    global_info : dict
        Global time reference frame information.
    """

    # Translate FITS deprecated scale into astropy scale, or else just convert
    # to lower case for further checks.
    global_info['TIMESYS'] = FITS_DEPRECATED_SCALES.get(global_info['TIMESYS'],
                                                      global_info['TIMESYS'].lower())

    # Verify global time scale
    if not global_info['TIMESYS'] in Time.SCALES:

        if global_info['TIMESYS'] == 'gps':
            warnings.warn(
                'Global time scale (TIMESYS) has a FITS recognized time scale '
                'value "GPS". In Astropy, "GPS" is a time from epoch format '
                'which runs synchronously with TAI; GPS is approximately '
                'TAI − 19 s. Hence, this format will be used.', AstropyUserWarning)
            global_info['FORMAT'] = 'gps'
            global_info['TIMESYS'] = 'tai'

        if global_info['TIMESYS'] == 'local':
            warnings.warn(
                'Global time scale (TIMESYS) has a FITS recognized time scale '
                'value "LOCAL". However, the standard states that "LOCAL" should be '
                'tied to one of the existing scales because it is intrinsically '
                'unreliable and/or ill-defined. Astropy will thus use the default '
                'global time scale "UTC" instead of "LOCAL".', AstropyUserWarning)
            global_info['TIMESYS'] = 'utc'

        else:
            raise AssertionError(
                'Global time scale (TIMESYS) should have a FITS recognized '
                'time scale value (got {!r}). The FITS standard states that '
                'the use of local time scales should be restricted to alternate '
                'coordinates.'.format(global_info['TIMESYS']))


def _verify_column_info(column_info, global_info):
    """
    Given the column-specific time reference frame information, verify that
    each column-specific time coordinate attribute has a valid value.
    Return True if the coordinate column is time, or else return False.

    Parameters
    ----------
    global_info : dict
        Global time reference frame information.
    column_info : dict
        Column-specific time reference frame override information.
    """

    scale = column_info.get('TCTYP', None)
    unit = column_info.get('TCUNI', None)
    pos = column_info.get('TRPOS', None)

    if pos is not None:

        if pos == 'TOPOCENTER':
            column_info['TRPOS'] = EarthLocation(global_info['OBSGEO-X'],
                                                 global_info['OBSGEO-Y'],
                                                 global_info['OBSGEO-Z'],
                                                 unit='m')
        else:
            column_info['TRPOS'] = None

    if scale is not None:

        if scale.lower() in Time.SCALES:
            column_info['TCTYP'] = scale.lower()
            return True

        if scale in FITS_DEPRECATED_SCALES.keys():
            column_info['TCTYP'] = FITS_DEPRECATED_SCALES[scale]
            return True

        if scale == 'TIME':
            column_info['TCTYP'] = global_info['TIMESYS']
            column_info['FORMAT'] = global_info.get('FORMAT', None)
            return True

        if scale == 'GPS':
            warnings.warn(
                'Table column {} has a FITS recognized time scale value "GPS". '
                'In Astropy, "GPS" is a time from epoch format which runs '
                'synchronously with TAI, as the scale GPS ≈ TAI − 19 s. '
                'Hence, this format will be used.'.format(column_info),
                AstropyUserWarning)
            column_info['FORMAT'] = 'gps'
            column_info['TCTYP'] = 'tai'
            return True

        if scale == 'LOCAL':
            warnings.warn(
                'Table column {} has a FITS recognized time scale value "LOCAL". '
                'However, the standard states that "LOCAL" should be tied to one '
                'of the existing scales because it is intrinsically unreliable '
                'and/or ill-defined. Astropy will thus use the global time scale '
                '(TIMESYS) as the default.'. format(column_info),
                AstropyUserWarning)
            column_info['TCTYP'] = global_info['TIMESYS']
            column_info['FORMAT'] = global_info.get('FORMAT', None)

        # Non-linear coordinate types have "4-3" form and are not time coordinates
        if TCTYP_RE_TYPE.match(scale[:5]) and TCTYP_RE_ALGO.match(scale[5:]):
            return False

        # Coordinate type is either an unrecognized local time scale
        # or a linear coordinate type
        return False

    if global_info['TIMESYS'] is None:
        return False

    if (unit is not None and unit in FITS_TIME_UNIT) or pos is not None:
        column_info['TCTYP'] = global_info['TIMESYS']
        return True

    return False


def _convert_global_time(table, global_info):
    """
    Convert the table metadata for time informational keywords
    to astropy Time.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table whose time metadata is to be converted.
    global_info : dict
        Global time reference frame information.
    """

    # Read in Global Informational keywords as Time
    for key, value in global_info.items():
        # FITS uses a subset of ISO-8601 for DATE-xxx
        if key.startswith('DATE'):
            if key not in table.meta:
                scale = 'utc' if key == 'DATE' else global_info['TIMESYS']
                try:
                    precision = len(value.split('.')[-1]) if '.' in value else 0
                    value = Time(value, format='fits', scale=scale,
                                 precision=precision)
                except ValueError:
                    pass
                table.meta[key] = value

        # MJD-xxx in MJD according to TIMESYS
        elif key.startswith('MJD-'):
            if key not in table.meta:
                try:
                    value = Time(value, format='mjd',
                                 scale=global_info['TIMESYS'])
                except ValueError:
                    pass
                table.meta[key] = value


def _convert_time_columns(col, column_info):
    """
    Convert time columns to astropy Time columns.

    Parameters
    ----------
    col : `~astropy.table.Column`
        The time coordinate column to be converted to Time.
    column_info : dict
        Column-specific time reference frame override information.
    """

    # The code might fail while attempting to read FITS files not written by astropy.

    # Read in time coordinate column as Time
    if col.shape[-1] == 2 and col.ndim > 1:
        col = Time(col[..., 0], col[..., 1], format='jd',
                   scale=column_info['TCTYP'], location=column_info['TRPOS'])
    else:
        warnings.warn(
            'Time column {} is not in the astropy required (jd1, jd2) format. '
            'Hence, it will not be read as astropy Time'.format(col.info.name),
            AstropyUserWarning)
        return


def fits_to_time(hdr, table):
    """
    Read FITS binary table time columns as `~astropy.time.Time`.

    This method reads the metadata associated with time coordinates, as
    stored in a FITS binary table header, converts time columns into
    `~astropy.time.Time` columns and reads global reference times as
    `~astropy.time.Time` instances.

    Parameters
    ----------
    hdr : `~astropy.io.fits.header.Header`
        FITS Header
    table : `~astropy.table.Table`
        The table whose time columns are to be read as Time

    Returns
    -------
    hdr : `~astropy.io.fits.header.Header`
        Modified FITS Header (time metadata removed)
    """

    # Set defaults for global time scale, reference, etc.
    global_info = {'TIMESYS': 'UTC',
                   'MJDREF': None,
                   'TIMEUNIT': 's'}

    # Set default dictionary for time columns
    time_columns = defaultdict(OrderedDict)

    # Make a "copy" (not just a view) of the input header, since it
    # may get modified.  the data is still a "view" (for now)
    hcopy = hdr.copy(strip=True)

    for key, value, comment in hdr.cards:
        if key in TIME_KEYWORDS:

            global_info[key] = value
            hcopy.remove(key)

        elif is_time_column_keyword(key):

            base, idx = re.match(r'([A-Z]+)([0-9]+)', key).groups()
            time_columns[int(idx)][base] = value
            hcopy.remove(key)

    if len(hcopy) != len(hdr):
        _verify_global_info(global_info)
        _convert_global_time(table, global_info)
        if time_columns:
            for idx, column_info in time_columns.items():
                if _verify_column_info(column_info, global_info):
                    colname = table.colnames[idx - 1]
                    table[colname] = _convert_time_column(table[colname],
                                                          column_info)

    return hcopy


def time_to_fits(table):
    """
    Replace Time columns in a Table with non-mixin columns containing
    each element as a vector of two doubles (jd1, jd2) and return a FITS
    header with appropriate time coordinate keywords.
    jd = jd1 + jd2 represents time in the Julian Date format with
    high-precision.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table whose Time columns are to be replaced.

    Returns
    -------
    table : `~astropy.table.Table`
        The table with replaced Time columns
    hdr : `~astropy.io.fits.header.Header`
        Header containing Cards associated with the FITS time coordinate
    """

    # Shallow copy of the input table
    newtable = table.copy(copy_data=False)

    # Global time coordinate frame keywords
    hdr = Header([Card(keyword=key, value=val[0], comment=val[1])
                  for key, val in GLOBAL_TIME_INFO.items()])

    time_cols = table.columns.isinstance(Time)

    # Geocentric Position
    pos_geo = None

    for col in time_cols:
        # The following is necessary to deal with multi-dimensional ``Time`` objects
        # (i.e. where Time.shape is non-trivial).
        jd12 = np.array([col.jd1, col.jd2])
        # Roll the 0th (innermost) axis backwards, until it lies in the last position
        # (jd12.ndim)
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
            if col.scale in BARYCENTRIC_SCALES:
                warnings.warn(
                    'Earth Location "TOPOCENTER" for Time Column "{}" is incompatabile '
                    'with scale "{}".'.format(col.info.name, col.scale.upper()),
                    AstropyUserWarning)
            if col.location.size > 1:
                raise ValueError('Vectorized Location of Time Column "{}" cannot be written, '
                                 'as it is not supported.'.format(col.info.name))
            if pos_geo is None:
                hdr.extend([Card(keyword='OBSGEO-'+ dim.upper(),
                            value=getattr(col.location, dim).to_value(unit='m'))
                            for dim in ('x', 'y', 'z')])
                pos_geo = col.location
            elif pos_geo != col.location:
                raise ValueError('Multiple Time Columns with different geocentric observatory '
                                 'locations ({}, {}, {}) , ({}, {}, {}) are encountered. '
                                 'These are not supported by the FITS standard.'.format(pos_geo.x,
                                 pos_geo.y, pos_geo.z, col.location.x, col.location.y,
                                 col.location.z))

    return newtable, hdr
