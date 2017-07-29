# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings
from collections import defaultdict, OrderedDict

import numpy as np

from . import Header, Card

from .column import COLUMN_TIME_KEYWORDS, is_time_column_keyword

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
TIME_KEYWORDS = {'TIMESYS' : 'scale', 'MJDREF' : 'ref_mjd',
                 'JDREF' : 'ref_jd', 'DATEREF' : 'ref_date',
                 'TREFPOS' : 'pos', 'TREFDIR' : 'dir',
                 'TIMEUNIT' : 'unit', 'TIMEOFFS' : 'offs',
                 'OBSGEO-X' : 'loc_x', 'OBSGEO-Y' : 'loc_y',
                 'OBSGEO-Z' : 'loc_z', 'DATE' : 'date',
                 'DATE-OBS' : 'date-obs', 'DATE-AVG' : 'date-avg',
                 'DATE-BEG' : 'date-beg', 'DATE-END' : 'date-end',
                 'MJD-OBS' : 'mjd-obs', 'MJD-AVG' : 'mjd-avg',
                 'MJD-BEG' : 'mjd-beg', 'MJD-END' : 'mjd-end'}

# Set astropy time global information
GLOBAL_TIME_INFO = {'TIMESYS' : ('UTC','Default time scale'),
                    'JDREF' : (0.0,'Time columns are jd = jd1 + jd2'),
                    'TREFPOS' : ('TOPOCENTER','Time reference position')}


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
    global_info['scale'] = FITS_DEPRECATED_SCALES.get(global_info['scale'],
                                                      global_info['scale'].lower())

    # Verify global time scale
    if not global_info['scale'] in Time.SCALES:
        if global_info['scale'] in ['gps', 'local']:
            warnings.warn(
                'Global time scale (TIMESYS) has a FITS recognized time scale '
                'value {!r}. However, since it is not a valid astropy time scale, '
                'it will not be considered.'.format(global_info['scale']),
                AstropyUserWarning)
            global_info['scale'] = None
        else:
            raise AssertionError(
                'Global time scale (TIMESYS) should have a FITS recognized '
                'time scale value (got {!r}). The FITS standard states that '
                'the use of local time scales should be restricted to alternate '
                'coordinates.'.format(global_info['scale']))


def _verify_column_info(global_info, time_col):
    """
    Given the column-specific time reference frame information, verify that
    each column-specific time coordinate attribute has a valid value.
    Return True if the coordinate column is time, or else return False.

    Parameters
    ----------
    global_info : dict
        Global time reference frame information.
    time_col : dict
        Column time reference frame information.
    """

    scale = time_col.get('scale', None)
    unit = time_col.get('unit', None)

    if scale is not None:

        if scale.lower() in Time.SCALES:
            time_col['scale'] = scale.lower()
            return True

        if scale in FITS_DEPRECATED_SCALES.keys():
            time_col['scale'] = FITS_DEPRECATED_SCALES[scale]
            return True

        if scale == 'TIME':
            time_col['scale'] = global_info['scale']
            return True if time_col['scale'] else False

        if scale in ['GPS', 'LOCAL']:
            warnings.warn(
                'Table column {} has a FITS recognized time scale value {}. '
                'However, since it is not a valid astropy time scale, '
                'it will not be converted to astropy Time.'.format(idx, scale),
                AstropyUserWarning)
            return False

        # Non-linear coordinate types have "4-3" form and are not time coordinates
        if TCTYP_RE_TYPE.match(scale[:5]) and TCTYP_RE_ALGO.match(scale[5:]):
            return False

        # Coordinate type is either an unrecognized local time scale
        # or a linear coordinate type
        return False

    if global_info['scale'] is None:
        return False

    if (unit is not None and unit in FITS_TIME_UNIT) or 'pos' in time_col:
        return True

    # Other cases are yet to be checked and will be done when reading starts
    # This includes checking TTYPEn and TUNITn
    return False


def _convert_global_time(table, global_info):
    """
    Convert the table metadata for time informational keywords
    to astropy Time.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table whose time columns are to be converted.
    global_info : dict
        Global time reference frame information.
    """

    # Read in Global Informational keywords as Time
    for key, value in global_info.items():
        # FITS uses a subset of ISO-8601 for DATE-xxx
        if key.startswith('date'):
            if key not in table.meta:
                precision = (lambda x: len(x.split('.')[1])
                             if '.' in x else 0)
                scale = 'utc' if key == 'date' else global_info['scale']
                if scale is not None:
                    try:
                        table.meta[key.upper()] = Time(value, scale=scale,
                                                   precision=precision(value))
                        continue
                    except ValueError:
                        pass
                table.meta[key.upper()] = value

        # MJD-xxx in MJD according to TIMESYS
        elif key.startswith('mjd-'):
            if key not in table.meta:
                if global_info['scale'] is not None:
                    try:
                        table.meta[key.upper()] = Time(value, scale=global_info['scale'],
                                                       format='mjd')
                        continue
                    except ValueError:
                        pass
                table.meta[key.upper()] = value


def _convert_time_columns(table, global_info, time_columns):
    """
    Convert time columns to astropy Time columns.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The table whose time columns are to be converted.
    global_info : dict
        Global time reference frame information.
    time_columns : dict
        Column-specific time reference frame override information.
    """

    # The code might fail while attempting to read FITS files not written by astropy.

    # Read in time coordinate columns as Time
    for idx, time_col in time_columns.items():
        time_colname = table.colnames[idx - 1]
        if table[time_colname].shape[-1] == 2 and table[time_colname].ndim > 1:
            table[time_colname] = Time(table[time_colname][...,0], table[time_colname][...,1],
                                       format='jd', scale=time_col['scale'])
        else:
            warnings.warn(
                'Time column {} is not in the astropy required (jd1, jd2) format. '
                'Hence, it will not be read as astropy Time'.format(col.info.name),
                AstropyUserWarning)
            continue
        # Add location if attribute ``pos`` is set
        try:
            if time_col['pos'] == 'TOPOCENTER':
                table[time_colname].location = EarthLocation(global_info['loc_x'],
                                                             global_info['loc_y'],
                                                             global_info['loc_z'],
                                                             unit='m')
        except KeyError:
            pass


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
    global_info = {'scale' : 'UTC',
                   'mjdref' : None,
                   'unit' : 's'}

    # Set default dictionary for time columns
    time_columns = defaultdict(OrderedDict)

    # Make a "copy" (not just a view) of the input header, since it
    # may get modified.  the data is still a "view" (for now)
    hcopy = hdr.copy(strip=True)

    for key, value, comment in hdr.cards:
        if (key in TIME_KEYWORDS):

            global_info[TIME_KEYWORDS[key]] = value
            hcopy.remove(key)

        elif (is_time_column_keyword(key)):

            base, idx = re.match(r'([A-Z]+)([0-9]+)', key).groups()
            time_columns[int(idx)][COLUMN_TIME_KEYWORDS[base]] = value
            hcopy.remove(key)

    if len(hcopy) != len(hdr):
        _verify_global_info(global_info)
        _convert_global_time(table, global_info)
        if time_columns:
            # Select only time_columns that are valid and issue warnings otherwise
            time_columns = {idx: time_col for idx, time_col in time_columns.items()
                            if _verify_column_info(global_info, time_col)}
            _convert_time_columns(table, global_info, time_columns)

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
                hdr.extend([Card(keyword='OBSGEO-'+ dim.upper(), value=getattr(col.location, dim).value)
                            for dim in ('x', 'y', 'z')])
                pos_geo = col.location
            elif pos_geo != col.location:
                raise ValueError('Multiple Time Columns with different geocentric observatory '
                                 'locations ({}, {}, {}) , ({}, {}, {}) are encountered. '
                                 'These are not supported by the FITS standard.'.format(pos_geo.x,
                                 pos_geo.y, pos_geo.z, col.location.x, col.location.y,
                                 col.location.z))

    return newtable, hdr
