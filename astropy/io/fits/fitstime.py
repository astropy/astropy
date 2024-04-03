# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings
from collections import OrderedDict, defaultdict

import numpy as np

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.table import Column, MaskedColumn
from astropy.table.column import col_copy
from astropy.time import Time, TimeDelta
from astropy.time.core import BARYCENTRIC_SCALES
from astropy.time.formats import FITS_DEPRECATED_SCALES
from astropy.utils.exceptions import AstropyUserWarning

from . import Card, Header

# The following is based on the FITS WCS Paper IV, "Representations of time
# coordinates in FITS".
# https://ui.adsabs.harvard.edu/abs/2015A%26A...574A..36R


# FITS WCS standard specified "4-3" form for non-linear coordinate types
TCTYP_RE_TYPE = re.compile(r"(?P<type>[A-Z]+)[-]+")
TCTYP_RE_ALGO = re.compile(r"(?P<algo>[A-Z]+)\s*")


# FITS Time standard specified time units
FITS_TIME_UNIT = ["s", "d", "a", "cy", "min", "h", "yr", "ta", "Ba"]


# Global time reference coordinate keywords
OBSGEO_XYZ = ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")
OBSGEO_LBH = ("OBSGEO-L", "OBSGEO-B", "OBSGEO-H")
TIME_KEYWORDS = (
    (
        "DATE",
        "DATE-AVG",
        "DATE-BEG",
        "DATE-END",
        "DATE-OBS",
        "DATEREF",
        "JDREF",
        "MJD-AVG",
        "MJD-BEG",
        "MJD-END",
        "MJD-OBS",
        "MJDREF",
        "TIMEOFFS",
        "TIMESYS",
        "TIMEUNIT",
        "TREFDIR",
        "TREFPOS",
    )
    + OBSGEO_LBH
    + OBSGEO_XYZ
)


# Column-specific time override keywords
COLUMN_TIME_KEYWORDS = ("TCTYP", "TCUNI", "TRPOS")


# Column-specific keywords regex
COLUMN_TIME_KEYWORD_REGEXP = f"({'|'.join(COLUMN_TIME_KEYWORDS)})[0-9]+"


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
GLOBAL_TIME_INFO = {
    "TIMESYS": ("UTC", "Default time scale"),
    "JDREF": (0.0, "Time columns are jd = jd1 + jd2"),
    "TREFPOS": ("TOPOCENTER", "Time reference position"),
}


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
    global_info["scale"] = FITS_DEPRECATED_SCALES.get(
        global_info["TIMESYS"], global_info["TIMESYS"].lower()
    )

    # Verify global time scale
    if global_info["scale"] not in Time.SCALES:
        # 'GPS' and 'LOCAL' are FITS recognized time scale values
        # but are not supported by astropy.

        if global_info["scale"] == "gps":
            warnings.warn(
                "Global time scale (TIMESYS) has a FITS recognized time scale "
                'value "GPS". In Astropy, "GPS" is a time from epoch format '
                "which runs synchronously with TAI; GPS is approximately 19 s "
                "ahead of TAI. Hence, this format will be used.",
                AstropyUserWarning,
            )
            # Assume that the values are in GPS format
            global_info["scale"] = "tai"
            global_info["format"] = "gps"

        if global_info["scale"] == "local":
            warnings.warn(
                "Global time scale (TIMESYS) has a FITS recognized time scale "
                'value "LOCAL". However, the standard states that "LOCAL" should be '
                "tied to one of the existing scales because it is intrinsically "
                "unreliable and/or ill-defined. Astropy will thus use the default "
                'global time scale "UTC" instead of "LOCAL".',
                AstropyUserWarning,
            )
            # Default scale 'UTC'
            global_info["scale"] = "utc"
            global_info["format"] = None

        else:
            raise AssertionError(
                "Global time scale (TIMESYS) should have a FITS recognized "
                "time scale value (got {!r}). The FITS standard states that "
                "the use of local time scales should be restricted to alternate "
                "coordinates.".format(global_info["TIMESYS"])
            )
    else:
        # Scale is already set
        global_info["format"] = None

    # Check if geocentric global location is specified
    obs_geo = [global_info[attr] for attr in OBSGEO_XYZ if attr in global_info]

    # Location full specification is (X, Y, Z)
    if len(obs_geo) == 3:
        global_info["location"] = EarthLocation.from_geocentric(*obs_geo, unit=u.m)
    else:
        # Check if geodetic global location is specified (since geocentric failed)

        # First warn the user if geocentric location is partially specified
        if obs_geo:
            warnings.warn(
                f"The geocentric observatory location {obs_geo} is not completely "
                "specified (X, Y, Z) and will be ignored.",
                AstropyUserWarning,
            )

        # Check geodetic location
        obs_geo = [global_info[attr] for attr in OBSGEO_LBH if attr in global_info]

        if len(obs_geo) == 3:
            global_info["location"] = EarthLocation.from_geodetic(*obs_geo)
        else:
            # Since both geocentric and geodetic locations are not specified,
            # location will be None.

            # Warn the user if geodetic location is partially specified
            if obs_geo:
                warnings.warn(
                    f"The geodetic observatory location {obs_geo} is not completely "
                    "specified (lon, lat, alt) and will be ignored.",
                    AstropyUserWarning,
                )
            global_info["location"] = None

    # Get global time reference
    # Keywords are listed in order of precedence, as stated by the standard
    for key, format_ in (("MJDREF", "mjd"), ("JDREF", "jd"), ("DATEREF", "fits")):
        if key in global_info:
            global_info["ref_time"] = {"val": global_info[key], "format": format_}
            break
    else:
        # If none of the three keywords is present, MJDREF = 0.0 must be assumed
        global_info["ref_time"] = {"val": 0, "format": "mjd"}


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
    scale = column_info.get("TCTYP", None)
    unit = column_info.get("TCUNI", None)
    location = column_info.get("TRPOS", None)

    if scale is not None:
        # Non-linear coordinate types have "4-3" form and are not time coordinates
        if TCTYP_RE_TYPE.match(scale[:5]) and TCTYP_RE_ALGO.match(scale[5:]):
            return False

        elif scale.lower() in Time.SCALES:
            column_info["scale"] = scale.lower()
            column_info["format"] = None

        elif scale in FITS_DEPRECATED_SCALES.keys():
            column_info["scale"] = FITS_DEPRECATED_SCALES[scale]
            column_info["format"] = None

        # TCTYPn (scale) = 'TIME' indicates that the column scale is
        # controlled by the global scale.
        elif scale == "TIME":
            column_info["scale"] = global_info["scale"]
            column_info["format"] = global_info["format"]

        elif scale == "GPS":
            warnings.warn(
                (
                    f'Table column "{column_info}" has a FITS recognized time scale '
                    'value "GPS". In Astropy, "GPS" is a time from epoch format which '
                    "runs synchronously with TAI; GPS runs ahead of TAI approximately "
                    "by 19 s. Hence, this format will be used."
                ),
                AstropyUserWarning,
            )
            column_info["scale"] = "tai"
            column_info["format"] = "gps"

        elif scale == "LOCAL":
            warnings.warn(
                (
                    f'Table column "{column_info}" has a FITS recognized time scale '
                    'value "LOCAL". However, the standard states that "LOCAL" should '
                    "be tied to one of the existing scales because it is intrinsically "
                    "unreliable and/or ill-defined. Astropy will thus use the global "
                    "time scale (TIMESYS) as the default."
                ),
                AstropyUserWarning,
            )
            column_info["scale"] = global_info["scale"]
            column_info["format"] = global_info["format"]

        else:
            # Coordinate type is either an unrecognized local time scale
            # or a linear coordinate type
            return False

    # If TCUNIn is a time unit or TRPOSn is specified, the column is a time
    # coordinate. This has to be tested since TCTYP (scale) is not specified.
    elif (unit is not None and unit in FITS_TIME_UNIT) or location is not None:
        column_info["scale"] = global_info["scale"]
        column_info["format"] = global_info["format"]

    # None of the conditions for time coordinate columns is satisfied
    else:
        return False

    # Check if column-specific reference position TRPOSn is specified
    if location is not None:
        # Observatory position (location) needs to be specified only
        # for 'TOPOCENTER'.
        if location == "TOPOCENTER":
            column_info["location"] = global_info["location"]
            if column_info["location"] is None:
                warnings.warn(
                    'Time column reference position "TRPOSn" value is "TOPOCENTER". '
                    "However, the observatory position is not properly specified. "
                    "The FITS standard does not support this and hence reference "
                    "position will be ignored.",
                    AstropyUserWarning,
                )
        else:
            column_info["location"] = None

    # Warn user about ignoring global reference position when TRPOSn is
    # not specified
    elif global_info["TREFPOS"] == "TOPOCENTER":
        if global_info["location"] is not None:
            warnings.warn(
                'Time column reference position "TRPOSn" is not specified. The '
                'default value for it is "TOPOCENTER", and the observatory position '
                "has been specified. However, for supporting column-specific location, "
                "reference position will be ignored for this column.",
                AstropyUserWarning,
            )
        column_info["location"] = None
    else:
        column_info["location"] = None

    # Get reference time
    column_info["ref_time"] = global_info["ref_time"]

    return True


def _get_info_if_time_column(col, global_info):
    """
    Check if a column without corresponding time column keywords in the
    FITS header represents time or not. If yes, return the time column
    information needed for its conversion to Time.
    This is only applicable to the special-case where a column has the
    name 'TIME' and a time unit.
    """
    # Column with TTYPEn = 'TIME' and lacking any TC*n or time
    # specific keywords will be controlled by the global keywords.
    if col.info.name.upper() == "TIME" and col.info.unit in FITS_TIME_UNIT:
        column_info = {
            "scale": global_info["scale"],
            "format": global_info["format"],
            "ref_time": global_info["ref_time"],
            "location": None,
        }

        if global_info["TREFPOS"] == "TOPOCENTER":
            column_info["location"] = global_info["location"]
            if column_info["location"] is None:
                warnings.warn(
                    f'Time column "{col.info.name}" reference position will be ignored '
                    "due to unspecified observatory position.",
                    AstropyUserWarning,
                )

        return column_info

    return None


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
    for key in global_info:
        # FITS uses a subset of ISO-8601 for DATE-xxx
        if key not in table.meta:
            try:
                table.meta[key] = _convert_time_key(global_info, key)
            except ValueError:
                pass


def _convert_time_key(global_info, key):
    """
    Convert a time metadata key to a Time object.

    Parameters
    ----------
    global_info : dict
        Global time reference frame information.
    key : str
        Time key.

    Returns
    -------
    astropy.time.Time

    Raises
    ------
    ValueError
        If key is not a valid global time keyword.
    """
    value = global_info[key]
    if key.startswith("DATE"):
        scale = "utc" if key == "DATE" else global_info["scale"]
        precision = len(value.split(".")[-1]) if "." in value else 0
        return Time(value, format="fits", scale=scale, precision=precision)
    # MJD-xxx in MJD according to TIMESYS
    elif key.startswith("MJD-"):
        return Time(value, format="mjd", scale=global_info["scale"])
    else:
        raise ValueError("Key is not a valid global time keyword")


def _convert_time_column(col, column_info):
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
    try:
        # ISO-8601 is the only string representation of time in FITS
        if col.info.dtype.kind in ["S", "U"]:
            # [+/-C]CCYY-MM-DD[Thh:mm:ss[.s...]] where the number of characters
            # from index 20 to the end of string represents the precision
            precision = max(int(col.info.dtype.str[2:]) - 20, 0)
            return Time(
                col,
                format="fits",
                scale=column_info["scale"],
                precision=precision,
                location=column_info["location"],
            )

        if column_info["format"] == "gps":
            return Time(col, format="gps", location=column_info["location"])

        # If reference value is 0 for JD or MJD, the column values can be
        # directly converted to Time, as they are absolute (relative
        # to a globally accepted zero point).
        if column_info["ref_time"]["val"] == 0 and column_info["ref_time"][
            "format"
        ] in ["jd", "mjd"]:
            # (jd1, jd2) where jd = jd1 + jd2
            if col.shape[-1] == 2 and col.ndim > 1:
                return Time(
                    col[..., 0],
                    col[..., 1],
                    scale=column_info["scale"],
                    format=column_info["ref_time"]["format"],
                    location=column_info["location"],
                )
            else:
                return Time(
                    col,
                    scale=column_info["scale"],
                    format=column_info["ref_time"]["format"],
                    location=column_info["location"],
                )

        # Reference time
        ref_time = Time(
            column_info["ref_time"]["val"],
            scale=column_info["scale"],
            format=column_info["ref_time"]["format"],
            location=column_info["location"],
        )

        # Elapsed time since reference time
        if col.shape[-1] == 2 and col.ndim > 1:
            delta_time = TimeDelta(col[..., 0], col[..., 1])
        else:
            delta_time = TimeDelta(col)

        return ref_time + delta_time
    except Exception as err:
        warnings.warn(
            f'The exception "{err}" was encountered while trying to convert the time '
            f'column "{col.info.name}" to Astropy Time.',
            AstropyUserWarning,
        )
        return col


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
    global_info = {"TIMESYS": "UTC", "TREFPOS": "TOPOCENTER"}

    # Set default dictionary for time columns
    time_columns = defaultdict(OrderedDict)

    # Make a "copy" (not just a view) of the input header, since it
    # may get modified.  the data is still a "view" (for now)
    hcopy = hdr.copy(strip=True)

    # Scan the header for global and column-specific time keywords
    for key, value, comment in hdr.cards:
        if key in TIME_KEYWORDS:
            global_info[key] = value
            hcopy.remove(key)

        elif is_time_column_keyword(key):
            base, idx = re.match(r"([A-Z]+)([0-9]+)", key).groups()
            time_columns[int(idx)][base] = value
            hcopy.remove(key)

        elif value in OBSGEO_XYZ and re.match("TTYPE[0-9]+", key):
            global_info[value] = table[value]

    # Verify and get the global time reference frame information
    _verify_global_info(global_info)
    _convert_global_time(table, global_info)

    # Columns with column-specific time (coordinate) keywords
    if time_columns:
        for idx, column_info in time_columns.items():
            # Check if the column is time coordinate (not spatial)
            if _verify_column_info(column_info, global_info):
                colname = table.colnames[idx - 1]
                # Convert to Time
                table[colname] = _convert_time_column(table[colname], column_info)

    # Check for special-cases of time coordinate columns
    for idx, colname in enumerate(table.colnames):
        if (idx + 1) not in time_columns:
            column_info = _get_info_if_time_column(table[colname], global_info)
            if column_info:
                table[colname] = _convert_time_column(table[colname], column_info)

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
        Header containing global time reference frame FITS keywords
    """
    # Make a light copy of table (to the extent possible) and clear any indices along
    # the way. Indices are not serialized and cause problems later, but they are not
    # needed here so just drop.  For Column subclasses take advantage of copy() method,
    # but for others it is required to actually copy the data if there are attached
    # indices.  See #8077 and #9009 for further discussion.
    new_cols = []
    for col in table.itercols():
        if isinstance(col, Column):
            new_col = col.copy(copy_data=False)  # Also drops any indices
        else:
            new_col = col_copy(col, copy_indices=False) if col.info.indices else col
        new_cols.append(new_col)
    newtable = table.__class__(new_cols, copy=False)
    newtable.meta = table.meta

    # Global time coordinate frame keywords
    hdr = Header(
        [
            Card(keyword=key, value=val[0], comment=val[1])
            for key, val in GLOBAL_TIME_INFO.items()
        ]
    )

    # Store coordinate column-specific metadata
    newtable.meta["__coordinate_columns__"] = defaultdict(OrderedDict)
    coord_meta = newtable.meta["__coordinate_columns__"]

    time_cols = table.columns.isinstance(Time)

    # Geocentric location
    location = None

    for col in time_cols:
        # By default, Time objects are written in full precision, i.e. we store both
        # jd1 and jd2 (serialize_method['fits'] = 'jd1_jd2'). Formatted values for
        # Time can be stored if the user explicitly chooses to do so.
        col_cls = MaskedColumn if col.masked else Column
        if col.info.serialize_method["fits"] == "formatted_value":
            newtable.replace_column(col.info.name, col_cls(col.value))
            continue

        # The following is necessary to deal with multi-dimensional ``Time`` objects
        # (i.e. where Time.shape is non-trivial).
        # Note: easier would be np.stack([col.jd1, col.jd2], axis=-1), but that
        # fails for np.ma.MaskedArray, as it returns just the data, ignoring the mask.
        jd12 = np.empty_like(col.jd1, shape=col.jd1.shape + (2,))
        jd12[..., 0] = col.jd1
        jd12[..., 1] = col.jd2
        newtable.replace_column(col.info.name, col_cls(jd12, unit="d"))

        # Time column-specific override keywords
        coord_meta[col.info.name]["coord_type"] = col.scale.upper()
        coord_meta[col.info.name]["coord_unit"] = "d"

        # Time column reference position
        if col.location is None:
            coord_meta[col.info.name]["time_ref_pos"] = None
            if location is not None:
                warnings.warn(
                    (
                        f'Time Column "{col.info.name}" has no specified location, '
                        "but global Time Position is present, which will be the "
                        "default for this column in FITS specification."
                    ),
                    AstropyUserWarning,
                )
        else:
            coord_meta[col.info.name]["time_ref_pos"] = "TOPOCENTER"
            # Compatibility of Time Scales and Reference Positions
            if col.scale in BARYCENTRIC_SCALES:
                warnings.warn(
                    (
                        f'Earth Location "TOPOCENTER" for Time Column "{col.info.name}" '
                        f'is incompatible with scale "{col.scale.upper()}".'
                    ),
                    AstropyUserWarning,
                )

            if location is None:
                # Set global geocentric location
                location = col.location
                if location.size > 1:
                    for dim in ("x", "y", "z"):
                        newtable.add_column(
                            Column(getattr(location, dim).to_value(u.m)),
                            name=f"OBSGEO-{dim.upper()}",
                        )
                else:
                    hdr.extend(
                        [
                            Card(
                                keyword=f"OBSGEO-{dim.upper()}",
                                value=getattr(location, dim).to_value(u.m),
                            )
                            for dim in ("x", "y", "z")
                        ]
                    )
            elif np.any(location != col.location):
                raise ValueError(
                    "Multiple Time Columns with different geocentric "
                    f"observatory locations ({location}, {col.location}) encountered."
                    "This is not supported by the FITS standard."
                )

    return newtable, hdr
