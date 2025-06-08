# Licensed under a 3-clause BSD style license - see LICENSE.rst
import datetime
import fnmatch
import functools
import re
import time
import warnings
from decimal import Decimal

import erfa
import numpy as np

import astropy.units as u
from astropy.utils.compat.optional_deps import HAS_MATPLOTLIB
from astropy.utils.decorators import classproperty, lazyproperty
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning
from astropy.utils.masked import Masked

from . import _parse_times, conf, utils
from .utils import day_frac, quantity_day_frac, two_product, two_sum

__all__ = [
    "TIME_DELTA_FORMATS",
    "TIME_FORMATS",
    "AstropyDatetimeLeapSecondWarning",
    "TimeBesselianEpoch",
    "TimeBesselianEpochString",
    "TimeCxcSec",
    "TimeDatetime",
    "TimeDatetime64",
    "TimeDecimalYear",
    "TimeDeltaDatetime",
    "TimeDeltaFormat",
    "TimeDeltaJD",
    "TimeDeltaNumeric",
    "TimeDeltaQuantityString",
    "TimeDeltaSec",
    "TimeEpochDate",
    "TimeEpochDateString",
    "TimeFITS",
    "TimeFormat",
    "TimeFromEpoch",
    "TimeGPS",
    "TimeISO",
    "TimeISOT",
    "TimeJD",
    "TimeJulianEpoch",
    "TimeJulianEpochString",
    "TimeMJD",
    "TimeNumeric",
    "TimePlotDate",
    "TimeString",
    "TimeUnique",
    "TimeUnix",
    "TimeUnixTai",
    "TimeYMDHMS",
    "TimeYearDayTime",
    "TimezoneInfo",
]

__doctest_skip__ = ["TimePlotDate"]

# These both get filled in at end after TimeFormat subclasses defined.
# It is important that these get populated by insertion order.
# This ensures, e.g., that 'isot' gets tried before 'fits'.
TIME_FORMATS = {}
TIME_DELTA_FORMATS = {}

# Translations between deprecated FITS timescales defined by
# Rots et al. 2015, A&A 574:A36, and timescales used here.
FITS_DEPRECATED_SCALES = {
    "TDT": "tt",
    "ET": "tt",
    "GMT": "utc",
    "UT": "utc",
    "IAT": "tai",
}


class AstropyDatetimeLeapSecondWarning(AstropyUserWarning):
    """Warning for leap second when converting to datetime.datetime object."""


def _regexify_subfmts(subfmts):
    """
    Iterate through each of the sub-formats and try substituting simple
    regular expressions for the strptime codes for year, month, day-of-month,
    hour, minute, second.  If no % characters remain then turn the final string
    into a compiled regex.  This assumes time formats do not have a % in them.

    This is done both to speed up parsing of strings and to allow mixed formats
    where strptime does not quite work well enough.
    """
    new_subfmts = []
    for subfmt_tuple in subfmts:
        subfmt_in = subfmt_tuple[1]
        if isinstance(subfmt_in, str):
            for strptime_code, regex in (
                ("%Y", r"(?P<year>\d\d\d\d)"),
                ("%m", r"(?P<mon>\d{1,2})"),
                ("%d", r"(?P<mday>\d{1,2})"),
                ("%H", r"(?P<hour>\d{1,2})"),
                ("%M", r"(?P<min>\d{1,2})"),
                ("%S", r"(?P<sec>\d{1,2})"),
            ):
                subfmt_in = subfmt_in.replace(strptime_code, regex)

            if "%" not in subfmt_in:
                subfmt_tuple = (
                    subfmt_tuple[0],
                    re.compile(subfmt_in + "$"),
                    subfmt_tuple[2],
                )
        new_subfmts.append(subfmt_tuple)

    return tuple(new_subfmts)


class TimeFormat:
    """
    Base class for time representations.

    Parameters
    ----------
    val1 : numpy ndarray, list, number, str, or bytes
        Values to initialize the time or times.  Bytes are decoded as ascii.
        Quantities with time units are allowed for formats where the
        interpretation is unambiguous.
    val2 : numpy ndarray, list, or number; optional
        Value(s) to initialize the time or times.  Only used for numerical
        input, to help preserve precision.
    scale : str
        Time scale of input value(s)
    precision : int
        Precision for seconds as floating point
    in_subfmt : str
        Select subformat for inputting string times
    out_subfmt : str
        Select subformat for outputting string times
    from_jd : bool
        If true then val1, val2 are jd1, jd2
    """

    _default_scale = "utc"  # As of astropy 0.4
    _default_precision = 3
    _min_precision = 0
    _max_precision = 9

    subfmts = ()
    _registry = TIME_FORMATS

    # Check that numeric inputs are finite (not nan or inf). This is overridden in
    # subclasses in which nan and inf are valid inputs.
    _check_finite = True

    def __init__(
        self, val1, val2, scale, precision, in_subfmt, out_subfmt, from_jd=False
    ):
        self.scale = scale  # validation of scale done later with _check_scale
        self.precision = precision
        self.in_subfmt = in_subfmt
        self.out_subfmt = out_subfmt

        self._jd1, self._jd2 = None, None

        if from_jd:
            self.jd1 = val1
            self.jd2 = val2
        else:
            val1, val2 = self._check_val_type(val1, val2)
            self.set_jds(val1, val2)

    def __init_subclass__(cls, **kwargs):
        # Register time formats that define a name, but leave out astropy_time since
        # it is not a user-accessible format and is only used for initialization into
        # a different format.
        if "name" in cls.__dict__ and cls.name != "astropy_time":
            # FIXME: check here that we're not introducing a collision with
            # an existing method or attribute; problem is it could be either
            # astropy.time.Time or astropy.time.TimeDelta, and at the point
            # where this is run neither of those classes have necessarily been
            # constructed yet.
            if "value" in cls.__dict__ and not hasattr(cls.value, "fget"):
                raise ValueError("If defined, 'value' must be a property")

            cls._registry[cls.name] = cls

        # If this class defines its own subfmts, preprocess the definitions.
        if "subfmts" in cls.__dict__:
            cls.subfmts = _regexify_subfmts(cls.subfmts)

        return super().__init_subclass__(**kwargs)

    @classmethod
    def _get_allowed_subfmt(cls, subfmt):
        """Get an allowed subfmt for this class, either the input ``subfmt``
        if this is valid or '*' as a default.  This method gets used in situations
        where the format of an existing Time object is changing and so the
        out_ or in_subfmt may need to be coerced to the default '*' if that
        ``subfmt`` is no longer valid.
        """
        try:
            cls._select_subfmts(subfmt)
        except ValueError:
            subfmt = "*"
        return subfmt

    @property
    def in_subfmt(self):
        return self._in_subfmt

    @in_subfmt.setter
    def in_subfmt(self, subfmt):
        # Validate subfmt value for this class, raises ValueError if not.
        self._select_subfmts(subfmt)
        self._in_subfmt = subfmt

    @property
    def out_subfmt(self):
        return self._out_subfmt

    @out_subfmt.setter
    def out_subfmt(self, subfmt):
        # Validate subfmt value for this class, raises ValueError if not.
        self._select_subfmts(subfmt)
        self._out_subfmt = subfmt

    @property
    def jd1(self):
        return self._jd1

    @jd1.setter
    def jd1(self, jd1):
        self._jd1 = _validate_jd_for_storage(jd1)
        if self._jd2 is not None:
            self._jd1, self._jd2 = _broadcast_writeable(self._jd1, self._jd2)

    @property
    def jd2(self):
        return self._jd2

    @jd2.setter
    def jd2(self, jd2):
        self._jd2 = _validate_jd_for_storage(jd2)
        if self._jd1 is not None:
            self._jd1, self._jd2 = _broadcast_writeable(self._jd1, self._jd2)

    @classmethod
    @functools.cache
    def fill_value(cls, subfmt):
        """
        Return a value corresponding to J2000 (2000-01-01 12:00:00) in this format.

        This is used as a fill value for masked arrays to ensure that any ERFA
        operations on the masked array will not fail due to the masked value.
        """
        tm = Time(2451545.0, format="jd", scale="utc")
        return tm.to_value(format=cls.name, subfmt=subfmt)

    def __len__(self):
        return len(self.jd1)

    @property
    def scale(self):
        """Time scale."""
        self._scale = self._check_scale(self._scale)
        return self._scale

    @scale.setter
    def scale(self, val):
        self._scale = val

    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, val):
        if val is None:
            val = self._default_precision
        # Verify precision is 0-9 (inclusive)
        if not (
            isinstance(val, int) and self._min_precision <= val <= self._max_precision
        ):
            raise ValueError(
                "precision attribute must be an int between "
                f"{self._min_precision} and {self._max_precision}"
            )
        self._precision = val

    def _check_finite_vals(self, val1, val2):
        """A helper function to TimeFormat._check_val_type that's meant to be
        optionally bypassed in subclasses that have _check_finite=False
        """
        # val1 cannot contain nan, but val2 can contain nan
        isfinite1 = np.isfinite(val1)
        if val1.size > 1:  # Calling .all() on a scalar is surprisingly slow
            isfinite1 = (
                isfinite1.all()
            )  # Note: arr.all() about 3x faster than np.all(arr)
        elif val1.size == 0:
            isfinite1 = False
        ok1 = (
            val1.dtype.kind == "f" and val1.dtype.itemsize >= 8 and isfinite1
        ) or val1.size == 0
        ok2 = (
            val2 is None
            or (
                val2.dtype.kind == "f"
                and val2.dtype.itemsize >= 8
                and not np.any(np.isinf(val2))
            )
            or val2.size == 0
        )
        if not (ok1 and ok2):
            raise TypeError(
                f"Input values for {self.name} class must be finite doubles"
            )

    def _check_val_type(self, val1, val2):
        """Input value validation, typically overridden by derived classes."""
        if self.__class__._check_finite:
            self._check_finite_vals(val1, val2)

        if getattr(val1, "unit", None) is not None:
            # Convert any quantity-likes to days first, attempting to be
            # careful with the conversion, so that, e.g., large numbers of
            # seconds get converted without losing precision because
            # 1/86400 is not exactly representable as a float.
            val1 = u.Quantity(val1, copy=False)
            if val2 is not None:
                val2 = u.Quantity(val2, copy=False)

            try:
                val1, val2 = quantity_day_frac(val1, val2)
            except u.UnitsError:
                raise u.UnitConversionError(
                    "only quantities with time units can be "
                    "used to instantiate Time instances."
                )
            # We now have days, but the format may expect another unit.
            # On purpose, multiply with 1./day_unit because typically it is
            # 1./erfa.DAYSEC, and inverting it recovers the integer.
            # (This conversion will get undone in format's set_jds, hence
            # there may be room for optimizing this.)
            factor = 1.0 / getattr(self, "unit", 1.0)
            if factor != 1.0:
                val1, carry = two_product(val1, factor)
                carry += val2 * factor
                val1, val2 = two_sum(val1, carry)

        elif getattr(val2, "unit", None) is not None:
            raise TypeError("Cannot mix float and Quantity inputs")

        if val2 is None:
            val2 = np.array(0, dtype=val1.dtype)

        def asarray_or_scalar(val):
            """
            Remove ndarray subclasses since for jd1/jd2 we want a pure ndarray
            or a Python or numpy scalar.
            """
            return val.view(np.ndarray) if isinstance(val, np.ndarray) else val

        return asarray_or_scalar(val1), asarray_or_scalar(val2)

    def _check_scale(self, scale):
        """
        Return a validated scale value.

        If there is a class attribute 'scale' then that defines the default /
        required time scale for this format.  In this case if a scale value was
        provided that needs to match the class default, otherwise return
        the class default.

        Otherwise just make sure that scale is in the allowed list of
        scales.  Provide a different error message if `None` (no value) was
        supplied.
        """
        if scale is None:
            scale = self._default_scale

        if scale not in TIME_SCALES:
            raise ScaleValueError(
                f"Scale value '{scale}' not in allowed values {TIME_SCALES}"
            )

        return scale

    def set_jds(self, val1, val2):
        """
        Set internal jd1 and jd2 from val1 and val2.  Must be provided
        by derived classes.
        """
        raise NotImplementedError

    def to_value(self, parent=None, out_subfmt=None):
        """
        Return time representation from internal jd1 and jd2 in specified
        ``out_subfmt``.

        This is the base method that ignores ``parent`` and uses the ``value``
        property to compute the output. This is done by temporarily setting
        ``self.out_subfmt`` and calling ``self.value``. This is required for
        legacy Format subclasses prior to astropy 4.0  New code should instead
        implement the value functionality in ``to_value()`` and then make the
        ``value`` property be a simple call to ``self.to_value()``.

        Parameters
        ----------
        parent : object
            Parent `~astropy.time.Time` object associated with this
            `~astropy.time.TimeFormat` object
        out_subfmt : str or None
            Output subformt (use existing self.out_subfmt if `None`)

        Returns
        -------
        value : numpy.array, numpy.ma.array
            Array or masked array of formatted time representation values
        """
        # Get value via ``value`` property, overriding out_subfmt temporarily if needed.
        if out_subfmt is not None:
            out_subfmt_orig = self.out_subfmt
            try:
                self.out_subfmt = out_subfmt
                value = self.value
            finally:
                self.out_subfmt = out_subfmt_orig
        else:
            value = self.value

        return value

    @property
    def value(self):
        raise NotImplementedError

    @classmethod
    def _select_subfmts(cls, pattern):
        """
        Return a list of subformats where name matches ``pattern`` using
        fnmatch.

        If no subformat matches pattern then a ValueError is raised.  A special
        case is a format with no allowed subformats, i.e. subfmts=(), and
        pattern='*'.  This is OK and happens when this method is used for
        validation of an out_subfmt.
        """
        if not isinstance(pattern, str):
            raise ValueError("subfmt attribute must be a string")
        elif pattern == "*":
            return cls.subfmts

        subfmts = [x for x in cls.subfmts if fnmatch.fnmatchcase(x[0], pattern)]
        if len(subfmts) == 0:
            if len(cls.subfmts) == 0:
                raise ValueError(f"subformat not allowed for format {cls.name}")
            else:
                subfmt_names = [x[0] for x in cls.subfmts]
                raise ValueError(
                    f"subformat {pattern!r} must match one of "
                    f"{subfmt_names} for format {cls.name}"
                )

        return subfmts

    @classmethod
    def _fill_masked_values(cls, val, val2, mask, in_subfmt):
        """Fill masked values with the fill value for this format.

        This also takes care of broadcasting the outputs to the correct shape.

        Parameters
        ----------
        val : ndarray
            Array of values
        val2 : ndarray, None
            Array of second values (or None)
        mask : ndarray
            Mask array
        in_subfmt : str
            Input subformat

        Returns
        -------
        val, val2 : ndarray
            Arrays with masked values filled with the fill value for this format.
            These are copies of the originals.
        """
        if val2 is None:
            val, mask = np.broadcast_arrays(val, mask)
        else:
            val, val2, mask = np.broadcast_arrays(val, val2, mask)
            val2 = val2.copy()
            val2[mask] = np.zeros_like(val2, shape=())

        # Fill value needs to comply with the specified input subformat. Usually this
        # is "*" for any matching input, but for a custom subformat the fill value
        # needs to be compatible with the specified subformat.
        fill_value = cls.fill_value(in_subfmt)

        # For string types ensure that the numpy string length is long enough to
        # hold the fill value for the specified subformat.
        if (val_kind := val.dtype.kind) in ("U", "S") and (
            new_width := len(fill_value)
        ) > val.dtype.itemsize // (4 if val_kind == "U" else 1):
            val = val.astype(f"{val_kind}{new_width}")  # Makes copy.
        else:
            val = val.copy()
        val[mask] = fill_value
        return val, val2


class TimeNumeric(TimeFormat):
    subfmts = (
        ("float", np.float64, None, np.add),
        ("long", np.longdouble, utils.longdouble_to_twoval, utils.twoval_to_longdouble),
        ("decimal", np.object_, utils.decimal_to_twoval, utils.twoval_to_decimal),
        ("str", np.str_, utils.decimal_to_twoval, utils.twoval_to_string),
        ("bytes", np.bytes_, utils.bytes_to_twoval, utils.twoval_to_bytes),
    )

    def _check_val_type(self, val1, val2):
        """Input value validation, typically overridden by derived classes."""
        # Save original state of val2 because the super()._check_val_type below
        # may change val2 from None to np.array(0). The value is saved in order
        # to prevent a useless and slow call to np.result_type() below in the
        # most common use-case of providing only val1.
        orig_val2_is_none = val2 is None

        if val1.dtype.kind == "f":
            val1, val2 = super()._check_val_type(val1, val2)
        elif not orig_val2_is_none or not (
            val1.dtype.kind in "US"
            or (
                val1.dtype.kind == "O"
                and all(isinstance(v, Decimal) for v in val1.flat)
            )
        ):
            raise TypeError(
                f"for {self.name} class, input should be doubles, string, or Decimal, "
                "and second values are only allowed for doubles."
            )

        val_dtype = (
            val1.dtype if orig_val2_is_none else np.result_type(val1.dtype, val2.dtype)
        )
        subfmts = self._select_subfmts(self.in_subfmt)
        for subfmt, dtype, convert, _ in subfmts:
            if np.issubdtype(val_dtype, dtype):
                break
        else:
            raise ValueError("input type not among selected sub-formats.")

        if convert is not None:
            try:
                val1, val2 = convert(val1, val2)
            except Exception:
                raise TypeError(
                    f"for {self.name} class, input should be (long) doubles, string, "
                    "or Decimal, and second values are only allowed for "
                    "(long) doubles."
                )

        return val1, val2

    def to_value(self, jd1=None, jd2=None, parent=None, out_subfmt=None):
        """
        Return time representation from internal jd1 and jd2.
        Subclasses that require ``parent`` or to adjust the jds should
        override this method.
        """
        # TODO: do this in __init_subclass__?
        if self.__class__.value.fget is not self.__class__.to_value:
            return self.value

        if jd1 is None:
            jd1 = self.jd1
        if jd2 is None:
            jd2 = self.jd2
        if out_subfmt is None:
            out_subfmt = self.out_subfmt
        subfmt = self._select_subfmts(out_subfmt)[0]
        kwargs = {}
        if subfmt[0] in ("str", "bytes"):
            unit = getattr(self, "unit", 1)
            digits = int(np.ceil(np.log10(unit / np.finfo(float).eps)))
            # TODO: allow a way to override the format.
            kwargs["fmt"] = f".{digits}f"
        return subfmt[3](jd1, jd2, **kwargs)

    value = property(to_value)


class TimeJD(TimeNumeric):
    """
    Julian Date time format.

    This represents the number of days since the beginning of
    the Julian Period.
    For example, 2451544.5 in JD is midnight on January 1, 2000.
    """

    name = "jd"

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2)


class TimeMJD(TimeNumeric):
    """
    Modified Julian Date time format.

    This represents the number of days since midnight on November 17, 1858.
    For example, 51544.0 in MJD is midnight on January 1, 2000.
    """

    name = "mjd"

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        jd1, jd2 = day_frac(val1, val2)
        jd1 += erfa.DJM0  # erfa.DJM0=2400000.5 (from erfam.h).
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(self, **kwargs):
        jd1 = self.jd1 - erfa.DJM0  # This cannot lose precision.
        jd2 = self.jd2
        return super().to_value(jd1=jd1, jd2=jd2, **kwargs)

    value = property(to_value)


def _check_val_type_not_quantity(format_name, val1, val2):
    # If val2 is a Quantity, the super() call that follows this check
    # will raise a TypeError.
    if hasattr(val1, "to") and getattr(val1, "unit", None) is not None:
        raise ValueError(
            f"cannot use Quantities for {format_name!r} format, as the unit of year "
            "is defined as 365.25 days, while the length of year is variable "
            "in this format. Use float instead."
        )


class TimeDecimalYear(TimeNumeric):
    """
    Time as a decimal year, with integer values corresponding to midnight of the first
    day of each year.

    The fractional part represents the exact fraction of the year, considering the
    precise number of days in the year (365 or 366). The following example shows
    essentially how the decimal year is computed::

      >>> from astropy.time import Time
      >>> tm = Time("2024-04-05T12:34:00")
      >>> tm0 = Time("2024-01-01T00:00:00")
      >>> tm1 = Time("2025-01-01T00:00:00")
      >>> print(2024 + (tm.jd - tm0.jd) / (tm1.jd - tm0.jd))  # doctest: +FLOAT_CMP
      2024.2609934729812
      >>> print(tm.decimalyear)  # doctest: +FLOAT_CMP
      2024.2609934729812

    Since for this format the length of the year varies between 365 and 366 days, it is
    not possible to use Quantity input, in which a year is always 365.25 days.

    This format is convenient for low-precision applications or for plotting data.
    """

    name = "decimalyear"

    def _check_val_type(self, val1, val2):
        _check_val_type_not_quantity(self.name, val1, val2)
        # if val2 is a Quantity, super() will raise a TypeError.
        return super()._check_val_type(val1, val2)

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.

        sum12, err12 = two_sum(val1, val2)
        iy_start = np.trunc(sum12).astype(int)
        extra, y_frac = two_sum(sum12, -iy_start)
        y_frac += extra + err12

        val = (val1 + val2).astype(np.double)
        iy_start = np.trunc(val).astype(int)

        imon = np.ones_like(iy_start)
        iday = np.ones_like(iy_start)
        ihr = np.zeros_like(iy_start)
        imin = np.zeros_like(iy_start)
        isec = np.zeros_like(y_frac)

        # Possible enhancement: use np.unique to only compute start, stop
        # for unique values of iy_start.
        scale = self.scale.upper().encode("ascii")
        jd1_start, jd2_start = erfa.dtf2d(scale, iy_start, imon, iday, ihr, imin, isec)
        jd1_end, jd2_end = erfa.dtf2d(scale, iy_start + 1, imon, iday, ihr, imin, isec)

        t_start = Time(jd1_start, jd2_start, scale=self.scale, format="jd")
        t_end = Time(jd1_end, jd2_end, scale=self.scale, format="jd")
        t_frac = t_start + (t_end - t_start) * y_frac

        self.jd1, self.jd2 = day_frac(t_frac.jd1, t_frac.jd2)

    def to_value(self, **kwargs):
        scale = self.scale.upper().encode("ascii")
        # precision=0
        iy_start, ims, ids, ihmsfs = erfa.d2dtf(scale, 0, self.jd1, self.jd2)
        imon = np.ones_like(iy_start)
        iday = np.ones_like(iy_start)
        ihr = np.zeros_like(iy_start)
        imin = np.zeros_like(iy_start)
        isec = np.zeros_like(self.jd1)

        # Possible enhancement: use np.unique to only compute start, stop
        # for unique values of iy_start.
        scale = self.scale.upper().encode("ascii")
        jd1_start, jd2_start = erfa.dtf2d(scale, iy_start, imon, iday, ihr, imin, isec)
        jd1_end, jd2_end = erfa.dtf2d(scale, iy_start + 1, imon, iday, ihr, imin, isec)
        # Trying to be precise, but more than float64 not useful.
        dt = (self.jd1 - jd1_start) + (self.jd2 - jd2_start)
        dt_end = (jd1_end - jd1_start) + (jd2_end - jd2_start)
        decimalyear = iy_start + dt / dt_end

        return super().to_value(jd1=decimalyear, jd2=np.float64(0.0), **kwargs)

    value = property(to_value)


class TimeFromEpoch(TimeNumeric):
    """
    Base class for times that represent the interval from a particular
    epoch as a numerical multiple of a unit time interval (e.g. seconds
    or days).
    """

    @classproperty(lazy=True)
    def _epoch(cls):
        # Ideally we would use `def epoch(cls)` here and not have the instance
        # property below. However, this breaks the sphinx API docs generation
        # in a way that was not resolved. See #10406 for details.
        return Time(
            cls.epoch_val,
            cls.epoch_val2,
            scale=cls.epoch_scale,
            format=cls.epoch_format,
        )

    @property
    def epoch(self):
        """Reference epoch time from which the time interval is measured."""
        return self._epoch

    def set_jds(self, val1, val2):
        """
        Initialize the internal jd1 and jd2 attributes given val1 and val2.
        For an TimeFromEpoch subclass like TimeUnix these will be floats giving
        the effective seconds since an epoch time (e.g. 1970-01-01 00:00:00).
        """
        # Form new JDs based on epoch time + time from epoch (converted to JD).
        # One subtlety that might not be obvious is that 1.000 Julian days in
        # UTC can be 86400 or 86401 seconds.  For the TimeUnix format the
        # assumption is that every day is exactly 86400 seconds, so this is, in
        # principle, doing the math incorrectly, *except* that it matches the
        # definition of Unix time which does not include leap seconds.

        # note: use divisor=1./self.unit, since this is either 1 or 1/86400,
        # and 1/86400 is not exactly representable as a float64, so multiplying
        # by that will cause rounding errors. (But inverting it as a float64
        # recovers the exact number)
        day, frac = day_frac(val1, val2, divisor=1.0 / self.unit)

        jd1 = self.epoch.jd1 + day
        jd2 = self.epoch.jd2 + frac

        # For the usual case that scale is the same as epoch_scale, we only need
        # to ensure that abs(jd2) <= 0.5. Since abs(self.epoch.jd2) <= 0.5 and
        # abs(frac) <= 0.5, we can do simple (fast) checks and arithmetic here
        # without another call to day_frac(). Note also that `round(jd2.item())`
        # is about 10x faster than `np.round(jd2)`` for a scalar.
        if self.epoch.scale == self.scale:
            jd1_extra = np.round(jd2) if jd2.shape else round(jd2.item())
            jd1 += jd1_extra
            jd2 -= jd1_extra

            self.jd1, self.jd2 = jd1, jd2
            return

        # Create a temporary Time object corresponding to the new (jd1, jd2) in
        # the epoch scale (e.g. UTC for TimeUnix) then convert that to the
        # desired time scale for this object.
        #
        # A known limitation is that the transform from self.epoch_scale to
        # self.scale cannot involve any metadata like lat or lon.
        try:
            tm = getattr(
                Time(jd1, jd2, scale=self.epoch_scale, format="jd"), self.scale
            )
        except Exception as err:
            raise ScaleValueError(
                f"Cannot convert from '{self.name}' epoch scale '{self.epoch_scale}' "
                f"to specified scale '{self.scale}', got error:\n{err}"
            ) from err

        self.jd1, self.jd2 = day_frac(tm._time.jd1, tm._time.jd2)

    def to_value(self, parent=None, **kwargs):
        # Make sure that scale is the same as epoch scale so we can just
        # subtract the epoch and convert
        if self.scale != self.epoch_scale:
            if parent is None:
                raise ValueError("cannot compute value without parent Time object")
            try:
                tm = getattr(parent, self.epoch_scale)
            except Exception as err:
                raise ScaleValueError(
                    f"Cannot convert from '{self.name}' epoch scale "
                    f"'{self.epoch_scale}' to specified scale '{self.scale}', "
                    f"got error:\n{err}"
                ) from err

            jd1, jd2 = tm._time.jd1, tm._time.jd2
        else:
            jd1, jd2 = self.jd1, self.jd2

        # This factor is guaranteed to be exactly representable, which
        # means time_from_epoch1 is calculated exactly.
        factor = 1.0 / self.unit
        time_from_epoch1 = (jd1 - self.epoch.jd1) * factor
        time_from_epoch2 = (jd2 - self.epoch.jd2) * factor

        return super().to_value(jd1=time_from_epoch1, jd2=time_from_epoch2, **kwargs)

    value = property(to_value)

    @property
    def _default_scale(self):
        return self.epoch_scale


class TimeUnix(TimeFromEpoch):
    """
    Unix time (UTC): seconds from 1970-01-01 00:00:00 UTC, ignoring leap seconds.

    For example, 946684800.0 in Unix time is midnight on January 1, 2000.

    NOTE: this quantity is not exactly unix time and differs from the strict
    POSIX definition by up to 1 second on days with a leap second.  POSIX
    unix time actually jumps backward by 1 second at midnight on leap second
    days while this class value is monotonically increasing at 86400 seconds
    per UTC day.
    """

    name = "unix"
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = "1970-01-01 00:00:00"
    epoch_val2 = None
    epoch_scale = "utc"
    epoch_format = "iso"


class TimeUnixTai(TimeUnix):
    """
    Unix time (TAI): SI seconds elapsed since 1970-01-01 00:00:00 TAI (see caveats).

    This will generally differ from standard (UTC) Unix time by the cumulative
    integral number of leap seconds introduced into UTC since 1972-01-01 UTC
    plus the initial offset of 10 seconds at that date.

    This convention matches the definition of linux CLOCK_TAI
    (https://www.cl.cam.ac.uk/~mgk25/posix-clocks.html),
    and the Precision Time Protocol
    (https://en.wikipedia.org/wiki/Precision_Time_Protocol), which
    is also used by the White Rabbit protocol in High Energy Physics:
    https://white-rabbit.web.cern.ch.

    Caveats:

    - Before 1972, fractional adjustments to UTC were made, so the difference
      between ``unix`` and ``unix_tai`` time is no longer an integer.
    - Because of the fractional adjustments, to be very precise, ``unix_tai``
      is the number of seconds since ``1970-01-01 00:00:00 TAI`` or equivalently
      ``1969-12-31 23:59:51.999918 UTC``.  The difference between TAI and UTC
      at that epoch was 8.000082 sec.
    - On the day of a positive leap second the difference between ``unix`` and
      ``unix_tai`` times increases linearly through the day by 1.0. See also the
      documentation for the `~astropy.time.TimeUnix` class.
    - Negative leap seconds are possible, though none have been needed to date.

    Examples
    --------
      >>> # get the current offset between TAI and UTC
      >>> from astropy.time import Time
      >>> t = Time('2020-01-01', scale='utc')
      >>> t.unix_tai - t.unix
      np.float64(37.0)

      >>> # Before 1972, the offset between TAI and UTC was not integer
      >>> t = Time('1970-01-01', scale='utc')
      >>> t.unix_tai - t.unix  # doctest: +FLOAT_CMP
      np.float64(8.000082)

      >>> # Initial offset of 10 seconds in 1972
      >>> t = Time('1972-01-01', scale='utc')
      >>> t.unix_tai - t.unix
      np.float64(10.0)
    """

    name = "unix_tai"
    epoch_val = "1970-01-01 00:00:00"
    epoch_scale = "tai"


class TimeCxcSec(TimeFromEpoch):
    """
    Chandra X-ray Center seconds from 1998-01-01 00:00:00 TT.
    For example, 63072064.184 is midnight on January 1, 2000.
    """

    name = "cxcsec"
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = "1998-01-01 00:00:00"
    epoch_val2 = None
    epoch_scale = "tt"
    epoch_format = "iso"


class TimeGPS(TimeFromEpoch):
    """GPS time: seconds from 1980-01-06 00:00:00 UTC
    For example, 630720013.0 is midnight on January 1, 2000.

    Notes
    -----
    This implementation is strictly a representation of the number of seconds
    (including leap seconds) since midnight UTC on 1980-01-06.  GPS can also be
    considered as a time scale which is ahead of TAI by a fixed offset
    (to within about 100 nanoseconds).

    For details, see https://www.usno.navy.mil/USNO/time/gps/usno-gps-time-transfer
    """

    name = "gps"
    unit = 1.0 / erfa.DAYSEC  # in days (1 day == 86400 seconds)
    epoch_val = "1980-01-06 00:00:19"
    # above epoch is the same as Time('1980-01-06 00:00:00', scale='utc').tai
    epoch_val2 = None
    epoch_scale = "tai"
    epoch_format = "iso"


class TimePlotDate(TimeFromEpoch):
    """
    Input for a `~matplotlib.axes.Axes` object with ax.xaxis.axis_date():
    1 + number of days from 0001-01-01 00:00:00 UTC.

    This can be used as follow::

      >>> import matplotlib.pyplot as plt
      >>> jyear = np.linspace(2000, 2001, 20)
      >>> t = Time(jyear, format='jyear', scale='utc')
      >>> fig, ax = plt.subplots()
      >>> ax.xaxis.axis_date()
      >>> ax.scatter(t.plot_date, jyear)
      >>> fig.autofmt_xdate()  # orient date labels at a slant
      >>> fig.show()

    For example, 730120.0003703703 is midnight on January 1, 2000.
    """

    # This corresponds to the zero reference time for matplotlib plot_date().
    # Note that TAI and UTC are equivalent at the reference time.
    name = "plot_date"
    unit = 1.0
    epoch_val = 1721424.5  # Time('0001-01-01 00:00:00', scale='tai').jd - 1
    epoch_val2 = None
    epoch_scale = "utc"
    epoch_format = "jd"

    @lazyproperty
    def epoch(self):
        """Reference epoch time from which the time interval is measured."""
        if HAS_MATPLOTLIB:
            from matplotlib.dates import get_epoch

            # Get the matplotlib date epoch as an ISOT string in UTC
            epoch_utc = get_epoch()
            from erfa import ErfaWarning

            with warnings.catch_warnings():
                # Catch possible dubious year warnings from erfa
                warnings.filterwarnings("ignore", category=ErfaWarning)
                _epoch = Time(epoch_utc, scale="utc", format="isot")
            _epoch.format = "jd"
        else:
            # If matplotlib is not installed then the epoch is '0001-01-01'
            _epoch = self._epoch

        return _epoch


class TimeStardate(TimeFromEpoch):
    """
    Stardate: date units from 2318-07-05 12:00:00 UTC.
    For example, stardate 41153.7 is 00:52 on April 30, 2363.
    See http://trekguide.com/Stardates.htm#TNG for calculations and reference points.
    """

    name = "stardate"
    unit = 0.397766856  # Stardate units per day
    epoch_val = "2318-07-05 11:00:00"  # Date and time of stardate 00000.00
    epoch_val2 = None
    epoch_scale = "tai"
    epoch_format = "iso"


class TimeUnique(TimeFormat):
    """
    Base class for time formats that can uniquely create a time object
    without requiring an explicit format specifier.  This class does
    nothing but provide inheritance to identify a class as unique.
    """


class TimeAstropyTime(TimeUnique):
    """
    Instantiate date from an Astropy Time object (or list thereof).

    This is purely for instantiating from a Time object.  The output
    format is the same as the first time instance.
    """

    name = "astropy_time"

    def __new__(
        cls, val1, val2, scale, precision, in_subfmt, out_subfmt, from_jd=False
    ):
        """
        Use __new__ instead of __init__ to output a class instance that
        is the same as the class of the first Time object in the list.
        """
        val1_0 = val1.item(0)
        if not (
            isinstance(val1_0, Time)
            and all(type(val) is type(val1_0) for val in val1.flat)
        ):
            raise TypeError(
                f"Input values for {cls.name} class must all be the same "
                "astropy Time type."
            )

        if scale is None:
            scale = val1_0.scale

        if val1.shape:
            vals = [getattr(val, scale)._time for val in val1]
            jd1 = np.concatenate([np.atleast_1d(val.jd1) for val in vals])
            jd2 = np.concatenate([np.atleast_1d(val.jd2) for val in vals])

            # Collect individual location values and merge into a single location.
            if any(tm.location is not None for tm in val1):
                if any(tm.location is None for tm in val1):
                    raise ValueError(
                        "cannot concatenate times unless all locations "
                        "are set or no locations are set"
                    )
                locations = []
                for tm in val1:
                    location = np.broadcast_to(
                        tm.location, tm._time.jd1.shape, subok=True
                    )
                    locations.append(np.atleast_1d(location))

                location = np.concatenate(locations)

            else:
                location = None
        else:
            val = getattr(val1_0, scale)._time
            jd1, jd2 = val.jd1, val.jd2
            location = val1_0.location

        OutTimeFormat = val1_0._time.__class__
        self = OutTimeFormat(
            jd1, jd2, scale, precision, in_subfmt, out_subfmt, from_jd=True
        )

        # Make a temporary hidden attribute to transfer location back to the
        # parent Time object where it needs to live.
        self._location = location

        return self


class TimeDatetime(TimeUnique):
    """
    Represent date as Python standard library `~datetime.datetime` object.

    Example::

      >>> from astropy.time import Time
      >>> from datetime import datetime
      >>> t = Time(datetime(2000, 1, 2, 12, 0, 0), scale='utc')
      >>> t.iso
      '2000-01-02 12:00:00.000'
      >>> t.tt.datetime
      datetime.datetime(2000, 1, 2, 12, 1, 4, 184000)
    """

    name = "datetime"

    def _check_val_type(self, val1, val2):
        if not all(isinstance(val, datetime.datetime) for val in val1.flat):
            raise TypeError(
                f"Input values for {self.name} class must be datetime objects"
            )
        if val2 is not None:
            raise ValueError(
                f"{self.name} objects do not accept a val2 but you provided {val2}"
            )
        return val1, None

    def set_jds(self, val1, val2):
        """Convert datetime object contained in val1 to jd1, jd2."""
        # Iterate through the datetime objects, getting year, month, etc.
        iterator = np.nditer(
            [val1, None, None, None, None, None, None],
            flags=["refs_ok", "zerosize_ok"],
            op_dtypes=[None] + 5 * [np.intc] + [np.double],
        )
        for val, iy, im, id, ihr, imin, dsec in iterator:
            dt = val.item()

            if dt.tzinfo is not None:
                dt = (dt - dt.utcoffset()).replace(tzinfo=None)

            iy[...] = dt.year
            im[...] = dt.month
            id[...] = dt.day
            ihr[...] = dt.hour
            imin[...] = dt.minute
            dsec[...] = dt.second + dt.microsecond / 1e6

        jd1, jd2 = erfa.dtf2d(
            self.scale.upper().encode("ascii"), *iterator.operands[1:]
        )
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(
        self, timezone=None, leap_second_strict="raise", parent=None, out_subfmt=None
    ):
        """
        Convert to (potentially timezone-aware) `~datetime.datetime` object.

        If ``timezone`` is not ``None``, return a timezone-aware datetime object.

        Since the `~datetime.datetime` class does not natively handle leap seconds, the
        behavior when converting a time within a leap second is controlled by the
        ``leap_second_strict`` argument. For example::

          >>> from astropy.time import Time
          >>> t = Time("2015-06-30 23:59:60.500")
          >>> print(t.to_datetime(leap_second_strict='silent'))
          2015-07-01 00:00:00.500000

        Parameters
        ----------
        timezone : {`~datetime.tzinfo`, None}, optional
            If not `None`, return timezone-aware datetime.
        leap_second_strict : str, optional
            If ``raise`` (default), raise an exception if the time is within a leap
            second. If ``warn`` then issue a warning. If ``silent`` then silently
            handle the leap second.

        Returns
        -------
        `~datetime.datetime`
            If ``timezone`` is not ``None``, output will be timezone-aware.
        """
        if out_subfmt is not None:
            # Out_subfmt not allowed for this format, so raise the standard
            # exception by trying to validate the value.
            self._select_subfmts(out_subfmt)

        if timezone is not None:
            if self._scale != "utc":
                raise ScaleValueError(
                    f"scale is {self._scale}, must be 'utc' when timezone is supplied."
                )

        # Rather than define a value property directly, we have a function,
        # since we want to be able to pass in timezone information.
        scale = self.scale.upper().encode("ascii")
        # 6 for microsec
        iys, ims, ids, ihmsfs = erfa.d2dtf(scale, 6, self.jd1, self.jd2)
        ihrs = ihmsfs["h"]
        imins = ihmsfs["m"]
        isecs = ihmsfs["s"]
        ifracs = ihmsfs["f"]
        iterator = np.nditer(
            [iys, ims, ids, ihrs, imins, isecs, ifracs, None],
            flags=["refs_ok", "zerosize_ok"],
            op_dtypes=7 * [None] + [object],
        )

        for iy, im, id, ihr, imin, isec, ifracsec, out in iterator:
            if isec >= 60:
                isec = isec - 1
                in_leap_second = True
            else:
                in_leap_second = False

            if timezone is not None:
                dt = datetime.datetime(
                    iy, im, id, ihr, imin, isec, ifracsec, tzinfo=TimezoneInfo()
                ).astimezone(timezone)
            else:
                dt = datetime.datetime(iy, im, id, ihr, imin, isec, ifracsec)

            if in_leap_second:
                dt += datetime.timedelta(seconds=1)
                msg = (
                    f"Time {dt} is within a leap second but `datetime` does not "
                    "support leap seconds. Use the `leap_second_strict` argument "
                    "of the `Time.to_datetime()` method with value of 'raise', 'warn', "
                    "or 'silent' to control how leap seconds are handled."
                )
                if leap_second_strict == "raise":
                    raise ValueError(msg)
                elif leap_second_strict == "warn":
                    warnings.warn(msg, AstropyDatetimeLeapSecondWarning)
                elif leap_second_strict != "silent":
                    raise ValueError(
                        f"leap_second_strict must be 'raise', 'warn', or 'silent', "
                        f"not '{leap_second_strict}'"
                    )

            out[...] = dt

        return iterator.operands[-1]

    value = property(to_value)


class TimeYMDHMS(TimeUnique):
    """
    ymdhms: A Time format to represent Time as year, month, day, hour,
    minute, second (thus the name ymdhms).

    Acceptable inputs must have keys or column names in the "YMDHMS" set of
    ``year``, ``month``, ``day`` ``hour``, ``minute``, ``second``:

    - Dict with keys in the YMDHMS set
    - NumPy structured array, record array or astropy Table, or single row
      of those types, with column names in the YMDHMS set

    One can supply a subset of the YMDHMS values, for instance only 'year',
    'month', and 'day'.  Inputs have the following defaults::

      'month': 1, 'day': 1, 'hour': 0, 'minute': 0, 'second': 0

    When the input is supplied as a ``dict`` then each value can be either a
    scalar value or an array.  The values will be broadcast to a common shape.

    Example::

      >>> from astropy.time import Time
      >>> t = Time({'year': 2015, 'month': 2, 'day': 3,
      ...           'hour': 12, 'minute': 13, 'second': 14.567},
      ...           scale='utc')
      >>> t.iso
      '2015-02-03 12:13:14.567'
      >>> t.ymdhms.year
      np.int32(2015)
    """

    name = "ymdhms"

    def _check_val_type(self, val1, val2):
        """
        This checks inputs for the YMDHMS format.

        It is bit more complex than most format checkers because of the flexible
        input that is allowed.  Also, it actually coerces ``val1`` into an appropriate
        dict of ndarrays that can be used easily by ``set_jds()``.  This is useful
        because it makes it easy to get default values in that routine.

        Parameters
        ----------
        val1 : ndarray or None
        val2 : ndarray or None

        Returns
        -------
        val1_as_dict, val2 : val1 as dict or None, val2 is always None

        """
        if val2 is not None:
            raise ValueError("val2 must be None for ymdhms format")

        ymdhms = ["year", "month", "day", "hour", "minute", "second"]

        if val1.dtype.names:
            # Convert to a dict of ndarray
            val1_as_dict = {name: val1[name] for name in val1.dtype.names}

        elif val1.shape == (0,):
            # Input was empty list [], so set to None and set_jds will handle this
            return None, None

        elif (
            val1.dtype.kind == "O"
            and val1.shape == ()
            and isinstance(val1.item(), dict)
        ):
            # Code gets here for input as a dict.  The dict input
            # can be either scalar values or N-d arrays.

            # Extract the item (which is a dict) and broadcast values to the
            # same shape here.
            names = val1.item().keys()
            values = val1.item().values()
            val1_as_dict = dict(zip(names, np.broadcast_arrays(*values)))

        else:
            raise ValueError("input must be dict or table-like")

        # Check that the key names now are good.
        names = val1_as_dict.keys()
        required_names = ymdhms[: len(names)]

        def comma_repr(vals):
            return ", ".join(repr(val) for val in vals)

        bad_names = set(names) - set(ymdhms)
        if bad_names:
            raise ValueError(
                f"{comma_repr(bad_names)} not allowed as YMDHMS key name(s)"
            )

        if set(names) != set(required_names):
            raise ValueError(
                f"for {len(names)} input key names "
                f"you must supply {comma_repr(required_names)}"
            )

        return val1_as_dict, val2

    def set_jds(self, val1, val2):
        if val1 is None:
            # Input was empty list []
            jd1 = np.array([], dtype=np.float64)
            jd2 = np.array([], dtype=np.float64)

        else:
            jd1, jd2 = erfa.dtf2d(
                self.scale.upper().encode("ascii"),
                val1["year"],
                val1.get("month", 1),
                val1.get("day", 1),
                val1.get("hour", 0),
                val1.get("minute", 0),
                val1.get("second", 0),
            )

        self.jd1, self.jd2 = day_frac(jd1, jd2)

    @property
    def value(self):
        scale = self.scale.upper().encode("ascii")
        iys, ims, ids, ihmsfs = erfa.d2dtf(scale, 9, self.jd1, self.jd2)

        out = np.empty(
            self.jd1.shape,
            dtype=[
                ("year", "i4"),
                ("month", "i4"),
                ("day", "i4"),
                ("hour", "i4"),
                ("minute", "i4"),
                ("second", "f8"),
            ],
        )
        out["year"] = iys
        out["month"] = ims
        out["day"] = ids
        out["hour"] = ihmsfs["h"]
        out["minute"] = ihmsfs["m"]
        out["second"] = ihmsfs["s"] + ihmsfs["f"] * 10 ** (-9)
        return out.view(np.recarray)


class TimezoneInfo(datetime.tzinfo):
    """
    Subclass of the `~datetime.tzinfo` object, used in the
    to_datetime method to specify timezones.

    It may be safer in most cases to use a timezone database package like
    pytz rather than defining your own timezones - this class is mainly
    a workaround for users without pytz.
    """

    @u.quantity_input(utc_offset=u.day, dst=u.day)
    def __init__(self, utc_offset=0 * u.day, dst=0 * u.day, tzname=None):
        """
        Parameters
        ----------
        utc_offset : `~astropy.units.Quantity`, optional
            Offset from UTC in days. Defaults to zero.
        dst : `~astropy.units.Quantity`, optional
            Daylight Savings Time offset in days. Defaults to zero
            (no daylight savings).
        tzname : str or None, optional
            Name of timezone

        Examples
        --------
        >>> from datetime import datetime
        >>> from astropy.time import TimezoneInfo  # Specifies a timezone
        >>> import astropy.units as u
        >>> utc = TimezoneInfo()    # Defaults to UTC
        >>> utc_plus_one_hour = TimezoneInfo(utc_offset=1*u.hour)  # UTC+1
        >>> dt_aware = datetime(2000, 1, 1, 0, 0, 0, tzinfo=utc_plus_one_hour)
        >>> print(dt_aware)
        2000-01-01 00:00:00+01:00
        >>> print(dt_aware.astimezone(utc))
        1999-12-31 23:00:00+00:00
        """
        if utc_offset == 0 and dst == 0 and tzname is None:
            tzname = "UTC"
        self._utcoffset = datetime.timedelta(utc_offset.to_value(u.day))
        self._tzname = tzname
        self._dst = datetime.timedelta(dst.to_value(u.day))

    def utcoffset(self, dt):
        return self._utcoffset

    def tzname(self, dt):
        return str(self._tzname)

    def dst(self, dt):
        return self._dst


class TimeString(TimeUnique):
    """
    Base class for string-like time representations.

    This class assumes that anything following the last decimal point to the
    right is a fraction of a second.

    **Fast C-based parser**

    Time format classes can take advantage of a fast C-based parser if the times
    are represented as fixed-format strings with year, month, day-of-month,
    hour, minute, second, OR year, day-of-year, hour, minute, second. This can
    be a factor of 20 or more faster than the pure Python parser.

    Fixed format means that the components always have the same number of
    characters. The Python parser will accept ``2001-9-2`` as a date, but the C
    parser would require ``2001-09-02``.

    A subclass in this case must define a class attribute ``fast_parser_pars``
    which is a `dict` with all of the keys below. An inherited attribute is not
    checked, only an attribute in the class ``__dict__``.

    - ``delims`` (tuple of int): ASCII code for character at corresponding
      ``starts`` position (0 => no character)

    - ``starts`` (tuple of int): position where component starts (including
      delimiter if present). Use -1 for the month component for format that use
      day of year.

    - ``stops`` (tuple of int): position where component ends. Use -1 to
      continue to end of string, or for the month component for formats that use
      day of year.

    - ``break_allowed`` (tuple of int): if true (1) then the time string can
          legally end just before the corresponding component (e.g. "2000-01-01"
          is a valid time but "2000-01-01 12" is not).

    - ``has_day_of_year`` (int): 0 if dates have year, month, day; 1 if year,
      day-of-year
    """

    def __init_subclass__(cls, **kwargs):
        if "fast_parser_pars" in cls.__dict__:
            fpp = cls.fast_parser_pars
            fpp = np.array(
                list(
                    zip(
                        map(chr, fpp["delims"]),
                        fpp["starts"],
                        fpp["stops"],
                        fpp["break_allowed"],
                    )
                ),
                _parse_times.dt_pars,
            )
            if cls.fast_parser_pars["has_day_of_year"]:
                fpp["start"][1] = fpp["stop"][1] = -1
            cls._fast_parser = _parse_times.create_parser(fpp)

        super().__init_subclass__(**kwargs)

    def _check_val_type(self, val1, val2):
        if val1.dtype.kind not in ("S", "U") and val1.size:
            raise TypeError(f"Input values for {self.name} class must be strings")
        if val2 is not None:
            raise ValueError(
                f"{self.name} objects do not accept a val2 but you provided {val2}"
            )
        return val1, None

    def parse_string(self, timestr, subfmts):
        """Read time from a single string, using a set of possible formats."""
        # Datetime components required for conversion to JD by ERFA, along
        # with the default values.
        components = ("year", "mon", "mday", "hour", "min", "sec")
        defaults = (None, 1, 1, 0, 0, 0)
        # Assume that anything following "." on the right side is a
        # floating fraction of a second.
        try:
            idot = timestr.rindex(".")
        except Exception:
            timestr_has_fractional_digits = False
        else:
            timestr, fracsec = timestr[:idot], timestr[idot:]
            fracsec = float(fracsec)
            timestr_has_fractional_digits = True

        for _, strptime_fmt_or_regex, _ in subfmts:
            if isinstance(strptime_fmt_or_regex, str):
                subfmt_has_sec = "%S" in strptime_fmt_or_regex
                try:
                    tm = time.strptime(timestr, strptime_fmt_or_regex)
                except ValueError:
                    continue
                else:
                    vals = [getattr(tm, "tm_" + component) for component in components]

            else:
                tm = re.match(strptime_fmt_or_regex, timestr)
                if tm is None:
                    continue
                tm = tm.groupdict()
                vals = [
                    int(tm.get(component, default))
                    for component, default in zip(components, defaults)
                ]
                subfmt_has_sec = "sec" in tm

            # Add fractional seconds if they were in the original time string
            # and the subformat has seconds. A time like "2022-08-01.123" will
            # never pass this for a format like ISO and will raise a parsing
            # exception.
            if timestr_has_fractional_digits:
                if subfmt_has_sec:
                    vals[-1] = vals[-1] + fracsec
                else:
                    continue

            return vals
        raise ValueError(f"Time {timestr} does not match {self.name} format")

    def set_jds(self, val1, val2):
        """Parse the time strings contained in val1 and set jd1, jd2."""
        # If specific input subformat is required then use the Python parser.
        # Also do this if Time format class does not define `use_fast_parser` or
        # if the fast parser is entirely disabled. Note that `use_fast_parser`
        # is ignored for format classes that don't have a fast parser.
        if (
            self.in_subfmt != "*"
            or "_fast_parser" not in self.__class__.__dict__
            or conf.use_fast_parser == "False"
        ):
            jd1, jd2 = self.get_jds_python(val1, val2)
        else:
            try:
                jd1, jd2 = self.get_jds_fast(val1, val2)
            except Exception:
                # Fall through to the Python parser unless fast is forced.
                if conf.use_fast_parser == "force":
                    raise
                else:
                    jd1, jd2 = self.get_jds_python(val1, val2)

        self.jd1 = jd1
        self.jd2 = jd2

    def get_jds_python(self, val1, val2):
        """Parse the time strings contained in val1 and get jd1, jd2."""
        # Select subformats based on current self.in_subfmt
        subfmts = self._select_subfmts(self.in_subfmt)
        # Be liberal in what we accept: convert bytes to ascii.
        # Here .item() is needed for arrays with entries of unequal length,
        # to strip trailing 0 bytes.
        to_string = (
            str if val1.dtype.kind == "U" else lambda x: str(x.item(), encoding="ascii")
        )
        iterator = np.nditer(
            [val1, None, None, None, None, None, None],
            flags=["zerosize_ok"],
            op_dtypes=[None] + 5 * [np.intc] + [np.double],
        )
        for val, iy, im, id, ihr, imin, dsec in iterator:
            val = to_string(val)
            (
                iy[...],
                im[...],
                id[...],
                ihr[...],
                imin[...],
                dsec[...],
            ) = self.parse_string(val, subfmts)

        jd1, jd2 = erfa.dtf2d(
            self.scale.upper().encode("ascii"), *iterator.operands[1:]
        )
        # The iterator above eats the mask...
        if isinstance(val1, Masked):
            jd1 = Masked(jd1, mask=val1.mask.copy())
        jd1, jd2 = day_frac(jd1, jd2)

        return jd1, jd2

    def get_jds_fast(self, val1, val2):
        """Use fast C parser to parse time strings in val1 and get jd1, jd2."""
        # Handle bytes or str input and convert to uint8.  We need to the
        # dtype _parse_times.dt_u1 instead of uint8, since otherwise it is
        # not possible to create a gufunc with structured dtype output.
        # See note about ufunc type resolver in pyerfa/erfa/ufunc.c.templ.
        if val1.dtype.kind == "U":
            # Note: val1.astype('S') is *very* slow, so we check ourselves
            # that the input is pure ASCII.
            val1_uint32 = val1.view((np.uint32, val1.dtype.itemsize // 4))
            if np.any(val1_uint32 > 127):
                raise ValueError("input is not pure ASCII")

            # It might be possible to avoid making a copy via astype with
            # cleverness in parse_times.c but leave that for another day.
            chars = val1_uint32.astype(_parse_times.dt_u1)

        else:
            chars = val1.view((_parse_times.dt_u1, val1.dtype.itemsize))

        # Call the fast parsing ufunc.
        time_struct = self._fast_parser(chars)
        jd1, jd2 = erfa.dtf2d(
            self.scale.upper().encode("ascii"),
            time_struct["year"],
            time_struct["month"],
            time_struct["day"],
            time_struct["hour"],
            time_struct["minute"],
            time_struct["second"],
        )
        return day_frac(jd1, jd2)

    def str_kwargs(self):
        """
        Generator that yields a dict of values corresponding to the
        calendar date and time for the internal JD values.
        """
        scale = (self.scale.upper().encode("ascii"),)
        iys, ims, ids, ihmsfs = erfa.d2dtf(scale, self.precision, self.jd1, self.jd2)

        # Get the str_fmt element of the first allowed output subformat
        _, _, str_fmt = self._select_subfmts(self.out_subfmt)[0]

        yday = None
        has_yday = "{yday:" in str_fmt

        ihrs = ihmsfs["h"]
        imins = ihmsfs["m"]
        isecs = ihmsfs["s"]
        ifracs = ihmsfs["f"]
        for iy, im, id, ihr, imin, isec, ifracsec in np.nditer(
            [iys, ims, ids, ihrs, imins, isecs, ifracs], flags=["zerosize_ok"]
        ):
            if has_yday:
                yday = datetime.datetime(iy, im, id).timetuple().tm_yday

            yield {
                "year": int(iy),
                "mon": int(im),
                "day": int(id),
                "hour": int(ihr),
                "min": int(imin),
                "sec": int(isec),
                "fracsec": int(ifracsec),
                "yday": yday,
            }

    def format_string(self, str_fmt, **kwargs):
        """Write time to a string using a given format.

        By default, just interprets str_fmt as a format string,
        but subclasses can add to this.
        """
        return str_fmt.format(**kwargs)

    @property
    def value(self):
        # Select the first available subformat based on current
        # self.out_subfmt
        subfmts = self._select_subfmts(self.out_subfmt)
        _, _, str_fmt = subfmts[0]

        # TODO: fix this ugly hack
        if self.precision > 0 and str_fmt.endswith("{sec:02d}"):
            str_fmt += ".{fracsec:0" + str(self.precision) + "d}"

        # Try to optimize this later.  Can't pre-allocate because length of
        # output could change, e.g. year rolls from 999 to 1000.
        outs = []
        for kwargs in self.str_kwargs():
            outs.append(str(self.format_string(str_fmt, **kwargs)))

        return np.array(outs).reshape(self.jd1.shape)


class TimeISO(TimeString):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DD HH:MM:SS.sss...".
    For example, 2000-01-01 00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = "iso"
    subfmts = (
        (
            "date_hms",
            "%Y-%m-%d %H:%M:%S",
            # XXX To Do - use strftime for output ??
            "{year:d}-{mon:02d}-{day:02d} {hour:02d}:{min:02d}:{sec:02d}",
        ),
        (
            "date_hm",
            "%Y-%m-%d %H:%M",
            "{year:d}-{mon:02d}-{day:02d} {hour:02d}:{min:02d}",
        ),
        ("date", "%Y-%m-%d", "{year:d}-{mon:02d}-{day:02d}"),
    )

    # Define positions and starting delimiter for year, month, day, hour,
    # minute, seconds components of an ISO time. This is used by the fast
    # C-parser parse_ymdhms_times()
    #
    #  "2000-01-12 13:14:15.678"
    #   01234567890123456789012
    #   yyyy-mm-dd hh:mm:ss.fff
    # Parsed as ('yyyy', '-mm', '-dd', ' hh', ':mm', ':ss', '.fff')
    fast_parser_pars = dict(
        delims=(0, ord("-"), ord("-"), ord(" "), ord(":"), ord(":"), ord(".")),
        starts=(0, 4, 7, 10, 13, 16, 19),
        stops=(3, 6, 9, 12, 15, 18, -1),
        # Break allowed *before*
        #              y  m  d  h  m  s  f
        break_allowed=(0, 0, 0, 1, 0, 1, 1),
        has_day_of_year=0,
    )

    def parse_string(self, timestr, subfmts):
        # Handle trailing 'Z' for UTC time
        if timestr.endswith("Z"):
            if self.scale != "utc":
                raise ValueError("Time input terminating in 'Z' must have scale='UTC'")
            timestr = timestr[:-1]
        return super().parse_string(timestr, subfmts)


class TimeISOT(TimeISO):
    """
    ISO 8601 compliant date-time format "YYYY-MM-DDTHH:MM:SS.sss...".
    This is the same as TimeISO except for a "T" instead of space between
    the date and time.
    For example, 2000-01-01T00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = "isot"
    subfmts = (
        (
            "date_hms",
            "%Y-%m-%dT%H:%M:%S",
            "{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}",
        ),
        (
            "date_hm",
            "%Y-%m-%dT%H:%M",
            "{year:d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}",
        ),
        ("date", "%Y-%m-%d", "{year:d}-{mon:02d}-{day:02d}"),
    )

    # See TimeISO for explanation
    fast_parser_pars = dict(
        delims=(0, ord("-"), ord("-"), ord("T"), ord(":"), ord(":"), ord(".")),
        starts=(0, 4, 7, 10, 13, 16, 19),
        stops=(3, 6, 9, 12, 15, 18, -1),
        # Break allowed *before*
        #              y  m  d  h  m  s  f
        break_allowed=(0, 0, 0, 1, 0, 1, 1),
        has_day_of_year=0,
    )


class TimeYearDayTime(TimeISO):
    """
    Year, day-of-year and time as "YYYY:DOY:HH:MM:SS.sss...".
    The day-of-year (DOY) goes from 001 to 365 (366 in leap years).
    For example, 2000:001:00:00:00.000 is midnight on January 1, 2000.

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date_hm': date + hours, mins
    - 'date': date
    """

    name = "yday"
    subfmts = (
        (
            "date_hms",
            "%Y:%j:%H:%M:%S",
            "{year:d}:{yday:03d}:{hour:02d}:{min:02d}:{sec:02d}",
        ),
        ("date_hm", "%Y:%j:%H:%M", "{year:d}:{yday:03d}:{hour:02d}:{min:02d}"),
        ("date", "%Y:%j", "{year:d}:{yday:03d}"),
    )

    # Define positions and starting delimiter for year, month, day, hour,
    # minute, seconds components of an ISO time. This is used by the fast
    # C-parser parse_ymdhms_times()
    #
    #  "2000:123:13:14:15.678"
    #   012345678901234567890
    #   yyyy:ddd:hh:mm:ss.fff
    # Parsed as ('yyyy', ':ddd', ':hh', ':mm', ':ss', '.fff')
    #
    # delims: character at corresponding `starts` position (0 => no character)
    # starts: position where component starts (including delimiter if present)
    # stops: position where component ends (-1 => continue to end of string)

    fast_parser_pars = dict(
        delims=(0, 0, ord(":"), ord(":"), ord(":"), ord(":"), ord(".")),
        starts=(0, -1, 4, 8, 11, 14, 17),
        stops=(3, -1, 7, 10, 13, 16, -1),
        # Break allowed before:
        #              y  m  d  h  m  s  f
        break_allowed=(0, 0, 0, 1, 0, 1, 1),
        has_day_of_year=1,
    )


class TimeDatetime64(TimeISOT):
    name = "datetime64"

    def _check_val_type(self, val1, val2):
        if not val1.dtype.kind == "M":
            if val1.size > 0:
                raise TypeError(
                    f"Input values for {self.name} class must be datetime64 objects"
                )
            else:
                val1 = np.array([], "datetime64[D]")
        if val2 is not None:
            raise ValueError(
                f"{self.name} objects do not accept a val2 but you provided {val2}"
            )

        return val1, None

    def set_jds(self, val1, val2):
        # If there are any masked values in the ``val1`` datetime64 array
        # ('NaT') then stub them with a valid date so downstream parse_string
        # will work.  The value under the mask is arbitrary but a "modern" date
        # is good.
        mask = np.isnat(val1)
        masked = np.any(mask)
        if masked:
            val1 = val1.copy()
            val1[mask] = "2000"

        # Make sure M(onth) and Y(ear) dates will parse and convert to bytestring
        if val1.dtype.name in ["datetime64[M]", "datetime64[Y]"]:
            val1 = val1.astype("datetime64[D]")
        val1 = val1.astype("S")

        # Standard ISO string parsing now
        super().set_jds(val1, val2)

        # Finally apply mask if necessary
        if masked:
            self.jd1 = Masked(self.jd1, mask=mask)
            self.jd2 = Masked(self.jd2, mask=mask)

    @property
    def value(self):
        precision = self.precision
        self.precision = 9
        ret = super().value
        self.precision = precision
        return ret.astype("datetime64")


class TimeFITS(TimeString):
    """
    FITS format: "[±Y]YYYY-MM-DD[THH:MM:SS[.sss]]".

    ISOT but can give signed five-digit year (mostly for negative years);

    The allowed subformats are:

    - 'date_hms': date + hours, mins, secs (and optional fractional secs)
    - 'date': date
    - 'longdate_hms': as 'date_hms', but with signed 5-digit year
    - 'longdate': as 'date', but with signed 5-digit year

    See Rots et al., 2015, A&A 574:A36 (arXiv:1409.7583).
    """

    name = "fits"
    subfmts = (
        (
            "date_hms",
            (
                r"(?P<year>\d{4})-(?P<mon>\d\d)-(?P<mday>\d\d)T"
                r"(?P<hour>\d\d):(?P<min>\d\d):(?P<sec>\d\d(\.\d*)?)"
            ),
            "{year:04d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}",
        ),
        (
            "date",
            r"(?P<year>\d{4})-(?P<mon>\d\d)-(?P<mday>\d\d)",
            "{year:04d}-{mon:02d}-{day:02d}",
        ),
        (
            "longdate_hms",
            (
                r"(?P<year>[+-]\d{5})-(?P<mon>\d\d)-(?P<mday>\d\d)T"
                r"(?P<hour>\d\d):(?P<min>\d\d):(?P<sec>\d\d(\.\d*)?)"
            ),
            "{year:+06d}-{mon:02d}-{day:02d}T{hour:02d}:{min:02d}:{sec:02d}",
        ),
        (
            "longdate",
            r"(?P<year>[+-]\d{5})-(?P<mon>\d\d)-(?P<mday>\d\d)",
            "{year:+06d}-{mon:02d}-{day:02d}",
        ),
    )
    # Add the regex that parses the scale and possible realization.
    # Support for this is deprecated.  Read old style but no longer write
    # in this style.
    subfmts = tuple(
        (
            subfmt[0],
            subfmt[1] + r"(\((?P<scale>\w+)(\((?P<realization>\w+)\))?\))?",
            subfmt[2],
        )
        for subfmt in subfmts
    )

    def parse_string(self, timestr, subfmts):
        """Read time and deprecated scale if present."""
        # Try parsing with any of the allowed sub-formats.
        for _, regex, _ in subfmts:
            tm = re.match(regex, timestr)
            if tm:
                break
        else:
            raise ValueError(f"Time {timestr} does not match {self.name} format")
        tm = tm.groupdict()
        # Scale and realization are deprecated and strings in this form
        # are no longer created.  We issue a warning but still use the value.
        if tm["scale"] is not None:
            warnings.warn(
                "FITS time strings should no longer have embedded time scale.",
                AstropyDeprecationWarning,
            )
            # If a scale was given, translate from a possible deprecated
            # timescale identifier to the scale used by Time.
            fits_scale = tm["scale"].upper()
            scale = FITS_DEPRECATED_SCALES.get(fits_scale, fits_scale.lower())
            if scale not in TIME_SCALES:
                raise ValueError(
                    f"Scale {scale!r} is not in the allowed scales "
                    f"{sorted(TIME_SCALES)}"
                )
            # If no scale was given in the initialiser, set the scale to
            # that given in the string.  Realization is ignored
            # and is only supported to allow old-style strings to be
            # parsed.
            if self._scale is None:
                self._scale = scale
            if scale != self.scale:
                raise ValueError(
                    f"Input strings for {self.name} class must all "
                    "have consistent time scales."
                )
        return [
            int(tm["year"]),
            int(tm["mon"]),
            int(tm["mday"]),
            int(tm.get("hour", 0)),
            int(tm.get("min", 0)),
            float(tm.get("sec", 0.0)),
        ]

    @property
    def value(self):
        """Convert times to strings, using signed 5 digit if necessary."""
        if "long" not in self.out_subfmt:
            # If we have times before year 0 or after year 9999, we can
            # output only in a "long" format, using signed 5-digit years.
            jd = self.jd1 + self.jd2
            if jd.size and (jd.min() < 1721425.5 or jd.max() >= 5373484.5):
                self.out_subfmt = "long" + self.out_subfmt
        return super().value


class TimeEpochDate(TimeNumeric):
    """
    Base class for support of Besselian and Julian epoch dates.
    """

    _default_scale = "tt"  # As of astropy 3.2, this is no longer 'utc'.

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # validate scale.
        epoch_to_jd = getattr(erfa, self.epoch_to_jd)
        jd1, jd2 = epoch_to_jd(val1 + val2)
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(self, **kwargs):
        jd_to_epoch = getattr(erfa, self.jd_to_epoch)
        value = jd_to_epoch(self.jd1, self.jd2)
        return super().to_value(jd1=value, jd2=np.float64(0.0), **kwargs)

    value = property(to_value)


class TimeBesselianEpoch(TimeEpochDate):
    """Besselian Epoch year as decimal value(s) like 1950.0.

    For information about this epoch format, see:
    `<https://en.wikipedia.org/wiki/Epoch_(astronomy)#Besselian_years>`_.

    The astropy Time class uses the ERFA functions ``epb2jd`` and ``epb`` to convert
    between Besselian epoch years and Julian dates. This is roughly equivalent to the
    following formula (see the wikipedia page for the reference)::

      B = 1900.0 + (Julian date - 2415020.31352) / 365.242198781

    Since for this format the length of the year varies, input needs to be floating
    point; it is not possible to use Quantity input, for which a year always equals
    365.25 days.

    The Besselian epoch year is used for expressing the epoch or equinox in older source
    catalogs, but it has been largely replaced by the Julian epoch year.
    """

    name = "byear"
    epoch_to_jd = "epb2jd"
    jd_to_epoch = "epb"

    def _check_val_type(self, val1, val2):
        _check_val_type_not_quantity(self.name, val1, val2)
        # FIXME: is val2 really okay here?
        return super()._check_val_type(val1, val2)


class TimeJulianEpoch(TimeEpochDate):
    """Julian epoch year as decimal value(s) like 2000.0.

    This format is based the Julian year which is exactly 365.25 days/year and a day is
    exactly 86400 SI seconds.

    The Julian epoch year is defined so that 2000.0 is 12:00 TT on January 1, 2000.
    Using astropy this is expressed as::

      >>> from astropy.time import Time
      >>> import astropy.units as u
      >>> j2000_epoch = Time("2000-01-01T12:00:00", scale="tt")
      >>> print(j2000_epoch.jyear)  # doctest: +FLOAT_CMP
      2000.0
      >>> print((j2000_epoch + 365.25 * u.day).jyear)  # doctest: +FLOAT_CMP
      2001.0

    The Julian year is commonly used in astronomy for expressing the epoch of a source
    catalog or the time of an observation. The Julian epoch year is sometimes written as
    a string like "J2001.5" with a preceding "J". You can initialize a ``Time`` object with
    such a string::

      >>> print(Time("J2001.5").jyear)  # doctest: +FLOAT_CMP
      2001.5

    See also: `<https://en.wikipedia.org/wiki/Julian_year_(astronomy)>`_.
    """

    name = "jyear"
    unit = erfa.DJY  # 365.25, the Julian year, for conversion to quantities
    epoch_to_jd = "epj2jd"
    jd_to_epoch = "epj"


class TimeEpochDateString(TimeString):
    """
    Base class to support string Besselian and Julian epoch dates
    such as 'B1950.0' or 'J2000.0' respectively.
    """

    _default_scale = "tt"  # As of astropy 3.2, this is no longer 'utc'.

    def set_jds(self, val1, val2):
        epoch_prefix = self.epoch_prefix
        # Be liberal in what we accept: convert bytes to ascii.
        to_string = (
            str if val1.dtype.kind == "U" else lambda x: str(x.item(), encoding="ascii")
        )
        iterator = np.nditer(
            [val1, None], op_dtypes=[val1.dtype, np.double], flags=["zerosize_ok"]
        )
        for val, years in iterator:
            try:
                time_str = to_string(val)
                epoch_type, year_str = time_str[0], time_str[1:]
                year = float(year_str)
                if epoch_type.upper() != epoch_prefix:
                    raise ValueError
            except (IndexError, ValueError, UnicodeEncodeError):
                raise ValueError(f"Time {val} does not match {self.name} format")
            else:
                years[...] = year

        self._check_scale(self._scale)  # validate scale.
        epoch_to_jd = getattr(erfa, self.epoch_to_jd)
        jd1, jd2 = epoch_to_jd(iterator.operands[-1])
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    @property
    def value(self):
        jd_to_epoch = getattr(erfa, self.jd_to_epoch)
        years = jd_to_epoch(self.jd1, self.jd2)
        # Use old-style format since it is a factor of 2 faster
        str_fmt = self.epoch_prefix + "%." + str(self.precision) + "f"
        outs = [str_fmt % year for year in years.flat]
        return np.array(outs).reshape(self.jd1.shape)


class TimeBesselianEpochString(TimeEpochDateString):
    """Besselian Epoch year as string value(s) like 'B1950.0'."""

    name = "byear_str"
    epoch_to_jd = "epb2jd"
    jd_to_epoch = "epb"
    epoch_prefix = "B"


class TimeJulianEpochString(TimeEpochDateString):
    """Julian Epoch year as string value(s) like 'J2000.0'."""

    name = "jyear_str"
    epoch_to_jd = "epj2jd"
    jd_to_epoch = "epj"
    epoch_prefix = "J"


class TimeDeltaFormat(TimeFormat):
    """Base class for time delta representations."""

    _registry = TIME_DELTA_FORMATS
    _default_precision = 3
    # Somewhat arbitrary values that are effectively no limit for precision.
    _min_precision = -99
    _max_precision = 99

    def _check_scale(self, scale):
        """
        Check that the scale is in the allowed list of scales, or is `None`.
        """
        if scale is not None and scale not in TIME_DELTA_SCALES:
            raise ScaleValueError(
                f"Scale value '{scale}' not in allowed values {TIME_DELTA_SCALES}"
            )

        return scale


class TimeDeltaNumeric(TimeDeltaFormat, TimeNumeric):
    _check_finite = False

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        self.jd1, self.jd2 = day_frac(val1, val2, divisor=1.0 / self.unit)

    def to_value(self, **kwargs):
        # Note that 1/unit is always exactly representable, so the
        # following multiplications are exact.
        factor = 1.0 / self.unit
        jd1 = self.jd1 * factor
        jd2 = self.jd2 * factor
        return super().to_value(jd1=jd1, jd2=jd2, **kwargs)

    value = property(to_value)


class TimeDeltaSec(TimeDeltaNumeric):
    """Time delta in SI seconds."""

    name = "sec"
    unit = 1.0 / erfa.DAYSEC  # for quantity input


class TimeDeltaJD(TimeDeltaNumeric, TimeUnique):
    """Time delta in Julian days (86400 SI seconds)."""

    name = "jd"
    unit = 1.0


class TimeDeltaDatetime(TimeDeltaFormat, TimeUnique):
    """Time delta in datetime.timedelta."""

    name = "datetime"

    def _check_val_type(self, val1, val2):
        if not all(isinstance(val, datetime.timedelta) for val in val1.flat):
            raise TypeError(
                f"Input values for {self.name} class must be datetime.timedelta objects"
            )
        if val2 is not None:
            raise ValueError(
                f"{self.name} objects do not accept a val2 but you provided {val2}"
            )
        return val1, None

    def set_jds(self, val1, val2):
        self._check_scale(self._scale)  # Validate scale.
        iterator = np.nditer(
            [val1, None, None],
            flags=["refs_ok", "zerosize_ok"],
            op_dtypes=[None, np.double, np.double],
        )

        day = datetime.timedelta(days=1)
        for val, jd1, jd2 in iterator:
            jd1[...], other = divmod(val.item(), day)
            jd2[...] = other / day

        self.jd1, self.jd2 = day_frac(iterator.operands[-2], iterator.operands[-1])

    @property
    def value(self):
        iterator = np.nditer(
            [self.jd1, self.jd2, None],
            flags=["refs_ok", "zerosize_ok"],
            op_dtypes=[None, None, object],
        )

        for jd1, jd2, out in iterator:
            jd1_, jd2_ = day_frac(jd1, jd2)
            out[...] = datetime.timedelta(days=jd1_, microseconds=jd2_ * 86400 * 1e6)

        return iterator.operands[-1]


class TimeDeltaQuantityString(TimeDeltaFormat, TimeUnique):
    """Time delta as a string with one or more Quantity components.

    This format provides a human-readable multi-scale string representation of a time
    delta. It is convenient for applications like a configuration file or a command line
    option.

    The format is specified as follows:

    - The string is a sequence of one or more components.
    - Each component is a number followed by an astropy unit of time.
    - For input, whitespace within the string is allowed but optional.
    - For output, there is a single space between components.
    - The allowed components are listed below.
    - The order (yr, d, hr, min, s) is fixed but individual components are optional.

    The allowed input component units are shown below:

    - "yr": years (365.25 days)
    - "d": days (24 hours)
    - "hr": hours (60 minutes)
    - "min": minutes (60 seconds)
    - "s": seconds

    .. Note:: These definitions correspond to physical units of time and are NOT
       calendar date intervals. Thus adding "1yr" to "2000-01-01 00:00:00" will give
       "2000-12-31 06:00:00" instead of "2001-01-01 00:00:00".

    The ``out_subfmt`` attribute specifies the components to be included in the string
    output.  The default is ``"multi"`` which represents the time delta as
    ``"<days>d <hours>hr <minutes>min <seconds>s"``, where only non-zero components are
    included.

    - "multi": multiple components, e.g. "2d 3hr 15min 5.6s"
    - "yr": years
    - "d": days
    - "hr": hours
    - "min": minutes
    - "s": seconds

    Examples
    --------
    >>> from astropy.time import Time, TimeDelta
    >>> import astropy.units as u

    >>> print(TimeDelta("1yr"))
    365d 6hr

    >>> print(Time("2000-01-01") + TimeDelta("1yr"))
    2000-12-31 06:00:00.000
    >>> print(TimeDelta("+3.6d"))
    3d 14hr 24min
    >>> print(TimeDelta("-3.6d"))
    -3d 14hr 24min
    >>> print(TimeDelta("1yr 3.6d", out_subfmt="d"))
    368.85d

    >>> td = TimeDelta(40 * u.hr)
    >>> print(td.to_value(format="quantity_str"))
    1d 16hr
    >>> print(td.to_value(format="quantity_str", subfmt="d"))
    1.667d
    >>> td.precision = 9
    >>> print(td.to_value(format="quantity_str", subfmt="d"))
    1.666666667d
    """

    name = "quantity_str"

    subfmts = (
        ("multi", None, None),
        ("yr", None, None),
        ("d", None, None),
        ("hr", None, None),
        ("min", None, None),
        ("s", None, None),
    )

    # Regex to parse "1.02yr 2.2d 3.12hr 4.322min 5.6s" where each element is optional
    # but the order is fixed. Each element is a float with optional exponent. Each
    # element is named.
    re_float = r"(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    re_ydhms = re.compile(
        rf"""^ \s*
        (?P<sign>[-+])? \s*  # Optional sign
        (?=[^-+\s])  # At least one character which is not a sign or whitespace
        ((?P<yr>{re_float}) \s* yr \s*)?
        ((?P<d>{re_float}) \s* d \s*)?
        ((?P<hr>{re_float}) \s* hr \s*)?
        ((?P<min>{re_float}) \s* min \s*)?
        ((?P<s>{re_float}) \s* s)?
        \s* $
        """,
        re.VERBOSE,
    )

    def _check_val_type(self, val1, val2):
        if val1.dtype.kind not in ("S", "U") and val1.size:
            raise TypeError(f"Input values for {self.name} class must be strings")
        if val2 is not None:
            raise ValueError(
                f"{self.name} objects do not accept a val2 but you provided {val2}"
            )
        return val1, None

    def parse_string(self, timestr):
        """Read time from a single string"""
        components = ("yr", "d", "hr", "min", "s")

        if (match := self.re_ydhms.match(timestr)) is None:
            raise ValueError(
                f"Time delta '{timestr}' does not match {self.name} format"
            )

        tm = match.groupdict()
        vals = [float(tm[component] or 0.0) for component in components]
        if tm["sign"] == "-":
            vals = [-val for val in vals]

        return vals

    def set_jds(self, val1, val2):
        """Parse the time strings contained in val1 and get jd1, jd2."""
        # Be liberal in what we accept: convert bytes to ascii.
        # Here .item() is needed for arrays with entries of unequal length,
        # to strip trailing 0 bytes.
        to_string = (
            str if val1.dtype.kind == "U" else lambda x: str(x.item(), encoding="ascii")
        )
        iterator = np.nditer(
            [val1, None, None, None, None, None],
            flags=["zerosize_ok"],
            op_dtypes=[None] + 5 * [np.double],
        )
        for val, yr, day, hr, min, sec in iterator:
            val = to_string(val)
            (
                yr[...],
                day[...],
                hr[...],
                min[...],
                sec[...],
            ) = self.parse_string(val)

        yrs, days, hrs, mins, secs = iterator.operands[1:]

        jd1 = yrs * 365.25 + days  # Exact in the case that yrs and days are integer
        jd2 = hrs / 24.0 + mins / 1440.0 + secs / 86400.0  # Inexact
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(self, parent=None, out_subfmt=None):
        out_subfmt = out_subfmt or self.out_subfmt
        subfmt = self._get_allowed_subfmt(out_subfmt)

        iterator = np.nditer(
            [self.jd1, self.jd2, None],
            flags=["refs_ok", "zerosize_ok"],
            op_dtypes=[None, None, object],
        )

        for jd1, jd2, out in iterator:
            jd = jd1 + jd2
            if jd < 0:
                jd1, jd2, jd = -jd1, -jd2, -jd  # Flip all signs
                sign = "-"
            else:
                sign = ""

            if subfmt in ["*", "multi"]:
                comps = self.get_multi_comps(jd1, jd2)

            else:
                value = (jd * u.day).to_value(subfmt)
                value = np.round(value, self.precision)
                comps = [f"{value}{subfmt}"]

            out[...] = sign + " ".join(comps)

        return iterator.operands[-1]

    def get_multi_comps(self, jd1, jd2):
        jd, remainder = two_sum(jd1, jd2)
        days = int(np.floor(jd))
        jd -= days
        jd += remainder

        hours = int(np.floor(jd * 24.0))
        jd -= hours / 24.0
        mins = int(np.floor(jd * 1440.0))
        jd -= mins / 1440.0
        secs = np.round(jd * 86400.0, self.precision)

        comp_vals = [days, hours, mins, secs]
        if secs >= 60.0:
            self.fix_comp_vals_overflow(comp_vals)

        comps = [
            f"{comp_val}{name}"
            for comp_val, name in zip(comp_vals, ("d", "hr", "min", "s"))
            if comp_val != 0
        ]
        if not comps:
            comps = ["0.0s"]

        return comps

    @staticmethod
    def fix_comp_vals_overflow(comp_vals):
        comp_maxes = (None, 24, 60, 60.0)
        for ii in [3, 2, 1]:
            comp_val = comp_vals[ii]
            comp_max = comp_maxes[ii]
            if comp_val >= comp_max:
                comp_vals[ii] -= comp_max
                comp_vals[ii - 1] += 1

    @property
    def value(self):
        return self.to_value()


def _validate_jd_for_storage(jd):
    if isinstance(jd, (float, int)):
        return np.array(jd, dtype=float)
    if isinstance(jd, np.generic) and (
        (jd.dtype.kind == "f" and jd.dtype.itemsize <= 8) or jd.dtype.kind in "iu"
    ):
        return np.array(jd, dtype=float)
    elif isinstance(jd, np.ndarray) and jd.dtype.kind == "f" and jd.dtype.itemsize == 8:
        return jd
    else:
        raise TypeError(
            "JD values must be arrays (possibly zero-dimensional) "
            f"of floats but we got {jd!r} of type {type(jd)}"
        )


def _broadcast_writeable(jd1, jd2):
    if jd1.shape == jd2.shape:
        return jd1, jd2
    # When using broadcast_arrays, *both* are flagged with
    # warn-on-write, even the one that wasn't modified, and
    # require "C" only clears the flag if it actually copied
    # anything.
    shape = np.broadcast(jd1, jd2).shape
    if jd1.shape == shape:
        s_jd1 = jd1
    else:
        s_jd1 = np.require(np.broadcast_to(jd1, shape), requirements=["C", "W"])
    if jd2.shape == shape:
        s_jd2 = jd2
    else:
        s_jd2 = np.require(np.broadcast_to(jd2, shape), requirements=["C", "W"])
    return s_jd1, s_jd2


# Import symbols from core.py that are used in this module. This succeeds
# because __init__.py imports format.py just before core.py.
from .core import TIME_DELTA_SCALES, TIME_SCALES, ScaleValueError, Time  # noqa: E402
