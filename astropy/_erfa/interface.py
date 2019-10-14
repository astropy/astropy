# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Interfaces to the ERFA library, in particular for leap seconds.
"""
from datetime import datetime, timedelta
from warnings import warn

import numpy as np

from astropy.utils.decorators import classproperty
from astropy.utils.exceptions import AstropyWarning

from .ufunc import get_leap_seconds, set_leap_seconds, dt_eraLEAPSECOND


class leap_seconds:
    """Leap second management.

    This singleton class allows access to ERFA's leap second table,
    using the methods 'get', 'set', and 'update'.

    One can also check expiration with 'expires' and 'expired'.

    Note that usage of the class is similar to a ``ScienceState`` class,
    but it cannot be used as a context manager.
    """
    _expires = None
    """Explicit expiration date inferred from leap-second table."""
    _expiration_days = 180
    """Number of days beyond last leap second at which table expires."""

    def __init__(self):
        raise RuntimeError("This class is a singleton.  Do not instantiate.")

    @classmethod
    def get(cls):
        """Get the current leap-second table used internally."""
        return get_leap_seconds()

    @classmethod
    def validate(cls, table):
        """Validate a leap-second table.

        Raises
        ------
        ValueError
            If the leap seconds in the table are not on the 1st of January or
            July, or if the sorted TAI-UTC do not increase in increments of 1.
        """
        try:
            day = table['day']
        except Exception:
            day = 1

        # Take care of astropy Table.
        if hasattr(table, '__array__'):
            table = table.__array__()[list(dt_eraLEAPSECOND.names)]

        table = np.array(table, dtype=dt_eraLEAPSECOND, copy=False,
                         ndmin=1)

        # Simple sanity checks.
        if table.ndim > 1:
            raise ValueError("can only pass in one-dimensional tables.")

        if not np.all(((day == 1) &
                       (table['month'] == 1) | (table['month'] == 7)) |
                      (table['year'] < 1972)):
            raise ValueError("leap seconds inferred that are not on "
                             "1st of January or 1st of July.")

        if np.any((table['year'][:-1] > 1970) &
                  (np.diff(table['tai_utc']) != 1)):
            raise ValueError("jump in TAI-UTC by something else than one.")

        return table

    @classmethod
    def set(cls, table=None):
        """Set the ERFA leap second table.

        Note that it is generally safer to update the leap-second table than
        to set it directly, since most tables do not have the pre-1970 changes
        in TAI-UTC that are part of the built-in ERFA table.

        Parameters
        ----------
        table : array-like or `None`
            Leap-second table that should at least hold columns of 'year',
            'month', and 'tai_utc'.  Only simple validation is done before it
            is being used, so care need to be taken that entries are correct.
            If `None`, reset the ERFA table to its built-in values.

        Raises
        ------
        ValueError
            If the leap seconds in the table are not on the 1st of January or
            July, or if the sorted TAI-UTC do not increase in increments of 1.
        """
        set_leap_seconds(None if table is None else cls.validate(table))
        cls._expires = getattr(table, 'expires', None)

    @classproperty
    def expires(cls):
        """The expiration date of the current ERFA table.

        This is either a date inferred from the last table used to update or
        set the leap-second array, or a number of days beyond the last leap
        second.
        """
        if cls._expires is None:
            last = cls.get()[-1]
            return (datetime(last['year'], last['month'], 1) +
                    timedelta(cls._expiration_days))
        else:
            return cls._expires

    @classproperty
    def expired(cls):
        """Whether the leap second table is valid beyond the present."""
        return cls.expires < datetime.now()

    @classmethod
    def update(cls, table):
        """Add any leap seconds not already present to the ERFA table.

        This method matches leap seconds with those present in the ERFA table,
        and extends the latter as necessary.

        If the ERFA leap seconds file was corrupted, it will be reset.

        If the table is corrupted, the ERFA file will be unchanged.

        Parameters
        ----------
        table : array-like or `~astropy.utils.iers.LeapSeconds`
            Array or table with TAI-UTC from leap seconds.  Should have
            'year', 'month', and 'tai_utc' columns.

        Returns
        -------
        n_update : int
            Number of items updated.

        Raises
        ------
        ValueError
            If the leap seconds in the table are not on the 1st of January or
            July, or if the sorted TAI-UTC do not increase in increments of 1.
        """
        expires = getattr(table, 'expires', None)
        table = cls.validate(table)

        # Get erfa table and check it is OK; if not, reset it.
        try:
            erfa_ls = cls.validate(cls.get())
        except Exception:
            cls.set()
            erfa_ls = cls.get()

        # Create the combined array and use it (validating the combination).
        ls = np.union1d(erfa_ls, table)
        cls.set(ls)

        # If the update table has an expiration beyond that inferred from
        # the new leap second second array, use it (but, now that the new
        # array is set, do not allow exceptions due to misformed expires).
        try:
            if expires is not None and expires > cls.expires:
                cls._expires = expires

        except Exception as exc:
            warn("table 'expires' attribute ignored as comparing it "
                 "with a datetime raised an error:\n" + str(exc),
                 AstropyWarning)

        return len(ls) - len(erfa_ls)
