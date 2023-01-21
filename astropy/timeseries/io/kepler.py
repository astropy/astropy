# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import numpy as np

from astropy.io import fits, registry
from astropy.table import MaskedColumn, Table
from astropy.time import Time, TimeDelta
from astropy.timeseries.sampled import TimeSeries

__all__ = ["kepler_fits_reader"]


def kepler_fits_reader(filename, unit_parse_strict="warn"):
    """
    This serves as the FITS reader for KEPLER or TESS files within
    astropy-timeseries.

    This function should generally not be called directly, and instead this
    time series reader should be accessed with the
    :meth:`~astropy.timeseries.TimeSeries.read` method::

        >>> from astropy.timeseries import TimeSeries
        >>> ts = TimeSeries.read('kplr33122.fits', format='kepler.fits')  # doctest: +SKIP

    Parameters
    ----------
    filename : `str` or `pathlib.Path`
        File to load.
    unit_parse_strict : str, optional
        Behaviour when encountering invalid column units in the FITS header.
        Default is "warn", which will emit a ``UnitsWarning`` and create a
        :class:`~astropy.units.core.UnrecognizedUnit`.
        Values are the ones allowed by the ``parse_strict`` argument of
        :class:`~astropy.units.core.Unit`: ``raise``, ``warn`` and ``silent``.

    Returns
    -------
    ts : `~astropy.timeseries.TimeSeries`
        Data converted into a TimeSeries.
    """
    hdulist = fits.open(filename)
    # Get the lightcurve HDU
    telescope = hdulist[0].header["telescop"].lower()

    if telescope == "tess":
        hdu = hdulist["LIGHTCURVE"]
    elif telescope == "kepler":
        hdu = hdulist[1]
    else:
        raise NotImplementedError(
            f"{hdulist[0].header['telescop']} is not implemented, only KEPLER or TESS"
            " are supported through this reader"
        )

    if hdu.header["EXTVER"] > 1:
        raise NotImplementedError(
            f"Support for {hdu.header['TELESCOP']} v{hdu.header['EXTVER']} files not"
            " yet implemented"
        )

    # Check time scale
    if hdu.header["TIMESYS"] != "TDB":
        raise NotImplementedError(
            f"Support for {hdu.header['TIMESYS']} time scale not yet implemented in"
            f" {hdu.header['TELESCOP']} reader"
        )

    tab = Table.read(hdu, format="fits", unit_parse_strict=unit_parse_strict)

    # Some KEPLER files have a T column instead of TIME.
    if "T" in tab.colnames:
        tab.rename_column("T", "TIME")

    for colname in tab.colnames:
        unit = tab[colname].unit
        # Make masks nan for any column which will turn into a Quantity
        # later. TODO: remove once we support Masked Quantities properly?
        if unit and isinstance(tab[colname], MaskedColumn):
            tab[colname] = tab[colname].filled(np.nan)
        # Fix units
        if unit == "e-/s":
            tab[colname].unit = "electron/s"
        if unit == "pixels":
            tab[colname].unit = "pixel"

        # Rename columns to lowercase
        tab.rename_column(colname, colname.lower())

    # Filter out NaN rows
    nans = np.isnan(tab["time"].data)
    if np.any(nans):
        warnings.warn(f"Ignoring {np.sum(nans)} rows with NaN times")
    tab = tab[~nans]

    # Time column is dependent on source and we correct it here
    reference_date = Time(
        hdu.header["BJDREFI"],
        hdu.header["BJDREFF"],
        scale=hdu.header["TIMESYS"].lower(),
        format="jd",
    )
    time = reference_date + TimeDelta(tab["time"].data, format="jd")
    time.format = "isot"

    # Remove original time column
    tab.remove_column("time")

    hdulist.close()

    return TimeSeries(time=time, data=tab)


registry.register_reader("kepler.fits", TimeSeries, kepler_fits_reader)
registry.register_reader("tess.fits", TimeSeries, kepler_fits_reader)
