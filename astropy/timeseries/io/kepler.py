# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import numpy as np

from astropy.io import registry, fits
from astropy.table import Table
from astropy.time import Time, TimeDelta

from astropy.timeseries.sampled import TimeSeries

__all__ = ["kepler_read"]


def kepler_read(filename):
    """
    This serves as the FITS reader for KEPLER or TESS files within astropy-timeseries.

    This allows reading a supported FITS file using syntax such as::
    >>> from astropy.timeseries.io import kepler
    >>> timeseries = kepler.kepler_reader('<name of fits file>')  # doctest: +SKIP

    Parameters
    ----------
    filename: `str`, `pathlib.Path`
        File to load.

    Returns
    -------
    `astropy.timeseries.sampled.TimeSeries`
        Data converted into a TimeSeries.

    """
    hdulist = fits.open(filename)
    telescop = hdulist[0].header['telescop'].lower()
    if telescop not in ["kepler", "tess"]:
        raise NotImplementedError("{} is not implemented, only KEPLER or TESS are "
                                  "supported through this reader".format(hdulist[0].header['telescop']))

    # Get the lightcurve HDU
    if telescop == 'tess':
        hdu = hdulist['LIGHTCURVE']
    if telescop == 'kepler':
        hdu = hdulist[1]

    if hdu.header['EXTVER'] > 1:
        raise NotImplementedError("Support for {0} v{1} files not yet "
                                  "implemented".format(hdu.header['TELESCOP'], hdu.header['EXTVER']))

    # Check time scale
    if hdu.header['TIMESYS'] != 'TDB':
        raise NotImplementedError("Support for {0} time scale not yet "
                                  "implemented in {1} reader".format(hdu.header['TIMESYS'], hdu.header['TELESCOP']))

    # TODO: I wonder if we need these checks.
    if telescop == 'kepler':
        if "VERSION" in hdu.header.keys() and hdu.header["VERSION"] != '2.0':
            warnings.warn("This is not a EVEREST pipeline version 2 KEPLER file.")
        else:
            warnings.warn("This seems to be a non EVEREST reduced KEPLER file.")

    tab = Table.read(hdu, format='fits')

    # Some KEPLER files have a T column instead of TIME.
    if "T" in tab.colnames:
        tab.rename_column("T", "TIME")

    for colname in tab.colnames:
        # Fix units
        if tab[colname].unit == 'e-/s':
            tab[colname].unit = 'electron/s'
        if tab[colname].unit == 'pixels':
            tab[colname].unit = 'pixel'

        # Rename columns to lowercase
        tab.rename_column(colname, colname.lower())

    # Filter out NaN rows
    # TODO: Maybe change them to 0 instead?
    nans = np.isnan(tab['time'].data)
    if np.any(nans):
        warnings.warn('Ignoring {0} rows with NaN times'.format(np.sum(nans)))
    tab = tab[~nans]

    # Time column is dependent on source and we correct it here
    reference_date = Time(hdu.header['BJDREFI'], hdu.header['BJDREFF'],
                          scale=hdu.header['TIMESYS'].lower(), format='jd')
    time = reference_date + TimeDelta(tab['time'].data)
    time.format = 'isot'

    # Remove original time column
    tab.remove_column('time')

    return TimeSeries(time=time, data=tab)


registry.register_reader('kepler.fits', TimeSeries, kepler_read)
registry.register_reader('tess.fits', TimeSeries, kepler_read)
