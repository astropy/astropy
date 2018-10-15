import warnings
import numpy as np

from ...io import registry, fits
from ...table import Table
from ...time import Time, TimeDelta
from ... import units as u

from ..sampled import SampledTimeSeries

__all__ = ['kepler_fits_reader']


def tess_fits_reader(filename):

    # Open FITS file
    hdulist = fits.open(filename)

    # Get lightcurve HDU and check version
    hdu = hdulist['LIGHTCURVE']
    if hdu.header['EXTVER'] > 1:
        raise NotImplementedError("Support for TESS v{0} files not yet implemented".format(hdu.header['EXTVER']))

    # Get reference date for times and check time scale
    if hdu.header['TIMESYS'] != 'TDB':
        raise NotImplementedError("Support for {0} time scale not yet implemented in TESS reader".format(hdu.header['TIMESYS']))
    reference_date = Time(hdu.header['BJDREFI'], hdu.header['BJDREFF'], scale='tdb', format='jd')

    # Parse table with regular FITS reader
    tab = Table.read(hdu, format='fits')

    for colname in tab.colnames:

        # Fix units
        if tab[colname].unit == 'e-/s':
            tab[colname].unit = 'electron/s'

        # Rename columns to lowercase
        tab.rename_column(colname, colname.lower())

    # Filter out NaN rows
    nans = np.isnan(tab['time'].data)
    if np.any(nans):
        warnings.warn('Ignoring {0} rows with NaN times'.format(np.sum(nans)))
    tab = tab[~nans]

    # Compute Time object
    time = reference_date + TimeDelta(tab['time'].data)

    # Remove original time column
    tab.remove_column('time')

    # Create time series
    ts = SampledTimeSeries(time=time, data=tab)

    return ts


registry.register_reader('tess.fits', SampledTimeSeries, tess_fits_reader)
