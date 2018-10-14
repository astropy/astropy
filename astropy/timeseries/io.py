from ..io import registry
from ..table import Table
from ..time import Time

from .timeseries import SampledTimeSeries


def kepler_fits_reader(filename):

    # Parse Kepler FITS file with regular FITS reader
    tab = Table.read(filename, format='fits')

    # Rename columns to lowercase
    for colname in tab.colnames:
        tab.rename_column(colname, colname.lower())

    # Compute Time object
    time = Time(tab['time'].data + 2454833, scale='tcb', format='jd')

    # Remove original time column
    tab.remove_column('time')

    # Create time series
    ts = SampledTimeSeries(time=time, data=tab)

    return ts


registry.register_reader('kepler.fits', SampledTimeSeries, kepler_fits_reader)
