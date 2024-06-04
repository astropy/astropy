.. _timeseries-masking:

Masking Values in Time Series
*****************************

.. warning:: Note that masking does not yet work for columns that have units.

Masking values is done in the same way as for |Table| objects (see
:ref:`masking_and_missing_values`). The most convenient way to use masking is to
initialize a |TimeSeries| object using the ``masked=True`` option.

Example
-------

.. EXAMPLE START: Masking Values in TimeSeries Objects

We start by initializing a |TimeSeries| object with ``masked=True``::

    >>> from astropy import units as u
    >>> from astropy.timeseries import TimeSeries
    >>> ts = TimeSeries(time_start='2016-03-22T12:30:31',
    ...                 time_delta=3 * u.s,
    ...                 n_samples=5, masked=True)

We can now add some data to our time series::

    >>> ts['flux'] = [1., -2., 5., -1., 4.]

As you can see, some of the values are negative. We can mask these using::

    >>> ts['flux'].mask = ts['flux'] < 0
    >>> ts
    <TimeSeries masked=True length=5>
              time            flux
              Time          float64
    ----------------------- -------
    2016-03-22T12:30:31.000     1.0
    2016-03-22T12:30:34.000      --
    2016-03-22T12:30:37.000     5.0
    2016-03-22T12:30:40.000      --
    2016-03-22T12:30:43.000     4.0

We can also access the mask values::

    >>> ts['flux'].mask
    array([False,  True, False,  True, False]...)

Masks are column-based, so masking a single cell does not mask the whole row.
Having masked cells then allows functions that normally understand masked values
and operate on columns to ignore the masked entries::

    >>> import numpy as np
    >>> np.min(ts['flux'])
    np.float64(1.0)
    >>> np.ma.median(ts['flux'])
    np.float64(4.0)

.. EXAMPLE END
