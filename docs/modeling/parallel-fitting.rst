.. _parallel-fitting:

Fitting models in parallel with N-dimensional data
**************************************************

In some cases, you may want to fit a model many times to data. For example, you
may have a spectral cube (with two celestial axes and one spectral axis) and you
want to fit a 1D model (which could be either a simple Gaussian model or a
complex compound model with multiple lines and a continuum) to each individual
spectrum in the cube. Alternatively, you may have a cube with two celestial
axes, one spectral axis, and one time axis, and you want to fit a 2D model to
each 2D celestial plane in the cube. Provided each model fit can be treated as
independent, there are significant performance benefits to carrying out these
model fits in parallel.

The :func:`~astropy.modeling.fitting.parallel_fit_dask` function is ideally
suited to these use cases. It makes it simple to set up fitting of M-dimensional
models to N-dimensional datasets and leverages the power of the `dask
<https://www.dask.org/>`_ package to efficiently parallelize the problem,
running it either on multiple processes of a single machine or in a distributed
environment. You do not need to know how to use dask in order to use this function,
but you will need to make sure you have `dask <https://www.dask.org/>`_
installed.

Note that the approach here is different from *model sets* which are described
in :ref:`example-fitting-model-sets`, which are a way of fitting a linear model
with a vector of parameters to a data array, as in that specific case the
fitting can be truly vectorized, and will likely not benefit from the approach
described here.

Getting started
===============

To demonstrate the use of this function, we will work through a simple
example of fitting a 1D model to a small spectral cube (if you are
interested in accessing the file, you can find it at
:download:`l1448_13co.fits <http://www.astropy.org/astropy-data/l1448/l1448_13co.fits>`,
but the code below will automatically download it).

.. The following block is to make sure 'data' and 'wcs' are defined if we are not running with --remote-data

.. plot::
   :context: close-figs
   :nofigs:

    >>> import numpy as np
    >>> from astropy.wcs import WCS
    >>> wcs = WCS(naxis=3)
    >>> wcs.wcs.ctype =  ['RA---SFL', 'DEC--SFL', 'VOPT']
    >>> wcs.wcs.crval = [57.66, 0., -9959.44378305]
    >>> wcs.wcs.crpix =  [-799.0, -4741.913, -187.0]
    >>> wcs.wcs.cdelt = [-0.006388889, 0.006388889, 66.42361]
    >>> wcs.wcs.cunit = ['deg', 'deg', 'm s-1']
    >>> wcs._naxis = [105, 105, 53]
    >>> wcs.wcs.set()
    >>> data = np.broadcast_to(np.exp(-(np.arange(53) - 25)**2 / 6 ** 2).reshape((53, 1, 1)), (53, 105, 105))

We start by downloading the cube and extracting the data and WCS:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> from astropy.wcs import WCS
    >>> from astropy.io import fits
    >>> from astropy.utils.data import get_pkg_data_filename

    >>> filename = get_pkg_data_filename('l1448/l1448_13co.fits')  # doctest: +REMOTE_DATA
    >>> with fits.open(filename) as hdulist:
    ...     data = hdulist[0].data
    ...     wcs = WCS(hdulist[0].header)  # doctest: +REMOTE_DATA

We extract a sub-cube spatially for the purpose of demonstration:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> data = data[:, 25:75, 35:85]
    >>> wcs = wcs[:, 25:75, 35:85]

This is a cube of a star-formation region traced by the 13CO line. We can look
at one of the channels:

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(subplot_kw=dict(projection=wcs, slices=('x', 'y', 20)))
    >>> ax.imshow(data[20, :, :])  # doctest: +IGNORE_OUTPUT

We can also extract a spectrum for one of the celestial positions:

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> fig, ax = plt.subplots(subplot_kw=dict(projection=wcs, slices=(5, 5, 'x')))
    >>> ax.plot(data[:, 5, 5])  # doctest: +IGNORE_OUTPUT

We now set up a model to fit this; we will use a simple Gaussian model,
with some reasonable initial guesses for the parameters:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> from astropy import units as u
    >>> from astropy.modeling.models import Gaussian1D
    >>> model = Gaussian1D(amplitude=1 * u.one, mean=4000 * u.m / u.s, stddev=500 * u.m / u.s)

The data does not have any units in this case, so we use ``u.one`` as
the unit, which indicates it is dimensionless.

Before fitting this to all spectra in the cube, it’s a good idea to test
the model with at least one of the spectra manually. To do this, we need to extract the x-axis of the spectra:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> import numpy as np
    >>> x = wcs.pixel_to_world(0, 0, np.arange(data.shape[0]))[1]
    >>> x
    <SpectralCoord
       (target: <ICRS Coordinate: (ra, dec, distance) in (deg, deg, kpc)
                    (57.66, 0., 1000.)
                 (pm_ra_cosdec, pm_dec, radial_velocity) in (mas / yr, mas / yr, km / s)
                    (0., 0., 0.)>)
      [2528.19489695, 2594.61850695, 2661.04211695, 2727.46572695,
       2793.88933695, 2860.31294695, 2926.73655695, 2993.16016695,
       ...
       5716.52817695, 5782.95178695, 5849.37539695, 5915.79900695,
       5982.22261695] m / s>

We can now carry out the fit:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> from astropy.modeling.fitting import TRFLSQFitter
    >>> fitter = TRFLSQFitter()
    >>> model_fit_single = fitter(model, x, data[:, 5, 5])

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, data[:, 5, 5], '.', label='data')  # doctest: +IGNORE_OUTPUT
    >>> ax.plot(x, model(x), label='initial model')  # doctest: +IGNORE_OUTPUT
    >>> ax.plot(x, model_fit_single(x), label='fitted model')  # doctest: +IGNORE_OUTPUT
    >>> ax.legend()  # doctest: +IGNORE_OUTPUT

The model seems to work! We can now use the
:func:`~astropy.modeling.fitting.parallel_fit_dask` function
to fit all spectra in the cube:

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> from astropy.modeling.fitting import parallel_fit_dask
    >>> model_fit = parallel_fit_dask(model=model,
    ...                               fitter=fitter,
    ...                               data=data,
    ...                               world=wcs,
    ...                               fitting_axes=0,
    ...                               data_unit=u.one,
    ...                               scheduler='synchronous')

The arguments in this case are as follows:

*  ``model=`` is the initial model. While in our case the initial
   parameters were specified as scalars, it is possible to pass in a
   model that has array parameters if you want to have different initial
   parameters as a function of location in the dataset.
*  ``fitter=`` is the fitter instance.
*  ``data=`` is the N-dimensional dataset, in our case the 3D spectral
   cube.
*  ``world=`` provides information about the world coordinates for the
   fit, for example the spectral coordinates for a spectrum. This can be
   specified in different ways, but above we have chosen to pass in the
   WCS object for the dataset, from which the spectral axis coordinates
   will be extracted.
*  ``fitting_axes=`` specifies which axis or axes include the data to
   fit. In our example, we are fitting the spectra,
   which in NumPy notation is the first axis in the cube, so we specify
   ``fitting_axes=0``.
*  ``data_unit=`` specifies the unit to use for the data. In our case,
   the data has no unit, but because we are using units for the spectral
   axis, we need to specify ``u.one`` here.

We can now take a look at the parameter maps:

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> fig, axs = plt.subplots(figsize=(10, 5), ncols=3)
    >>> ax1 = axs[0]
    >>> ax1.set_title('Amplitude')  # doctest: +IGNORE_OUTPUT
    >>> ax1.imshow(model_fit.amplitude.value, vmin=0, vmax=5, origin='lower')  # doctest: +IGNORE_OUTPUT
    >>> ax2 = axs[1]
    >>> ax2.set_title('Mean')  # doctest: +IGNORE_OUTPUT
    >>> ax2.imshow(model_fit.mean.value, vmin=2500, vmax=6000, origin='lower')  # doctest: +IGNORE_OUTPUT
    >>> ax3 = axs[2]
    >>> ax3.set_title('Standard deviation')  # doctest: +IGNORE_OUTPUT
    >>> ax3.imshow(model_fit.stddev.value, vmin=0, vmax=2000, origin='lower')  # doctest: +IGNORE_OUTPUT

There are a number of pixels that appear to have issues. Inspecting the
histogram of means, we can see that a lot of values are not at all in
the spectral range we are fitting:

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> fig, ax = plt.subplots()
    >>> ax.hist(model_fit.mean.value.ravel(), bins=100)  # doctest: +IGNORE_OUTPUT
    >>> ax.set(yscale='log', xlabel='mean', ylabel='number')  # doctest: +IGNORE_OUTPUT

We can set the bounds on the mean and try the fit again

.. plot::
   :context: close-figs
   :include-source:
   :nofigs:

    >>> model.mean.bounds = (3000, 6000) * u.km / u.s
    >>> model_fit = parallel_fit_dask(model=model,
    ...                               fitter=fitter,
    ...                               data=data,
    ...                               world=wcs,
    ...                               fitting_axes=0,
    ...                               data_unit=u.one,
    ...                               scheduler='synchronous')

and we can visualize the results:

.. plot::
   :context: close-figs
   :include-source:
   :align: center

    >>> fig, axs = plt.subplots(figsize=(10, 5), ncols=3)
    >>> ax1 = axs[0]
    >>> ax1.set_title('Amplitude')  # doctest: +IGNORE_OUTPUT
    >>> ax1.imshow(model_fit.amplitude.value, vmin=0, vmax=5, origin='lower')  # doctest: +IGNORE_OUTPUT
    >>> ax2 = axs[1]
    >>> ax2.set_title('Mean')  # doctest: +IGNORE_OUTPUT
    >>> ax2.imshow(model_fit.mean.value, vmin=2500, vmax=6000, origin='lower')  # doctest: +IGNORE_OUTPUT
    >>> ax3 = axs[2]
    >>> ax3.set_title('Standard deviation')  # doctest: +IGNORE_OUTPUT
    >>> ax3.imshow(model_fit.stddev.value, vmin=0, vmax=2000, origin='lower')  # doctest: +IGNORE_OUTPUT

The amplitude map no longer contains any problematic pixels.

World input
===========

The example above demonstrated that it is possible to pass in a
:class:`astropy.wcs.WCS` object to the ``world=`` argument in order to determine
the world coordinates for the fit (e.g. the spectral axis values for a spectral
fit). It is also possible to pass in a tuple of arrays - if you do this, the
tuple should have one item per fitting axis. It is most efficient to pass in a
tuple of 1D arrays, but if the world coordinates vary over the axes being
iterated over, you can also pass in a tuple of N-d arrays, giving the
coordinates of each individual pixel (it is also possible to pass in arrays that
are not 1D but also not fully N-d as long as they can be broadcasted to the data
shape).

Multiprocessing
===============

By default, :func:`~astropy.modeling.fitting.parallel_fit_dask` will make use
of multi-processing to parallelize the fitting. If you write a script to
carry out the fitting, you will likely need to move your code inside a::

    if __name__ == "__main__":

        ...

clause as otherwise Python will execute the whole code in the script many times,
and potentially recursively, rather than just parallelizing the fitting.

Performance
===========

The :func:`~astropy.modeling.fitting.parallel_fit_dask` function splits the data
into chunks, each of which is then sent to a different process. The size of
these chunks is critical to obtaining good performance. If we split the data
into one chunk per fit, the process would be inefficient due to significant
overhead from inter-process communication. Conversely, if we split the data into
fewer chunks than there are available processes, we will not utilize all the
available computational power. If we split the data into slightly more chunks
than there are processes, inefficiencies can arise as well. For example,
splitting the data into five chunks with four available processes means the four
processes will first fit four chunks, and then a single process will be held up
fitting the remaining chunk. Therefore, it is important to carefully consider
how the data is split.

To control the splitting of the data, use the ``chunk_n_max=`` keyword argument.
This determines how many individual fits will be carried out in each chunk. For
example, when fitting a model to individual spectra in a spectral cube, setting
``chunk_n_max=100`` means each chunk will contain 100 spectra. As a general
guide, you will likely want to set this to be roughly the number of fits to be
carried out in the data divided by several times the number of available
processes. For example, if you need to fit 100,000 spectra and have 8 processes
available, setting ``chunk_n_max=1000`` would be reasonable. This configuration
would break the data into 100 chunks, meaning each process will need to handle
approximately a dozen chunks. Additionally, fitting 1,000 spectra per chunk will
take enough time to avoid being dominated by communication overhead.

The default value for ``chunk_n_max`` is 500.

Diagnostics
===========

One of the challenges of fitting a model many different times is understanding
what went wrong when issues arise. By default, if a fit fails with a warning or
an exception, the parameters for that fit will be set to NaN, and no warning or
exception will be shown to the user. However, it can be helpful to have more
information, such as the specific error or exception that occurred.

You can control this by setting the ``diagnostics=`` argument. This allows you
to choose whether to output information about:

* Failed fits with errors (``diagnostics='error'``),
* Fits with errors or warnings (``diagnostics='error+warn'``), or
* All fits (``diagnostics='all'``).

If the ``diagnostics`` option is specified, you will also need to specify
``diagnostics_path``, which should be the path to a folder that will contain all
the output. Each fit that needs to be output will be assigned a sub-folder named
after the indices along the axes of the data (excluding the fitting axes). The
output will include (if appropriate):

* ``error.log``, containing details of any exceptions that occurred
* ``warn.log``, containing any warnings

You may also want to automatically create a plot of the fit, inspect the data
being fit, or examine the model. To do this, you can pass a function to
``diagnostics_callable``. See :func:`~astropy.modeling.fitting.parallel_fit_dask`
for more information about the arguments this function should accept.

Schedulers
==========

By default, :func:`~astropy.modeling.fitting.parallel_fit_dask` will make use of
the ``'processes'`` scheduler, which means that multiple processes on your local
machine can be used. You can override the scheduler being used with the
``scheduler=`` keyword argument. You can either set this to the name of a
scheduler (such as ``'synchronous'``), or you can set it to ``'default'`` in order
to make use of whatever is the currently active dask scheduler, which allows
you for example to set up a `dask.distributed
<https://distributed.dask.org/en/stable/>`_ scheduler.
