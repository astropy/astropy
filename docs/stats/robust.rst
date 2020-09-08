.. _stats-robust:

*****************************
Robust Statistical Estimators
*****************************

Robust statistics provides reliable estimates of basic statistics for complex
distributions. The statistics package includes several robust statistical
functions that are commonly used in astronomy. This includes methods for
rejecting outliers as well as statistical description of the underlying
distributions.

In addition to the functions mentioned here, models can be fit with outlier
rejection using :func:`~astropy.modeling.fitting.FittingWithOutlierRemoval`.

Sigma Clipping
==============

Sigma clipping provides a fast method for identifying outliers in a
distribution. For a distribution of points, a center and a standard
deviation are calculated. Values which are less or more than a
specified number of standard deviations from a center value are
rejected. The process can be iterated to further reject outliers.

The `astropy.stats` package provides both a functional and
object-oriented interface for sigma clipping. The function is called
:func:`~astropy.stats.sigma_clip` and the class is called
:class:`~astropy.stats.SigmaClip`. By default, they both return a
masked array where the rejected points are masked.

Examples
--------

..
  EXAMPLE START
  Functional Sigma Clipping with astropy.stats.sigma_clip

We can start by generating some data that has a mean of 0 and standard
deviation of 0.2, but with outliers:

.. doctest-requires:: scipy

     >>> import numpy as np
     >>> import scipy.stats as stats
     >>> np.random.seed(0)
     >>> x = np.arange(200)
     >>> y = np.zeros(200)
     >>> c = stats.bernoulli.rvs(0.35, size=x.shape)
     >>> y += (np.random.normal(0., 0.2, x.shape) +
     ...       c*np.random.normal(3.0, 5.0, x.shape))

Now we can use :func:`~astropy.stats.sigma_clip` to perform sigma
clipping on the data:

.. doctest-requires:: scipy

     >>> from astropy.stats import sigma_clip
     >>> filtered_data = sigma_clip(y, sigma=3, maxiters=10)

The output masked array then can be used to calculate statistics on
the data, fit models to the data, or otherwise explore the data.

..
  EXAMPLE END

..
  EXAMPLE START
  Object-Oriented Sigma Clipping with the astropy.stats.SigmaClip Class

To perform the same sigma clipping with the
:class:`~astropy.stats.SigmaClip` class:

.. doctest-requires:: scipy

     >>> from astropy.stats import SigmaClip
     >>> sigclip = SigmaClip(sigma=3, maxiters=10)
     >>> print(sigclip)  # doctest: +SKIP
     <SigmaClip>
        sigma: 3
        sigma_lower: None
        sigma_upper: None
        maxiters: 10
        cenfunc: <function median at 0x108dbde18>
        stdfunc: <function std at 0x103ab52f0>
     >>> filtered_data = sigclip(y)

Note that once the ``sigclip`` instance is defined above, it can be
applied to other data using the same already defined sigma-clipping
parameters.

..
  EXAMPLE END

For basic statistics, :func:`~astropy.stats.sigma_clipped_stats` is a
convenience function to calculate the sigma-clipped mean, median, and
standard deviation of an array. As can be seen, rejecting the
outliers returns accurate values for the underlying distribution.

..
  EXAMPLE START
  Calculating the Sigma-Clipped Mean, Median, and Standard Deviation of an Array

To use :func:`~astropy.stats.sigma_clipped_stats` for sigma-clipped statistics
calculation:

.. doctest-requires:: scipy

     >>> from astropy.stats import sigma_clipped_stats
     >>> y.mean(), np.median(y), y.std()  # doctest: +FLOAT_CMP
     (0.86586417693378226, 0.03265864495523732, 3.2913811977676444)
     >>> sigma_clipped_stats(y, sigma=3, maxiters=10)  # doctest: +FLOAT_CMP
     (-0.0020337793767186197, -0.023632809025713953, 0.19514652532636906)

:func:`~astropy.stats.sigma_clip` and
:class:`~astropy.stats.SigmaClip` can be combined with other robust
statistics to provide improved outlier rejection as well.

.. plot::
    :include-source:

    import numpy as np
    import scipy.stats as stats
    from matplotlib import pyplot as plt
    from astropy.stats import sigma_clip, mad_std

    # Generate fake data that has a mean of 0 and standard deviation of 0.2 with outliers
    np.random.seed(0)
    x = np.arange(200)
    y = np.zeros(200)
    c = stats.bernoulli.rvs(0.35, size=x.shape)
    y += (np.random.normal(0., 0.2, x.shape) +
          c*np.random.normal(3.0, 5.0, x.shape))

    filtered_data = sigma_clip(y, sigma=3, maxiters=1, stdfunc=mad_std)

    # plot the original and rejected data
    plt.figure(figsize=(8,5))
    plt.plot(x, y, '+', color='#1f77b4', label="original data")
    plt.plot(x[filtered_data.mask], y[filtered_data.mask], 'x',
             color='#d62728', label="rejected data")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc=2, numpoints=1)

.. automodapi:: astropy.stats.sigma_clipping

..
  EXAMPLE END

Median Absolute Deviation
=========================

The median absolute deviation (MAD) is a measure of the spread of a
distribution and is defined as ``median(abs(a - median(a)))``. The
MAD can be calculated using `~astropy.stats.median_absolute_deviation`. For a
normal distribution, the MAD is related to the standard deviation by a factor
of 1.4826, and a convenience function, `~astropy.stats.mad_std`, is
available to apply the conversion.

.. note::

   A function can be supplied to the
   `~astropy.stats.median_absolute_deviation` to specify the median
   function to be used in the calculation. Depending on the version
   of NumPy and whether the array is masked or contains irregular
   values, significant performance increases can be had by
   preselecting the median function. If the median function is not
   specified, `~astropy.stats.median_absolute_deviation` will attempt
   to select the most relevant function according to the input data.


Biweight Estimators
===================

A set of functions are included in the `astropy.stats` package that use the
biweight formalism. These functions have long been used in astronomy,
particularly to calculate the velocity dispersion of galaxy clusters [1]_. The
following set of tasks are available for biweight measurements:

.. automodapi:: astropy.stats.biweight


References
----------

.. [1] Beers, Flynn, and Gebhardt (1990; AJ 100, 32) (https://ui.adsabs.harvard.edu/abs/1990AJ....100...32B)
