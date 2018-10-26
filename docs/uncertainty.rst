.. |quantity| replace:: :class:`~astropy.units.Quantity`
.. |distribution| replace:: :class:`~astropy.uncertainty.Distribution`

.. _astropy-uncertainty:

*******************************************************
Uncertainties and Distributions (`astropy.uncertainty`)
*******************************************************

.. warning::

    `astropy.uncertainty` is currently a work-in-progress, and thus it is quite
    possible there will be API changes in later versions of Astropy. If you have
    specific ideas for how it might be improved, please  let us know on the
    `astropy-dev mailing list`_ or at http://feedback.astropy.org .

Introduction
============

In addition to |quantity|, astropy provides a |distribution| object to represent
statistical distributions in a form that acts as a drop-in replacement for
|quantity| or regulay Numpy arrays. Used in this manner, |distribution| provides
uncertainty propogation at the cost of additional computation.  It can also
more generally represent sampled distributions for e.g., Monte Carlo calculation
techniques.

The core object for this feature is the |distribution|.  Currently, all
such distributions are Monte Carlo sampled.  While this means each distribution
may take more memory, it allows arbitrarily complex operations to be performed
on distributions while maintaining their correlation structure. Some specific
well-behaved distributions (e.g., the Normal distribution) have
analytic forms which may eventually enable a more compact/efficient
representation, but this is not currently implemented.

Getting Started
---------------

To demonstrate a simple use case for distributions, consider the problem of
uncertainty propogation of normal distributions.  Assume there are two
measurements you wish to add.  We start with some initial imports/setup::

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy import uncertainty as unc
  >>> np.random.seed(12345)  # ensures reproducible example numbers

Now we create two |distribution| objects to represent our distributions::

  >>> a = unc.normal(1*u.kpc, std=30*u.pc)
  >>> b = unc.normal(2*u.kpc, std=40*u.pc)

For normal distributions, the centers should add as expected, and the standard
deviations add in quadrature.  We can check these results (to the limits of our
Monte Carlo sampling) trivially with |distribution| arithmetic and attributes::

  >>> c = a + b
  >>> c # doctest: +ELLIPSIS
  <QuantityDistribution [..., 3.06148029,...] kpc with n_samples=1000>
  >>> c.pdf_mean # doctest: +FLOAT_CMP
  <Quantity 3.0005627 kpc>
  >>> c.pdf_std.to(u.pc) # doctest: +FLOAT_CMP
  <Quantity 51.4783738 pc>

Indeed these are close to the expectations. While this may seem unnecessary for
the simple Gaussian case, for more complex distributions or airthmetic
operations where error analysis becomes untenable, |distribution| still powers
through:

  >>> d = unc.poisson(3*u.kpc)
  >>> e = unc.uniform_center_width(3*u.kpc, 800*u.pc)
  >>> f = (c * d * e) ** (1/3)
  >>> f.pdf_mean # doctest: +FLOAT_CMP
  <Quantity 2.83327524 kpc>
  >>> f.pdf_std # doctest: +FLOAT_CMP
  <Quantity 0.79948566 kpc>
  >>> f.hist(bins=50) # doctest: +SKIP

..plot::
  :align: center

  import numpy as np
  from astropy import units as u
  from astropy import uncertainty as unc
  np.random.seed(12345)
  a = unc.normal(1*u.kpc, std=30*u.pc)
  b = unc.normal(2*u.kpc, std=40*u.pc)
  c = a + b
  d = unc.poisson(3*u.kpc)
  e = unc.uniform_center_width(3*u.kpc, 800*u.pc)
  f = (c * d * e) ** (1/3)
  f.hist(bins=50)



Using `astropy.uncertainty`
===========================

Creating distributions
----------------------

The most direct way to create a distribution is to use an array or quantity
that carries the samples in the *last* dimension:

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy import uncertainty as unc
  >>> np.random.seed(123456)  # ensures "random" numbers match examples below

  >>> unc.Distribution(np.random.poisson(12, (1000)))  # doctest: +ELLIPSIS
  ndarrayDistribution([..., 12,...]) with n_samples=1000
  >>> pq = np.random.poisson([1, 5, 30, 400], (1000, 4)).T * u.kpc # note the transpose, required to get the sampling on the *last* axis
  >>> distr = unc.Distribution(pq)
  >>> distr # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
             [...],
             [...],
             [...]] kpc with n_samples=1000>
  >>> distr.unit
  Unit("kpc")

Note the distinction for these two distributions: the first is built from an
array and therfore does not have |quantity| attributes like ``unit``, while the
latter does.  This is reflected in how they interact with other objects - for
example the ``ndarrayDistribution`` will not combine with unitful |quantity|
objects.

For commonly-used distributions, helper functions exist to make creating them
easier. Below demonstrates several equivalent ways to create a normal/Gaussian
distribution:

  >>> centerq = [1, 5, 30, 400]*u.kpc
  >>> n_distr = unc.normal(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
  >>> n_distr = unc.normal(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
  >>> n_distr = unc.normal(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)
  >>> n_distr.distribution.shape
  (4, 1000)
  >>> unc.normal(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100).distribution.shape
  (4, 100)
  >>> unc.normal(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=20000).distribution.shape
  (4, 20000)


Additionally, Poisson and uniform |distribution| creation functions exist:

  >>> unc.poisson(centerq) # doctest: +FLOAT_CMP
  <QuantityDistribution [[  2.,   2.,   2., ...,   1.,   1.,   3.],
               [  7.,   3.,   2., ...,   1.,   3.,   5.],
               [ 32.,  27.,  30., ...,  31.,  39.,  29.],
               [425., 387., 379., ..., 404., 384., 377.]] kpc with n_samples=1000>
  >>> uwidth = [10, 20, 10, 55]*u.pc
  >>> unc.uniform_center_width(centerq, uwidth) # doctest: +FLOAT_CMP
  <QuantityDistribution [[  1.00264599,   1.00203114,   1.00228494, ...,   1.00280875,
                  1.00077583,   1.00213164],
               [  5.00861993,   4.99102886,   5.00254059, ...,   4.99831796,
                  5.00073192,   4.99709326],
               [ 29.99540657,  29.99856712,  30.00200227, ...,  29.99766923,
                 30.0014714 ,  29.99719035],
               [400.01991204, 399.98233216, 399.98231204, ..., 399.99071525,
                400.01797099, 400.01086062]] kpc with n_samples=1000>
  >>> unc.uniform(lower=centerq-uwidth/2,  upper=centerq+uwidth/2)  # doctest: +FLOAT_CMP
  <QuantityDistribution [[  0.99578122,   1.00056763,   0.99984519, ...,   0.99531055,
                  1.00396487,   1.00018797],
               [  5.0076216 ,   5.00672847,   5.00009831, ...,   5.00216818,
                  4.99087752,   5.00813196],
               [ 30.0031185 ,  29.99844025,  29.99764698, ...,  29.99581011,
                 30.00328488,  29.99903985],
               [399.97772381, 400.00168039, 399.98895828, ..., 399.97706371,
                400.00833078, 399.99777102]] kpc with n_samples=1000>

Users are free to create their own distribution classes following similar
patterns.


Using Distributions
-------------------

This object now acts much like a |quantity| for all but the non-sampled
dimension, but with additional statistical operations that work on the sampled
distributions:

  >>> distr.shape
  (4,)
  >>> distr.size
  4
  >>> distr.unit
  Unit("kpc")
  >>> distr.n_samples
  1000
  >>> distr.pdf_mean # doctest: +FLOAT_CMP
  <Quantity [  0.998,   5.017,  30.085, 400.345] kpc>
  >>> distr.pdf_std # doctest: +FLOAT_CMP
  <Quantity [ 0.97262326,  2.32222114,  5.47629208, 20.6328373 ] kpc>
  >>> distr.pdf_var # doctest: +FLOAT_CMP
  <Quantity [  0.945996,   5.392711,  29.989775, 425.713975] kpc2>
  >>> distr.pdf_median
  <Quantity [   1.,   5.,  30., 400.] kpc>
  >>> distr.pdf_mad  # Median absolute deviation # doctest: +FLOAT_CMP
  <Quantity [ 1.,  2.,  4., 14.] kpc>
  >>> distr.pdf_smad  # Median absolute deviation, rescaled to match std for normal # doctest: +FLOAT_CMP
  <Quantity [ 1.48260222,  2.96520444,  5.93040887, 20.75643106] kpc>
  >>> distr.pdf_percentiles([10, 50, 90])
  <Quantity [[  0. ,   2. ,  23. , 374. ],
             [  1. ,   5. ,  30. , 400. ],
             [  2. ,   8. ,  37.1, 427. ]] kpc>

If need be, the underlying array can then be accessed from the ``distribution``
attribute:

  >>> distr.distribution
  <Quantity [[  0.,   0.,   1., ...,   1.,   0.,   1.],
             [  7.,   3.,   4., ...,   3.,   2.,   5.],
             [ 27.,  32.,  35., ...,  37.,  21.,  40.],
             [421., 373., 389., ..., 405., 391., 369.]] kpc>
  >>> distr.distribution.shape
  (4, 1000)

A |quantity| distribution interact naturally with non-|distribution| quantities,
essentially assuming the |quantity| is a dirac delta distribution:

  >>> distrplus = distr + [2000,0,0,500]*u.pc
  >>> distrplus.pdf_median
  <Quantity [   3. ,   5. ,  30. , 400.5] kpc>
  >>> distrplus.pdf_var
  <Quantity [  0.945996,   5.392711,  29.989775, 425.713975] kpc2>


It also operates as expected with other distributions  (But see below for a
discussion of covariances):

  >>> another_distr = unc.Distribution((np.random.randn(1000,4)*[1000,.01 , 3000, 10] + [2000, 0, 0, 500]).T * u.pc)
  >>> combined_distr = distr + another_distr
  >>> combined_distr.pdf_median
  <Quantity [  3.01847755,   4.99999576,  29.60559788, 400.49176321] kpc>
  >>> combined_distr.pdf_var
  <Quantity [  1.8427705 ,   5.39271147,  39.5343726 , 425.71324244] kpc2>


Covariance in distributions
---------------------------

One of the main applications for distributions is unceratinty propogation, which
critically requires proper treatment of covariance. This comes naturally in the
Monte Carlo sampling approach used by the |distribution| class, as long as
proper care is taken with sampling error.

To start with a simple example, two un-correlated distributions should produce
an un-correlated joint distribution plot:

.. plot::
  :context: close-figs
  :include-source:
  :align: center

  >>> import numpy as np
  >>> np.random.seed(12345)  # produce repeatable plots
  >>> from astropy import units as u
  >>> from astropy import uncertainty as unc
  >>> from matplotlib import pyplot as plt # doctest: +SKIP
  >>> n1 = unc.normal(center=0., std=1, n_samples=10000)
  >>> n2 = unc.normal(center=0., std=2, n_samples=10000)
  >>> plt.scatter(n1.distribution, n2.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP

Indeed, the distributions are independent.  If we instead construct a covariant
pair of gaussians, it is immediately apparent:

.. plot::
  :context: close-figs
  :include-source:
  :align: center

  >>> ncov = np.random.multivariate_normal([0, 0], [[1, .5], [.5, 2]], size=10000)
  >>> n1 = unc.Distribution(ncov[:, 0])
  >>> n2 = unc.Distribution(ncov[:, 1])
  >>> plt.scatter(n1.distribution, n2.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP


Most importantly, the proper correlated structure is preserved or generated as
expected by appropriate arithmetic operations. For example, ratios of
uncorrelated normal distribution gain covariances if the axes are not
independent, as in this simulation of iron, hydrogen, and oxygen abundances in
a hypothetical collection of stars:

.. plot::
  :context: close-figs
  :include-source:
  :align: center

  >>> fe_abund = unc.normal(center=-2, std=.25, n_samples=10000)
  >>> o_abund = unc.normal(center=-6., std=.5, n_samples=10000)
  >>> h_abund = unc.normal(center=-0.7, std=.1, n_samples=10000)
  >>> feh = fe_abund - h_abund
  >>> ofe = o_abund - fe_abund
  >>> plt.scatter(ofe.distribution, feh.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlabel('[Fe/H]') # doctest: +SKIP
  >>> plt.ylabel('[O/Fe]') # doctest: +SKIP

This demonstrates that the correlations naturally arise from the variables, but
there is no need to explicitly account for it: the sampling process naturally
recovers correlations that are present.

An important note of warning, however, is that the covariance is only preserved
if the sampling axes are exactly matched sample-by-sample.  If they are not, all
covariance information is (silently) lost:

.. plot::
  :context: close-figs
  :include-source:
  :align: center

  >>> n2_wrong = unc.Distribution(ncov[::-1, 1])  #reverse the sampling axis order
  >>> plt.scatter(n1.distribution, n2_wrong.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP


Moreover, an insufficiently-sampled distribution may give poor estimates or
hide correlations.  The example below is the same as the covariant gaussian
example above, but with 200x fewer samples:


.. plot::
  :context: close-figs
  :include-source:
  :align: center

  >>> ncov = np.random.multivariate_normal([0, 0], [[1, .5], [.5, 2]], size=50)
  >>> n1 = unc.Distribution(ncov[:, 0])
  >>> n2 = unc.Distribution(ncov[:, 1])
  >>> plt.scatter(n1.distribution, n2.distribution, s=5, lw=0) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP
  >>> np.cov(n1.distribution, n2.distribution) # doctest: +FLOAT_CMP
  array([[1.04667972, 0.19391617],
         [0.19391617, 1.50899902]])


The covaraiance structure is much less apparent by eye, and this is reflected
in significant discrepencies between the input and output covariance matrix.
In general this is an intrinsic trade-off using sampled distributions: a smaller
number of samples is computationally more efficient, but leads to larger
uncertainties in any of  the relevant quantities.  These tend to be of order
:math:`\sqrt{n_{\rm samples}}` in any derived quantity, but that depends on the
complexity of the distribution in question.  You have been warned!


Reference/API
=============

.. automodapi:: astropy.uncertainty
