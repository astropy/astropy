.. |quantity| replace:: :class:`~astropy.units.Quantity`
.. |distribution| replace:: :class:`~astropy.units.Distribution`

.. _unit_distributions:

Distributions
*************

In addition to |quantity|, astropy provides a |distribution| object to represent
statistical distributions in a form that acts as a drop-in replacement for
|quantity| or regulay Numpy arrays. USed in this manner, |distribution| provides
uncertainty propogation at the cost of additional computation.  It can also
more generally represent sampled distributions for e.g., Monte Carlo calculation
techniques.

The core object for this feature is the |distribution|.  Currently, all
such distributions are Monte Carlo sampled.  While this means each distribution
may take more memory, it allows arbitrarily complex operations to be performed
on distributions while maintaining their correlation structure. Some specific
well-behaved distributions (e.g., `~astropy.units.NormalDistribution`) have
analytic forms which may eventually enable a more compact/efficient
representation, but this is not currently implemented.


Creating distributions
======================

The most conceptually straightforward way to create a distribution is to use an
array or quantity that carries the samples in the *last* dimension:

  >>> import numpy as np
  >>> from astropy import units as u
  >>> np.random.seed(123456)  # ensures "random" numbers match examples below

  >>> u.Distribution(np.random.poisson(12, (1000))) # doctest: +ELLIPSIS
  ndarrayDistribution([..., 12,...]) with n_samples=1000
  >>> pq = np.random.poisson([1, 5, 30, 400], (1000, 4)).T * u.kpc # note the transpose, required to get the sampling on the *last* axis
  >>> distr = u.Distribution(pq)
  >>> distr # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
             [...],
             [...],
             [...]] kpc with n_samples=1000>

For commonly-used distributions, helper classes exist  to make creating them
easier. Below demonstrates several equivalent ways to create a normal/Gaussian
distribution:

  >>> centerq = [1, 5, 30, 400]*u.kpc
  >>> n_distr = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
  >>> n_distr = u.NormalDistribution(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
  >>> n_distr = u.NormalDistribution(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)
  >>> n_distr.distribution.shape
  (4, 1000)
  >>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100).distribution.shape
  (4, 100)
  >>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=20000).distribution.shape
  (4, 20000)


Additionally, Poisson and uniform |distribution| classes exist:

  >>> p_dist = u.PoissonDistribution(centerq)
  >>> uwidth = [10, 20, 10, 55]*u.pc
  >>> u_dist = u.UniformDistribution(lower=centerq-uwidth/2,  upper=centerq+uwidth/2)

Users are free to create their own distribution classes following similar
patterns.


Using Distributions
===================

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
  >>> distr.percentiles([10, 50, 90])
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

  >>> another_distr = u.Distribution((np.random.randn(1000,4)*[1000,.01 , 3000, 10] + [2000, 0, 0, 500]).T, unit=u.pc)
  >>> combined_distr = distr + another_distr
  >>> combined_distr.pdf_median
  <Quantity [  2.90856297,   4.99999764,  30.09385367, 400.50056651] kpc>
  >>> combined_distr.pdf_var
  <Quantity [  2.0051053 ,   5.39271159,  38.24442151, 425.70428603] kpc2>


Covariance in distributions
===========================

One of the main applications for distributions is unceratinty propogation, which
critically requires proper treatment of covariance. This comes naturally in the
Monte Carlo sampling approach used by the |distribution| class, as long as
proper care is taken with sampling error.

To start with a simple example, two un-correlated distributions should produce
an un-correlated joint distribution plot:

.. plot::
  :context:
  :include-source:
  :align: center

  >>> from matplotlib import pyplot as plt # doctest: +SKIP
  >>> n1 = u.NormalDistribution(center=0., std=1, n_samples=10000)
  >>> n2 = u.NormalDistribution(center=0., std=2, n_samples=10000)
  >>> plt.scatter(n1.distribution, n2.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP

Indeed, the distributions are independent.  If we instead construct a covariant
pair of gaussians, it is immediately apparent:

.. plot::
  :context:
  :include-source:
  :align: center

  >>> ncov = np.random.multivariate_normal([0, 0], [[1, .5], [.5, 2]], size=10000)
  >>> n1 = u.Distribution(ncov[:, 0])
  >>> n2 = u.Distribution(ncov[:, 1])
  >>> plt.scatter(n1.distribution, n2.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP


Most importantly, the proper correlated structure is preserved or generated as
expected by appropriate arithmetic operations. For example, ratios of
uncorrelated normal distribution gain covariances if the axes are not
independent, as in this simulation of iron, hydrogen, and oxygen abundances in
a hypothetical collection of stars:

.. plot::
  :context:
  :include-source:
  :align: center

  >>> fe_abund = u.NormalDistribution(center=-2, std=.25, n_samples=10000)
  >>> o_abund = u.NormalDistribution(center=-6., std=.5, n_samples=10000)
  >>> h_abund = u.NormalDistribution(center=-0.7, std=.1, n_samples=10000)
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
  :context:
  :include-source:
  :align: center

  >>> n2_wrong = u.Distribution(ncov[::-1, 1])  #reverse the sampling axis order
  >>> plt.scatter(n1.distribution, n2_wrong.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP


Moreover, an insufficiently-sampled distribution may give poor estimates or
hide correlations.  The example below is the same as the covariant gaussian
example above, but with 100x fewer samples:


.. plot::
  :context:
  :include-source:
  :align: center

  >>> ncov = np.random.multivariate_normal([0, 0], [[1, .5], [.5, 2]], size=100)
  >>> n1 = u.Distribution(ncov[:, 0])
  >>> n2 = u.Distribution(ncov[:, 1])
  >>> plt.scatter(n1.distribution, n2.distribution, s=5, lw=0) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP

In general this is an intrinsic trade-off using sampled distributions: a smaller
number of samples is computationally more efficient, but leads to larger
uncertainties in any of  the relevant quantities.  These tend to be of order
:math:`\sqrt(n_samples)` in any derived quantity, but that depends on the
complexity of the distribution in question.  You have been warned!
