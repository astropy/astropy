.. |Quantity| replace:: :class:`~astropy.units.Quantity`
.. |Distribution| replace:: :class:`~astropy.uncertainty.Distribution`
.. |ndarray| replace:: :class:`numpy.ndarray`

.. _astropy-uncertainty:

*******************************************************
Uncertainties and Distributions (`astropy.uncertainty`)
*******************************************************

.. note::

    `astropy.uncertainty` is relatively new (``astropy`` v3.1), and thus it is
    possible there will be API changes in upcoming versions of ``astropy``. If
    you have specific ideas for how it might be improved, please let us know on
    the `astropy-dev mailing list`_ or at http://feedback.astropy.org.

Introduction
============

``astropy`` provides a |Distribution| object to represent statistical
distributions in a form that acts as a drop-in replacement for a |Quantity|
object or a regular |ndarray|. Used in this manner, |Distribution| provides
uncertainty propagation at the cost of additional computation. It can also more
generally represent sampled distributions for Monte Carlo calculation
techniques, for instance.

The core object for this feature is the |Distribution|. Currently, all
such distributions are Monte Carlo sampled. While this means each distribution
may take more memory, it allows arbitrarily complex operations to be performed
on distributions while maintaining their correlation structure. Some specific
well-behaved distributions (e.g., the normal distribution) have
analytic forms which may eventually enable a more compact and efficient
representation. In the future, these may provide a coherent uncertainty
propagation mechanism to work with `~astropy.nddata.NDData`. However, this is
not currently implemented. Hence, details of storing uncertainties for
`~astropy.nddata.NDData` objects can be found in the :ref:`astropy_nddata`
section.

Getting Started
===============

To demonstrate a basic use case for distributions, consider the problem of
uncertainty propagation of normal distributions. Assume there are two
measurements you wish to add, each with normal uncertainties. We start
with some initial imports and setup::

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy import uncertainty as unc
  >>> np.random.seed(12345)  # ensures reproducible example numbers

Now we create two |Distribution| objects to represent our distributions::

  >>> a = unc.normal(1*u.kpc, std=30*u.pc, n_samples=10000)
  >>> b = unc.normal(2*u.kpc, std=40*u.pc, n_samples=10000)

For normal distributions, the centers should add as expected, and the standard
deviations add in quadrature. We can check these results (to the limits of our
Monte Carlo sampling) trivially with |Distribution| arithmetic and attributes::

  >>> c = a + b
  >>> c # doctest: +ELLIPSIS
  <QuantityDistribution [...] kpc with n_samples=10000>
  >>> c.pdf_mean() # doctest: +FLOAT_CMP
  <Quantity 2.99970555 kpc>
  >>> c.pdf_std().to(u.pc) # doctest: +FLOAT_CMP
  <Quantity 50.07120457 pc>

Indeed these are close to the expectations. While this may seem unnecessary for
the basic Gaussian case, for more complex distributions or arithmetic
operations where error analysis becomes untenable, |Distribution| still powers
through::

  >>> d = unc.uniform(center=3*u.kpc, width=800*u.pc, n_samples=10000)
  >>> e = unc.Distribution(((np.random.beta(2,5, 10000)-(2/7))/2 + 3)*u.kpc)
  >>> f = (c * d * e) ** (1/3)
  >>> f.pdf_mean() # doctest: +FLOAT_CMP
  <Quantity 2.99786227 kpc>
  >>> f.pdf_std() # doctest: +FLOAT_CMP
  <Quantity 0.08330476 kpc>
  >>> from matplotlib import pyplot as plt # doctest: +SKIP
  >>> from astropy.visualization import quantity_support # doctest: +SKIP
  >>> with quantity_support():
  ...     plt.hist(f.distribution, bins=50) # doctest: +SKIP

.. plot::

  import numpy as np
  from astropy import units as u
  from astropy import uncertainty as unc
  from astropy.visualization import quantity_support
  from matplotlib import pyplot as plt
  np.random.seed(12345)
  a = unc.normal(1*u.kpc, std=30*u.pc, n_samples=10000)
  b = unc.normal(2*u.kpc, std=40*u.pc, n_samples=10000)
  c = a + b
  d = unc.uniform(center=3*u.kpc, width=800*u.pc, n_samples=10000)
  e = unc.Distribution(((np.random.beta(2,5, 10000)-(2/7))/2 + 3)*u.kpc)
  f = (c * d * e) ** (1/3)
  with quantity_support():
      plt.hist(f.distribution, bins=50)


Using `astropy.uncertainty`
===========================

Creating Distributions
----------------------

.. EXAMPLE START: Creating Distributions Using Arrays or Quantities

The most direct way to create a distribution is to use an array or |Quantity|
that carries the samples in the *last* dimension::

  >>> import numpy as np
  >>> from astropy import units as u
  >>> from astropy import uncertainty as unc
  >>> np.random.seed(123456)  # ensures "random" numbers match examples below
  >>> unc.Distribution(np.random.poisson(12, (1000)))  # doctest: +ELLIPSIS
  NdarrayDistribution([..., 12,...]) with n_samples=1000
  >>> pq = np.random.poisson([1, 5, 30, 400], (1000, 4)).T * u.ct # note the transpose, required to get the sampling on the *last* axis
  >>> distr = unc.Distribution(pq)
  >>> distr # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
             [...],
             [...],
             [...]] ct with n_samples=1000>

Note the distinction for these two distributions: the first is built from an
array and therefore does not have |Quantity| attributes like ``unit``, while the
latter does have these attributes. This is reflected in how they interact with
other objects, for example, the ``NdarrayDistribution`` will not combine with
|Quantity| objects containing units.

.. EXAMPLE END

.. EXAMPLE START: Creating Distributions Using Helper Functions

For commonly used distributions, helper functions exist to make creating them
more convenient. The examples below demonstrate several equivalent ways to
create a normal/Gaussian distribution::

  >>> center = [1, 5, 30, 400]
  >>> n_distr = unc.normal(center*u.kpc, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=1000)
  >>> n_distr = unc.normal(center*u.kpc, var=[0.04, 2.25, 16, 1]*u.kpc**2, n_samples=1000)
  >>> n_distr = unc.normal(center*u.kpc, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2, n_samples=1000)
  >>> n_distr.distribution.shape
  (4, 1000)
  >>> unc.normal(center*u.kpc, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100).distribution.shape
  (4, 100)
  >>> unc.normal(center*u.kpc, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=20000).distribution.shape
  (4, 20000)

Additionally, Poisson and uniform |Distribution| creation functions exist::

  >>> unc.poisson(center*u.count, n_samples=1000) # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
               [...],
               [...],
               [...]] ct with n_samples=1000>
  >>> uwidth = [10, 20, 10, 55]*u.pc
  >>> unc.uniform(center=center*u.kpc, width=uwidth, n_samples=1000) # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
               [...],
               [...],
               [...]] kpc with n_samples=1000>
  >>> unc.uniform(lower=center*u.kpc - uwidth/2,  upper=center*u.kpc + uwidth/2, n_samples=1000)  # doctest: +ELLIPSIS
  <QuantityDistribution [[...],
               [...],
               [...],
               [...]] kpc with n_samples=1000>

.. EXAMPLE END

Users are free to create their own distribution classes following similar
patterns.

Using Distributions
-------------------

.. EXAMPLE START: Accessing Properties of Distributions

This object now acts much like a |Quantity| or |ndarray| for all but the
non-sampled dimension, but with additional statistical operations that work on
the sampled distributions::

  >>> distr.shape
  (4,)
  >>> distr.size
  4
  >>> distr.unit
  Unit("ct")
  >>> distr.n_samples
  1000
  >>> distr.pdf_mean() # doctest: +FLOAT_CMP
  <Quantity [  0.998,   5.017,  30.085, 400.345] ct>
  >>> distr.pdf_std() # doctest: +FLOAT_CMP
  <Quantity [ 0.97262326,  2.32222114,  5.47629208, 20.6328373 ] ct>
  >>> distr.pdf_var() # doctest: +FLOAT_CMP
  <Quantity [  0.945996,   5.392711,  29.989775, 425.713975] ct2>
  >>> distr.pdf_median()
  <Quantity [   1.,   5.,  30., 400.] ct>
  >>> distr.pdf_mad()  # Median absolute deviation # doctest: +FLOAT_CMP
  <Quantity [ 1.,  2.,  4., 14.] ct>
  >>> distr.pdf_smad()  # Median absolute deviation, rescaled to match std for normal # doctest: +FLOAT_CMP
  <Quantity [ 1.48260222,  2.96520444,  5.93040887, 20.75643106] ct>
  >>> distr.pdf_percentiles([10, 50, 90])
  <Quantity [[  0. ,   2. ,  23. , 374. ],
             [  1. ,   5. ,  30. , 400. ],
             [  2. ,   8. ,  37.1, 427. ]] ct>
  >>> distr.pdf_percentiles([.1, .5, .9]*u.dimensionless_unscaled)
  <Quantity [[  0. ,   2. ,  23. , 374. ],
            [  1. ,   5. ,  30. , 400. ],
            [  2. ,   8. ,  37.1, 427. ]] ct>

If need be, the underlying array can then be accessed from the ``distribution``
attribute::

  >>> distr.distribution  # doctest: +ELLIPSIS
  <Quantity [[...1...],
             [...5...],
             [...27...],
             [...405...]] ct>
  >>> distr.distribution.shape
  (4, 1000)

.. EXAMPLE END

.. EXAMPLE START: Interaction Between Quantity Objects and Distributions

A |Quantity| distribution interacts naturally with non-|Distribution|
|Quantity| objects, assuming the |Quantity| is a Dirac delta distribution::

  >>> distr_in_kpc = distr * u.kpc/u.count  # for the sake of round numbers in examples
  >>> distrplus = distr_in_kpc + [2000,0,0,500]*u.pc
  >>> distrplus.pdf_median()
  <Quantity [   3. ,   5. ,  30. , 400.5] kpc>
  >>> distrplus.pdf_var() # doctest: +FLOAT_CMP
  <Quantity [  0.945996,   5.392711,  29.989775, 425.713975] kpc2>

It also operates as expected with other distributions (but see below for a
discussion of covariances)::

  >>> another_distr = unc.Distribution((np.random.randn(1000,4)*[1000,.01 , 3000, 10] + [2000, 0, 0, 500]).T * u.pc)
  >>> combined_distr = distr_in_kpc + another_distr
  >>> combined_distr.pdf_median()  # doctest: +FLOAT_CMP
  <Quantity [  3.01847755,   4.99999576,  29.60559788, 400.49176321] kpc>
  >>> combined_distr.pdf_var()  # doctest: +FLOAT_CMP
  <Quantity [  1.8427705 ,   5.39271147,  39.5343726 , 425.71324244] kpc2>

.. EXAMPLE END

Covariance in Distributions and Discrete Sampling Effects
---------------------------------------------------------

One of the main applications for distributions is uncertainty propagation, which
critically requires proper treatment of covariance. This comes naturally in the
Monte Carlo sampling approach used by the |Distribution| class, as long as
proper care is taken with sampling error.

.. EXAMPLE START: Covariance in Distributions

To start with a basic example, two un-correlated distributions should produce
an un-correlated joint distribution plot:

.. plot::
  :context: close-figs
  :include-source:

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

Indeed, the distributions are independent. If we instead construct a covariant
pair of Gaussians, it is immediately apparent:

.. plot::
  :context: close-figs
  :include-source:

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

.. EXAMPLE END

.. EXAMPLE START: Preserving Covariance in Distributions

An important note of warning, however, is that the covariance is only preserved
if the sampling axes are exactly matched sample by sample. If they are not, all
covariance information is (silently) lost:

.. plot::
  :context: close-figs
  :include-source:

  >>> n2_wrong = unc.Distribution(ncov[::-1, 1])  #reverse the sampling axis order
  >>> plt.scatter(n1.distribution, n2_wrong.distribution, s=2, lw=0, alpha=.5) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP

Moreover, an insufficiently sampled distribution may give poor estimates or
hide correlations. The example below is the same as the covariant Gaussian
example above, but with 200x fewer samples:


.. plot::
  :context: close-figs
  :include-source:

  >>> ncov = np.random.multivariate_normal([0, 0], [[1, .5], [.5, 2]], size=50)
  >>> n1 = unc.Distribution(ncov[:, 0])
  >>> n2 = unc.Distribution(ncov[:, 1])
  >>> plt.scatter(n1.distribution, n2.distribution, s=5, lw=0) # doctest: +SKIP
  >>> plt.xlim(-4, 4) # doctest: +SKIP
  >>> plt.ylim(-4, 4) # doctest: +SKIP
  >>> np.cov(n1.distribution, n2.distribution) # doctest: +FLOAT_CMP
  array([[1.04667972, 0.19391617],
         [0.19391617, 1.50899902]])

The covariance structure is much less apparent by eye, and this is reflected
in significant discrepancies between the input and output covariance matrix.
In general this is an intrinsic trade-off using sampled distributions: a smaller
number of samples is computationally more efficient, but leads to larger
uncertainties in any of the relevant quantities. These tend to be of order
:math:`\sqrt{n_{\rm samples}}` in any derived quantity, but that depends on the
complexity of the distribution in question.

.. EXAMPLE END

.. note that if this section gets too long, it should be moved to a separate
   doc page - see the top of performance.inc.rst for the instructions on how to do
   that
.. include:: performance.inc.rst

Reference/API
=============

.. automodapi:: astropy.uncertainty
