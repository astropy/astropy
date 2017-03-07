>>> import numpy as np
>>> from astropy import units as u


Create a quantity from a pre-existing sampled distribution.  The MC direction
is the *first* dimension.

>>> parr = np.random.poisson([1, 5, 30, 400],(1000, 4))
>>> distr = u.Distribution(parr, u.kpc)

Or equivalently:

>>> pq = np.random.poisson([1, 5, 30, 400],(1000, 4))*u.kpc
>>> distr = u.Distribution(pq)


You can ask it about itself and it acts like an array including the dimension,
but can be "quantityified" as needed with standard statistics

>>> distr.shape
(1000, 4)
>>> distr.size
4000
>>> distr.nmc
1000
>>> distr.pdf_shape
(4, )
>>> distr.pdf_mean
<Quantity [   1.059,   4.939,  30.006, 400.915] kpc>
>>> distr.pdf_std
<Quantity [  1.03610762,  2.28719894,  5.32916166, 19.88888572] kpc>
>>> distr.pdf_var
<Quantity [   1.073519,   5.231279,  28.399964, 395.567775] kpc2>
>>> distr.pdf_median
<Quantity [   1.,   5.,  30., 400.] kpc>
>>> distr.pdf_mad  # Median absolute deviation
<Quantity [  1. ,  2. ,  4. , 13.5] kpc>
>>> distr.pdf_smad  # Median absolute deviation, rescaled to match std for normal
<Quantity [  1.4826,  2.9652,  5.9304, 20.0151] kpc>


Can also ask for more complex statistical summaries:

>>> distr.percentiles([.1, .5, .9])
<Quantity [[   0.,   2.,  23., 376.],
           [   1.,   5.,  30., 400.],
           [   2.,   8.,  37., 426.]] kpc>


The distribution should interact correctly with non-Distribution quantities:

>>> distrplus = distr + [2000,0,0,500]*u.pc
>>> distrplus.pdf_median
<Quantity [   3. ,   5. ,  30. , 400.5] kpc>
>>> distrplus.pdf_var
<Quantity [   1.073519,   5.231279,  28.399964, 395.567775] kpc2>


It should *also* combine reasonably with othe distributions:

>>> another_distr = u.Distribution(np.random.randn(1000,4)*[1000,.01 , 3000, 10] + [2000, 0, 0, 500], u.pc)
>>> combined_distr = distr + another_distr
>>> combined_distr.pdf_median
<Quantity [   2.94735032,   4.99999564,  29.67887559, 400.5067096 ] kpc>
>>> combined_distr.pdf_var
<Quantity [   2.06640821,   5.23127931,  38.66344651, 395.5653249 ] kpc2>


For commonly-used distributions, provide helpers to make creating them easier.
Note that all of the normal ones are equivalent, but note that the units must
make sense:

>>> centerq = [1, 5, 30, 400]*u.kpc
>>> n_distr = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
>>> n_distr = u.NormalDistribution(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
>>> n_distr = u.NormalDistribution(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)
>>> p_dist = u.PoissonDistribution(centerq)
>>> u_dist = u.UniformDistribution(centerq, width=[10, 20, 10, 55]*u.pc)

Can specify how many samples are desired:

>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, nmc=100).shape
(100, 4)
>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, nmc=20000).shape
(20000, 4)



"Stretch goal" (i.e., in the plan but could be done later):

Have the named/analytic distributions actually propagate analytically when it makes sense, instead of Monte Carlo:

>>> real_distr1 = u.NormalDistribution(10*u.kpc, std=3*u.kpc, nmc=None)
>>> real_distr2 = u.NormalDistribution(10*u.kpc, std=4*u.kpc, nmc=None)
>>> distr12_plus = real_distr1 + real_distr2
>>> distr12_plus.q
<Quantity 20.0 kpc>
>>> distr12_plus.unc  # analytically should be the quadrature sum of the components
<Quantity 5.0 kpc>

>>> distr12_mul = real_distr1 * real_distr2
>>> distr12_mul.pdf_
<Quantity 19.79716280988112 kpc>
>>> distr12_plus.pdf_std  # analytically should be the quadrature sum of the *pecentage* uncertainty of the components
<Quantity 4.196109517371489 kpc>

But they should still combine with MC distributions and sample in that case

>>> sampled_distr2 = u.Distribution(np.random.randn(100)*4+10, u.kpc)
>>> sample12_plus = real_distr1 + sampled_distr2
>>> sample12_plus.pdf_mean
<Quantity 20.343156238130238 kpc>
>>> sample12_plus.pdf_std  # analytically should be the quadrature sum of the components
<Quantity 5.137387685596549 kpc>


MORE DEBATABLE ITEMS:

Allow stateful choice of what "expected" and "uncertainty" mean:

>>> u.Distribution.q_type = 'median'
>>> distr.q
<Quantity [   3. ,   5. ,  30. , 400.5] kpc>
>>> u.Distribution.unc_type = 'std'
>>> distr.unc
<Quantity [  1.03610762,  2.28719894,  5.32916166, 19.88888572] kpc>

Allow stateful definition of default nmc:

>>> u.Distribution.nmc = 10000
>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc).shape
(10000, 4)

Allow easy-access to histograms:

>>> distr.hist(bins='knuth')
<matplotlib hist plot comes up>

Add stricter assumptions on allowable units with discrete distributions:

>>> u.PoissonDistribution([1, 5, 30, 400]*u.kpc)
UnitsError("Poisson distribution is discrete, so it doesn't make sense to use it with a unit like kpc which is a continuous physical value")
