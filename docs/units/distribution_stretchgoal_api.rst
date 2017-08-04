.. doctest-skip-all
>>> # code requires imports and init from distribution_api.rst

"Stretch goals" (i.e., in the plan but could be done later):

Have the named/analytic distributions actually propagate analytically when it makes sense, instead of Monte Carlo:

>>> real_distr1 = u.NormalDistribution(10*u.kpc, std=3*u.kpc, n_samples=None)
>>> real_distr2 = u.NormalDistribution(10*u.kpc, std=4*u.kpc, n_samples=None)
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

>>> u.Distribution.n_samples = 10000
>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc).shape
(10000, 4)

Allow easy-access to histograms:

>>> distr.hist(bins='knuth')
<matplotlib hist plot comes up>

Add stricter assumptions on allowable units with discrete distributions:

>>> u.PoissonDistribution([1, 5, 30, 400]*u.kpc)
UnitsError("Poisson distribution is discrete, so it doesn't make sense to use it with a unit like kpc which is a continuous physical value")
