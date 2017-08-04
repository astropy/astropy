>>> import numpy as np
>>> from astropy import units as u
>>> np.random.seed(12345)  # ensures "random" numbers match examples

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
>>> distr.n_samples
1000
>>> distr.distr_shape
(4,)
>>> distr.pdf_mean
<Quantity [   1.029,   4.985,  30.132, 400.463] kpc>
>>> distr.pdf_std
<Quantity [  1.01200741,  2.14866819,  5.42093866, 20.24402705] kpc>
>>> distr.pdf_var
<Quantity [   1.024159,   4.616775,  29.386576, 409.820631] kpc2>
>>> distr.pdf_median
<Quantity [   1.,   5.,  30., 400.] kpc>
>>> distr.pdf_mad  # Median absolute deviation
<Quantity [  1.,  1.,  4., 13.] kpc>
>>> distr.pdf_smad  # Median absolute deviation, rescaled to match std for normal
<Quantity [  1.48260222,  1.48260222,  5.93040887, 19.27382884] kpc>


Can also ask for more complex statistical summaries:

>>> distr.percentiles([10, 50, 90])
<Quantity [[   0.,   2.,  23., 375.],
           [   1.,   5.,  30., 400.],
           [   2.,   8.,  37., 426.]] kpc>


The distribution should interact correctly with non-Distribution quantities:

>>> distrplus = distr + [2000,0,0,500]*u.pc
>>> distrplus.pdf_median
<Quantity [   3. ,   5. ,  30. , 400.5] kpc>
>>> distrplus.pdf_var
<Quantity [   1.024159,   4.616775,  29.386576, 409.820631] kpc2>


It should *also* combine reasonably with othe distributions:

>>> another_distr = u.Distribution(np.random.randn(1000,4)*[1000,.01 , 3000, 10] + [2000, 0, 0, 500], u.pc)
>>> combined_distr = distr + another_distr
>>> combined_distr.pdf_median
<Quantity [   2.9548952 ,   4.99999855,  29.93483557, 400.50685423] kpc>
>>> combined_distr.pdf_var
<Quantity [   2.17250083,   4.6167747 ,  37.46238268, 409.82738255] kpc2>


For commonly-used distributions, provide helpers to make creating them easier.
Note that all of the normal ones are equivalent, but note that the units must
make sense:

>>> centerq = [1, 5, 30, 400]*u.kpc
>>> n_distr = u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc)
>>> n_distr = u.NormalDistribution(centerq, var=[0.04, 2.25, 16, 1]*u.kpc**2)
>>> n_distr = u.NormalDistribution(centerq, ivar=[25, 0.44444444, 0.625, 1]*u.kpc**-2)
>>> p_dist = u.PoissonDistribution(centerq)
>>> uwidth = [10, 20, 10, 55]*u.pc
>>> u_dist = u.UniformDistribution(lower=centerq-uwidth/2,  upper=centerq+uwidth/2)

Can specify how many samples are desired:

>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=100).shape
(100, 4)
>>> u.NormalDistribution(centerq, std=[0.2, 1.5, 4, 1]*u.kpc, n_samples=20000).shape
(20000, 4)
