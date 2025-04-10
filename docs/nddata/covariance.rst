
.. _nddata-covariance:

Covariance
**********

Overview
========

For a data vector, :math:`{\mathbf x} = \{x_0, x_1, ...\}`, the covariance
between any two elements :math:`x_i` and :math:`x_j` define the elements of the
*covariance matrix*

.. math::

    \Sigma_{ij} = \rho_{ij} \sigma_i \sigma_j,

where :math:`\rho_{ij}` are the elements of the *correlation matrix* and
:math:`V_i \equiv \sigma^2_i` is the variance in :math:`x_i`.  The covariance
matrix is, by definition, symmetric and positive-semidefinite (all eigenvalues
are non-negative).

The `~astropy.nddata.covariance.Covariance` object is a general utility for
constructing, visualizing, and storing two-dimensional covariance matrices.  To
minimize its memory footprint, the class uses sparse matrices (i.e., the module
requires `scipy.sparse`) and only stores the variance vector and the upper
triangle of the correlation matrix.

The class provides two convenient *static* methods for swapping between a full
covariance matrix (`~astropy.nddata.covariance.Covariance.revert_correlation`)
and the combination of a variance vector and correlation matrix
(`~astropy.nddata.covariance.Covariance.to_correlation`).

.. _nddata-covariance-intro:

Introductory Examples
---------------------

As a general introduction, let :math:`{\mathbf x}` contain 10 measurements.  Let
us set the correlation coefficient between adjacent measurements to 0.5
(:math:`\rho_{ij} = 0.5\ {\rm for}\ |j-i| = 1`), to 0.2 for next but one
measurements (:math:`\rho_{ij} = 0.2\ {\rm for}\ |j-i| = 2`), and to 0
otherwise.  If we adopt a uniform unity variance for all elements of
:math:`{\mathbf x}`, we can directly construct this (banded) covariance matrix
in python as follows:

>>> import numpy as np
>>>
>>> # Create the covariance matrix as a dense array
>>> npts = 10
>>> c = (np.diag(np.full(npts-2, 0.2, dtype=float), k=-2)
...         + np.diag(np.full(npts-1, 0.5, dtype=float), k=-1)
...         + np.diag(np.full(npts, 1.0, dtype=float), k=0)
...         + np.diag(np.full(npts-1, 0.5, dtype=float), k=1)
...         + np.diag(np.full(npts-2, 0.2, dtype=float), k=2))
>>> c
array([[1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. ],
       [0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. , 0. ],
       [0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. ],
       [0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. ],
       [0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. ],
       [0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. ],
       [0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. ],
       [0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2],
       [0. , 0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5],
       [0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. ]])

In this case, the correlation matrix and the covariance matrix are *identical*
because the elements of the variance vector are all unity.

With a correlation matrix, we can construct the covariance matrix with any
arbitrary variance vector.  Continuing the example above, the following creates
a new covariance matrix with a variance vector with all elements equal to 4:

>>> new_sig = np.full(npts, 2., dtype=float)
>>> new_c = c * new_sig[:,None] * new_sig[None,:]
>>> new_c
array([[4. , 2. , 0.8, 0. , 0. , 0. , 0. , 0. , 0. , 0. ],
       [2. , 4. , 2. , 0.8, 0. , 0. , 0. , 0. , 0. , 0. ],
       [0.8, 2. , 4. , 2. , 0.8, 0. , 0. , 0. , 0. , 0. ],
       [0. , 0.8, 2. , 4. , 2. , 0.8, 0. , 0. , 0. , 0. ],
       [0. , 0. , 0.8, 2. , 4. , 2. , 0.8, 0. , 0. , 0. ],
       [0. , 0. , 0. , 0.8, 2. , 4. , 2. , 0.8, 0. , 0. ],
       [0. , 0. , 0. , 0. , 0.8, 2. , 4. , 2. , 0.8, 0. ],
       [0. , 0. , 0. , 0. , 0. , 0.8, 2. , 4. , 2. , 0.8],
       [0. , 0. , 0. , 0. , 0. , 0. , 0.8, 2. , 4. , 2. ],
       [0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.8, 2. , 4. ]])

Or likewise for heteroscedastic data:

>>> new_sig = 1. + np.absolute(np.arange(npts) - npts//2).astype(float)
>>> new_sig
array([6., 5., 4., 3., 2., 1., 2., 3., 4., 5.])
>>> new_c = c * new_sig[:,None] * new_sig[None,:]
>>> new_c
array([[36. , 15. ,  4.8,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
       [15. , 25. , 10. ,  3. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
       [ 4.8, 10. , 16. ,  6. ,  1.6,  0. ,  0. ,  0. ,  0. ,  0. ],
       [ 0. ,  3. ,  6. ,  9. ,  3. ,  0.6,  0. ,  0. ,  0. ,  0. ],
       [ 0. ,  0. ,  1.6,  3. ,  4. ,  1. ,  0.8,  0. ,  0. ,  0. ],
       [ 0. ,  0. ,  0. ,  0.6,  1. ,  1. ,  1. ,  0.6,  0. ,  0. ],
       [ 0. ,  0. ,  0. ,  0. ,  0.8,  1. ,  4. ,  3. ,  1.6,  0. ],
       [ 0. ,  0. ,  0. ,  0. ,  0. ,  0.6,  3. ,  9. ,  6. ,  3. ],
       [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  1.6,  6. , 16. , 10. ],
       [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  3. , 10. , 25. ]])

Finally, note that reordering the elements of :math:`\mathbf{x}` also requires
reordering its covariance matrix (just as you would need to reorder the error
vector for an uncorrelated dataset).  The example below first generates a vector
with a shuffled set of indices.  The reordered :math:`\mathbf{x}` vector would
be constructed by setting ``reordered_x = x[i]`` and the covariance matrix would
be reordered using `numpy.ix_`, as follows:

>>> rng = np.random.default_rng(99)
>>> i = np.arange(npts)
>>> rng.shuffle(i)
>>> i
array([4, 6, 0, 3, 8, 1, 2, 5, 7, 9])
>>> reordered_c = c[np.ix_(i,i)]
>>> reordered_c
array([[1. , 0.2, 0. , 0.5, 0. , 0. , 0.2, 0.5, 0. , 0. ],
       [0.2, 1. , 0. , 0. , 0.2, 0. , 0. , 0.5, 0.5, 0. ],
       [0. , 0. , 1. , 0. , 0. , 0.5, 0.2, 0. , 0. , 0. ],
       [0.5, 0. , 0. , 1. , 0. , 0.2, 0.5, 0.2, 0. , 0. ],
       [0. , 0.2, 0. , 0. , 1. , 0. , 0. , 0. , 0.5, 0.5],
       [0. , 0. , 0.5, 0.2, 0. , 1. , 0.5, 0. , 0. , 0. ],
       [0.2, 0. , 0.2, 0.5, 0. , 0.5, 1. , 0. , 0. , 0. ],
       [0.5, 0.5, 0. , 0.2, 0. , 0. , 0. , 1. , 0.2, 0. ],
       [0. , 0.5, 0. , 0. , 0.5, 0. , 0. , 0.2, 1. , 0.2],
       [0. , 0. , 0. , 0. , 0.5, 0. , 0. , 0. , 0.2, 1. ]])

Note that the diagonal of ``reordered_c`` is still unity (all elements of
:math:`\mathbf{x}` are perfectly correlated with themselves), but the
off-diagonal terms have been rearranged to maintain the pre-existing
correlations.

In N-dimensions
---------------

    Covariance matrices of higher dimensional arrays are always assumed to be
    stored following row-major indexing.  That is, the covariance value
    :math:`\Sigma_{ij}` for an image of size :math:`(nx,ny)` is the covariance
    between image pixels :math:`I_{x_i,y_i}` and :math:`I_{x_j,y_j}`, where
    :math:`i = x_i + nx y_i` and, similarly, :math:`j = x_j + nx y_j`.



.. _nddata-covariance-construction:

Construction
============

There are numerous methods that can be used to construct
`~astropy.nddata.covariance.Covariance` objects, as well as methods that can be
used to store and reload them.

In *all* of the following examples, the object ``c`` is the banded covariance
array created at the beginning of the :ref:`nddata-covariance-covariance-access`
section.

Instantiating from pre-existing arrays
--------------------------------------

The simplest instantiation methods are based on using data that are already
available.

To create a `~astropy.nddata.covariance.Covariance` object from a
variance vector:

.. doctest-requires:: scipy

    >>> # Create from a variance vector
    >>> var = np.ones(3, dtype=float)
    >>> # Create from the Covariance object
    >>> covar = Covariance.from_variance(var)
    >>> # Test its contents
    >>> bool(np.array_equal(covar.to_dense(), np.identity(3)))
    True

In this case, the variance is unity for all elements of the data array such that
the covariance matrix is diagonal and identical to the identity matrix.

.. note::
    
    Wrapping the result of `~numpy.array_equal` with the ``bool`` operator above
    is done just to be sure that the returned value is ``True``, regardless of
    the version of numpy installed.

To create a `~astropy.nddata.covariance.Covariance` object from a "dense" (i.e.,
fully populated) covariance matrix:

.. doctest-requires:: scipy

    >>> # Instantiate from a covariance array
    >>> covar = Covariance(array=c)
    >>> bool(np.array_equal(covar.to_dense(), c))
    True
    >>> covar.to_dense()
    array([[1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. , 0. , 0. ],
           [0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. , 0. ],
           [0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. , 0. ],
           [0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. , 0. ],
           [0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. , 0. ],
           [0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. , 0. ],
           [0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2, 0. ],
           [0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5, 0.2],
           [0. , 0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. , 0.5],
           [0. , 0. , 0. , 0. , 0. , 0. , 0. , 0.2, 0.5, 1. ]])

Here, we instantiated the object using a pre-created matrix.

.. important::
    
    The last statement uses `~astropy.nddata.covariance.Covariance.to_dense` to
    access the array; see :ref:`nddata-covariance-covariance-access`.

Instantiating from random samples
---------------------------------

You can construct a covariance matrix based on samples from a distribution:

.. doctest-requires:: scipy

    >>> # Set the mean to 0 for all elements
    >>> m = np.zeros(npts, dtype=float)
    >>>
    >>> # Sample the multivariate normal distribution with the provided
    >>> # mean and covariance.
    >>> s = rng.multivariate_normal(m, c, size=100000)
    >>>
    >>> # Construct the covariance matrix from the random samples
    >>> covar = Covariance.from_samples(s.T, cov_tol=0.1)
    >>>
    >>> # Test that the known input covariance matrix is close to the
    >>> # measured covariance from the random samples
    >>> bool(np.all(np.absolute(c - covar.to_dense()) < 0.02))
    True

Here, we have drawn samples from a known multivariate normal distribution with a
mean of zero (``m``) and a known covariance matrix (``c``), defined for the 10
(``npts``) elements in the dataset (e.g., 10 pixels in a spectrum).  The code
checks the reconstruction of the known covariance matrix against the result
built from these random samples using
`~astropy.nddata.covariance.Covariance.from_samples`.

Instantiating from a matrix multiplication
------------------------------------------

Linear operations on a dataset (e.g., binning or smoothing) can be written as
matrix multiplications of the form

.. math::

    {\mathbf y} = {\mathbf T}\ {\mathbf x},

where :math:`{\mathbf T}` is a transfer matrix of size :math:`N_y\times N_x`,
:math:`{\mathbf x}` is a vector of size :math:`N_x`, and :math:`{\mathbf y}` is
a vector of length :math:`{N_y}` that results from the multiplication.  If
:math:`{\mathbf \Sigma}_x` is the covariance matrix for :math:`{\mathbf x}`, then
the covariance matrix for :math:`{\mathbf Y}` is

.. math::

    {\mathbf \Sigma}_y = {\mathbf T}\ {\mathbf \Sigma}_x\ {\mathbf T}^\top.

The example below shows how to build a covariance matrix from a matrix
multiplication using
`~astropy.nddata.covariance.Covariance.from_matrix_multiplication`:

.. doctest-requires:: scipy

    >>> # Construct a dataset
    >>> x = np.arange(npts, dtype=float)
    >>>
    >>> # Construct a transfer matrix that simply selects the elements at
    >>> # indices 0, 2, and 4
    >>> t = np.zeros((3,npts), dtype=float)
    >>> t[0,0] = 1.0
    >>> t[1,2] = 1.0
    >>> t[2,4] = 1.0
    >>>
    >>> # Get y
    >>> y = np.dot(t, x)
    >>> y
    array([0., 2., 4.])
    >>>
    >>> # Construct the covariance matrix
    >>> covar = Covariance.from_matrix_multiplication(t, c)
    >>> 
    >>> # Test the result
    >>> _c = (np.diag(np.full(3-1, 0.2, dtype=float), k=-1)
    ...         + np.diag(np.full(3, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(3-1, 0.2, dtype=float), k=1))
    >>> _c
    array([[1. , 0.2, 0. ],
           [0.2, 1. , 0.2],
           [0. , 0.2, 1. ]])
    >>> bool(np.array_equal(covar.to_dense(), _c))
    True

.. _nddata-covariance-data-access:

Accessing the data
==================

.. _nddata-covariance-covariance-access:

The `~astropy.nddata.covariance.Covariance` object is primarily a storage and IO
utility. Internally, the object stores the covariance matrix as a variance
vector and the upper triangle of the correlation matrix.  **This means that you
cannot directly access a covariance value within the object itself**; you must use
the functions described below.

Covariance Matrix
-----------------

There are two ways to access the full covariance matrix: Use 
`~astropy.nddata.covariance.Covariance.to_sparse` to produce a sparse matrix and
`~astropy.nddata.covariance.Covariance.to_dense` for a dense matrix.  The output
of these two methods can be used as you would use any `scipy.sparse.csr_matrix`
or `numpy.ndarray` object, respectively.

.. _nddata-covariance-correl-access:

Variance Vector and Correlation Matrix
--------------------------------------

The variance vector is stored as an accessible property
(`~astropy.nddata.covariance.Covariance.variance`), but note that the property
is immutable.

Access to the full correlation matrix is provided using
`~astropy.nddata.covariance.Covariance.to_sparse` to produce a sparse matrix or
`~astropy.nddata.covariance.Covariance.to_dense` and setting ``correlation =
True``.   

.. _nddata-covariance-fitsio:

FITS file I/O methods
=====================

`~astropy.nddata.covariance.Covariance` objects can be saved as a binary table
in a FITS file using the `~astropy.nddata.covariance.Covariance.write` method.
To reload the covariance matrix, use the
`~astropy.nddata.covariance.Covariance.from_fits` instantiation method:

.. doctest-requires:: scipy

    >>> import numpy as np
    >>> from astropy.nddata.covariance import Covariance
    >>> ofile = 'test_covar_io.fits'
    >>> m = np.zeros(10, dtype=float)
    >>> c = (np.diag(np.full(10-2, 0.2, dtype=float), k=-2)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=-1)
    ...         + np.diag(np.full(10, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=1)
    ...         + np.diag(np.full(10-2, 0.2, dtype=float), k=2))
    >>> s = np.random.multivariate_normal(m, c, size=100000)
    >>> covar = Covariance.from_samples(s.T, cov_tol=0.1)
    >>> covar.write(ofile)
    >>> from astropy.io import fits
    >>> with fits.open(ofile) as hdu:
    ...     hdu.info()
    Filename: test_covar_io.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       7   ()
      1  VAR           1 ImageHDU         9   (10,)   float64
      2  CORREL        1 BinTableHDU     18   27R x 3C   [K, K, D]
    >>> _covar = Covariance.from_fits(ofile)
    >>> bool(np.allclose(covar.to_dense(), _covar.to_dense()))
    True

The details of how the covariance data are stored are described by the
`~astropy.nddata.covariance.Covariance.write` method.
