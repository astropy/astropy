
.. _nddata-covariance:

Covariance
**********

Overview
========

The `~astropy.nddata.covariance.Covariance` object is a general utility for
constructing, visualizing, and storing two-dimensional covariance matrices.

.. _nddata-covariance-construction:

Construction
============

Beyond the nominal instantiation method, there are numerous methods used to
construct a `~astropy.nddata.covariance.Covariance` object.

Instantiating from pre-existing arrays
--------------------------------------

The simplest approaches are if you have a variance vector or a pre-calculated
covariance matrix as a dense array and you want to construct a
`~astropy.nddata.covariance.Covariance` object for further calculations::

.. doctest-requires:: scipy

    >>> import numpy as np
    >>> from astropy.nddata.covariance import Covariance
    >>> # Create from a variance vector
    >>> var = np.ones(3, dtype=float)
    >>> covar = Covariance.from_variance(var)
    >>> np.array_equal(covar.toarray(), np.identity(3))
    True
    >>> # Instantiate from a covariance array
    >>> c = (np.diag(np.full(10-2, 0.2, dtype=float), k=-2)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=-1)
    ...         + np.diag(np.full(10, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=1)
    ...         + np.diag(np.full(10-2, 0.2, dtype=float), k=2))
    >>> covar = Covariance(array=c)
    >>> np.array_equal(covar.toarray(), c)
    True
    >>> covar.toarray()
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

Note the use of the `~astropy.nddata.covariance.Covariance.toarray` above to
access the array; see :ref:`nddata-covariance-access`.

Instantiating from random samples
---------------------------------

You can construct a covariance matrix based on samples from a distribution::

.. doctest-requires:: scipy

    >>> import numpy as np
    >>> from astropy.nddata.covariance import Covariance
    >>> m = np.zeros(10, dtype=float)
    >>> c = (np.diag(np.full(10-2, 0.2, dtype=float), k=-2)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=-1)
    ...         + np.diag(np.full(10, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=1)
    ...         + np.diag(np.full(10-2, 0.2, dtype=float), k=2))
    >>> s = np.random.multivariate_normal(m, c, size=100000)
    >>> covar = Covariance.from_samples(s.T, cov_tol=0.1)
    >>> np.all(np.absolute(c - covar.toarray()) < 0.02)
    np.True_

Here, we have drawn samples from a known multivariate normal distribution with a
given covariance between its 10 axes and checked the reconstruction of the
covariance matrix from those samples using
`~astropy.nddata.covariance.Covariance.from_samples`. This construction method
is a simple wrapper for `numpy.cov` and
`~astropy.nddata.covariance.Covariance.from_array`.

Instantiating from a matrix multiplication
------------------------------------------

You can construct the covariance matrix that results from a matrix
multiplication::

.. doctest-requires:: scipy

    >>> import numpy as np
    >>> from astropy.nddata.covariance import Covariance
    >>> c = (np.diag(np.full(10-2, 0.2, dtype=float), k=-2)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=-1)
    ...         + np.diag(np.full(10, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(10-1, 0.5, dtype=float), k=1)
    ...         + np.diag(np.full(10-2, 0.2, dtype=float), k=2))
    >>> x = np.ones(10, dtype=float)
    >>> t = np.zeros((3,10), dtype=float)
    >>> t[0,0] = 1.0
    >>> t[1,2] = 1.0
    >>> t[2,4] = 1.0
    >>> y = np.dot(t, x)
    >>> _c = (np.diag(np.full(3-1, 0.2, dtype=float), k=-1)
    ...         + np.diag(np.full(3, 1.0, dtype=float), k=0)
    ...         + np.diag(np.full(3-1, 0.2, dtype=float), k=1))
    >>> _c
    array([[1. , 0.2, 0. ],
            [0.2, 1. , 0.2],
            [0. , 0.2, 1. ]])
    >>> covar = Covariance.from_matrix_multiplication(t, c)
    >>> np.array_equal(covar.toarray(), _c)
    True

Finally, you can construct the covariance matrix from a previous instance that
was saved to a FITS file using the :ref:`nddata-covariance-fitsio`.

.. _nddata-covariance-access:

Accessing the covariance data
=============================

The `~astropy.nddata.covariance.Covariance` object is primarily a storage and IO
utility. Internally, the object only keeps the upper triangle of the matrix,
which means that use of the ``cov`` attribute is *not* recommended unless you
know what you're doing.

There are two ways to access the full covariance matrix: Use 
`~astropy.nddata.covariance.Covariance.full` to produce a sparse matrix and
`~astropy.nddata.covariance.Covariance.toarray` for a dense matrix.  The output
of these two methods can be used as you would use any `scipy.sparse.csr_matrix`
or `numpy.ndarray` object, respectively.

.. _nddata-covariance-correl:

Toggling between covariance and correlation matrices
====================================================

The `~astropy.nddata.covariance.Covariance` object allows you to toggle between
the full covariance matrix, :math:`{\mathbf C}` and a correlation matrix,
:math:`{\mathbf \rho}`, where

.. math::

    \rho_{ij} = \frac{C_{ij}}{(V_i V_j)^{1/2}}

and :math:`{\mathbf V}` is the variance vector (the diagonal elements of
:math:`{\mathbf C}`). To convert a `~astropy.nddata.covariance.Covariance`
object to a correlation matrix (or ensure that it already is one), use
`~astropy.nddata.covariance.Covariance.to_correlation`. To revert back to a
covariance matrix, use
`~astropy.nddata.covariance.Covariance.revert_correlation`.

.. _nddata-covariance-fitsio:

FITS file I/O methods
=====================

`~astropy.nddata.covariance.Covariance` objects can be saved as a binary table
in a FITS file using the `~astropy.nddata.covariance.Covariance.write` method.
To reload the covariance matrix, use the
`~astropy.nddata.covariance.Covariance.from_fits` instantiation method::

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
    >>> np.allclose(covar.toarray(), _covar.toarray())
    True

The details of how the covariance data are stored are described by the
`~astropy.nddata.covariance.Covariance.write` method.
