
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
requires `scipy.sparse`) and only stores the upper triangle of the covariance
matrix.

The class provides two convenient *static* methods for swapping between a full
covariance matrix (`~astropy.nddata.covariance.Covariance.revert_correlation`)
and the combination of a variance vector and correlation matrix
(`~astropy.nddata.covariance.Covariance.to_correlation`).

.. _nddata-covariance-intro:

Introductory Examples
---------------------

As a general introduction to covariance matrices, let :math:`{\mathbf x}`
contain 10 measurements.  Let the correlation coefficient between adjacent
measurements be 0.5 (:math:`\rho_{ij} = 0.5\ {\rm for}\ |j-i| = 1`), 0.2 for
next but one measurements (:math:`\rho_{ij} = 0.2\ {\rm for}\ |j-i| = 2`), and 0
otherwise.  If we adopt unity variance for all elements of :math:`{\mathbf x}`,
we can directly construct the (banded) covariance matrix in python as follows:

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

>>> new_var = np.full(npts, 4., dtype=float)
>>> new_c = c * np.sqrt(new_var[:,None] * new_var[None,:])
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

>>> new_var = (1. + np.absolute(np.arange(npts) - npts//2).astype(float))**2
>>> new_var
array([36., 25., 16.,  9.,  4.,  1.,  4.,  9., 16., 25.])
>>> new_c = c * np.sqrt(new_var[:,None] * new_var[None,:])
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

.. note::

    The `~astropy.nddata.covariance.Covariance` class provides a convenience
    function for creating a new `~astropy.nddata.covariance.Covariance` instance
    with the same correlation matrix but a new variance vector; see
    :ref:`here<covariance-apply-new-variance>`.

Reordering and Subsets
----------------------

When reordering or down-selecting subsets of the elements of :math:`\mathbf{x}`,
these changes must be propagated to the associated covariance matrix, just as
would be needed for the error vector for an uncorrelated dataset.

The example below first generates a vector with a shuffled set of indices.  The
reordered :math:`\mathbf{x}` vector would be constructed by setting
``reordered_x = x[i]`` and the covariance matrix would be reordered using
`numpy.ix_`, as follows:

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

Creating a covariance matrix for a subset of data is a very similar operation.
If we want the covariance matrix for the first 3 elements of the data vector, we
can do the following:

>>> i = np.arange(3)
>>> sub_c = c[np.ix_(i,i)]
>>> sub_c
array([[1. , 0.5, 0.2],
       [0.5, 1. , 0.5],
       [0.2, 0.5, 1. ]])

.. note::

    The `~astropy.nddata.covariance.Covariance` class provides a convenience
    function for matching the covariance data to a slice of its parent data array;
    see :ref:`here<covariance-match-to-data-slice>`.

In N-dimensions
---------------

Covariance matrices can be constructed for arrays of higher dimensionality by
flattening the data arrays.  For a row-major array flattening order, one can
adopt the convention that :math:`\Sigma_{ij}` for an image of size
:math:`(N_x,N_y)` is the covariance between image pixels :math:`I_{x_i,y_i}` and
:math:`I_{x_j,y_j}`, where :math:`i = x_i + N_x\ y_i` and :math:`j = x_j + N_x\
y_j`.

As an example, let the covariance matrix ``c``, used throughout this section, be
the covariance matrix for a :math:`5 \times 2` array, instead of a 10-element
vector.  The complication is determining the mapping from the data array to the
relevant covariance element; we can do this using `numpy` functions as follows.
To determine the covariance between elements ``data[1,0]`` and ``data[2,0]``, we
convert the indices from the ``data`` to find a covariance of 0.2:

>>> data_array_shape = (5,2)
>>> i_data = (np.array([1]), np.array([0]))
>>> j_data = (np.array([2]), np.array([0]))
>>> i_cov = np.ravel_multi_index(i_data, data_array_shape)
>>> j_cov = np.ravel_multi_index(j_data, data_array_shape)
>>> i_cov, j_cov
(array([2]), array([4]))
>>> c[i_cov, j_cov]
array([0.2])

The inverse operation (determining the indices of the data array given the
indices in the covariance matrix) uses `~numpy.unravel_index` (cf. ``i_data``):

>>> np.unravel_index(i_cov, data_array_shape)
(array([1]), array([0]))

.. note::

    The `~astropy.nddata.covariance.Covariance` class provides convenience
    functions for switching between the data array and covariance matrix
    indexing when working with higher dimensionality data arrays;
    see :ref:`here<covariance-nd-indexing>`.

.. _nddata-covariance-construction:

Construction
============

Many methods are provided to construct a `~astropy.nddata.covariance.Covariance`
object.  In *all* of the following examples, the object ``c`` is the banded
covariance array created at the beginning of the :ref:`nddata-covariance-intro`
section.

Instantiating from pre-existing arrays
--------------------------------------

The simplest instantiation methods are based on using data that are already
available.

To create a `~astropy.nddata.covariance.Covariance` object from a
variance vector:

.. doctest-requires:: scipy

    >>> from astropy.nddata.covariance import Covariance
    >>> # Create from a variance vector
    >>> var = np.ones(3, dtype=float)
    >>> # Create from the Covariance object
    >>> covar = Covariance.from_variance(var)
    >>> # Test its contents
    >>> print(np.array_equal(covar.to_dense(), np.identity(3)))
    True

In this case, the variance is unity for all elements of the data array such that
the covariance matrix is diagonal and identical to the identity matrix.

To create a `~astropy.nddata.covariance.Covariance` object from a "dense" (i.e.,
fully populated) covariance matrix:

.. doctest-requires:: scipy

    >>> # Instantiate from a covariance array
    >>> covar = Covariance(array=c)
    >>> print(np.array_equal(covar.to_dense(), c))
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

.. important::

    The last statement uses `~astropy.nddata.covariance.Covariance.to_dense` to
    access the array; see :ref:`nddata-covariance-data-access`.

Above, the base instantiation method is used; however, the
`~astropy.nddata.covariance.Covariance.from_array` method is also provided.  The
primary difference is that the latter allows limits to be imposed on the
(absolute value of the) correlation or covariance values.

Finally, note that, by default, all instantiations of a
`~astropy.nddata.covariance.Covariance` object check that the input matrix is
symmetric.  If it is not, a warning is issued.  To skip the check and the
warning, set ``assume_symmetric=True``.  Regardless of whether or not the check
is performed, the object *only stores the upper triangle of the input matrix*
effectively meaning that any asymmetry in the matrix is lost when it is
ingested.

Instantiating from random samples
---------------------------------

You can construct a covariance matrix based on samples from a distribution using
`~astropy.nddata.covariance.Covariance.from_samples`:

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
    >>> print(np.all(np.absolute(c - covar.to_dense()) < 0.02))
    True

Here, we have drawn samples from a known multivariate normal distribution with a
mean of zero (``m``) and a known covariance matrix (``c``), defined for the 10
(``npts``) elements in the dataset (e.g., 10 pixels in a spectrum).  The code
checks the reconstruction of the known covariance matrix against the result
built from these random samples.

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
the covariance matrix for :math:`{\mathbf y}` is

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
    >>> print(np.array_equal(covar.to_dense(), _c))
    True

In N-dimensions
---------------

All of the instantiation methods above allow you to define the "data shape" of
the data array for the associated covariance matrix.  Following the previous
N-dimensional example, let ``c`` be the covariance matrix for a :math:`5 \times
2` array, instead of a 10-element vector.

.. doctest-requires:: scipy

    >>> data_array_shape
    (5, 2)
    >>> covar = Covariance(array=c, data_shape=data_array_shape)
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

The covariance matrix looks identical, but the higher dimensionality will affect
its :ref:`nddata-covariance-coord-access`.

.. _nddata-covariance-data-access:

Accessing the data
==================

The `~astropy.nddata.covariance.Covariance` object is primarily a storage
utility. Internally, the object only stores the upper triangle of the covariance
matrix.  **This means that you should not directly access a covariance value
within the object itself**; you must use the functions described below.

.. _nddata-covariance-covariance-access:

Covariance Matrix
-----------------

There are two ways to access the full covariance matrix:

- Use `~astropy.nddata.covariance.Covariance.to_sparse` to produce a sparse matrix or

- Use `~astropy.nddata.covariance.Covariance.to_dense` for a dense matrix.

The output of these two methods can be used as you would use any
`scipy.sparse.csr_matrix` or `numpy.ndarray` object, respectively.

.. _nddata-covariance-correl-access:

Variance Vector and Correlation Matrix
--------------------------------------

The variance vector is stored as an accessible property
(`~astropy.nddata.covariance.Covariance.variance`), but note that the property
is immutable.

Access to the full correlation matrix is provided using
`~astropy.nddata.covariance.Covariance.to_sparse` to produce a sparse matrix or
`~astropy.nddata.covariance.Covariance.to_dense` for a dense matrix by setting
the keyword argument ``correlation = True``.

.. _nddata-covariance-coord-access:

Coordinate Data
---------------

Although more useful as preparation for storage, the covariance data can also be
accessed in coordinate format:

.. doctest-requires:: scipy

    >>> covar = Covariance(array=c)
    >>> i, j, cij = covar.coordinate_data()
    >>> print(np.array_equal(covar.to_dense()[i,j], cij))
    True

The arrays returned by `~astropy.nddata.covariance.Covariance.coordinate_data`
provide the matrix coordinates (``i`` and ``j``) for the non-zero covariance
values (``cij``).

.. _nddata-covariance-table:

File IO
=======

The primary way to write/read `~astropy.nddata.covariance.Covariance` objects is
by first parsing the data into a `~astropy.table.Table` using the
`~astropy.nddata.covariance.Covariance.to_table` method:

.. doctest-requires:: scipy

    >>> covar = Covariance(array=c)
    >>> tbl = covar.to_table()
    >>> tbl.meta
    {'COVSHAPE': '(10, 10)'}
    >>> tbl[:3]
    <Table length=3>
    INDXI INDXJ COVARIJ
    int64 int64 float64
    ----- ----- -------
        0     0     1.0
        0     1     0.5
        0     2     0.2

The output above just shows the first 3 rows of the table to demonstrate that
the non-zero elements of the covariance matrix are stored in "coordinate
format."  Specifically, the data is provided in three columns:

- ``'INDXI'``: The row index in the covariance matrix (:math:`i`).

- ``'INDXJ'``: The column index in the covariance matrix (:math:`j`).

- ``'COVARIJ'``: The covariance value (:math:`\Sigma_{ij}`).

The table also contains the following metadata:

- ``'COVSHAPE'``: The shape of the covariance matrix.

- ``'BUNIT'``: (If defined) The string representation of the covariance units.

- ``'COVDSHP'``: (If the dimensionality is greater than 1) The shape of the
  associated data array.

For higher dimensional arrays, the coordinate data are automatically reshaped so
that the indices correspond to the data array.  For example,

.. doctest-requires:: scipy

    >>> data_array_shape
    (5, 2)
    >>> covar = Covariance(array=c, data_shape=data_array_shape)
    >>> tbl = covar.to_table()
    >>> tbl.meta
    {'COVSHAPE': '(10, 10)', 'COVDSHP': '(5, 2)'}
    >>> tbl[:3]
    <Table length=3>
     INDXI    INDXJ   COVARIJ
    int64[2] int64[2] float64
    -------- -------- -------
      0 .. 0   0 .. 0     1.0
      0 .. 0   0 .. 1     0.5
      0 .. 0   1 .. 0     0.2
    >>> tbl['INDXI'][0]
    array([0, 0])

.. warning::

    Recall that the storage of covariance matrices for higher
    dimensional data always assumes a row-major storage order.

The inverse operation is also provided to instantiate a
`~astropy.nddata.covariance.Covariance` object from a table.  Continuing the
N-dimensional example above:

.. doctest-requires:: scipy

    >>> _covar = Covariance.from_table(tbl)
    >>> _covar.data_shape
    (5, 2)
    >>> _covar.to_dense()
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

Use of the `~astropy.nddata.covariance.Covariance.to_table` and
`~astropy.nddata.covariance.Covariance.from_table` methods can be used with
Astropy's unified file I/O system to read and write the covariance matrices.

For example, to write the covariance matrix to table and reload it:

.. doctest-requires:: scipy

    >>> ofile = 'test_covar_io.fits'
    >>> covar = Covariance(array=c)
    >>> tbl = covar.to_table()
    >>> tbl.write(ofile, format='fits')
    >>> from astropy.io import fits
    >>> with fits.open(ofile) as hdu:
    ...     hdu.info()
    ...
    Filename: test_covar_io.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       4   ()
      1                1 BinTableHDU     15   27R x 3C   [K, K, D]
    >>> from astropy.table import Table
    >>> _tbl = Table.read(ofile, format='fits')
    >>> _covar = Covariance.from_table(_tbl)
    >>> print(np.array_equal(covar.to_dense(), _covar.to_dense()))
    True

Utility Functions
=================

.. _covariance-apply-new-variance:

Renormalizing the variance
--------------------------

To create a new covariance matrix that maintains the same correlations as an
existing matrix but a different variance, you can apply a new variance
normalization (following the examples in the :ref:`introductory section
<nddata-covariance-intro>`).  The `~astropy.nddata.covariance.Covariance` object
provides a convenience function for this.

.. doctest-requires:: scipy

    >>> covar_var1 = Covariance(array=c)
    >>> covar_var1.to_dense()
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
    >>> var4 = np.full(c.shape[0], 4.0, dtype=float)
    >>> covar_var4 = covar_var1.apply_new_variance(var4)
    >>> covar_var4.to_dense()
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

.. _covariance-match-to-data-slice:

Matching the covariance data to a slice of its parent data array
----------------------------------------------------------------

To adjust a `~astropy.nddata.covariance.Covariance` object so that it is
appropriate for a slice of its parent data array, use
`~astropy.nddata.covariance.Covariance.match_to_data_slice`.  For example, to
create a matrix with every other entry:

.. doctest-requires:: scipy

    >>> covar = Covariance(array=c)
    >>> sub_covar = covar.match_to_data_slice(np.s_[::2])
    >>> sub_covar
    <Covariance; shape = (5, 5)>
    >>> sub_covar.to_dense()
    array([[1. , 0.2, 0. , 0. , 0. ],
           [0.2, 1. , 0.2, 0. , 0. ],
           [0. , 0.2, 1. , 0.2, 0. ],
           [0. , 0. , 0.2, 1. , 0.2],
           [0. , 0. , 0. , 0.2, 1. ]])

or to adjust for a reordering of the parent data array:

.. doctest-requires:: scipy

    >>> covar = Covariance(array=c)
    >>> rng = np.random.default_rng(99)
    >>> reorder = np.arange(covar.shape[0])
    >>> rng.shuffle(reorder)
    >>> reorder
    array([4, 6, 0, 3, 8, 1, 2, 5, 7, 9])
    >>> reorder_covar = covar.match_to_data_slice(reorder)
    >>> reorder_covar.to_dense()
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

.. _covariance-nd-indexing:

Data-to-covariance Indexing Transformations
-------------------------------------------

For higher dimensional arrays, two methods are provided to ease conversion
between data array and covariance matrix indexing.  Following examples above,
define the ten elements in the covariance matrix as coming from a :math:`5
\times 2` array, then find the indices in the data array for the covariance
values at indices covariance values at matrix locations ``(0,3)``, ``(1,4)``,
and ``(2,3)``:

.. doctest-requires:: scipy

    >>> covar = Covariance(array=c, data_shape=data_array_shape)
    >>> i_data, j_data = covar.covariance_to_data_indices([0,1,2], [3,4,3])
    >>> i_data
    (array([0, 0, 1]), array([0, 1, 0]))
    >>> j_data
    (array([1, 2, 1]), array([1, 0, 1]))

This shows that the covariance elements provide the covariance between
``data[0,0]`` and ``data[1,1]``, elements ``data[0,1]`` and ``data[2,0]``, and
elements ``data[1,0]`` and ``data[1,1]``.

The inverse operation gives the covariance indices for a specified set of
data-array indices.  Keeping the indices we defined above:

.. doctest-requires:: scipy

    >>> i_cov, j_cov = covar.data_to_covariance_indices(i_data, j_data)
    >>> i_cov, j_cov
    (array([0, 1, 2]), array([3, 4, 3]))

