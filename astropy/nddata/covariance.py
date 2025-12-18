# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Implements a class used to store and manipulate covariance matrices.
"""

import warnings

import numpy as np

from astropy import table
from astropy.units import Quantity
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

from .nddata import NDUncertainty

# Required scipy.sparse imports
if HAS_SCIPY:
    from scipy.sparse import coo_matrix, csr_matrix, find, isspmatrix_csr, triu
else:
    err = "Use of 'astropy.nddata.Covariance' requires 'scipy.sparse' module!"

    def find(*args, **kwargs):
        raise ModuleNotFoundError(err)

    def triu(*args, **kwargs):
        raise ModuleNotFoundError(err)

    def csr_matrix(*args, **kwargs):
        raise ModuleNotFoundError(err)

    def coo_matrix(*args, **kwargs):
        raise ModuleNotFoundError(err)

    def isspmatrix_csr(*args, **kwargs):
        raise ModuleNotFoundError(err)


__all__ = ["Covariance"]


# Disabling doctests when scipy isn't present.
__doctest_requires__ = {"*": ["scipy"]}


def _get_csr(arr):
    """
    Helper method used to confirm the array is a `~scipy.sparse.csr_matrix` or
    try to convert it to one.

    Parameters
    ----------
    arr : array-like
        An array that either is a `~scipy.sparse.csr_matrix` or can be converted
        to one.

    Returns
    -------
    `~scipy.sparse.csr_matrix`
        Converted or original matrix.
    """
    # Check the input type
    if isspmatrix_csr(arr):
        return arr
    try:
        return csr_matrix(arr)
    except ValueError:
        raise TypeError(
            "Input matrix is not a scipy.sparse.csr_matrix and could "
            "not be converted to one."
        )


def _impose_sparse_value_threshold(arr, threshold):
    """
    Remove values from a sparse matrix if their absolute value is below the
    provided threshold.

    Parameters
    ----------
    arr : `~scipy.sparse.csr_matrix`
        Array to manipulate.

    threshold : :obj:`float`
        Threshold value

    Returns
    -------
    `~scipy.sparse.csr_matrix`
        Manipulated or original matrix.
    """
    i, j, aij = find(arr)
    index = np.logical_not(np.absolute(aij) < threshold)
    if all(index):
        return arr
    return coo_matrix((aij[index], (i[index], j[index])), shape=arr.shape).tocsr()


def _parse_shape(shape):
    """
    Parse a string representation of an array shape into a tuple.

    Parameters
    ----------
    shape : str
        String representation of the tuple.  It should only contain
        comma-separated integers and parentheses.

    Returns
    -------
    tuple
        Tuple with the shape.
    """
    return tuple([int(n) for n in shape.strip("()").split(",") if len(n) > 0])


class Covariance(NDUncertainty):
    r"""
    A general utility for storing, manipulating, and I/O of covariance matrices.

    Covariance matrices are symmetric by definition, :math:`\Sigma_{ij} =
    \Sigma_{ji}`.  The object therefore only stores the upper triangle of the
    matrix using a `scipy.sparse.csr_matrix`.  By default, instantiation will
    check for symmetry and issue a warning if the matrix is not symmetric.  This
    check can be skipped using the ``assume_symmetric`` keyword.  However, by
    virtue of how the data is stored, symmetry is *always imposed* on the
    matrix.  That is, if a non-symmetric matrix is used to instantiate a
    `Covariance` object, the stored data will yield a matrix that is different
    from the original input.

    Covariance matrices of higher dimensional arrays are always assumed to be
    stored following row-major indexing.  For example, the covariance value
    :math:`\Sigma_{ij}` for an image of size :math:`(N_x,N_y)` is the covariance
    between image pixels :math:`I_{x_i,y_i}` and :math:`I_{x_j,y_j}`, where
    :math:`i = x_i + N_x y_i` and, similarly, :math:`j = x_j + N_x y_j`.

    See :ref:`nddata-covariance` for additional documentation and examples.

    Parameters
    ----------
    array : array-like, `~scipy.sparse.csr_matrix`
        Covariance matrix to store.  If the array is not a
        `~scipy.sparse.csr_matrix` instance, it must be convertible to one.  To
        match the calling sequence for `NDUncertainty`, ``array`` has a default
        value of None, but the array *must* be provided for this `Covariance`
        object.

    data_shape : :obj:`tuple`, optional
        The covariance data is for a higher dimensional array with this shape.
        For example, if the covariance data is for a 2D image with shape
        ``(nx,ny)``, set ``data_shape=(nx,ny)``; the shape of the covariance
        array must then be ``(nx*ny, nx*ny)``.  If None, any higher
        dimensionality is ignored.

    assume_symmetric : bool, optional
        Assume the matrix is symmetric.  This means that a check for symmetry is
        not performed, and the user is not warned if the matrix is not
        symmetric.

    unit : unit-like, optional
        Unit for the covariance values.

    Raises
    ------
    TypeError
        Raised if the input array not a `~scipy.sparse.csr_matrix` object and
        cannot be converted to one.

    ValueError
        Raised if ``data_shape`` is provided and the input covariance matrix
        ``array`` does not have the expected shape or if ``array`` is None.
    """

    def __init__(self, array=None, data_shape=None, assume_symmetric=False, unit=None):
        if array is None:
            raise ValueError("Covariance object cannot be instantiated with None.")

        # Ingest the matrix
        self._cov = triu(
            Covariance._ingest_matrix(array, assume_symmetric=assume_symmetric)
        )
        # Save the diagonal as a variance array for convenience
        self._var = self._cov.diagonal()

        # Set the shape and check it; note self._cov must be defined so that
        # call to self.shape below is valid.
        self._data_shape = data_shape
        if self._data_shape is not None and np.prod(self._data_shape) != self.shape[0]:
            raise ValueError(
                "Product of ``data_shape`` must match the covariance axis length."
            )

        # Workspace for index mapping from flattened to original data arrays
        self._data_index_map = None

        super().__init__(array=self._cov, copy=False, unit=unit)

    @staticmethod
    def _ingest_matrix(arr, assume_symmetric=False):
        """
        Helper method to ingest a covariance or correlation matrix.

        This function converts the input to a '~scipy.sparse.csr_matrix` using
        :func:`_get_csr`, and checks that the array is 2D, square, and symmetric.

        Parameters
        ----------
        arr : array-like
            An array that either is a `~scipy.sparse.csr_matrix` or can be converted
            to one.

        assume_symmetric : bool, optional
            Assume the matrix is symmetric.  This means that a check for
            symmetry is not performed, and the user is not warned if the matrix
            is not symmetric.

        Returns
        -------
        `~scipy.sparse.csr_matrix`
            Converted or original matrix.
        """
        # Make sure it's a sparse matrix or can be converted to one.
        _arr = _get_csr(arr)

        # Check that it's 2D
        if _arr.ndim != 2:
            raise ValueError("Covariance arrays must be 2-dimensional.")
        # Check that it's square
        if _arr.shape[0] != _arr.shape[1]:
            raise ValueError("Covariance matrices must be square.")

        # Skip the symmetry check, if requested
        if assume_symmetric:
            return _arr

        # Check that it's symmetric
        flip_diff = _arr - _arr.T
        if not np.allclose(flip_diff.data, np.zeros_like(flip_diff.data)):
            warnings.warn(
                "Asymmetry detected in covariance/correlation matrix.  Matrix will be modified "
                "to be symmetric using its upper triangle.",
                AstropyUserWarning,
            )
            _arr = triu(_arr) + triu(_arr, 1).T
        return _arr

    @property
    def shape(self):
        """Tuple with the shape of the covariance matrix"""
        return self._cov.shape

    @property
    def nnz(self):
        """
        The number of non-zero (NNZ) elements in the full covariance matrix,
        *including* both the upper and lower triangles.
        """
        return self.stored_nnz * 2 - self._cov.shape[0]

    @property
    def stored_nnz(self):
        """
        The number of non-zero elements stored by the object, which only
        counts the non-zero elements in the upper triangle.
        """
        return self._cov.nnz

    @property
    def variance(self):
        """
        The diagonal of the covariance matrix.
        """
        return self._var

    @variance.setter
    def variance(self, value):
        raise NotImplementedError(
            "Directly setting variance values is not allowed for Covariance objects."
        )

    @property
    def uncertainty_type(self):
        """``"cov"``: `Covariance` implements a covariance matrix."""
        return "cov"

    @property
    def quantity(self):
        """
        The covariance matrix as an dense `~astropy.units.Quantity` object.
        """
        return Quantity(self.to_dense(), self.unit, copy=False, dtype=self._cov.dtype)

    def _data_unit_to_uncertainty_unit(self, value):
        """
        Return the uncertainty unit for covariances given the data unit.
        """
        return value**2

    def __repr__(self):
        return f"<{self.__class__.__name__}; shape = {self.shape}>"

    # Skip error propagation for now
    def _propagate_add(self, other_uncert, result_data, correlation):
        return None

    def _propagate_subtract(self, other_uncert, result_data, correlation):
        return None

    def _propagate_multiply(self, other_uncert, result_data, correlation):
        return None

    def _propagate_divide(self, other_uncert, result_data, correlation):
        return None

    @classmethod
    def from_samples(cls, samples, cov_tol=None, rho_tol=None, **kwargs):
        r"""
        Build a covariance object using discrete samples.

        The covariance is generated using `~numpy.cov` for a set of discretely
        sampled data for an :math:`N`-dimensional parameter space.

        Parameters
        ----------
        samples : `~numpy.ndarray`
            Array with samples drawn from an :math:`N`-dimensional parameter
            space. The shape of the input array must be :math:`N_{\rm par}\times
            N_{\rm samples}`.

        cov_tol : :obj:`float`, optional
            The absolute value of any *covariance matrix* entry less than this
            is assumed to be equivalent to (and set to) 0.

        rho_tol : :obj:`float`, optional
            The absolute value of any *correlation coefficient* less than this
            is assumed to be equivalent to (and set to) 0.

        **kwargs : dict, optional
            Passed directly to main instantiation method.

        Returns
        -------
        `Covariance`
            An :math:`N_{\rm par}\times N_{\rm par}` covariance matrix built
            using the provided samples.

        Raises
        ------
        ValueError
            Raised if the input array is not 2D or if the number of samples (length
            of the second axis) is less than 2.
        """
        if samples.ndim != 2:
            raise ValueError("Input samples for covariance matrix must be a 2D array!")
        if samples.shape[1] < 2:
            raise ValueError("Fewer than two samples provided!")
        return Covariance.from_array(
            np.cov(samples), cov_tol=cov_tol, rho_tol=rho_tol, **kwargs
        )

    @classmethod
    def from_array(cls, covar, cov_tol=None, rho_tol=None, **kwargs):
        r"""
        Define a covariance object from an array.

        .. note::

            The only difference between this method and the direct instantiation
            method (i.e., ``Covariance(array=covar)``) is that it can be used to
            impose tolerances on the covariance value and/or correlation
            coefficients.

        Parameters
        ----------
        covar : array-like
            Array with the covariance data. The object must be either a
            `~scipy.sparse.csr_matrix` or an object that can be converted to
            one.  It must also be 2-dimensional and square.

        cov_tol : :obj:`float`, optional
            The absolute value of any *covariance* matrix entry less than this
            is assumed to be equivalent to (and set to) 0.

        rho_tol : :obj:`float`, optional
            The absolute value of any *correlation coefficient* less than this
            is assumed to be equivalent to (and set to) 0.

        **kwargs : dict, optional
            Passed directly to main instantiation method.

        Returns
        -------
        `Covariance`
            The covariance matrix built using the provided array.
        """
        # Get the assume_symmetric flag, either from kwargs or as the default
        # value
        assume_symmetric = kwargs.pop("assume_symmetric", False)
        # Convert the covariance to a correlation matrix.  If rho_tol is None,
        # this just serves to symmetrize the matrix if it's not already.  Set
        # assume_symmetric to True hereafter
        var, rho = Covariance.to_correlation(covar, assume_symmetric=assume_symmetric)
        if rho_tol is not None:
            rho = _impose_sparse_value_threshold(rho, rho_tol)
        _covar = Covariance.revert_correlation(var, rho, assume_symmetric=True)
        if cov_tol is not None:
            _covar = _impose_sparse_value_threshold(_covar, cov_tol)
        return cls(array=_covar, assume_symmetric=True, **kwargs)

    @classmethod
    def from_table(cls, triu_covar):
        r"""
        Construct the covariance matrix from a table with the non-zero elements
        of the upper triangle of the covariance matrix in coordinate format.

        This is the inverse operation of :func:`to_table`.  The class can read
        covariance data written by other programs *as long as they have a
        commensurate format*; see :func:`to_table`.

        Parameters
        ----------
        triu_covar : `~astropy.table.Table`
            The non-zero elements of the upper triangle of the covariance matrix
            in coordinate format; see :func:`to_table`.

        Returns
        -------
        `Covariance`
            The covariance matrix constructed from the tabulated data.

        Raises
        ------
        ValueError
            Raised if ``triu_covar.meta`` is ``None``, if the provided variance array
            does not have the correct size, or if the data is multidimensional
            and the table columns do not have the right shape.
        """
        # Read shapes
        if "COVSHAPE" not in triu_covar.meta:
            raise ValueError("Table meta dictionary *must* contain COVSHAPE")

        shape = _parse_shape(triu_covar.meta["COVSHAPE"])
        data_shape = (
            _parse_shape(triu_covar.meta["COVDSHP"])
            if "COVDSHP" in triu_covar.meta
            else None
        )

        # Number of non-zero elements
        nnz = len(triu_covar)

        # Read coordinate data
        # WARNING: If the data is written correctly, it should always be true that i<=j
        if data_shape is None:
            i = triu_covar["INDXI"].data
            j = triu_covar["INDXJ"].data
        else:
            ndim = triu_covar["INDXI"].shape[1]
            if len(data_shape) != ndim:
                raise ValueError("Mismatch between COVDSHP keyword and tabulated data.")
            i = np.ravel_multi_index(triu_covar["INDXI"].data.T, data_shape)
            j = np.ravel_multi_index(triu_covar["INDXJ"].data.T, data_shape)

        # Units
        unit = triu_covar.meta.get("BUNIT", None)

        # Set covariance data
        cij = triu_covar["COVARIJ"].data
        # NOTE: the astype conversion of cij when instantiating the matrix below
        # is because of how scipy.sparse restricts instantiation of sparse
        # arrays.  It doesn't like big-endian byte order.  To reproduce the
        # underlying error:
        #   >>> import numpy as np
        #   >>> from scipy.sparse._sputils import getdtype
        #   >>> getdtype(np.dtype('>f8'))
        #   Traceback (most recent call last):
        #   File "<stdin>", line 1, in <module>
        #   File "/Users/westfall/.virtualenvs/astropy/lib/python3.12/site-packages/scipy/sparse/_sputils.py", line 137, in getdtype
        #       raise ValueError(f"scipy.sparse does not support dtype {newdtype.name}. "
        #   ValueError: scipy.sparse does not support dtype float64. The only supported types are: bool, int8, uint8, int16, uint16, int32, uint32, int64, uint64, longlong, ulonglong, float32, float64, longdouble, complex64, complex128, clongdouble.
        #   >>> getdtype(np.dtype('<f8'))
        #   dtype('float64')
        cov = coo_matrix((cij.astype(cij.dtype.type), (i, j)), shape=shape).tocsr()

        # Instantiate.  Set assume_symmetric to true to avoid the warning from
        # the _ingest_matrix method
        return cls(array=cov, data_shape=data_shape, unit=unit, assume_symmetric=True)

    @classmethod
    def from_matrix_multiplication(cls, T, covar, **kwargs):
        r"""
        Construct the covariance matrix that results from a matrix
        multiplication.

        Linear operations on a dataset (e.g., binning or smoothing) can be
        written as matrix multiplications of the form

        .. math::

            {\mathbf y} = {\mathbf T}\ {\mathbf x},

        where :math:`{\mathbf T}` is a transfer matrix of size :math:`N_y\times
        N_x`, :math:`{\mathbf x}` is a vector of size :math:`N_x`, and
        :math:`{\mathbf y}` is a vector of length :math:`{N_y}` that results
        from the multiplication.  If :math:`{\mathbf \Sigma}_x` is the
        covariance matrix for :math:`{\mathbf x}`, then the covariance matrix
        for :math:`{\mathbf Y}` is

        .. math::

            {\mathbf \Sigma}_y = {\mathbf T}\ {\mathbf \Sigma}_x\
            {\mathbf T}^\top.

        If ``covar`` is provided as a vector of length :math:`N_x`, it is
        assumed that the elements of :math:`{\mathbf X}` are independent and the
        provided vector gives the *variance* in each element; i.e., the provided
        data represent the diagonal of :math:`{\mathbf \Sigma}`.

        Parameters
        ----------
        T : `~scipy.sparse.csr_matrix`, `~numpy.ndarray`
            Transfer matrix.  See above.

        covar : `~scipy.sparse.csr_matrix`, `~numpy.ndarray`
            Covariance matrix.  See above.

        **kwargs : dict, optional
            Passed directly to main instantiation method.

        Returns
        -------
        `Covariance`
            The covariance matrix resulting from the matrix multiplication.

        Raises
        ------
        ValueError
            Raised if the provided arrays are not two dimensional or if there is
            a shape mismatch.
        """
        if T.ndim != 2:
            raise ValueError("Input transfer matrix must be two-dimensional.")
        nx = T.shape[1]
        if covar.shape != (nx, nx) and covar.shape != (nx,):
            raise ValueError(
                f"Shape of input variance matrix must be either ({nx}, {nx}) or ({nx},)."
            )
        # If it isn't already, convert T to a csr_matrix
        _T = T if isinstance(T, csr_matrix) else csr_matrix(T)
        # Set the covariance matrix in X
        _covar = (
            coo_matrix((covar, (np.arange(nx), np.arange(nx))), shape=(nx, nx)).tocsr()
            if covar.ndim == 1
            else (covar if isinstance(covar, csr_matrix) else csr_matrix(covar))
        )
        # Construct the covariance matrix
        return cls(_T.dot(_covar.dot(_T.transpose())).tocsr(), **kwargs)

    @classmethod
    def from_variance(cls, variance, **kwargs):
        """
        Construct a diagonal covariance matrix using the provided variance.

        Parameters
        ----------
        variance : `~numpy.ndarray`
            The variance vector.

        **kwargs : dict, optional
            Passed directly to main instantiation method.

        Returns
        -------
        `Covariance`
            The diagonal covariance matrix.
        """
        return cls(csr_matrix(np.diagflat(variance)), **kwargs)

    def to_sparse(self, correlation=False):
        """
        Return the full covariance matrix as a `~scipy.sparse.csr_matrix`
        object.

        This method is essentially equivalent to `to_dense` except that it
        returns a sparse array.

        Parameters
        ----------
        correlation : :obj:`bool`, optional
            Return the *correlation* matrix.  If False, return the covariance
            matrix.

        Returns
        -------
        `~scipy.sparse.csr_matrix`
            The sparse matrix with both the upper and lower triangles filled
            (with symmetric information).
        """
        cov = triu(self._cov) + triu(self._cov, 1).T
        if not correlation:
            return cov
        return Covariance.to_correlation(cov, assume_symmetric=True)[1]

    def apply_new_variance(self, var):
        """
        Using the same correlation coefficients, return a new `Covariance`
        object with the provided variance.

        Parameters
        ----------
        var : array-like
            Variance vector. Must have a length that matches this `Covariance`
            instance; e.g., if this instance is ``cov``, the length of ``var``
            must be ``cov.shape[0]``).  Note that, if the covariance is for
            higher dimensional data, this variance array *must* be flattened to
            1D.

        Returns
        -------
        `Covariance`
            A covariance matrix with the same shape and correlation coefficients
            and this object, but with the provided variance.

        Raises
        ------
        ValueError
            Raised if the length of the variance vector is incorrect.
        """
        _var = np.asarray(var)
        if _var.shape != self._var.shape:
            raise ValueError(
                f"Provided variance has incorrect shape.  Expected {self._var.shape}, "
                f"found {_var.shape}."
            )
        i, j, cij = find(self._cov)
        _cov = coo_matrix(
            (cij * np.sqrt(_var[i] / self._var[i] * _var[j] / self._var[j]), (i, j)),
            shape=self.shape,
        ).tocsr()
        return Covariance(
            array=_cov, data_shape=self._data_shape, assume_symmetric=True
        )

    def copy(self):
        """
        Return a copy of this Covariance object.

        Returns
        -------
        `Covariance`
            A copy of the current covariance matrix.
        """
        # Create the new Covariance instance with a copy of the data
        return Covariance(
            array=self._cov.copy(),
            data_shape=self._data_shape,
            assume_symmetric=True,
            unit=self.unit,
        )

    def to_dense(self, correlation=False):
        """
        Return the full covariance matrix as a `numpy.ndarray` object (a "dense"
        array).

        Parameters
        ----------
        correlation : bool, optional
            Flag to return the correlation matrix, instead of the covariance
            matrix.  Note that setting this to ``True`` does *not* also return the
            variance vector.

        Returns
        -------
        `~numpy.ndarray`
            Dense array with the full covariance matrix.
        """
        return self.to_sparse(correlation=correlation).toarray()

    def find(self, correlation=False):
        """
        Find the non-zero values in the **full** covariance matrix (not just the
        upper triangle).

        This is a simple wrapper for `to_sparse` and `~scipy.sparse.find`.

        Parameters
        ----------
        correlation : bool, optional
            Flag to return the correlation data, instead of the covariance data.
            Note that setting this to ``True`` does *not* also return the variance
            vector.

        Returns
        -------
        i, j : `numpy.ndarray`
            Arrays containing the index coordinates of the non-zero values in
            the covariance (or correlation) matrix.
        c : `numpy.ndarray`
            The non-zero covariance (or correlation) matrix values located at
            the provided ``i,j`` coordinates.
        """
        return find(self.to_sparse(correlation=correlation))

    def covariance_to_data_indices(self, i, j):
        r"""
        Given indices along the two axes of the covariance matrix, return the
        relevant indices in the data array.  This is the inverse of
        :func:`data_to_covariance_indices`.

        Parameters
        ----------
        i : `~numpy.ndarray`
            1D array with the index along the first axis of the covariance
            matrix.  Must be in the range :math:`0...n-1`, where :math:`n` is
            the length of the covariance-matrix axes.

        j : `~numpy.ndarray`
            1D array with the index along the second axis of the covariance
            matrix.  Must be in the range :math:`0...n-1`, where :math:`n` is
            the length of the covariance-matrix axes.

        Returns
        -------
        i_data, i_data : tuple, `numpy.ndarray`
            If `data_shape` is not defined, the input arrays are simply returned
            (and not copied).  Otherwise, the code uses `~numpy.unravel_index`
            to calculate the relevant data-array indices; each element in the
            two-tuple is itself a tuple of :math:`N_{\rm dim}` arrays, one array
            per dimension of the data array.

        Raises
        ------
        ValueError
            Raised if the provided indices fall outside the range of covariance
            matrix.

        """
        if self._data_shape is None:
            if np.any(
                (i < 0) | (i > self.shape[0] - 1) | (j < 0) | (j > self.shape[1] - 1)
            ):
                raise ValueError(
                    "Some indices not valid for covariance matrix with shape "
                    f"{self.shape}."
                )
            return i, j
        return np.unravel_index(
            np.atleast_1d(i).ravel(), self._data_shape
        ), np.unravel_index(np.atleast_1d(j).ravel(), self._data_shape)

    def data_to_covariance_indices(self, i, j):
        r"""
        Given indices of elements in the source data array, return the matrix
        coordinates with the associated covariance.  This is the inverse of
        :func:`covariance_to_data_indices`.

        Parameters
        ----------
        i : array-like, tuple
            A tuple of :math:`N_{\rm dim}` array-like objects providing the
            indices of elements in the N-dimensional data array.  This can be an
            array-like object if ``data_shape`` is undefined, in which case the
            values must be in the range :math:`0...n-1`, where :math:`n` is the
            length of the data array.

        j : array-like, tuple
            The same as ``i``, but providing a second set of coordinates at which
            to access the covariance.

        Returns
        -------
        i_covar, j_covar : `numpy.ndarray`
            Arrays providing the indices in the covariance matrix associated
            with the provided data array coordinates.  If ``data_shape`` is not
            defined, the input arrays are simply returned (and not copied).
            Otherwise, the code uses `~numpy.ravel_multi_index` to calculate the
            relevant covariance indices.

        Raises
        ------
        ValueError
            Raised if the provided indices fall outside the range of data array,
            or if the length of the ``i`` or ``j`` tuples is not :math:`N_{\rm
            dim}`.

        """
        if self._data_shape is None:
            if np.any(
                (i < 0) | (i > self.shape[0] - 1) | (j < 0) | (j > self.shape[1] - 1)
            ):
                raise ValueError(
                    "Some indices not valid for covariance matrix with shape "
                    f"{self.shape}."
                )
            return i, j
        if len(i) != len(self.data_shape):
            raise ValueError(
                "Length of input coordinate list (i) is incorrect; expected "
                f"{len(self.data_shape)}, found {len(i)}"
            )
        if len(j) != len(self.data_shape):
            raise ValueError(
                "Length of input coordinate list (j) is incorrect; expected "
                f"{len(self.data_shape)}, found {len(i)}"
            )
        return np.ravel_multi_index(i, self.data_shape), np.ravel_multi_index(
            j, self.data_shape
        )

    def coordinate_data(self, reshape=False):
        r"""
        Construct data arrays with the non-zero covariance components in
        coordinate format.

        Coordinate format means that the covariance matrix data is provided in
        three columns providing :math:`\Sigma_{ij}` and the (0-indexed) matrix
        coordinates :math:`i,j`.

        This procedure is primarily used when constructing the data arrays for
        storage.  Matching the class convention, the returned data only includes
        the upper triangle.

        Parameters
        ----------
        reshape : :obj:`bool`, optional
            If ``reshape`` is ``True`` and `data_shape` is defined, the :math:`i,j`
            indices are converted to the expected data-array indices; see
            :func:`covariance_to_data_indices`.  These can be reverted to the
            coordinates in the covariance matrix using
            :func:`data_to_covariance_indices`.

        Returns
        -------
        i, j : tuple, `numpy.ndarray`
            The row and column indices, :math:`i,j`: of the covariance matrix.
            If reshaping, these are tuples with the index arrays along each of
            the reshaped axes.
        cij : `numpy.ndarray`
            The covariance, :math:`\Sigma_{ij}`, between array elements at
            indices :math:`i` and :math:`j`.

        Raises
        ------
        ValueError
            Raised if ``reshape`` is True but `data_shape` is undefined.
        """
        if reshape and self._data_shape is None:
            raise ValueError(
                "If reshaping, the shape of the data before flattening to the "
                "covariance array (``data_shape``) must be defined."
            )

        # Get the data (only stores the upper triangle!)
        i, j, cij = find(self._cov)

        # Return the data.
        if reshape:
            # Reshape the indices and the variance array.
            return (
                np.unravel_index(i, self._data_shape),
                np.unravel_index(j, self._data_shape),
                cij,
            )
        return i, j, cij

    def to_table(self):
        r"""
        Return the covariance data in a `~astropy.table.Table` using coordinate
        format.

        Coordinate format means that the covariance matrix data is provided in
        three columns providing :math:`\Sigma_{ij}` and the (0-indexed) matrix
        coordinates :math:`i,j`.

        The output table has three columns:

        - ``'INDXI'``: The row index in the covariance matrix.

        - ``'INDXJ'``: The column index in the covariance matrix.

        - ``'COVARIJ'``: The covariance at the relevant :math:`i,j` coordinate.

        The table also contains the following metadata:

        - ``'COVSHAPE'``: The shape of the covariance matrix

        - ``'BUNIT'``: (If ``unit`` is defined) The string representation of the
          covariance units.

        - ``'COVDSHP'``: (If ``data_shape`` is defined) The shape of the
          associated data array.

        If ``data_shape`` is set, the covariance matrix indices are reformatted
        to match the coordinates in the N-dimensional array.

        .. warning::

            Recall that the storage of covariance matrices for higher
            dimensional data always assumes a row-major storage order.

        Objects instantiated by this method can be used to re-instantiate the
        `Covariance` object using `from_table`.

        Returns
        -------
        `~astropy.table.Table`
            Table with the covoariance matrix in coordinate format and the
            relevant metadata.
        """
        meta = {}
        meta["COVSHAPE"] = str(self.shape)
        if self.unit is not None:
            meta["BUNIT"] = self.unit.to_string()
        reshape = self._data_shape is not None
        i, j, cij = self.coordinate_data(reshape=reshape)
        triu_nnz = cij.size
        if reshape:
            meta["COVDSHP"] = str(self._data_shape)
            i = np.column_stack(i)
            j = np.column_stack(j)
            coo_shape = (i.shape[1],)
        else:
            coo_shape = None
        return table.Table(
            [
                table.Column(
                    data=i, name="INDXI", dtype=int, length=triu_nnz, shape=coo_shape
                ),
                table.Column(
                    data=j, name="INDXJ", dtype=int, length=triu_nnz, shape=coo_shape
                ),
                table.Column(data=cij, name="COVARIJ", dtype=float, length=triu_nnz),
            ],
            meta=meta,
        )

    @property
    def data_shape(self):
        """
        The expected shape of the data array associated with this covariance array.
        """
        return (self.shape[0],) if self._data_shape is None else self._data_shape

    @property
    def data_index_map(self):
        """
        An array mapping the index along each axis of the covariance matrix to
        the shape of the associated data array.
        """
        if self._data_index_map is None:
            self._data_index_map = np.arange(self.shape[0])
            if self._data_shape is not None:
                self._data_index_map = self._data_index_map.reshape(self._data_shape)
        return self._data_index_map

    def match_to_data_slice(self, data_slice):
        """
        Return a new `Covariance` instance that is matched to a slice of its
        parent data array.

        Parameters
        ----------
        data_slice : slice, array-like
            Anything that can be used to slice a `numpy.ndarray`.  To generate a
            slice using syntax that mimics accessing numpy array elements, use
            `numpy.s_`; see examples
            :ref:`here<covariance-match-to-data-slice>`.

        Returns
        -------
        `Covariance`
            A new covariance object for the sliced data array.
        """
        remap = self.data_index_map[data_slice]
        index = remap.ravel()
        return Covariance(
            self.to_sparse()[np.ix_(index, index)],
            data_shape=None if len(remap.shape) == 1 else remap.shape,
        )

    @staticmethod
    def to_correlation(cov, assume_symmetric=False):
        r"""
        Convert a covariance matrix into a correlation matrix by dividing each
        element by the variances.

        Specifically, extract ``var`` (:math:`V_i = C_{ii} \equiv \sigma^2_i`)
        and convert ``cov`` from a covariance matrix with elements
        :math:`C_{ij}` to a correlation matrix with :math:`\rho_{ij}` such that

        .. math::

            C_{ij} \equiv \rho_{ij} \sigma_i \sigma_j.

        To revert a variance vector and correlation matrix back to a covariance
        matrix, use :func:`revert_correlation`.

        Parameters
        ----------
        cov : array-like
            Covariance matrix to convert.  Must be a `~scipy.sparse.csr_matrix`
            instance or convertible to one.

        assume_symmetric : bool, optional
            Assume the matrix is symmetric.  This means that a check for
            symmetry is not performed, and the user is not warned if the matrix
            is not symmetric.

        Returns
        -------
        var : `numpy.ndarray`
            Variance vector
        rho : `~scipy.sparse.csr_matrix`
            Correlation matrix

        Raises
        ------
        ValueError
            Raised if the input array is not 2D and square.
        """
        # Ingest the matrix
        _cov = Covariance._ingest_matrix(cov, assume_symmetric=assume_symmetric)
        # Save the diagonal
        var = _cov.diagonal()
        # Find all the non-zero elements
        i, j, cij = find(_cov)
        rho = coo_matrix(
            (cij / np.sqrt(var[i] * var[j]), (i, j)), shape=_cov.shape
        ).tocsr()
        return var, rho

    @staticmethod
    def revert_correlation(var, rho, assume_symmetric=False):
        r"""
        Revert a variance vector and correlation matrix into a covariance matrix.

        This is the reverse operation of `to_correlation`.

        Parameters
        ----------
        var : `~numpy.ndarray`
            Variance vector.  Length must match the diagonal of ``rho``.

        rho : `~numpy.ndarray`, `~scipy.sparse.csr_matrix`
            Correlation matrix.  Diagonal must have the same length as ``var``.

        assume_symmetric : bool, optional
            Assume the matrix is symmetric.  This means that a check for
            symmetry is not performed, and the user is not warned if the matrix
            is not symmetric.

        Returns
        -------
        `~scipy.sparse.csr_matrix`
            Covariance matrix.
        """
        i, j, rhoij = find(
            Covariance._ingest_matrix(rho, assume_symmetric=assume_symmetric)
        )
        return coo_matrix(
            (rhoij * np.sqrt(var[i] * var[j]), (i, j)), shape=rho.shape
        ).tocsr()
