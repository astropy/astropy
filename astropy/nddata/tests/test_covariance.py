from pathlib import Path

import numpy as np
import pytest

from astropy import units
from astropy.nddata import covariance
from astropy.table import Table
from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

scipy_required = pytest.mark.skipif(not HAS_SCIPY, reason="scipy not installed")


def mock_cov():
    # Build a bogus covariance matrix
    return (
        np.diag(np.full(10 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(10, 1.0, dtype=float), k=0)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(10 - 2, 0.2, dtype=float), k=2)
    )


def mock_cov_2d():
    # Build a bogus covariance matrix that is for a 3,2 2D array
    c = (
        np.diag(np.full(6 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(6 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(6, 1.0, dtype=float), k=0)
        + np.diag(np.full(6 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(6 - 2, 0.2, dtype=float), k=2)
    )
    return (3, 2), c


def mock_cov_3d():
    # Build a bogus covariance matrix that is for a 3,2,3 3D array
    c = (
        np.diag(np.full(18 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(18 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(18, 1.0, dtype=float), k=0)
        + np.diag(np.full(18 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(18 - 2, 0.2, dtype=float), k=2)
    )
    return (3, 2, 3), c


def test_scipy_funcs():
    # Functions should always be callable, regardless of whether or not scipy is installed
    assert callable(covariance.find), "find should be callable"
    assert callable(covariance.triu), "triu should be callable"
    assert callable(covariance.csr_matrix), "csr_matrix should be callable"
    assert callable(covariance.coo_matrix), "coo_matrix should be callable"
    assert callable(covariance.isspmatrix_csr), "isspmatrix_csr should be callable"

    if not HAS_SCIPY:
        # Should raise an error if scipy is not available
        with pytest.raises(ModuleNotFoundError):
            covariance.find()
        with pytest.raises(ModuleNotFoundError):
            covariance.triu()
        with pytest.raises(ModuleNotFoundError):
            covariance.csr_matrix()
        with pytest.raises(ModuleNotFoundError):
            covariance.coo_matrix()
        with pytest.raises(ModuleNotFoundError):
            covariance.isspmatrix_csr()


@scipy_required
def test_get_csr():
    c = mock_cov()
    # Convert a numpy array
    carr = covariance._get_csr(c)
    assert np.array_equal(carr.toarray(), c), "np.array not converted properly"
    # Convert a list
    clist = covariance._get_csr(c.tolist())
    assert np.array_equal(clist.toarray(), c), "List not converted properly"
    # Return input array if already a csr_matrix
    ccsr = covariance.csr_matrix(c)
    _ccsr = covariance._get_csr(ccsr)
    assert _ccsr is ccsr, "Should return reference to the input"
    # Fail if cannot be converted
    with pytest.raises(TypeError):
        cov = covariance.Covariance(array="test")


@scipy_required
def test_impose_sparse_value_threshold():
    c = covariance.csr_matrix(mock_cov())

    _c = covariance._impose_sparse_value_threshold(c, 0.0)
    assert c is _c, "All elements are above 0, so array should be returned exactly."

    _c = covariance._impose_sparse_value_threshold(c, 0.21)
    assert c is not _c, "Should remove elements"
    assert _c.nnz + 16 == c.nnz, "Incorrect number of elements removed."

    _c = covariance._impose_sparse_value_threshold(c, 0.51)
    assert c is not _c, "Should remove elements"
    assert _c.nnz == _c.shape[0], "Should remove all but the diagonal"


def test_parse_shape():
    shape = "(2,1)"
    assert covariance._parse_shape(shape) == (2, 1), f"Bad parse for {shape}"
    shape = "(2,1,5)"
    assert covariance._parse_shape(shape) == (2, 1, 5), f"Bad parse for {shape}"
    shape = "(2,1,5,)"
    assert covariance._parse_shape(shape) == (2, 1, 5), f"Bad parse for {shape}"
    shape = "(2,  1,  5,)"
    assert covariance._parse_shape(shape) == (2, 1, 5), f"Bad parse for {shape}"


@scipy_required
def test_ingest():
    # It must be possible to convert the input array into a csr_matrix
    with pytest.raises(TypeError):
        cov = covariance.Covariance._ingest_matrix("test")

    # Must be 2D.  If 1D, the array is successfully converted to a csr_matrix,
    # but fails the 2D check.  So the error is ValueError.
    arr1d = np.ones(27, dtype=float)
    with pytest.raises(ValueError):
        cov = covariance.Covariance._ingest_matrix(arr1d)

    # Must be 2D.  If 3D (or higher), the array fails to be converted to a
    # csr_matrix, so the error is TypeError.
    arr3d = np.ones(27, dtype=float).reshape(3, 3, 3)
    with pytest.raises(TypeError):
        cov = covariance.Covariance._ingest_matrix(arr3d)

    # Must be square
    with pytest.raises(ValueError):
        cov = covariance.Covariance._ingest_matrix(arr3d.reshape(9, 3))

    # Should be symmetric, but just throw a warning for now
    c = np.array([[1, 0], [1, 1]]).astype(float)
    with pytest.warns(AstropyUserWarning):
        cov = covariance.Covariance._ingest_matrix(c)

    # Can skip the symmetry check, the warning won't be thrown
    cov = covariance.Covariance._ingest_matrix(c, assume_symmetric=True)
    assert np.array_equal(cov.toarray(), c), (
        "_ingest_matrix should not correct symmetry issue if the check is skipped"
    )


@scipy_required
def test_bad_init_type():
    # Cannot be None
    with pytest.raises(ValueError):
        cov = covariance.Covariance()

    # Should be symmetric, but just throw a warning for now
    c = np.array([[1, 0], [1, 1]]).astype(float)
    with pytest.warns(AstropyUserWarning):
        cov = covariance.Covariance(array=c)

    # If the symmetry check is skipped, the warning won't be thrown but the
    # symmetry will be imposed anyway.
    cov = covariance.Covariance(array=c, assume_symmetric=True)
    assert np.array_equal(cov.to_dense(), np.identity(2)), (
        "Instantiation should correct symmetry"
    )


@scipy_required
def test_shape_mismatch():
    data_shape, c = mock_cov_2d()
    bad_shape = (2, 2)
    assert bad_shape != data_shape, (
        "Shapes should not match for the test to work properly"
    )
    with pytest.raises(ValueError):
        cov = covariance.Covariance(array=c, data_shape=bad_shape)


@scipy_required
def test_uncertainty_string():
    cov = covariance.Covariance(array=mock_cov())
    assert cov.uncertainty_type == "cov", "Uncertainty type changed"


@scipy_required
def test_quantity():
    data_shape, c = mock_cov_2d()
    cov = covariance.Covariance(array=c, data_shape=data_shape, unit="Jy^2")
    q = cov.quantity

    assert isinstance(q, units.Quantity), "Wrong type"
    assert q.unit == cov.unit, "Unit mismatch"
    assert np.array_equal(q.value, cov.to_dense()), "Array mismatch"


@scipy_required
def test_init():
    c = mock_cov()

    # Should not fault if directly instantiated from a numpy.array
    cov = covariance.Covariance(array=c)

    # Convert to a CSR sparse matrix
    c_csr = covariance.csr_matrix(c)
    # And instantiate
    cov = covariance.Covariance(array=c_csr)

    # Directly access the covariance array, convert it to a dense matrix, and
    # check it against the original array.
    assert np.array_equal(cov.to_dense(), c), "Ingested array does not match input"


@scipy_required
def test_stored_nnz():
    c = mock_cov()
    cov = covariance.Covariance(array=c)
    assert cov.stored_nnz == np.sum(np.triu(c) > 0), (
        "Number of stored non-zero elements differ"
    )


@scipy_required
def test_nnz():
    c_csr = covariance.csr_matrix(mock_cov())
    cov = covariance.Covariance(array=c_csr)
    assert c_csr.nnz == cov.nnz, "Number of non-zero elements differ"


@scipy_required
def test_indices():
    # Test when data_shape is not defined
    cov = covariance.Covariance(array=mock_cov())

    # Test out of bounds indices
    i = np.array([10, 1, 2])
    j = np.array([13, 4, 3])
    with pytest.raises(ValueError):
        cov.covariance_to_data_indices(i, j)
    with pytest.raises(ValueError):
        cov.data_to_covariance_indices(i, j)

    # Test in bounds
    i_cov = np.array([0, 1, 2])
    j_cov = np.array([3, 4, 3])
    i_data, j_data = cov.covariance_to_data_indices(i_cov, j_cov)
    assert np.array_equal(i_cov, i_data) and np.array_equal(j_cov, j_data), (
        "Should return output"
    )
    _i_cov, _j_cov = cov.data_to_covariance_indices(i_data, j_data)
    assert np.array_equal(_i_cov, i_data) and np.array_equal(_j_cov, j_data), (
        "Should return output"
    )

    # Test multi-dimensional data
    data_shape, c = mock_cov_2d()
    cov = covariance.Covariance(array=c, data_shape=data_shape)
    i_cov = np.array([0, 1, 2])
    j_cov = np.array([3, 4, 3])
    i_data, j_data = cov.covariance_to_data_indices(i_cov, j_cov)
    assert len(i_data) == 2, "Should return indices for each of 2 dimensions."
    assert (i_data[0][0], i_data[1][0]) == (0, 0), "Coordinate conversion error"

    _i_cov, _j_cov = cov.data_to_covariance_indices(i_data, j_data)
    assert np.array_equal(i_cov, _i_cov), "Inverse operation error"
    assert np.array_equal(j_cov, _j_cov), "Inverse operation error"

    # Test shape mismatch
    with pytest.raises(ValueError):
        cov.data_to_covariance_indices(i_data + (np.array([0, 0, 0]),), j_data)
    with pytest.raises(ValueError):
        cov.data_to_covariance_indices(i_data, j_data + (np.array([0, 0, 0]),))


@scipy_required
def test_coo():
    # 1D
    cov = covariance.Covariance(array=covariance.csr_matrix(mock_cov()))
    i, j, cij = cov.coordinate_data()
    assert i.size == cov.stored_nnz, "Coordinate data length is the incorrect size"

    # Cannot reshape when data_shape is not defined
    with pytest.raises(ValueError):
        cov.coordinate_data(reshape=True)

    # 2D
    data_shape, c = mock_cov_2d()
    c_csr = covariance.csr_matrix(c)
    cov = covariance.Covariance(array=c_csr, data_shape=data_shape)
    # Try without reshaping
    ic, jc, cij = cov.coordinate_data(reshape=False)
    assert isinstance(ic, np.ndarray), (
        "Index object should be an array if not reshaping"
    )
    assert ic.size == cov.stored_nnz, "Incorrect number of non-zero elements"
    # Try with reshaping
    i, j, cij = cov.coordinate_data(reshape=True)
    assert len(i) == len(data_shape), "Dimensionality does not match"
    assert i[0].size == cov.stored_nnz, "Incorrect number of non-zero elements"

    # Make sure we recover the same covariance matrix indices
    assert np.array_equal(ic, np.ravel_multi_index(i, cov.data_shape)), (
        "Bad covariance index mapping"
    )

    # 3D
    data_shape, c = mock_cov_3d()
    c_csr = covariance.csr_matrix(c)
    cov = covariance.Covariance(array=c_csr, data_shape=data_shape)
    i, j, cij = cov.coordinate_data(reshape=True)
    assert len(i) == len(data_shape), "Dimensionality does not match"
    assert i[0].size == cov.stored_nnz, "Incorrect number of non-zero elements"


@scipy_required
def test_copy():
    cov = covariance.Covariance(array=covariance.csr_matrix(mock_cov()))
    _cov = cov.copy()
    assert cov is not _cov, "Objects have the same reference"
    assert cov._cov is not _cov._cov, "Object arrays have the same reference"
    assert np.array_equal(cov.to_dense(), _cov.to_dense()), "Arrays should be equal"


@scipy_required
def test_tbls():
    cov = covariance.Covariance(array=covariance.csr_matrix(mock_cov()))
    covar = cov.to_table()
    assert isinstance(covar, Table), "correlation data should be output as a table"
    assert len(covar) == np.sum(np.triu(mock_cov()) > 0), (
        "Incorrect number of table entries"
    )
    assert len(covar.colnames) == 3, "Incorrect number of columns"
    assert covar["INDXI"].ndim == 1, "Incorrect shape for index array"

    _cov = covariance.Covariance.from_table(covar)
    assert np.array_equal(cov.to_dense(), _cov.to_dense()), (
        "Bad convert/revert from tables"
    )

    # Test failure when meta is None
    with pytest.raises(ValueError):
        _covar = covar.copy()
        _covar.meta = {}
        _cov = covariance.Covariance.from_table(_covar)

    data_shape, c = mock_cov_3d()
    cov = covariance.Covariance(array=covariance.csr_matrix(c), data_shape=data_shape)
    covar = cov.to_table()
    assert len(covar) == np.sum(np.triu(c) > 0), "Incorrect number of table entries"
    assert len(covar.colnames) == 3, "Incorrect number of columns"
    assert covar["INDXI"].ndim == 2, "Incorrect shape for index array"

    _cov = covariance.Covariance.from_table(covar)
    assert np.array_equal(cov.to_dense(), _cov.to_dense()), (
        "Bad convert/revert from tables"
    )

    # Introduce a shape mismatch
    covar.meta["COVDSHP"] = "(3,6)"
    with pytest.raises(ValueError):
        _cov = covariance.Covariance.from_table(covar)


@scipy_required
def test_samples():
    # Fixed-seed RNG for repeatability
    rng = np.random.default_rng(seed=8001)

    m = np.zeros(10, dtype=float)
    c = mock_cov()

    # Shape must be 2D
    s = rng.multivariate_normal(m, c)
    with pytest.raises(ValueError):
        covar = covariance.Covariance.from_samples(s, cov_tol=0.1)

    # Shape must be at least two samples
    s = np.expand_dims(s, 1).T
    with pytest.raises(ValueError):
        covar = covariance.Covariance.from_samples(s.T, cov_tol=0.1)

    # Draw samples
    s = rng.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = covariance.Covariance.from_samples(s.T, cov_tol=0.1)

    # Check the values are very nearly the same as the input
    assert np.all(np.absolute(c - covar.to_dense()) < 0.02), (
        "Covariances are too different"
    )

    # Check that `find` returns the same indices
    # NOTE: For some reason the ordering of np.where and
    # scipy.sparse.find are different.  So I run lexsort before checking
    # if the arrays are the same.

    coo = np.array([[i, j] for i, j in zip(*np.where(c > 0))])
    srt = np.lexsort((coo[:, 0], coo[:, 1]))
    coo = coo[srt, :]

    _coo = np.array([[i, j] for i, j, k in zip(*covar.find())])
    srt = np.lexsort((_coo[:, 0], _coo[:, 1]))
    _coo = _coo[srt, :]

    assert np.array_equal(coo, _coo), "Did not find the same indices"


@scipy_required
def test_mult():
    c = mock_cov()
    x = np.ones(10, dtype=float)

    # Transfer matrix must be 2D
    t = np.ones(3, dtype=float)
    with pytest.raises(ValueError):
        covar = covariance.Covariance.from_matrix_multiplication(t, c)

    # Uncorrelated
    t = np.zeros((3, 10), dtype=float)
    t[0, 0] = 1.0
    t[1, 3] = 1.0
    t[2, 6] = 1.0

    # Sigma has the wrong shape
    with pytest.raises(ValueError):
        covar = covariance.Covariance.from_matrix_multiplication(t, c[1:, 1:])
    with pytest.raises(ValueError):
        covar = covariance.Covariance.from_matrix_multiplication(t, np.diag(c)[1:])

    y = np.dot(t, x)

    covar = covariance.Covariance.from_matrix_multiplication(t, c)
    assert np.array_equal(covar.to_dense(), np.identity(3)), (
        "Result should be uncorrelated."
    )

    # Correlated by 0.2
    t = np.zeros((3, 10), dtype=float)
    t[0, 0] = 1.0
    t[1, 2] = 1.0
    t[2, 4] = 1.0
    _c = (
        np.diag(np.full(3 - 1, 0.2, dtype=float), k=-1)
        + np.diag(np.full(3, 1.0, dtype=float), k=0)
        + np.diag(np.full(3 - 1, 0.2, dtype=float), k=1)
    )
    y = np.dot(t, x)
    covar = covariance.Covariance.from_matrix_multiplication(t, c)
    assert np.array_equal(covar.to_dense(), _c), (
        "Result should have off-diagonals = 0.2"
    )

    # Correlated by 0.5 and 0.2
    t = np.zeros((3, 10), dtype=float)
    t[0, 0] = 1.0
    t[1, 1] = 1.0
    t[2, 2] = 1.0
    _c = (
        np.diag(np.full(3 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(3 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(3, 1.0, dtype=float), k=0)
        + np.diag(np.full(3 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(3 - 2, 0.2, dtype=float), k=2)
    )
    y = np.dot(t, x)
    covar = covariance.Covariance.from_matrix_multiplication(t, c)
    assert np.array_equal(covar.to_dense(), _c), (
        "Result should have off-diagonals = 0.5,0.2"
    )


@scipy_required
def test_var():
    var = np.ones(3, dtype=float)
    covar = covariance.Covariance.from_variance(var)
    assert np.array_equal(covar.to_dense(), np.identity(3)), (
        "Result should be an identity matrix"
    )


@scipy_required
def test_shape():
    cov = covariance.Covariance(array=mock_cov())
    assert len(cov.data_shape) == 1, "Incorrect dimensionality"
    assert cov.data_shape == (mock_cov().shape[0],), "Bad data shape"

    # 2D
    data_shape, c = mock_cov_2d()
    cov = covariance.Covariance(array=c, data_shape=data_shape)
    assert len(cov.data_shape) == 2, "Incorrect dimensionality"
    assert cov.data_shape == data_shape, "Bad data shape"

    # 3D
    data_shape, c = mock_cov_3d()
    cov = covariance.Covariance(array=c, data_shape=data_shape)
    assert len(cov.data_shape) == 3, "Incorrect dimensionality"
    assert cov.data_shape == data_shape, "Bad data shape"


@scipy_required
def test_match_to_data_slice():
    c = mock_cov()
    cov = covariance.Covariance(array=c)

    # Check the index map
    assert np.array_equal(cov.data_index_map, np.arange(10)), "Bad index map"

    # 1D slice
    sub_cov = cov.match_to_data_slice(np.s_[:5])
    assert isinstance(sub_cov, covariance.Covariance), (
        "Submatrix should be a Covariance instance"
    )
    assert np.array_equal(sub_cov.to_dense(), c[:5, :5]), "Bad submatrix"

    # 1D reorder
    rng = np.random.default_rng(99)
    reorder = np.arange(cov.shape[0])
    rng.shuffle(reorder)
    sub_cov = cov.match_to_data_slice(reorder)
    assert isinstance(sub_cov, covariance.Covariance), (
        "Submatrix should be a Covariance instance"
    )
    assert np.array_equal(sub_cov.to_dense(), c[np.ix_(reorder, reorder)]), (
        "Bad matrix reorder"
    )

    # 2D
    data_shape, c = mock_cov_2d()
    cov = covariance.Covariance(array=c, data_shape=data_shape)

    # Check the index map
    assert cov.data_index_map.shape == data_shape, "Bad index map shape"
    assert np.array_equal(cov.data_index_map, np.arange(6).reshape(data_shape)), (
        "Bad index map"
    )

    # Reduce dimensionality
    sub_cov = cov.match_to_data_slice(np.s_[:, 0])
    assert len(sub_cov.data_shape) == 1, "Submatrix should have reduced dimensionality"
    assert sub_cov.shape == (3, 3), "Incorrect covariance submatrix shape"
    # Still 2D
    sub_cov = cov.match_to_data_slice(np.s_[1:, :])
    assert sub_cov.data_shape == (2, 2), "Submatrix does not have the correct shape"

    # 3D
    data_shape, c = mock_cov_3d()
    cov = covariance.Covariance(array=c, data_shape=data_shape)
    # Reduce dimensionality to 1D
    sub_cov = cov.match_to_data_slice(np.s_[:, 0, 0])
    assert len(sub_cov.data_shape) == 1, "Submatrix should have reduced dimensionality"
    assert sub_cov.shape == (3, 3), "Incorrect covariance submatrix shape"
    # Reduce dimensionality to 2D
    sub_cov = cov.match_to_data_slice(np.s_[1:, :, 0])
    assert sub_cov.data_shape == (2, 2), "Submatrix does not have the correct shape"
    # Still 3D
    sub_cov = cov.match_to_data_slice(np.s_[1:, :, :-1])
    assert sub_cov.data_shape == (2, 2, 2), "Submatrix does not have the correct shape"
    assert sub_cov.shape[0] == np.prod(sub_cov.data_shape), "Shape mismatch"


@scipy_required
def test_correl():
    # Since the diagonal is 1, this is actually a correlation matrix
    c = mock_cov()
    # Set two covariance arrays with different variances but the same correlations
    cov1 = covariance.Covariance(array=c * 4.0)
    cov2 = covariance.Covariance(array=c * 2.0)
    rho1 = cov1.to_sparse(correlation=True)
    # Should be the same as the input
    assert np.allclose(rho1.toarray(), c), "Correlation matrix changed"

    # Should match cov1
    rho2 = cov2.to_sparse(correlation=True)
    assert np.allclose(rho1.toarray(), rho2.toarray()), (
        "Correlation matrices should be identical"
    )
    assert np.allclose(cov1.variance / 2, cov2.variance), "Variances incorrect"


@scipy_required
def test_to_from_correl():
    c = mock_cov()
    cov1 = covariance.Covariance(array=c * 4.0)

    var, rho1 = covariance.Covariance.to_correlation(cov1.to_sparse())
    assert np.array_equal(rho1.toarray(), c), "Bad correlation matrix recovery"

    _cov1 = covariance.Covariance.revert_correlation(var, rho1)
    assert np.array_equal(cov1.to_dense(), _cov1.toarray())


@scipy_required
def test_newvar():
    c = mock_cov()
    cov1 = covariance.Covariance(array=c)
    var = np.full(c.shape[0], 4.0, dtype=float)
    with pytest.raises(NotImplementedError):
        cov1.variance = var
    cov2 = cov1.apply_new_variance(var)
    assert np.allclose(var, cov2.variance), "Variance does not match request"
    var2, rho2 = covariance.Covariance.to_correlation(cov2.to_dense())
    assert np.allclose(cov1.to_dense(), rho2.toarray()), (
        "Correlation matrices do not match"
    )


@scipy_required
def test_array():
    # Construct the Covariance matrix from a pre-calculated array
    c = mock_cov()
    covar = covariance.Covariance.from_array(c)
    # Should be the same as the identity matrix.
    assert np.array_equal(covar.to_dense(), c), "Arrays should be identical"

    # Copy
    _c = c.copy()
    # Introduce an asymmetry in the upper triangle
    _c[0, 9] = 1.0
    with pytest.warns(AstropyUserWarning):
        cov = covariance.Covariance.from_array(_c)

    # Instantiate again but assume symmetry
    cov = covariance.Covariance.from_array(_c, assume_symmetric=True)
    assert not np.array_equal(cov.to_dense(), c), (
        "Asymmetries in the upper triangle should be kept"
    )
    _c = c.copy()
    # Introduce an asymmetry in the lower triangle
    _c[9, 0] = 1.0
    cov = covariance.Covariance.from_array(_c, assume_symmetric=True)
    assert np.array_equal(cov.to_dense(), c), (
        "Asymmetries in the lower triangle should be ignored"
    )

    # Test rho tolerance (cov tolerance is test elsewhere)
    rho_tol = 0.3
    covar = covariance.Covariance.from_array(c, rho_tol=rho_tol)
    assert not np.array_equal(covar.to_dense(), c), (
        "Tolerance should have removed values"
    )
    _c = covar.to_dense()
    assert not np.any((_c > 0) & (_c < rho_tol)), (
        "Array includes elements below tolerance"
    )


@scipy_required
def test_io():
    # Set the file name
    ofile = Path("test_covar_io.fits")
    # Erase it if it already exists
    if ofile.is_file():
        ofile.unlink()

    # Create the covariance object
    covar = covariance.Covariance(array=mock_cov())
    # Write it to a table
    tbl = covar.to_table()
    # Write the table to a fits file
    tbl.write(ofile, format="fits")
    # Read it back in as a Table
    _tbl = Table.read(ofile, format="fits")
    # Use the table to instantiate a Covariance object
    _covar = covariance.Covariance.from_table(_tbl)
    # Check that the IO was successful
    assert np.array_equal(covar.to_dense(), _covar.to_dense()), "Array changed"
    # Delete the file
    ofile.unlink()
