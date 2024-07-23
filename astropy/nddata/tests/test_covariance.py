from pathlib import Path

import numpy as np
import pytest

from astropy import units
from astropy.io import fits
from astropy.nddata import covariance
from astropy.table import Table
from astropy.utils.compat.optional_deps import HAS_SCIPY

if HAS_SCIPY:
    from scipy.sparse import csr_matrix

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


# NOTE: At least for now, this does not require the @scipy_required decorator
# because an empty initialization doesn't use any scipy files.
# def test_empty_init():
#     cov = covariance.Covariance()
#     assert cov.array is None, "Array should not be defined"


@scipy_required
def test_bad_init_type():
    # It must be possible to convert the input array into a csr_matrix
    with pytest.raises(TypeError):
        cov = covariance.Covariance("test")


@scipy_required
def test_shape_mismatch():
    raw_shape, c = mock_cov_2d()
    bad_shape = (2, 2)
    assert (
        bad_shape != raw_shape
    ), "Shapes should not match for the test to work properly"
    with pytest.raises(ValueError):
        cov = covariance.Covariance(c, raw_shape=bad_shape)


def test_uncertainty_string():
    cov = covariance.Covariance(mock_cov())
    assert cov.uncertainty_type == "cov", "Uncertainty type changed"


@scipy_required
def test_quantity():
    raw_shape, c = mock_cov_2d()
    cov = covariance.Covariance(c, raw_shape=raw_shape, unit="Jy^2")
    q = cov.quantity

    assert isinstance(q, units.Quantity), "Wrong type"
    assert q.unit == cov.unit, "Unit mismatch"
    assert np.array_equal(q.value, cov.toarray()), "Array mismatch"


def test_bad_init_type():
    # It must be possible to convert the input array into a csr_matrix
    with pytest.raises(TypeError):
        cov = covariance.Covariance(array="test")


def test_shape_mismatch():
    raw_shape, c = mock_cov_2d()
    bad_shape = (2, 2)
    assert (
        bad_shape != raw_shape
    ), "Shapes should not match for the test to work properly"
    with pytest.raises(ValueError):
        cov = covariance.Covariance(array=c, raw_shape=bad_shape)


def test_uncertainty_string():
    cov = covariance.Covariance()
    assert cov.uncertainty_type == "cov", "Uncertainty type changed"


def test_quantity():
    raw_shape, c = mock_cov_2d()
    cov = covariance.Covariance(array=c, raw_shape=raw_shape, unit="Jy^2")
    q = cov.quantity

    assert isinstance(q, units.Quantity), "Wrong type"
    assert q.unit == cov.unit, "Unit mismatch"
    assert np.array_equal(q.value, cov.toarray()), "Array mismatch"


@scipy_required
def test_init():
    c = mock_cov()

    # Should not fault if directly instantiated from a numpy.array
    cov = covariance.Covariance(c)

    # Convert to a CSR sparse matrix
    c_csr = csr_matrix(c)
    # And instantiate
    cov = covariance.Covariance(c_csr)

    # Directly access the covariance array, convert it to a dense matrix, and
    # check it against the original array.
    assert np.array_equal(cov.toarray(), c), "Ingested array does not match input"


@scipy_required
def test_stored_nnz():
    c = mock_cov()
    cov = covariance.Covariance(c)
    assert cov.stored_nnz == np.sum(
        np.triu(c) > 0
    ), "Number of stored non-zero elements differ"


@scipy_required
def test_nnz():
    c_csr = csr_matrix(mock_cov())
    cov = covariance.Covariance(c_csr)
    assert c_csr.nnz == cov.nnz, "Number of non-zero elements differ"


@scipy_required
def test_indices():
    # Test when raw_shape is not defined
    cov = covariance.Covariance(mock_cov())

    # Test out of bounds indices
    i = np.array([10, 1, 2])
    j = np.array([13, 4, 3])
    with pytest.raises(ValueError):
        cov.cov2raw_indices(i, j)
    with pytest.raises(ValueError):
        cov.raw2cov_indices(i, j)

    # Test in bounds
    i_cov = np.array([0, 1, 2])
    j_cov = np.array([3, 4, 3])
    i_data, j_data = cov.cov2raw_indices(i_cov, j_cov)
    assert np.array_equal(i_cov, i_data) and np.array_equal(
        j_cov, j_data
    ), "Should return output"
    _i_cov, _j_cov = cov.raw2cov_indices(i_data, j_data)
    assert np.array_equal(_i_cov, i_data) and np.array_equal(
        _j_cov, j_data
    ), "Should return output"

    # Test multi-dimensional data
    raw_shape, c = mock_cov_2d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    i_cov = np.array([0, 1, 2])
    j_cov = np.array([3, 4, 3])
    i_data, j_data = cov.cov2raw_indices(i_cov, j_cov)
    assert len(i_data) == 2, "Should return indices for each of 2 dimensions."
    assert (i_data[0][0], i_data[1][0]) == (0, 0), "Coordinate conversion error"

    _i_cov, _j_cov = cov.raw2cov_indices(i_data, j_data)
    assert np.array_equal(i_cov, _i_cov), "Inverse operation error"
    assert np.array_equal(j_cov, _j_cov), "Inverse operation error"

    # Test shape mismatch
    with pytest.raises(ValueError):
        cov.raw2cov_indices(i_data + (np.array([0, 0, 0]),), j_data)
    with pytest.raises(ValueError):
        cov.raw2cov_indices(i_data, j_data + (np.array([0, 0, 0]),))


@scipy_required
def test_coo():
    # 1D
    cov = covariance.Covariance(csr_matrix(mock_cov()))
    i, j, rhoij, var = cov.coordinate_data()
    assert i.size == cov._rho.nnz, "Coordinate data length is the incorrect size"
    assert var.ndim == 1, "Incorrect dimensionality"

    # Cannot reshape when raw_shape is not defined
    with pytest.raises(ValueError):
        cov.coordinate_data(reshape=True)

    # 2D
    raw_shape, c = mock_cov_2d()
    c_csr = csr_matrix(c)
    cov = covariance.Covariance(c_csr, raw_shape=raw_shape)
    # Try without reshaping
    ic, jc, rhoij, var = cov.coordinate_data(reshape=False)
    assert isinstance(
        ic, np.ndarray
    ), "Index object should be an array if not reshaping"
    assert ic.size == cov._rho.nnz, "Incorrect number of non-zero elements"
    assert var.shape == (np.prod(raw_shape),), "Variance array has incorrect shape"
    # Try with reshaping
    i, j, rhoij, var = cov.coordinate_data(reshape=True)
    assert len(i) == len(raw_shape), "Dimensionality does not match"
    assert i[0].size == cov._rho.nnz, "Incorrect number of non-zero elements"
    assert var.shape == raw_shape, "Variance array has incorrect shape"

    # Make sure we recover the same covariance matrix indices
    assert np.array_equal(
        ic, np.ravel_multi_index(i, cov.raw_shape)
    ), "Bad covariance index mapping"

    # 3D
    raw_shape, c = mock_cov_3d()
    c_csr = csr_matrix(c)
    cov = covariance.Covariance(c_csr, raw_shape=raw_shape)
    i, j, rhoij, var = cov.coordinate_data(reshape=True)
    assert len(i) == len(raw_shape), "Dimensionality does not match"
    assert i[0].size == cov._rho.nnz, "Incorrect number of non-zero elements"
    assert var.shape == raw_shape, "Variance array has incorrect shape"


@scipy_required
def test_copy():
    cov = covariance.Covariance(csr_matrix(mock_cov()))
    _cov = cov.copy()
    assert cov is not _cov, "Objects have the same reference"
    assert cov._rho is not _cov._rho, "Object arrays have the same reference"
    assert np.array_equal(cov.toarray(), _cov.toarray()), "Arrays should be equal"

    # Convert to correlation
    cov.to_correlation()
    _cov = cov.copy()
    assert cov.cov is not _cov.cov, "Object arrays have the same reference"
    assert _cov.is_correlation, "Should still be a correlation matrix"
    cov.revert_correlation()
    _cov.revert_correlation()
    assert np.array_equal(cov.toarray(), _cov.toarray()), "Arrays should be equal"


@scipy_required
def test_tbls():
    cov = covariance.Covariance(csr_matrix(mock_cov()))
    var, correl = cov.to_tables()
    assert isinstance(var, np.ndarray), "variance should be output as an array"
    assert isinstance(correl, Table), "correlation data should be output as a table"
    assert len(correl) == np.sum(
        np.triu(mock_cov()) > 0
    ), "Incorrect number of table entries"
    assert len(correl.colnames) == 3, "Incorrect number of columns"
    assert correl["INDXI"].ndim == 1, "Incorrect shape for index array"

    _cov = covariance.Covariance.from_tables(var, correl)
    assert np.array_equal(
        cov.toarray(), _cov.toarray()
    ), "Bad convert/revert from tables"

    raw_shape, c = mock_cov_3d()
    cov = covariance.Covariance(csr_matrix(c), raw_shape=raw_shape)
    var, correl = cov.to_tables()
    assert len(correl) == np.sum(np.triu(c) > 0), "Incorrect number of table entries"
    assert len(correl.colnames) == 3, "Incorrect number of columns"
    assert correl["INDXI"].ndim == 2, "Incorrect shape for index array"
    assert (
        correl["INDXI"].shape[1] == var.ndim
    ), "Dimensionality mismatch between var and indices"

    _cov = covariance.Covariance.from_tables(var, correl)
    assert np.array_equal(
        cov.toarray(), _cov.toarray()
    ), "Bad convert/revert from tables"


@scipy_required
def test_samples():
    # Fixed-seed RNG for repeatability
    rng = np.random.default_rng(seed=8001)

    m = np.zeros(10, dtype=float)
    c = mock_cov()

    # Draw samples
    s = rng.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = covariance.Covariance.from_samples(s.T, cov_tol=0.1)

    # Check the values are very nearly the same as the input
    assert np.all(
        np.absolute(c - covar.toarray()) < 0.02
    ), "Covariances are too different"

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

    # Uncorrelated
    t = np.zeros((3, 10), dtype=float)
    t[0, 0] = 1.0
    t[1, 3] = 1.0
    t[2, 6] = 1.0

    y = np.dot(t, x)

    covar = covariance.Covariance.from_matrix_multiplication(t, c)
    assert np.array_equal(
        covar.toarray(), np.identity(3)
    ), "Result should be uncorrelated."

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
    assert np.array_equal(covar.toarray(), _c), "Result should have off-diagonals = 0.2"

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
    assert np.array_equal(
        covar.toarray(), _c
    ), "Result should have off-diagonals = 0.5,0.2"


@scipy_required
def test_var():
    var = np.ones(3, dtype=float)
    covar = covariance.Covariance.from_variance(var)
    assert np.array_equal(
        covar.toarray(), np.identity(3)
    ), "Result should be an identity matrix"


@scipy_required
def test_shape():
    cov = covariance.Covariance(mock_cov())
    assert len(cov.data_shape) == 1, "Incorrect dimensionality"
    assert cov.data_shape == (mock_cov().shape[0],), "Bad data shape"

    # 2D
    raw_shape, c = mock_cov_2d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    assert len(cov.data_shape) == 2, "Incorrect dimensionality"
    assert cov.data_shape == raw_shape, "Bad data shape"

    # 3D
    raw_shape, c = mock_cov_3d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    assert len(cov.data_shape) == 3, "Incorrect dimensionality"
    assert cov.data_shape == raw_shape, "Bad data shape"


@scipy_required
def test_sub_matrix():
    c = mock_cov()
    cov = covariance.Covariance(c)

    # 1D
    sub_cov = cov.sub_matrix(np.s_[:5])
    assert isinstance(
        sub_cov, covariance.Covariance
    ), "Submatrix should be a Covariance instance"
    assert np.array_equal(sub_cov.toarray(), c[:5, :5]), "Bad submatrix"

    # 2D
    raw_shape, c = mock_cov_2d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    # Reduce dimensionality
    sub_cov = cov.sub_matrix(np.s_[:, 0])
    assert sub_cov.raw_shape is None, "Submatrix should have reduced dimensionality"
    assert sub_cov.shape == (3, 3), "Incorrect covariance submatrix shape"
    # Still 2D
    sub_cov = cov.sub_matrix(np.s_[1:, :])
    assert sub_cov.raw_shape == (2, 2), "Submatrix does not have the correct shape"

    # 3D
    raw_shape, c = mock_cov_3d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    # Reduce dimensionality to 1D
    sub_cov = cov.sub_matrix(np.s_[:, 0, 0])
    assert sub_cov.raw_shape is None, "Submatrix should have reduced dimensionality"
    assert sub_cov.shape == (3, 3), "Incorrect covariance submatrix shape"
    # Reduce dimensionality to 2D
    sub_cov = cov.sub_matrix(np.s_[1:, :, 0])
    assert sub_cov.raw_shape == (2, 2), "Submatrix does not have the correct shape"
    # Still 3D
    sub_cov = cov.sub_matrix(np.s_[1:, :, :-1])
    assert sub_cov.raw_shape == (2, 2, 2), "Submatrix does not have the correct shape"
    assert sub_cov.shape[0] == np.prod(sub_cov.raw_shape), "Shape mismatch"


@scipy_required
def test_correl():
    # Since the diagonal is 1, this is actually a correlation matrix
    c = mock_cov()
    # Set two covariance arrays with different variances but the same correlations
    cov1 = covariance.Covariance(c * 4.0)
    cov2 = covariance.Covariance(c * 2.0)
    rho1 = cov1.full(correlation=True)
    # Should be the same as the input
    assert np.allclose(rho1.toarray(), c), "Correlation matrix changed"

    # Should match cov1
    rho2 = cov2.full(correlation=True)
    assert np.allclose(
        rho1.toarray(), rho2.toarray()
    ), "Correlation matrices should be identical"
    assert np.allclose(cov1._var / 2, cov2._var), "Variances incorrect"


@scipy_required
def test_newvar():
    c = mock_cov()
    cov1 = covariance.Covariance(c)
    var = np.full(c.shape[0], 4.0, dtype=float)
    cov2 = cov1.apply_new_variance(var)
    assert np.allclose(var, cov2._var), "Variance does not match request"
    var2, rho2 = covariance.Covariance.to_correlation(cov2.toarray())
    assert np.allclose(
        cov1.toarray(), rho2.toarray()
    ), "Correlation matrices do not match"


@scipy_required
def test_array():
    # Construct the Covariance matrix from a pre-calculated array
    c = mock_cov()
    covar = covariance.Covariance.from_array(c)
    # Should be the same as the identity matrix.
    assert np.array_equal(covar.toarray(), c), "Arrays should be identical"

    # Test rho tolerance (cov tolerance is test elsewhere)
    rho_tol = 0.3
    covar = covariance.Covariance.from_array(c, rho_tol=rho_tol)
    assert not np.array_equal(
        covar.toarray(), c
    ), "Tolerance should have removed values"
    _c = covar.toarray()
    assert not np.any(
        (_c > 0) & (_c < rho_tol)
    ), "Array includes elements below tolerance"

@scipy_required
def test_io():
    # Clean up in case of a failure
    ofile = Path("test_covar_io.fits")
    if ofile.is_file():
        ofile.unlink()

    # 1D
    cov = covariance.Covariance(mock_cov(), unit="km^2")
    # Rescale the variance
    cov = cov.apply_new_variance(np.full(cov.shape[0], 2.0))
    # Write
    cov.write(ofile)

    # Test file exists
    with pytest.raises(FileExistsError):
        cov.write(ofile)

    # Test contents of fits file
    with fits.open(ofile) as hdu:
        assert len(hdu) == 3, "Incorrect number of extensions written"
        assert hdu[1].name == "VAR", "Default name changed"
        assert hdu[2].name == "CORREL", "Default name changed"
        assert len(hdu["CORREL"].data.columns) == 3, "Incorrect number of table columns"
        assert hdu["CORREL"].data.columns.names == [
            "INDXI",
            "INDXJ",
            "RHOIJ",
        ], "Column names changed"
        assert len(hdu["CORREL"].data["INDXI"].shape) == 1, "Column should only be 1D"
        assert hdu["CORREL"].header["BUNIT"].strip() == "km2", "Unit wrong"

    with pytest.raises(ValueError):
        # Must define extension with correlation matrix
        _cov = covariance.Covariance.from_fits(ofile, covar_ext=None)

    # Read
    _cov = covariance.Covariance.from_fits(ofile)
    # Arrays should be the same
    assert np.allclose(cov.toarray(), _cov.toarray()), "Bad 1D I/O"
    # Units should be the same
    assert cov.unit == _cov.unit, "Units changed"

    # Read but ignore the variance extension
    _cov = covariance.Covariance.from_fits(ofile, var_ext=None)
    # This sets the variance to unity, so this should be the same as the
    # original covariance returned by mock_cov()
    assert np.array_equal(_cov.toarray(), mock_cov()), "Bad read"

    # Clean-up
    ofile.unlink()

    # ND
    raw_shape, c = mock_cov_3d()
    cov = covariance.Covariance(c, raw_shape=raw_shape)
    cov.write(ofile)
    # Test contents of fits file
    with fits.open(ofile) as hdu:
        assert len(hdu) == 3, "Incorrect number of extensions written"
        assert hdu[1].name == "VAR", "Default name changed"
        assert hdu[2].name == "CORREL", "Default name changed"
        assert len(hdu["CORREL"].data.columns) == 3, "Incorrect number of table columns"
        assert len(hdu["CORREL"].data["INDXI"].shape) == 2, "Column should be ND"
        assert hdu["CORREL"].data["INDXI"].shape[1] == 3, "Data is 3D"
    # Read
    _cov = covariance.Covariance.from_fits(ofile)
    # Arrays should be the same
    assert np.allclose(cov.toarray(), _cov.toarray()), "Bad ND I/O"
    # Clean-up
    ofile.unlink()
