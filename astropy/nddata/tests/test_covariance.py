from pathlib import Path

import numpy as np
from scipy import sparse

from astropy.io import fits
from astropy.nddata import Covariance
from astropy.table import Table


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


def test_empty_init():
    cov = Covariance()
    assert cov.array is None, "Array should not be defined"


def test_init():
    c = mock_cov()

    # Should not if directly instantiated from a numpy.array
    cov = Covariance(array=c)

    # Convert to a CSR sparse matrix
    c_csr = sparse.csr_matrix(c)
    # And instantiate
    cov = Covariance(array=c_csr)

    # Directly access the covariance array, convert it to a dense matrix, and
    # check it against the original array.
    assert np.array_equal(cov.cov.toarray(), c), "Ingested array does not match input"

    # Recreate, but force it to only keep the upper triangle
    cov = Covariance(array=c_csr, impose_triu=True)

    # The arrays should no longer be equal because the lower triangle is missing
    # when directly accessing the cov attribute
    assert not np.array_equal(cov.cov.toarray(), c), "Should remove lower triangle"

    # They should be identical when using the Covariance.toarray() function
    # instead.
    assert np.array_equal(cov.toarray(), c), "Array should be filled"


def test_nnz():
    c_csr = sparse.csr_matrix(mock_cov())
    cov = Covariance(array=c_csr)

    assert c_csr.nnz == cov.nnz, "Number of non-zero elements differ"

    cov = Covariance(array=c_csr, impose_triu=True)
    ndiag = c_csr.shape[0]
    triu_nnz = (c_csr.nnz - ndiag) // 2 + ndiag

    assert (
        cov.nnz == triu_nnz
    ), "Number of non-zero elements in the upper triangle differ"
    assert (
        c_csr.nnz == cov.full().nnz
    ), "Number of non-zero elements in full array differ"


def test_indices():
    raw_shape, c = mock_cov_2d()
    cov = Covariance(array=c, impose_triu=True, raw_shape=raw_shape)
    i_cov = np.array([0, 1, 2])
    j_cov = np.array([3, 4, 3])
    i_data, j_data = cov.cov2raw_indices(i_cov, j_cov)
    assert len(i_data) == 2, "Should return indices for each of 2 dimensions."
    assert (i_data[0][0], i_data[1][0]) == (0, 0), "Coordinate conversion error"

    _i_cov, _j_cov = cov.raw2cov_indices(i_data, j_data)
    assert np.array_equal(i_cov, _i_cov), "Inverse operation error"
    assert np.array_equal(j_cov, _j_cov), "Inverse operation error"


def test_coo():
    # 1D
    cov = Covariance(array=sparse.csr_matrix(mock_cov()))
    i, j, rhoij, var = cov.coordinate_data()
    assert i.size == cov.nnz, "Coordinate data length is the incorrect size"
    assert var.ndim == 1, "Incorrect dimensionality"

    # 2D
    raw_shape, c = mock_cov_2d()
    c_csr = sparse.csr_matrix(c)
    cov = Covariance(array=c_csr, impose_triu=True, raw_shape=raw_shape)
    # Try without reshaping
    ic, jc, rhoij, var = cov.coordinate_data(reshape=False)
    assert isinstance(
        ic, np.ndarray
    ), "Index object should be an array if not reshaping"
    assert ic.size == cov.nnz, "Incorrect number of non-zero elements"
    assert var.shape == (np.prod(raw_shape),), "Variance array has incorrect shape"
    # Try with reshaping
    i, j, rhoij, var = cov.coordinate_data(reshape=True)
    assert len(i) == len(raw_shape), "Dimensionality does not match"
    assert i[0].size == cov.nnz, "Incorrect number of non-zero elements"
    assert var.shape == raw_shape, "Variance array has incorrect shape"

    # Make sure we recover the same covariance matrix indices
    assert np.array_equal(
        ic, np.ravel_multi_index(i, cov.raw_shape)
    ), "Bad covariance index mapping"

    # 3D
    raw_shape, c = mock_cov_3d()
    c_csr = sparse.csr_matrix(c)
    cov = Covariance(array=c_csr, impose_triu=True, raw_shape=raw_shape)
    i, j, rhoij, var = cov.coordinate_data(reshape=True)
    assert len(i) == len(raw_shape), "Dimensionality does not match"
    assert i[0].size == cov.nnz, "Incorrect number of non-zero elements"
    assert var.shape == raw_shape, "Variance array has incorrect shape"


def test_copy():
    cov = Covariance(array=sparse.csr_matrix(mock_cov()))
    _cov = cov.copy()
    assert cov is not _cov, "Objects have the same reference"
    assert cov.cov is not _cov.cov, "Object arrays have the same reference"
    assert np.array_equal(cov.toarray(), _cov.toarray()), "Arrays should be equal"


def test_tbls():
    cov = Covariance(array=sparse.csr_matrix(mock_cov()))
    var, correl = cov.output_tables()
    assert isinstance(var, np.ndarray), "variance should be output as an array"
    assert isinstance(correl, Table), "correlation data should be output as a table"
    assert len(correl) == 44, "Incorrect number of table entries"
    assert len(correl.colnames) == 3, "Incorrect number of columns"
    assert correl["INDXI"].ndim == 1, "Incorrect shape for index array"

    raw_shape, c = mock_cov_3d()
    cov = Covariance(array=sparse.csr_matrix(c), impose_triu=True, raw_shape=raw_shape)
    var, correl = cov.output_tables()
    assert len(correl) == 51, "Incorrect number of table entries"
    assert len(correl.colnames) == 3, "Incorrect number of columns"
    assert correl["INDXI"].ndim == 2, "Incorrect shape for index array"
    assert (
        correl["INDXI"].shape[1] == var.ndim
    ), "Dimensionality mismatch between var and indices"


def test_samples():
    # Fixed-seed RNG for repeatability
    rng = np.random.default_rng(seed=8001)

    m = np.zeros(10, dtype=float)
    c = mock_cov()

    # Draw samples
    s = rng.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = Covariance.from_samples(s.T, cov_tol=0.1)

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


def test_mult():
    c = mock_cov()
    x = np.ones(10, dtype=float)

    # Uncorrelated
    t = np.zeros((3, 10), dtype=float)
    t[0, 0] = 1.0
    t[1, 3] = 1.0
    t[2, 6] = 1.0

    y = np.dot(t, x)

    covar = Covariance.from_matrix_multiplication(t, c)
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
    covar = Covariance.from_matrix_multiplication(t, c)
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
    covar = Covariance.from_matrix_multiplication(t, c)
    assert np.array_equal(
        covar.toarray(), _c
    ), "Result should have off-diagonals = 0.5,0.2"


def test_var():
    var = np.ones(3, dtype=float)
    covar = Covariance.from_variance(var)
    assert np.array_equal(
        covar.toarray(), np.identity(3)
    ), "Result should be an identity matrix"


def test_correl():
    # Since the diagonal is 1, this is actually a correlation matrix
    c = mock_cov()
    # Set two covariance arrays with different variances but the same correlations
    cov1 = Covariance(array=c * 4.0, impose_triu=True)
    cov2 = Covariance(array=c * 2.0, impose_triu=True)
    # Convert them both to correlation matrices
    assert (
        not cov1.is_correlation
    ), "Should start as a covariance, not correlation, matrix"
    assert cov1.var is None, "Variance is not yet defined"
    cov1.to_correlation()
    assert cov1.is_correlation, "Should have been flagged as a correlation matrix"
    # Should be the same as the input
    assert np.allclose(cov1.toarray(), c), "Correlation matrix changed"

    # Save the current cov matrix
    _cov2matrix = cov2.toarray()
    # Convert to correlation
    cov2.to_correlation()
    # Should match cov1
    assert np.allclose(
        cov1.toarray(), cov2.toarray()
    ), "Correlation matrices should be identical"
    assert np.allclose(cov1.var / 2, cov2.var), "Variances incorrect"

    # Revert the correlation
    cov2.revert_correlation()
    # Should match original matrix
    assert np.allclose(_cov2matrix, cov2.toarray())


def test_newvar():
    c = mock_cov()
    cov1 = Covariance(array=c, impose_triu=True)
    var = np.full(c.shape[0], 4.0, dtype=float)
    cov2 = cov1.apply_new_variance(var)
    assert np.allclose(var, cov2.variance()), "Variance does not match request"
    cov2.to_correlation()
    assert np.allclose(
        cov1.toarray(), cov2.toarray()
    ), "Correlation matrices do not match"


def test_array():
    # Construct the Covariance matrix from a pre-calculated array
    c = mock_cov()
    covar = Covariance.from_array(c)
    # Should be the same as the identity matrix.
    assert np.array_equal(covar.toarray(), c), "Arrays should be identical"


def test_io():
    # Clean up in case of a failure
    ofile = Path("test_covar_io.fits")
    if ofile.is_file():
        ofile.unlink()

    # 1D
    cov = Covariance(array=mock_cov(), unit="km^2")
    # Write
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
    # Read
    _cov = Covariance.from_fits(ofile, quiet=True)
    # Arrays should be the same
    assert np.allclose(cov.toarray(), _cov.toarray()), "Bad 1D I/O"
    # Units should be the same
    assert cov.unit == _cov.unit, "Units changed"
    # Clean-up
    ofile.unlink()

    # ND
    raw_shape, c = mock_cov_3d()
    cov = Covariance(array=c, impose_triu=True, raw_shape=raw_shape)
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
    _cov = Covariance.from_fits(ofile, quiet=True)
    # Arrays should be the same
    assert np.allclose(cov.toarray(), _cov.toarray()), "Bad ND I/O"
    # Clean-up
    ofile.unlink()
