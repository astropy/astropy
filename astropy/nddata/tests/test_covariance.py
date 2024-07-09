from pathlib import Path

import numpy as np

from astropy.nddata.covariance import Covariance


def test_samples():
    # Fixed-seed RNG for repeatability
    rng = np.random.default_rng(seed=8001)

    # Build a bogus covariance matrix
    m = np.zeros(10, dtype=float)
    c = (
        np.diag(np.full(10 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(10, 1.0, dtype=float), k=0)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(10 - 2, 0.2, dtype=float), k=2)
    )

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
    # Build a bogus covariance matrix
    c = (
        np.diag(np.full(10 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(10, 1.0, dtype=float), k=0)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(10 - 2, 0.2, dtype=float), k=2)
    )

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


def test_array():
    # Construct the Covariance matrix from a pre-calculated array
    c = (
        np.diag(np.full(10 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(10, 1.0, dtype=float), k=0)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(10 - 2, 0.2, dtype=float), k=2)
    )
    covar = Covariance.from_array(c)
    # Should be the same as the identity matrix.
    assert np.array_equal(covar.toarray(), c), "Arrays should be identical"


def test_io():
    rng = np.random.default_rng(seed=8001)

    # Clean up in case of a failure
    ofile = Path("test_covar_io.fits")
    if ofile.is_file():
        ofile.unlink()

    # Build a bogus covariance matrix
    m = np.zeros(10, dtype=float)
    c = (
        np.diag(np.full(10 - 2, 0.2, dtype=float), k=-2)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=-1)
        + np.diag(np.full(10, 1.0, dtype=float), k=0)
        + np.diag(np.full(10 - 1, 0.5, dtype=float), k=1)
        + np.diag(np.full(10 - 2, 0.2, dtype=float), k=2)
    )

    # Draw samples
    s = rng.multivariate_normal(m, c, size=100000)

    # Instantiate
    covar = Covariance.from_samples(s.T, cov_tol=0.1)
    # Write
    covar.write(ofile)
    # Read
    _covar = Covariance.from_fits(ofile)
    # Should be the same
    assert np.allclose(covar.toarray(), _covar.toarray()), "Bad I/O"
    # Clean-up
    ofile.unlink()
