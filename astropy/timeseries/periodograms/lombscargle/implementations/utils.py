import warnings
from math import factorial

import numpy as np

from astropy.utils.compat.optional_deps import HAS_SCIPY
from astropy.utils.exceptions import AstropyUserWarning

if HAS_SCIPY:
    from scipy.fft import next_fast_len
    from scipy.special import jv, lambertw


def bitceil(N):
    """
    Find the bit (i.e. power of 2) immediately greater than or equal to N
    Note: this works for numbers up to 2 ** 64.
    Roughly equivalent to int(2 ** np.ceil(np.log2(N))).
    """
    return 1 << int(N - 1).bit_length()


def extirpolate(x, y, N=None, M=4):
    """
    Extirpolate the values (x, y) onto an integer grid range(N),
    using lagrange polynomial weights on the M nearest points.

    Parameters
    ----------
    x : array-like
        array of abscissas
    y : array-like
        array of ordinates
    N : int
        number of integer bins to use. For best performance, N should be larger
        than the maximum of x
    M : int
        number of adjoining points on which to extirpolate.

    Returns
    -------
    yN : ndarray
        N extirpolated values associated with range(N)

    Examples
    --------
    >>> rng = np.random.default_rng(0)
    >>> x = 100 * rng.random(20)
    >>> y = np.sin(x)
    >>> y_hat = extirpolate(x, y)
    >>> x_hat = np.arange(len(y_hat))
    >>> f = lambda x: np.sin(x / 10)
    >>> np.allclose(np.sum(y * f(x)), np.sum(y_hat * f(x_hat)))
    True

    Notes
    -----
    This code is based on the C implementation of spread() presented in
    Numerical Recipes in C, Second Edition (Press et al. 1989; p.583).
    """
    x, y = map(np.ravel, np.broadcast_arrays(x, y))

    if N is None:
        N = int(np.max(x) + 0.5 * M + 1)

    # Now use legendre polynomial weights to populate the results array;
    # This is an efficient recursive implementation (See Press et al. 1989)
    result = np.zeros(N, dtype=y.dtype)

    # first take care of the easy cases where x is an integer
    integers = x % 1 == 0
    np.add.at(result, x[integers].astype(int), y[integers])
    x, y = x[~integers], y[~integers]

    # For each remaining x, find the index describing the extirpolation range.
    # i.e. ilo[i] < x[i] < ilo[i] + M with x[i] in the center,
    # adjusted so that the limits are within the range 0...N
    ilo = np.clip((x - M // 2).astype(int), 0, N - M)
    numerator = y * np.prod(x - ilo - np.arange(M)[:, np.newaxis], 0)
    denominator = factorial(M - 1)

    for j in range(M):
        if j > 0:
            denominator *= j / (j - M)
        ind = ilo + (M - 1 - j)
        np.add.at(result, ind, numerator / (denominator * (x - ind)))
    return result


if HAS_SCIPY:

    def get_lra_params(x, N, eps):
        """
        Compute parameters for the Low Rank Approximation (LRA) NUFFT.

        Parameters
        ----------
        x : array-like
            Non-uniform sample locations, shape (n,).
        N : int
            Required number of grid points in the frequency domain.
        eps : float
            Desired relative error tolerance.

        Returns
        -------
        s : ndarray, int
            Wrapped grid indices, shape (n,).
        gamma : float
            Maximum absolute deviation from nearest uniform grid point.
        K : int
            Number of FFTs required to achieve desired precision (eps) when round-off errors are ignored.
        """
        s = np.mod(np.rint(N * x), N).astype(int)
        er = (N * x - np.rint(N * x) + 0.5) % 1 - 0.5
        gamma = np.max(np.abs(er))
        if gamma <= eps:
            K = 1
        else:
            val = lambertw(np.log(10 / eps) / (7 * gamma)).real
            K = int(np.ceil(5 * gamma * np.exp(val)))
        return s, gamma, K

    def chebyshev_polynomials(n, x):
        """
        Compute Chebyshev polynomials of the first kind T_0(x) through T_n(x).

        Parameters
        ----------
        n : int
            Maximum polynomial degree to compute.
        x : array-like
            Points at which to evaluate polynomials, shape (m,).

        Returns
        -------
        T : ndarray
            Polynomial values, shape (m, n+1). T[:, k] contains T_k(x).
        """
        x = np.asarray(x)
        N = len(x)
        T = np.empty((N, n + 1), dtype=np.float64)
        T[:, 0] = 1.0
        if n > 0:
            T[:, 1] = x
        for k in range(1, n):
            T[:, k + 1] = 2 * x * T[:, k] - T[:, k - 1]
        return T

    def bessel_coefficients(K, gamma):
        """
        Compute Bessel function coefficients for NUFFT low-rank approximation.

        Parameters
        ----------
        K : int
            Truncation parameter for series expansion.
        gamma : float
            Maximum grid deviation parameter from get_lra_params.

        Returns
        -------
        cfs : ndarray
            Complex coefficient matrix, shape (K, K).

        Notes
        -----
        Only elements where (p+q) is even are non-zero.
        """
        cfs = np.zeros((K, K), dtype=np.complex128)
        arg = -gamma * np.pi * 0.5
        # Create full K x K arrays for p and q using meshgrid.
        p, q = np.meshgrid(np.arange(K), np.arange(K), indexing="ij")
        mask = (q - p) % 2 == 0
        cfs[mask] = (
            4
            * (1j) ** q[mask]
            * jv((p[mask] + q[mask]) * 0.5, arg)
            * jv((q[mask] - p[mask]) * 0.5, arg)
        )
        cfs[0, :] *= 0.5
        cfs[:, 0] *= 0.5
        return cfs

    def construct_UV(x, gamma, K, N):
        """
        Construct low-rank approximation matrices for NUFFT.

        Parameters
        ----------
        x : array-like
            Non-uniform sample locations, shape (N,).
        gamma : float
            Maximum grid deviation from get_lra_params.
        K : int
            Truncation parameter from get_lra_params.
        N : int
            Number of frequency grid points.

        Returns
        -------
        U : ndarray
            Sample-space factor matrix, shape (M, K).
        V : ndarray
            Frequency-space factor matrix, shape (N, K).

        Notes
        -----
        Constructs matrices such that the NUFFT kernel exp(-2πi x ω) ≈ U @ V^H.
        U combines Chebyshev polynomials of (x - grid)/gamma with Bessel coefficients.
        V contains Chebyshev polynomials of frequency grid points.
        """
        # Compute the periodic deviation for each sample point.
        er = (N * x - np.rint(N * x) + 0.5) % 1 - 0.5

        # Compute Chebyshev polynomials on a scaled variable for U.
        if gamma == 0:
            Tcheb_U = np.ones((len(x), K))
        else:
            Tcheb_U = chebyshev_polynomials(K - 1, er / gamma)
        B = bessel_coefficients(K, gamma)
        U = np.exp(-1j * np.pi * er)[:, None] * (Tcheb_U @ B)

        # For V: compute Chebyshev polynomials on the frequency variable.
        omega = np.arange(N)
        X = 2.0 * omega / N - 1
        Tcheb_V = chebyshev_polynomials(K - 1, X)
        V = Tcheb_V.astype(np.complex128)

        return U, V

    def nufft1_lra(x, y, N, eps):
        """
        Adjoint NUFFT  (Type-I) using Low Rank Approximation.

        Computes: s[k] = sum_{j=0}^{M-1} y[j] * exp(2πi x[j] k) for k = 0, ..., N-1

        Parameters
        ----------
        x : array-like
            Non-uniform sample locations in shape (M,).
        y : array-like
            Complex input values at non-uniform points, shape (M,).
        N : int
            Number of desired frequency bins.
        eps : float
            Desired relative error tolerance.

        Returns
        -------
        y : ndarray
            Complex Fourier coefficients, shape (N,).

        Notes
        -----
        Implements the algorithm from Ruiz-Antolin & Townsend (2018) using
        a low-rank approximation of the NUFFT kernel.

        This implementation based on the MIT-licensed FastTransforms.jl Julia library.
        """
        # Not required, although it makes the computation noticeably faster
        Nfft = next_fast_len(max(N, len(x)))

        # Compute grid indices and parameters for the low-rank approximation.
        T, gamma, K = get_lra_params(x, Nfft, eps)
        U, V = construct_UV(x, gamma, K, Nfft)

        # Scattering step: accumulate contributions via bin-counting.
        temp = np.zeros((Nfft, K), dtype=np.complex128)
        product = np.conjugate(U) * y[:, None]
        for k in range(K):
            real_part = np.bincount(T, weights=np.real(product[:, k]), minlength=Nfft)
            imag_part = np.bincount(T, weights=np.imag(product[:, k]), minlength=Nfft)
            temp[:, k] = real_part + 1j * imag_part

        # Compute the inverse FFT (with normalization by Nfft) and combine with V.
        temp_ifft = np.fft.ifft(temp, axis=0) * Nfft
        s = np.sum(np.conjugate(V) * temp_ifft, axis=1)

        # Return the result truncated to the desired frequency grid length.
        return s[:N]


def trig_sum(
    t,
    h,
    df,
    N,
    f0=0,
    freq_factor=1,
    oversampling=5,
    use_fft=True,
    Mfft=4,
    prefer_lra=True,
    eps=5e-13,
):
    """Compute (approximate) trigonometric sums for a number of frequencies.

    This routine computes weighted sine and cosine sums::

        S_j = sum_i { h_i * sin(2 pi * f_j * t_i) }
        C_j = sum_i { h_i * cos(2 pi * f_j * t_i) }

    Where f_j = freq_factor * (f0 + j * df) for the values j in 1 ... N.
    The sums can be computed either by a brute force O[N^2] method, or
    by an FFT-based O[Nlog(N)] method.

    Parameters
    ----------
    t : array-like
        array of input times
    h : array-like
        array weights for the sum
    df : float
        frequency spacing
    N : int
        number of frequency bins to return
    f0 : float, optional
        The low frequency to use
    freq_factor : float, optional
        Factor which multiplies the frequency
    use_fft : bool
        if True, use the approximate FFT algorithm to compute the result.
    prefer_lra : bool
        if True and SciPy is available, use the more accurate Low Rank Approximation by Ruiz-Antolin and Townsend.
        If False or SciPy is not available use Press & Rybicki's Lagrangian extirpolation instead.
    oversampling : int (default = 5)
        oversampling freq_factor for the approximation; roughly the number of
        time samples across the highest-frequency sinusoid. This parameter
        contains the trade-off between accuracy and speed. Not referenced
        if use_fft is False or prefer_lra is True and SciPy is available.
    Mfft : int
        The number of adjacent points to use in the FFT approximation.
        Not referenced if use_fft is False  or prefer_lra is True and SciPy is available.
    eps : float
        Desired relative error tolerance. Not referenced
        if use_fft is False or SciPy is not available.

    Returns
    -------
    S, C : ndarray
        summation arrays for frequencies f = df * np.arange(1, N + 1)
    """
    if df <= 0:
        raise ValueError("df must be positive")

    # Ensure t and h are 1D arrays.
    t, h = map(np.ravel, np.broadcast_arrays(t, h))

    if prefer_lra and not HAS_SCIPY:
        warnings.warn(
            "SciPy is not installed. The low-rank approximation (LRA) method for NUFFT "
            "requires SciPy and will be disabled. Astropy will use the less accurate "
            "Lagrangian interpolation-based method instead.\n"
            "To improve accuracy, install SciPy.\n"
            "To silence this warning, explicitly set `prefer_lra=False`.",
            AstropyUserWarning,
        )
        prefer_lra = False
        print("\n")

    df *= freq_factor
    f0 *= freq_factor

    if use_fft:
        t0 = t.min()
        if prefer_lra:
            if f0 != 0:
                h = h * np.exp(2j * np.pi * f0 * (t - t0))

            tnorm = (t - t0) * df

            fftgrid = nufft1_lra(tnorm, h, N, eps)

            if t0 > 0:
                f = f0 + df * np.arange(N)
                fftgrid = fftgrid * np.exp(2j * np.pi * t0 * f)

            C = fftgrid.real
            S = fftgrid.imag

        else:
            # required size of fft is the power of 2 above the oversampling rate
            Nfft = bitceil(int(N * oversampling))

            if f0 > 0:
                h = h * np.exp(2j * np.pi * f0 * (t - t0))

            tnorm = ((t - t0) * Nfft * df) % Nfft
            grid = extirpolate(tnorm, h, Nfft, Mfft)

            fftgrid = np.fft.ifft(grid)[:N]
            if t0 != 0:
                f = f0 + df * np.arange(N)
                fftgrid *= np.exp(2j * np.pi * t0 * f)

            C = Nfft * fftgrid.real
            S = Nfft * fftgrid.imag

    else:
        f = f0 + df * np.arange(N)
        C = np.dot(h, np.cos(2 * np.pi * f * t[:, np.newaxis]))
        S = np.dot(h, np.sin(2 * np.pi * f * t[:, np.newaxis]))

    return S, C
