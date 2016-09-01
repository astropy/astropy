from __future__ import print_function, division

import warnings
from math import factorial
import numpy as np
from ....extern.six.moves import range, map


def add_at(arr, ind, vals):
    """Utility that computes np.add.at()

    The fast version is available only in Numpy 1.8+; for older versions of
    numpy this defaults to a slower computation.
    """
    if hasattr(np.ufunc, 'at'):
        return np.add.at(arr, ind, vals)
    else:
        warnings.warn("Using slow replacement for numpy.add.at(). "
                      "For ~100x faster results update to numpy 1.8+")
        arr = np.asarray(arr)
        ind, vals = np.broadcast_arrays(ind, vals)
        unq = np.unique(ind)
        arr[unq] += [vals[ind == i].sum() for i in unq]


def bitceil(N):
    """
    Find the bit (i.e. power of 2) immediately greater than or equal to N
    Note: this works for numbers up to 2 ** 64.
    Roughly equivalent to int(2 ** np.ceil(np.log2(N)))
    """
    if hasattr(int, 'bit_length'):
        # Python 2.7 and 3.x
        return 1 << int(N - 1).bit_length()
    else:
        N = int(N) - 1
        for i in [1, 2, 4, 8, 16, 32]:
            N |= N >> i
        return N + 1


def extirpolate(x, y, N=None, M=4):
    """
    Extirpolate the values (x, y) onto an integer grid range(N),
    using lagrange polynomial weights on the M nearest points.
    Parameters
    ----------
    x : array_like
        array of abscissas
    y : array_like
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

    Example
    -------
    >>> rng = np.random.RandomState(0)
    >>> x = 100 * rng.rand(20)
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
    integers = (x % 1 == 0)
    add_at(result, x[integers].astype(int), y[integers])
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
        add_at(result, ind, numerator / (denominator * (x - ind)))
    return result


def trig_sum(t, h, df, N, f0=0, freq_factor=1,
             oversampling=5, use_fft=True, Mfft=4):
    """Compute (approximate) trigonometric sums for a number of frequencies
    This routine computes weighted sine and cosine sums:
        S_j = sum_i { h_i * sin(2 pi * f_j * t_i) }
        C_j = sum_i { h_i * cos(2 pi * f_j * t_i) }
    Where f_j = freq_factor * (f0 + j * df) for the values j in 1 ... N.
    The sums can be computed either by a brute force O[N^2] method, or
    by an FFT-based O[Nlog(N)] method.

    Parameters
    ----------
    t : array_like
        array of input times
    h : array_like
        array weights for the sum
    df : float
        frequency spacing
    N : int
        number of frequency bins to return
    f0 : float (optional, default=0)
        The low frequency to use
    freq_factor : float (optional, default=1)
        Factor which multiplies the frequency
    use_fft : bool
        if True, use the approximate FFT algorithm to compute the result.
        This uses the FFT with Press & Rybicki's Lagrangian extirpolation.
    oversampling : int (default = 5)
        oversampling freq_factor for the approximation; roughly the number of
        time samples across the highest-frequency sinusoid. This parameter
        contains the tradeoff between accuracy and speed. Not referenced
        if use_fft is False.
    Mfft : int
        The number of adjacent points to use in the FFT approximation.
        Not referenced if use_fft is False.

    Returns
    -------
    S, C : ndarrays
        summation arrays for frequencies f = df * np.arange(1, N + 1)
    """
    df *= freq_factor
    f0 *= freq_factor

    if df <= 0:
        raise ValueError("df must be positive")
    t, h = map(np.ravel, np.broadcast_arrays(t, h))

    if use_fft:
        Mfft = int(Mfft)
        if Mfft <= 0:
            raise ValueError("Mfft must be positive")

        # required size of fft is the power of 2 above the oversampling rate
        Nfft = bitceil(N * oversampling)
        t0 = t.min()

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
