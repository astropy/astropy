Convolution Performance
=======================

Convolution can be costly in terms of memory use or computation.  These
considerations are usually in some tension with one another.

For convolution with a small kernel, `~astropy.convolution.convolve` is
generally faster than FFT-based convolution
(`~astropy.convolution.convolve_fft`), while for large kernels, i.e.,
those where the kernel size approaches the array size, FFT-based convolution
is often faster.  However, the FFT-based approach requires having two full-sized
arrays in memory, and so even when it is faster, it can sometimes impose
a prohibitive memory requirement.

FFT-based convolution is also strongly dependent on array size.  Images with
sizes factorable into multiples of 2, 3, and 5 are more efficient than others,
and images with prime sizes are much more expensive.  For example, convolving
two arrays with size [59,59] or [61,61] is about 3-4x more expensive than
convolving two arrays with size [60,60] on a modern Intel processor.

Because of the cost associated with awkward array sizes,
`~astropy.convolution.convolve_fft` has an ``fft_pad`` option that will
increase the size of the array to the nearest "fast" size using scipy's
`next_fast_len
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.next_fast_len.html>`_.
This option increases the memory usage but will often speed up the execution of
a convolution operation.
