# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import time

import numpy as np
from astropy import units as u

from ..modeling.core import FittableModel, custom_model
from ..extern.six.moves import range

__all__ = ['discretize_model']


class DiscretizationError(Exception):
    """
    Called when discretization of models goes wrong.
    """


class KernelSizeError(Exception):
    """
    Called when size of kernels is even.
    """


def add_kernel_arrays_1D(array_1, array_2):
    """
    Add two 1D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        new_array = array_1.copy()
        center = array_1.size // 2
        slice_ = slice(center - array_2.size // 2,
                       center + array_2.size // 2 + 1)
        new_array[slice_] += array_2
        return new_array
    elif array_2.size > array_1.size:
        new_array = array_2.copy()
        center = array_2.size // 2
        slice_ = slice(center - array_1.size // 2,
                       center + array_1.size // 2 + 1)
        new_array[slice_] += array_1
        return new_array
    return array_2 + array_1


def add_kernel_arrays_2D(array_1, array_2):
    """
    Add two 2D kernel arrays of different size.

    The arrays are added with the centers lying upon each other.
    """
    if array_1.size > array_2.size:
        new_array = array_1.copy()
        center = [axes_size // 2 for axes_size in array_1.shape]
        slice_x = slice(center[1] - array_2.shape[1] // 2,
                        center[1] + array_2.shape[1] // 2 + 1)
        slice_y = slice(center[0] - array_2.shape[0] // 2,
                        center[0] + array_2.shape[0] // 2 + 1)
        new_array[slice_y, slice_x] += array_2
        return new_array
    elif array_2.size > array_1.size:
        new_array = array_2.copy()
        center = [axes_size // 2 for axes_size in array_2.shape]
        slice_x = slice(center[1] - array_1.shape[1] // 2,
                        center[1] + array_1.shape[1] // 2 + 1)
        slice_y = slice(center[0] - array_1.shape[0] // 2,
                        center[0] + array_1.shape[0] // 2 + 1)
        new_array[slice_y, slice_x] += array_1
        return new_array
    return array_2 + array_1


def discretize_model(model, x_range, y_range=None, mode='center', factor=10):
    """
    Function to evaluate analytical model functions on a grid.

    So far the function can only deal with pixel coordinates.

    Parameters
    ----------
    model : `~astropy.modeling.FittableModel` or callable.
        Analytic model function to be discretized. Callables, which are not an
        instances of `~astropy.modeling.FittableModel` are passed to
        `~astropy.modeling.custom_model` and then evaluated.
    x_range : tuple
        x range in which the model is evaluated. The difference between the
        upper an lower limit must be a whole number, so that the output array
        size is well defined.
    y_range : tuple, optional
        y range in which the model is evaluated. The difference between the
        upper an lower limit must be a whole number, so that the output array
        size is well defined. Necessary only for 2D models.
    mode : str, optional
        One of the following modes:
            * ``'center'`` (default)
                Discretize model by taking the value
                at the center of the bin.
            * ``'linear_interp'``
                Discretize model by linearly interpolating
                between the values at the corners of the bin.
                For 2D models interpolation is bilinear.
            * ``'oversample'``
                Discretize model by taking the average
                on an oversampled grid.
            * ``'integrate'``
                Discretize model by integrating the model
                over the bin using `scipy.integrate.quad`.
                Very slow.
    factor : float or int
        Factor of oversampling. Default = 10.

    Returns
    -------
    array : `numpy.array`
        Model value array

    Notes
    -----
    The ``oversample`` mode allows to conserve the integral on a subpixel
    scale. Here is the example of a normalized Gaussian1D:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        from astropy.modeling.models import Gaussian1D
        from astropy.convolution.utils import discretize_model
        gauss_1D = Gaussian1D(1 / (0.5 * np.sqrt(2 * np.pi)), 0, 0.5)
        y_center = discretize_model(gauss_1D, (-2, 3), mode='center')
        y_corner = discretize_model(gauss_1D, (-2, 3), mode='linear_interp')
        y_oversample = discretize_model(gauss_1D, (-2, 3), mode='oversample')
        plt.plot(y_center, label='center sum = {0:3f}'.format(y_center.sum()))
        plt.plot(y_corner, label='linear_interp sum = {0:3f}'.format(y_corner.sum()))
        plt.plot(y_oversample, label='oversample sum = {0:3f}'.format(y_oversample.sum()))
        plt.xlabel('pixels')
        plt.ylabel('value')
        plt.legend()
        plt.show()


    """
    if not callable(model):
        raise TypeError('Model must be callable.')
    if not isinstance(model, FittableModel):
        model = custom_model(model)()
    ndim = model.n_inputs
    if ndim > 2:
        raise ValueError('discretize_model only supports 1-d and 2-d models.')

    if not float(np.diff(x_range)).is_integer():
        raise ValueError("The difference between the upper an lower limit of"
                         " 'x_range' must be a whole number.")

    if y_range:
        if not float(np.diff(y_range)).is_integer():
            raise ValueError("The difference between the upper an lower limit of"
                             " 'y_range' must be a whole number.")

    if ndim == 2 and y_range is None:
        raise ValueError("y range not specified, but model is 2-d")
    if ndim == 1 and y_range is not None:
        raise ValueError("y range specified, but model is only 1-d.")
    if mode == "center":
        if ndim == 1:
            return discretize_center_1D(model, x_range)
        elif ndim == 2:
            return discretize_center_2D(model, x_range, y_range)
    elif mode == "linear_interp":
        if ndim == 1:
            return discretize_linear_1D(model, x_range)
        if ndim == 2:
            return discretize_bilinear_2D(model, x_range, y_range)
    elif mode == "oversample":
        if ndim == 1:
            return discretize_oversample_1D(model, x_range, factor)
        if ndim == 2:
            return discretize_oversample_2D(model, x_range, y_range, factor)
    elif mode == "integrate":
        if ndim == 1:
            return discretize_integrate_1D(model, x_range)
        if ndim == 2:
            return discretize_integrate_2D(model, x_range, y_range)
    else:
        raise DiscretizationError('Invalid mode.')


def discretize_center_1D(model, x_range):
    """
    Discretize model by taking the value at the center of the bin.
    """
    x = np.arange(*x_range)
    return model(x)


def discretize_center_2D(model, x_range, y_range):
    """
    Discretize model by taking the value at the center of the pixel.
    """
    x = np.arange(*x_range)
    y = np.arange(*y_range)
    x, y = np.meshgrid(x, y)
    return model(x, y)


def discretize_linear_1D(model, x_range):
    """
    Discretize model by performing a linear interpolation.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    values_intermediate_grid = model(x)
    return 0.5 * (values_intermediate_grid[1:] + values_intermediate_grid[:-1])


def discretize_bilinear_2D(model, x_range, y_range):
    """
    Discretize model by performing a bilinear interpolation.
    """
    # Evaluate model 0.5 pixel outside the boundaries
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    x, y = np.meshgrid(x, y)
    values_intermediate_grid = model(x, y)

    # Mean in y direction
    values = 0.5 * (values_intermediate_grid[1:, :]
                    + values_intermediate_grid[:-1, :])
    # Mean in x direction
    values = 0.5 * (values[:, 1:]
                    + values[:, :-1])
    return values


def discretize_oversample_1D(model, x_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)

    values = model(x)

    # Reshape and compute mean
    values = np.reshape(values, (x.size // factor, factor))
    return values.mean(axis=1)[:-1]


def discretize_oversample_2D(model, x_range, y_range, factor=10):
    """
    Discretize model by taking the average on an oversampled grid.
    """
    # Evaluate model on oversampled grid
    x = np.arange(x_range[0] - 0.5 * (1 - 1 / factor),
                  x_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)

    y = np.arange(y_range[0] - 0.5 * (1 - 1 / factor),
                  y_range[1] + 0.5 * (1 + 1 / factor), 1. / factor)
    x_grid, y_grid = np.meshgrid(x, y)
    values = model(x_grid, y_grid)

    # Reshape and compute mean
    shape = (y.size // factor, factor, x.size // factor, factor)
    values = np.reshape(values, shape)
    return values.mean(axis=3).mean(axis=1)[:-1, :-1]


def discretize_integrate_1D(model, x_range):
    """
    Discretize model by integrating numerically the model over the bin.
    """
    from scipy.integrate import quad
    # Set up grid
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    values = np.array([])

    # Integrate over all bins
    for i in range(x.size - 1):
        values = np.append(values, quad(model, x[i], x[i + 1])[0])
    return values


def discretize_integrate_2D(model, x_range, y_range):
    """
    Discretize model by integrating the model over the pixel.
    """
    from scipy.integrate import dblquad
    # Set up grid
    x = np.arange(x_range[0] - 0.5, x_range[1] + 0.5)
    y = np.arange(y_range[0] - 0.5, y_range[1] + 0.5)
    values = np.empty((y.size - 1, x.size - 1))

    # Integrate over all pixels
    for i in range(x.size - 1):
        for j in range(y.size - 1):
            values[j, i] = dblquad(lambda y, x: model(x, y), x[i], x[i + 1],
                                   lambda x: y[j], lambda x: y[j + 1])[0]
    return values

def profile_memory_usage_convolve_fft(sizes=[10,50,100,128,200,256,300,512,600,768],
                                      kernelsizefractions=[1./8., 1./4., 1./2., 1],
                                      dtype='float64',
                                      complex_dtype='complex',
                                     ):
    """
    Utility to profile the memory usage of convolve_fft.  Will loop over a
    range of parameters and image sizes and keep track of how much memory is
    used and how long the convolutions take.

    Requires memory_profiler and timeit

    Parameters
    ----------
    sizes : list
        List of symmetric 2D image sizes.  i.e., a size of 256 implies a
        256x256 image.
    kernelsizefractions : list
        List of floats.  Each image will be convolved with a kernel of size
        ``(size * kernelsizefraction)``.  These should be <1 in general.
    dtype : str
        The dtype of the arrays to be convolved.  The default is 'float64',
        as with numpy, but changing this will change the initial amount of
        memory located.  Using a smaller ``dtype`` (e.g., ``float32`) may
        result in a *larger* memory use within convolve_fft, unless
        ``complex_dtype`` is also set to be lower.
    complex_dtype : str
        The dtype to use during the convolution.  Options range from
        ``complex32`` to ``complex256``.
    """
    from memory_profiler import memory_usage
    import timeit

    # this can't be imported at module level because... I don't know why.
    from .convolve import convolve_fft, convolve

    results = []

    for size in sizes:
        kernelsizes = [int(np.ceil(size * ksv)) for ksv in kernelsizefractions]
        for kernelsize,kernelsizefraction in zip(kernelsizes,kernelsizefractions):
            for psf_pad in (True,False):
                for fft_pad in (True,False):
                    for convolver in ('convolve', 'convolve_fft'):
                        
                        if kernelsize % 2 == 0:
                            # required for direct convolution
                            kernelsize = kernelsize - 1
                            if kernelsize == 0:
                                kernelsize = 1

                        x=np.ones([size,size], dtype=dtype)
                        y=np.ones([kernelsize,kernelsize], dtype=dtype)

                        # time the memory profiling test to determine how many
                        # times to repeat the timing profiling later
                        t0 = time.time()
                        # this checks the memory usage many times over the course of the function call
                        try:
                            if convolver == 'convolve_fft':
                                memuse = memory_usage((convolve_fft, (x,y,),
                                                       dict(psf_pad=psf_pad,
                                                            fft_pad=fft_pad,
                                                            complex_dtype=complex_dtype,
                                                            boundary=None)),
                                                      interval=0.01)
                            else:
                                memuse = memory_usage((convolve, (x,y,),
                                                       dict(boundary=None)),
                                                      interval=0.01)
                        except ValueError:
                            # don't crash if the files get to be too large
                            continue
                        peak_memuse = (max(memuse)-min(memuse))*u.MB

                        if time.time() - t0 < 1:
                            # for tests that go fast, repeat them a few times for robustness' sake
                            number = 10
                        else:
                            # for tests that go slow, do them once - we're accurate enough already
                            number = 1

                        timersetup = ('from astropy.convolution import {convolver}; import numpy as np;'
                                      'x=np.ones([{size},{size}], dtype="{dtype}");'
                                      'y=np.ones([{kernelsize},{kernelsize}], dtype="{dtype}");')
                        if convolver == 'convolve_fft':
                            timercall = ('convolve_fft(x,y,psf_pad={psf_pad},fft_pad={fft_pad},'
                                         ' complex_dtype="{complex_dtype}", boundary=None)')
                        elif convolver == 'convolve':
                            timercall = ('convolve(x,y, boundary=None)')

                        timeuse = timeit.repeat(timercall.format(psf_pad=psf_pad,
                                                                 fft_pad=fft_pad,
                                                                 convolver=convolver,
                                                                 complex_dtype=complex_dtype),
                                                setup=(timersetup.format(size=size,
                                                                         convolver=convolver,
                                                                         kernelsize=kernelsize,
                                                                         dtype=dtype)),
                                                repeat=3, number=number)

                        these_results = {'peak_memuse':peak_memuse,
                                         'psf_pad': psf_pad,
                                         'fft_pad': fft_pad,
                                         'size': size,
                                         'kernelsize': kernelsize,
                                         'size_mb': (x.nbytes*u.byte).to(u.MB),
                                         'kernelsize_mb': (y.nbytes*u.byte).to(u.MB),
                                         'meantime': np.mean(timeuse)*u.s,
                                         'kernelsizefraction': kernelsizefraction,
                                         'convolver': convolver,
                                         }


                        print("convolver={convolver:12} psf_pad={psf_pad:2}, fft_pad={fft_pad:2}, size={size:5d}  kernelsize={kernelsize:5d}"
                              " Mem usage max={peak_memuse:8.1f}.  Array, kernel size are {size_mb:5.1f} and {kernelsize_mb:5.1f}."
                              " Time use was {meantime:6.3f}"
                              .format(**these_results))

                        results.append(these_results)

    return results

def plot_convolve_fft_profiling_results(results):
    """
    Create a plot of the results from profile_memory_usage_convolve_fft
    """
    import pylab as pl

    # fig = pl.figure(1, figsize=(14,8))
    # fig.clf()
    # ax1 = fig.add_subplot(2,1,1)
    # ax2 = fig.add_subplot(2,1,2)

    kernelsizefractions = list(set([x['kernelsizefraction'] for x in results]))

    phasefig = pl.figure(5)
    phasefig.clf()
    phaseax = phasefig.gca()


    labeltemplate1 = "psf={psf:1} fft={fft:1} KSV={ksv:0.3f} {convolver}"
    labeltemplate = "KSV={ksv:0.3f} {convolver}"
    # kept for possible future reuse, but this is still too crowded
    styles = {(True,False): '--',
              (True,True): '-',
              (False,True): ":",
              (False,False): "-."}
    figure = {(True,False): pl.figure(1),
              (True,True): pl.figure(2),
              (False,True): pl.figure(3),
              (False,False): pl.figure(4),
             }
    axes = {}
    for ind,fig in figure.items():
        fig.clf()
        fig.suptitle("psf_pad={0} fft_pad={1}".format(str(ind[0]), str(ind[1])))
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        axes[ind] = (ax1,ax2)

    for convolver, linestyle in [('convolve','--'), ('convolve_fft','-')]:
        for kernelsizefraction in kernelsizefractions:
            for psf_pad in (True,False):
                for fft_pad in (True,False):

                    sizes = np.array([x['size']
                                      for x in results
                                      if (x['psf_pad'] == psf_pad and
                                          x['fft_pad'] == fft_pad and
                                          x['convolver'] == convolver and
                                          x['kernelsizefraction'] == kernelsizefraction)])
                    times = np.array([x['meantime'].value
                                      for x in results
                                      if (x['psf_pad'] == psf_pad and
                                          x['fft_pad'] == fft_pad and
                                          x['convolver'] == convolver and
                                          x['kernelsizefraction'] == kernelsizefraction)])
                    memus = np.array([x['peak_memuse'].value
                                      for x in results
                                      if (x['psf_pad'] == psf_pad and
                                          x['fft_pad'] == fft_pad and
                                          x['convolver'] == convolver and
                                          x['kernelsizefraction'] == kernelsizefraction)])
                    if len(sizes) == 0:
                        continue

                    ax1,ax2 = axes[(psf_pad, fft_pad)]
                    ax1.plot(sizes, times, label=labeltemplate.format(psf=str(psf_pad)[0],
                                                                      fft=str(fft_pad)[0],
                                                                      convolver=convolver,
                                                                      ksv=kernelsizefraction),
                             #linestyle=styles[(psf_pad,fft_pad)],
                             linestyle=linestyle,
                             linewidth=2.0,
                             alpha=0.5,
                            )
                    ax2.plot(sizes, memus, label=labeltemplate.format(psf=str(psf_pad)[0],
                                                                      fft=str(fft_pad)[0],
                                                                      convolver=convolver,
                                                                      ksv=kernelsizefraction),
                             #linestyle=styles[(psf_pad,fft_pad)],
                             linestyle=linestyle,
                             linewidth=2.0,
                             alpha=0.5,
                            )

                    phaseax.plot(memus, times, 'o',
                                 label=labeltemplate1.format(psf=str(psf_pad)[0],
                                                             fft=str(fft_pad)[0],
                                                             convolver=convolver,
                                                             ksv=kernelsizefraction),
                                 markeredgecolor='none', alpha=0.5)


    for ind,fig in figure.items():
        ax1,ax2 = axes[ind]
        ax1.set_ylabel("Execution time (s)")
        ax2.set_ylabel("Peak memory use (MB)")
        ax2.set_xlabel("1D size of image")

        for ax in (ax1,ax2):
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            xmax = np.max([L.get_data()[0].max() for L in ax.get_lines()])
            ax.set_xlim(0, xmax)
            ymax = np.max([L.get_data()[1].max() for L in ax.get_lines()])
            ax.set_ylim(0, ymax*1.1)

        ax1.legend(loc='center left', bbox_to_anchor=(1.05, 0.50))
        #ax2.legend(loc='best', bbox_to_anchor)

    box = phaseax.get_position()
    xmax = np.max([L.get_data()[0].max() for L in phaseax.get_lines()])
    phaseax.set_xlim(0, xmax)
    ymax = np.max([L.get_data()[1].max() for L in phaseax.get_lines()])
    phaseax.set_ylim(0, ymax*1.1)
    phaseax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

    phaseax.set_ylabel("Execution time (s)")
    phaseax.set_xlabel("Peak memory use (MB)")
    phaseax.legend(loc='center left', bbox_to_anchor=(1.05, 0.50))

def parse_profiling_printed_results(filename, valid_ksfs=None):
    """
    In case the profiler crashes, you can copy and paste the results into a
    file and parse it with this tool

    Parameters
    ----------
    valid_ksvs : list
        Valid kernel size fractions.  Meant to fix roundoff errors.
    """
    with open(filename, 'r') as f:
        results = [{'psf_pad': bool(int(row[12])),
                    'fft_pad': bool(int(row[27])),
                    'size': int(row[36:40]),
                    'kernelsize': int(row[54:59]),
                    'peak_memuse': float(row[74:81])*u.MB,
                    'meantime': float(row[154:][:-2])*u.s}
                   for row in f.readlines()]
    for r in results:
        ksf = float(r['kernelsize'])/r['size']
        if valid_ksfs is not None:
            closest = np.argmin(np.abs(ksf - np.array(valid_ksfs)))
            ksf = valid_ksfs[closest]
        r['kernelsizefraction'] = ksf

    return results
