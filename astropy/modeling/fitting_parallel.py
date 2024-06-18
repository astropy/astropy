def parallel_fit_model_nd(*, model, fitter, data, fitting_axes, world, method='dask'):
    """
    Fit a model in parallel to an N-dimensional dataset.

    Axes in the N-dimensional dataset are considered to be either 'fitting
    axes' or 'iterating axes'. As a specific example, if fitting a
    spectral cube with two celestial and one spectral axis, then if fitting a
    1D model to each spectrum in the cube, the spectral axis would be a fitting
    axis and the celestial axes would be iterating axes.

    Parameters
    ----------
    model : :class:`astropy.modeling.Model`
        The model to fit, specifying the initial parameter values. The shape
        of the parameters should be broadcastable to the shape of the iterating
        axes.
    fitter : :class:`astropy.modeling.Fitter`
        The fitter to use in the fitting process.
    data : `numpy.ndarray` or `dask.array.core.Array`
        The N-dimensional data to fit.
    fitting_axes : int or tuple
        The axes to keep for the fitting (other axes will be sliced/iterated over)
    world : dict or APE-14-WCS
        This can be specified either as a dictionary mapping fitting axes to
        world axis values, or as a WCS for the whole cube. If the former, then
        the values in the dictionary can be either 1D arrays, or can be given
        as N-dimensional arrays with shape broadcastable to the data shape.
    method : { 'dask' }
        The framework to use for the parallelization.
        """
    pass
