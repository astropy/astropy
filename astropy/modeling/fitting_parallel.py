from math import prod

import numpy as np
from dask import array as da

__all__ = ["parallel_fit_model_nd"]


def fit_models_to_chunk(
    data,
    *parameters,
    block_info=None,
    model=None,
    fitter=None,
    world=None,
):
    """
    Function that gets passed to map_blocks and will fit models to a specific
    chunk of the data.
    """
    if data.ndim == 0 or data.size == 0 or block_info is None or block_info == []:
        return np.array(parameters)

    # BUG: fitter doesn't work correctly if pickled/unpickled
    from astropy.modeling.fitting import LMLSQFitter

    fitter = LMLSQFitter()

    # Iterate over the first axis and fit each model in turn
    for i in range(data.shape[0]):
        # Make a copy of the reference model and inject parameters
        model_i = model.copy()
        for ipar, name in enumerate(model_i.param_names):
            setattr(model_i, name, parameters[ipar][i])

        # Do the actual fitting
        model_fit = fitter(model_i, *world, data[i])

        # Put fitted parameters back into parameters arrays. These arrays are
        # created in-memory by dask and are local to this process so should be
        # safe to modify in-place
        for ipar, name in enumerate(model_fit.param_names):
            parameters[ipar][i] = getattr(model_fit, name).value

    return np.array(parameters)


def parallel_fit_model_nd(
    *,
    model,
    fitter,
    data,
    fitting_axes,
    world,
    chunk_n_max=None,
):
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
    chunk_n_max : int
        Maximum number of fits to include in a chunk. If this is made too large,
        then the workload will not be split properly over processes, and if it is
        too small it may be inefficient.
    """
    # Sanitize fitting_axes and determine iterating_axes
    ndim = data.ndim
    if not isinstance(fitting_axes, tuple):
        fitting_axes = (fitting_axes,)
    fitting_axes = tuple([(fi if fi > 0 else ndim - fi) for fi in fitting_axes])
    iterating_axes = tuple([i for i in range(ndim) if i not in fitting_axes])

    # Determine the shape along the fitting dimensions and the iterating dimensions
    fitting_shape = tuple([data.shape[i] for i in fitting_axes])
    iterating_shape = tuple([data.shape[i] for i in iterating_axes])

    # Make sure the input array is a dask array
    data = da.asarray(data)

    # Move all iterating dimensions to the front and flatten. We do this so
    # that fit_models_to_chunk can be agnostic of the complexity of the
    # iterating dimensions.
    original_axes = iterating_axes + fitting_axes
    new_axes = tuple(range(ndim))
    data = da.moveaxis(data, original_axes, new_axes)
    data = data.reshape((prod(iterating_shape),) + fitting_shape)

    # Re-index world if a dict, make it a tuple in order of fitting_axes
    if isinstance(world, dict):
        world = tuple([world[axis] for axis in fitting_axes])

    # Rechunk the array so that it is not chunked along the fitting axes
    chunk_shape = ("auto",) + (-1,) * len(fitting_axes)
    if chunk_n_max:
        data = data.rechunk(
            chunk_shape, block_size_limit=chunk_n_max * data.dtype.itemsize
        )
    else:
        data = data.rechunk(chunk_shape)

    # Extract the parameters arrays from the model, in the order in which they
    # appear in param_names, convert to dask arrays, and broadcast to shape of
    # iterable axes. We need to rechunk to the same chunk shape as the data
    # so that chunks line up and for map_blocks to work properly.
    parameter_arrays = []
    for name in model.param_names:
        values = getattr(model, name).value
        array = (
            da.broadcast_to(da.from_array(values), iterating_shape)
            .reshape(iterating_shape)
            .ravel()
            .rechunk(data.chunksize)
        )
        parameter_arrays.append(array)

    # Define a model with default parameters to pass in to fit_models_to_chunk without copying all the parameter data

    constraints = {}
    for constraint in model.parameter_constraints:
        constraints[constraint] = getattr(model, constraint)

    simple_model = model.__class__(**constraints)

    result = da.map_blocks(
        fit_models_to_chunk,
        data,
        *parameter_arrays,
        chunks=(3,) + data.chunksize,
        enforce_ndim=True,
        dtype=float,
        new_axis=0,
        model=simple_model,
        fitter=fitter,
        world=world,
    )

    parameter_arrays_fitted = result.compute(scheduler="processes")

    # Set up new parameter arrays with fitted values
    parameters = {}
    for i, name in enumerate(model.param_names):
        parameters[name] = parameter_arrays_fitted[i].reshape(iterating_shape)

    # Instantiate new fitted model
    model_fitted = model.__class__(**parameters, **constraints)

    return model_fitted
