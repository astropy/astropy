import os
import traceback
import warnings
from math import ceil, log10, prod

import numpy as np
from dask import array as da

from astropy.modeling import CompoundModel, models

__all__ = ["parallel_fit_model_nd"]


def _copy_with_new_parameters(model, parameters, shape=None):
    if isinstance(model, CompoundModel):
        if shape is None:
            new_model = model.copy()
        else:
            new_model = _compound_model_with_array_parameters(model, shape)
        for name, value in parameters.items():
            setattr(new_model, name, value)
    else:
        constraints = {}
        for constraint in model.parameter_constraints:
            constraints[constraint] = getattr(model, constraint)
        # HACK: we need a more general way, probably in astropy.modeling,
        # to do this kind of copy.
        if isinstance(model, models.Polynomial1D):
            args = (model.degree,)
        else:
            args = ()
        new_model = model.__class__(*args, **parameters, **constraints)
    return new_model


def _compound_model_with_array_parameters(model, shape):
    if isinstance(model, CompoundModel):
        return CompoundModel(
            model.op,
            _compound_model_with_array_parameters(model.left, shape),
            _compound_model_with_array_parameters(model.right, shape),
        )
    else:
        parameters = {name: np.zeros(shape) for name in model.param_names}
        return _copy_with_new_parameters(model, parameters)


def fit_models_to_chunk(
    data,
    *parameters,
    block_info=None,
    model=None,
    fitter=None,
    world=None,
    diagnostics=None,
    diagnostics_path=None,
    iterating_shape=None,
):
    """
    Function that gets passed to map_blocks and will fit models to a specific
    chunk of the data.
    """
    if data.ndim == 0 or data.size == 0 or block_info is None or block_info == []:
        return np.array(parameters)

    parameters = np.array(parameters)

    index = tuple([slice(None), slice(None)] + [0] * (parameters.ndim - 2))
    parameters = parameters[index]

    # Iterate over the first axis and fit each model in turn
    for i in range(data.shape[0]):
        # Make a copy of the reference model and inject parameters
        model_i = model.copy()
        for ipar, name in enumerate(model_i.param_names):
            setattr(model_i, name, parameters[ipar, i])

        output = diagnostics == "all"
        error = ""
        all_warnings = []

        # Do the actual fitting
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                model_fit = fitter(model_i, *world, data[i])
                all_warnings.extend(w)
        except Exception as exc:
            model_fit = None
            if diagnostics == "failed":
                output = True
            error = traceback.format_exc()
            for ipar, name in enumerate(model_i.param_names):
                parameters[ipar, i] = np.nan
        else:
            # Put fitted parameters back into parameters arrays. These arrays are
            # created in-memory by dask and are local to this process so should be
            # safe to modify in-place
            for ipar, name in enumerate(model_fit.param_names):
                parameters[ipar, i] = getattr(model_fit, name).value

        if diagnostics == "failed+warn" and len(all_warnings) > 0:
            output = True

        if output:
            # Construct a folder name based on the iterating index. Currently i
            # i a 1-d index but we need to re-convert it back to an N-dimensional
            # index.

            index = tuple(int(idx) for idx in np.unravel_index(i, iterating_shape))
            maxlen = int(ceil(log10(max(iterating_shape))))
            fmt = "{0:0" + str(maxlen) + "d}"
            index_folder = os.path.join(
                diagnostics_path, "_".join(fmt.format(idx) for idx in index)
            )
            os.makedirs(index_folder, exist_ok=True)

            # Output error, if any
            if error:
                with open(os.path.join(index_folder, "error.log"), "w") as f:
                    f.write(error)

            if all_warnings:
                with open(os.path.join(index_folder, "warn.log"), "w") as f:
                    for warning in all_warnings:
                        f.write(f"{warning}\n")

            # Make a plot, if model is 1D
            if data.ndim == 2:  # 2 here because extra iterating dimension
                import matplotlib.pyplot as plt

                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title(str(index))
                ax.plot(world[0], data[i], "k.")
                if model_fit is None:
                    ax.text(0.1, 0.9, "Fit failed!", color="r", transform=ax.transAxes)
                else:
                    xmodel = np.linspace(*ax.get_xlim(), 100) * world[0].unit
                    ax.plot(xmodel, model_fit(xmodel), color="r")
                fig.savefig(os.path.join(index_folder, "fit.png"))
                plt.close(fig)

    parameters = parameters.reshape((parameters.shape[0], parameters.shape[1], 1))

    return np.broadcast_to(parameters, (parameters.shape[0],) + data.shape)


def parallel_fit_model_nd(
    *,
    model,
    fitter,
    data,
    fitting_axes,
    world=None,
    chunk_n_max=None,
    diagnostics=None,
    diagnostics_path=None,
    use_default_scheduler=False,
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
    world : `None` or dict or APE-14-WCS
        This can be specified either as a dictionary mapping fitting axes to
        world axis values, or as a WCS for the whole cube. If the former, then
        the values in the dictionary can be either 1D arrays, or can be given
        as N-dimensional arrays with shape broadcastable to the data shape. If
        not specified, the fitting is carried out in pixel coordinates.
    method : { 'dask' }
        The framework to use for the parallelization.
    chunk_n_max : int
        Maximum number of fits to include in a chunk. If this is made too large,
        then the workload will not be split properly over processes, and if it is
        too small it may be inefficient.
    diagnostics : { None | 'failed' | 'failed+warn' | 'all' }, optional
        Whether to output diagnostic information for fits. This can be either
        `None` (nothing), ``'failed'`` (output information for failed fits), or
        ``'all'`` (output information for all fits).
    diagnostics_path : str, optional
        If `diagnostics` is not `None`, this should be the path to a folder in
        which a folder will be made for each fit that is output.
    use_default_scheduler : bool, optional
        If `False` (the default), a local multi-processing scheduler will be
        used. If `True`, whatever is the current default scheduler will be
        used. Set this to `True` if using e.g. dask.distributed.
    """
    if diagnostics in (None, "failed", "failed+warn", "all"):
        if diagnostics is not None:
            if diagnostics_path is None:
                raise ValueError("diagnostics_path should be set")
            else:
                os.makedirs(diagnostics_path, exist_ok=True)
    else:
        raise ValueError(
            "diagnostics should be None, 'failed', 'failed+warn', or 'all'"
        )

    # Sanitize fitting_axes and determine iterating_axes
    ndim = data.ndim
    if not isinstance(fitting_axes, tuple):
        fitting_axes = (fitting_axes,)
    fitting_axes = tuple([(fi if fi >= 0 else ndim - fi) for fi in fitting_axes])
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
    elif world is None:
        world = tuple([np.arange(size) for size in fitting_shape])

    # NOTE: dask tends to choose chunks that are too large for this kind of
    # problem, so if chunk_n_max is not specified, we try and determine a
    # sensible value of chunk_n_max.

    # Determine a reasonable chunk_n_max if not specified
    if not chunk_n_max:
        chunk_n_max = max(10, data.shape[0] // 100)

    # Rechunk the array so that it is not chunked along the fitting axes
    chunk_shape = (chunk_n_max,) + (-1,) * len(fitting_axes)
    data = data.rechunk(chunk_shape)

    # Extract the parameters arrays from the model, in the order in which they
    # appear in param_names, convert to dask arrays, and broadcast to shape of
    # iterable axes. We need to rechunk to the same chunk shape as the data
    # so that chunks line up and for map_blocks to work properly.
    # https://github.com/dask/dask/issues/11188
    parameter_arrays = []
    for name in model.param_names:
        values = getattr(model, name).value
        array = (
            da.broadcast_to(da.from_array(values), iterating_shape)
            .reshape(iterating_shape)
            .ravel()
        )
        array = array.reshape(array.shape + (1,) * len(fitting_shape))
        array = da.broadcast_to(array, array.shape[:1] + fitting_shape).rechunk(
            data.chunksize
        )
        parameter_arrays.append(array)

    # Define a model with default parameters to pass in to fit_models_to_chunk without copying all the parameter data

    simple_model = _copy_with_new_parameters(model, {})

    result = da.map_blocks(
        fit_models_to_chunk,
        data,
        *parameter_arrays,
        # chunks=(len(parameter_arrays),) + data.chunksize,
        enforce_ndim=True,
        dtype=float,
        new_axis=0,
        model=simple_model,
        fitter=fitter,
        world=world,
        diagnostics=diagnostics,
        diagnostics_path=diagnostics_path,
        iterating_shape=iterating_shape,
    )

    if use_default_scheduler:
        compute_kwargs = {}
    else:
        compute_kwargs = {"scheduler": "processes"}

    parameter_arrays_fitted = result.compute(**compute_kwargs)

    parameter_arrays_fitted = parameter_arrays_fitted[:, :, 0]

    # Set up new parameter arrays with fitted values
    parameters = {}
    for i, name in enumerate(model.param_names):
        parameters[name] = parameter_arrays_fitted[i].reshape(iterating_shape)

    # Instantiate new fitted model
    model_fitted = _copy_with_new_parameters(model, parameters, shape=iterating_shape)

    return model_fitted
