import os
import traceback
import warnings
from copy import deepcopy
from math import ceil, log10, prod

import numpy as np
from dask import array as da

from astropy.modeling import CompoundModel, models
from astropy.wcs.wcsapi import BaseLowLevelWCS

__all__ = ["parallel_fit_model_nd"]


def _pixel_to_world_values_block(*pixel, wcs=None):
    world = wcs.pixel_to_world_values(*pixel[::-1])[::-1]
    world = np.array(world)
    return world


def _wcs_to_world_dask(wcs, data):
    # Given a WCS and a data shape, return an iterable of dask arrays
    # representing the world coordinates of the array.
    pixel = tuple([np.arange(size) for size in data.shape])
    pixel_nd = da.meshgrid(*pixel, indexing="ij")
    world = da.map_blocks(
        _pixel_to_world_values_block,
        *pixel_nd,
        wcs=deepcopy(wcs),
        new_axis=0,
        chunks=(3,) + data.chunksize,
    )
    return tuple([world[idx] for idx in range(len(data.shape))])


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
    combined,
    block_info=None,
    model=None,
    fitter=None,
    world=None,
    diagnostics=None,
    diagnostics_path=None,
    iterating_shape=None,
    fitter_kwargs=None,
    iterating_axes=None,
    fitting_axes=None,
):
    """
    Function that gets passed to map_blocks and will fit models to a specific
    chunk of the data.
    """
    if fitter_kwargs is None:
        fitter_kwargs = {}

    # Start off by re-ordering axes so that iterating axes come first followed
    # by fitting axes
    original_axes = tuple([0] + [idx + 1 for idx in (iterating_axes + fitting_axes)])
    new_axes = tuple(range(combined.ndim))
    combined = da.moveaxis(combined, original_axes, new_axes)

    data = combined[0]
    if world == "arrays":
        parameters = combined[1 : -model.n_inputs]
        world_arrays = combined[-model.n_inputs :]
    else:
        parameters = combined[1:]

    if (
        combined.ndim == 0
        or combined.size == 0
        or block_info is None
        or block_info == []
    ):
        return parameters

    # Because of the way map_blocks works, we need to have all arrays passed
    # to map_blocks have the same shape, even though for the parameters this
    # means there are extra unneeded dimensions. We slice these out here.
    index = tuple([slice(None)] * (1 + len(iterating_axes)) + [0] * len(fitting_axes))
    parameters = parameters[index]

    # The world argument is used to pass through 1D arrays of world coordinates
    # (otherwise world_arrays is used) so if the model has more than one
    # dimension we need to make these arrays N-dimensional.
    if world != "arrays":
        if model.n_inputs > 1:
            world_values = np.meshgrid(*world, indexing="ij")
        else:
            world_values = world

    iterating_shape_chunk = data.shape[: len(iterating_axes)]

    for index in np.ndindex(iterating_shape_chunk):
        # If all data values are NaN, just set parameters to NaN and move on
        if np.all(np.isnan(data[index])):
            for ipar, name in enumerate(model.param_names):
                parameters[(ipar,) + index] = np.nan
            continue

        # Make a copy of the reference model and inject parameters
        model_i = model.copy()
        for ipar, name in enumerate(model_i.param_names):
            setattr(model_i, name, parameters[(ipar,) + index])

        output = diagnostics == "all"
        error = ""
        all_warnings = []

        if world == "arrays":
            world_values = tuple([w[index] for w in world_arrays])

        # Do the actual fitting
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                model_fit = fitter(model_i, *world_values, data[index], **fitter_kwargs)
                all_warnings.extend(w)
        except Exception as exc:
            model_fit = None
            if diagnostics == "failed":
                output = True
            error = traceback.format_exc()
            for ipar, name in enumerate(model_i.param_names):
                parameters[(ipar,) + index] = np.nan
        else:
            # Put fitted parameters back into parameters arrays. These arrays are
            # created in-memory by dask and are local to this process so should be
            # safe to modify in-place
            for ipar, name in enumerate(model_fit.param_names):
                parameters[(ipar,) + index] = getattr(model_fit, name).value

        if diagnostics == "failed+warn" and len(all_warnings) > 0:
            output = True

        if output:
            # Construct a folder name based on the iterating index. Currently i
            # i a 1-d index but we need to re-convert it back to an N-dimensional
            # index.

            index_abs = np.array(index) + np.array(
                [block_info[0]["array-location"][idx + 1][0] for idx in iterating_axes]
            )
            # index = tuple(int(idx) for idx in np.unravel_index(i_abs, iterating_shape))
            maxlen = int(ceil(log10(max(iterating_shape))))
            fmt = "{0:0" + str(maxlen) + "d}"
            index_folder = os.path.join(
                diagnostics_path, "_".join(fmt.format(idx) for idx in index_abs)
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
            if len(fitting_axes) == 1:  # 2 here because extra iterating dimension
                import matplotlib.pyplot as plt

                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.set_title(str(index))
                ax.plot(world_values[0], data[index], "k.")
                if model_fit is None:
                    ax.text(0.1, 0.9, "Fit failed!", color="r", transform=ax.transAxes)
                else:
                    xmodel = np.linspace(*ax.get_xlim(), 100)
                    if hasattr(world[0], "unit"):
                        xmodel = xmodel * world[0].unit
                    ax.plot(xmodel, model_fit(xmodel), color="r")
                fig.savefig(os.path.join(index_folder, "fit.png"))
                plt.close(fig)

    return parameters


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
    scheduler=None,
    fitter_kwargs=None,
    preserve_native_chunks=False,
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
    scheduler : str, optional
        If not specified, a local multi-processing scheduler will be
        used. If ``'default'``, whatever is the current default scheduler will be
        used. You can also set this to anything that would be passed to
        ``array.compute(scheduler=...)``
    fitter_kwargs : None or dict
        Keyword arguments to pass to the fitting when it is called.
    preserve_native_chunks : bool, optional
        If `True`, the native data chunks will be used, although an error will
        be raised if this chunk size does not include the whole fitting axes.
    """
    if scheduler is None:
        scheduler = "processes"

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

    if not isinstance(fitting_axes, tuple):
        fitting_axes = (fitting_axes,)

    # Check dimensionality
    if model.n_inputs != len(fitting_axes):
        raise ValueError(
            f"Model is {model.n_inputs}-dimensional, but got "
            f"{len(fitting_axes)} value(s) in fitting_axes="
        )

    for fi in fitting_axes:
        if fi <= -data.ndim or fi > data.ndim - 1:
            raise ValueError(
                f"Fitting index {fi} out of range for {data.ndim}-dimensional data"
            )

    if preserve_native_chunks and not isinstance(data, da.core.Array):
        raise ValueError(
            "Can only set preserve_native_chunks=True if input data is a dask array"
        )

    # Sanitize fitting_axes and determine iterating_axes
    ndim = data.ndim
    fitting_axes = tuple([(fi if fi >= 0 else ndim - fi) for fi in fitting_axes])
    iterating_axes = tuple([i for i in range(ndim) if i not in fitting_axes])

    # Determine the shape along the fitting dimensions and the iterating dimensions
    fitting_shape = tuple([data.shape[i] for i in fitting_axes])
    iterating_shape = tuple([data.shape[i] for i in iterating_axes])

    # Make sure the input array is a dask array
    if not isinstance(data, da.core.Array):
        data = da.asarray(data)

    world_arrays = []
    if isinstance(world, BaseLowLevelWCS):
        if world.pixel_n_dim != data.ndim:
            raise ValueError(
                f"The WCS pixel_n_dim ({world.pixel_n_dim}) does not match the data dimensionality ({data.ndim})"
            )

        # Note that in future we could in principle consider supporting cases
        # where the number of world dimensions does not match the number of
        # data dimensions, provided the model returns a different number of
        # outputs than it takes inputs. For example, one could consider fitting
        # a 1D dataset with a model that takes two inputs and returns one
        # output if the WCS provides two coordinates for each 1D pixel.
        # However, this is a very advanced and unusual use case, so we don't
        # cater for this for now.
        if world.world_n_dim != data.ndim:
            raise ValueError(
                "The WCS world_n_dim ({world.world_n_dim}) does not match the data dimensionality ({data.ndim})"
            )

        # Construct dask arrays of world coordinates for every pixel in the cube.
        # We will then iterate over this in map_blocks.
        world_arrays = _wcs_to_world_dask(world, data)

        # Extract world arrays for just fitting dimensions
        world_arrays = [world_arrays[idx] for idx in fitting_axes]

        world = "arrays"
    elif isinstance(world, tuple):
        for w in world:
            if w.shape != data.shape:
                raise ValueError(
                    f"arrays in world tuple should have same shape as data (expected {data.shape}, got {w.shape})"
                )
        # Extract world arrays for just fitting dimensions
        world_arrays = [da.asarray(world[idx]) for idx in fitting_axes]
        world = "arrays"
    elif isinstance(world, dict):
        # Re-index world if a dict, make it a tuple in order of fitting_axes
        for axis in fitting_axes:
            if axis not in world:
                raise KeyError(f"Values for axis {axis} missing from world")
            elif len(world[axis]) != data.shape[axis]:
                raise ValueError(
                    f"world[{axis}] has length {len(world[axis])} but data along dimension {axis} has length {data.shape[axis]}"
                )
        world = tuple([world[axis] for axis in fitting_axes])
    elif world is None:
        world = tuple([np.arange(size) for size in fitting_shape])

    if preserve_native_chunks:
        for idx in fitting_axes:
            if data.chunksize[idx] != data.shape[idx]:
                raise ValueError(
                    "When using preserve_native_chunks=True, the chunk size should match the data size along the fitting axes"
                )
    else:
        # Rechunk the array so that it is not chunked along the fitting axes
        chunk_shape = tuple(
            "auto" if idx in iterating_axes else -1 for idx in range(ndim)
        )
        if chunk_n_max:
            block_size_limit = chunk_n_max * prod(fitting_shape) * data.dtype.itemsize
        else:
            block_size_limit = None
        data = data.rechunk(chunk_shape, block_size_limit=block_size_limit)

    # Extract the parameters arrays from the model, in the order in which they
    # appear in param_names, convert to dask arrays, and broadcast to shape of
    # iterable axes. We need to rechunk to the same chunk shape as the data
    # so that chunks line up and for map_blocks to work properly.
    # https://github.com/dask/dask/issues/11188
    parameter_arrays = []
    for name in model.param_names:
        values = getattr(model, name).value
        array = da.broadcast_to(da.from_array(values), iterating_shape).reshape(
            iterating_shape
        )
        array = array.reshape(
            tuple(
                data.shape[idx] if idx in iterating_axes else 1 for idx in range(ndim)
            )
        )
        array = da.broadcast_to(array, data.shape)
        parameter_arrays.append(array)

    # Define a model with default parameters to pass in to fit_models_to_chunk without copying all the parameter data

    simple_model = _copy_with_new_parameters(model, {})

    # Define a single combined array that contains data, parameter arrays, and
    # optionally world arrays

    combined_array = da.stack([data] + parameter_arrays + world_arrays)

    # Make sure that the combined array is not chunked along the first axis. At
    # this point we are also assuming/hoping that the remainder of the chunk
    # size is given by that of ``data``.

    combined_array = combined_array.rechunk(
        (combined_array.shape[0],) + combined_array.chunksize[1:]
    )

    result = da.map_blocks(
        fit_models_to_chunk,
        combined_array,
        enforce_ndim=True,
        dtype=float,
        new_axis=0,
        drop_axis=fitting_axes,
        model=simple_model,
        fitter=fitter,
        world=world,
        diagnostics=diagnostics,
        diagnostics_path=diagnostics_path,
        iterating_shape=iterating_shape,
        iterating_axes=iterating_axes,
        fitting_axes=fitting_axes,
        fitter_kwargs=fitter_kwargs,
    )

    if scheduler == "default":
        compute_kwargs = {}
    else:
        compute_kwargs = {"scheduler": scheduler}

    parameter_arrays_fitted = result.compute(**compute_kwargs)

    # Set up new parameter arrays with fitted values
    parameters = {}
    for i, name in enumerate(model.param_names):
        parameters[name] = parameter_arrays_fitted[i].reshape(iterating_shape)

    # Instantiate new fitted model
    model_fitted = _copy_with_new_parameters(model, parameters, shape=iterating_shape)

    return model_fitted
