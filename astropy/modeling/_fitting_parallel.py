import traceback
import warnings
from copy import deepcopy
from math import ceil, log10, prod
from pathlib import Path

import numpy as np

import astropy.units as u
from astropy.modeling.utils import _combine_equivalency_dict
from astropy.nddata import NDUncertainty, StdDevUncertainty, support_nddata
from astropy.wcs.wcsapi import BaseHighLevelWCS, BaseLowLevelWCS
from astropy.wcs.wcsapi.wrappers import SlicedLowLevelWCS

__all__ = ["FitInfoArrayContainer", "parallel_fit_dask"]


def _pixel_to_world_values_block(*pixel, wcs=None):
    """
    Convert a block of pixel values to world values using a WCS. This is for
    use in map_blocks.
    """
    world = wcs.array_index_to_world_values(*pixel)
    world = np.array(world)
    return world


def _wcs_to_world_dask(wcs, shape, chunks=None):
    """
    Given a WCS and a data shape, return an iterable of dask arrays
    representing the world coordinates of the array.
    """
    import dask.array as da

    pixel = tuple([np.arange(size) for size in shape])
    pixel_nd = da.meshgrid(*pixel, indexing="ij")
    world = da.map_blocks(
        _pixel_to_world_values_block,
        *pixel_nd,
        wcs=deepcopy(wcs),
        new_axis=0,
        chunks=(len(shape),) + chunks,
    )
    return tuple([world[idx] for idx in range(len(shape))])


def _copy_with_new_parameters(model, parameters):
    # Make a copy of the model, setting the parameters to new values
    model_new = model.copy()
    model_new._reset_parameters(**parameters)
    return model_new


class FitInfoArrayContainer:
    """
    This class is intended to contain the object array of all fit_info values
    and provide a convenience method to access specific items from fit_info
    as arrays.
    """

    def __init__(self, fit_info_array):
        self._fit_info_array = fit_info_array

    @property
    def shape(self):
        return self._fit_info_array.shape

    @property
    def ndim(self):
        return self._fit_info_array.ndim

    def __getitem__(self, item):
        result = self._fit_info_array[item]
        if hasattr(result, "ndim"):
            return FitInfoArrayContainer(result)
        else:
            return result

    def get_property_as_array(self, name):
        """
        Return an array of one of the fit information properties

        Parameters
        ----------
        name : str
            The name of a property present on the individual fit information
            objects.
        """
        array = None
        for index in np.ndindex(self._fit_info_array.shape):
            fit_info = self._fit_info_array[index]
            if fit_info is not None:
                value = np.array(getattr(fit_info, name))
                if array is None:
                    array = np.zeros(self.shape + value.shape, dtype=value.dtype)
                if value.shape != array.shape[self.ndim :]:
                    raise ValueError(
                        "Property {name} does not have consistent shape in all fit_info"
                    )
                array[index] = value
        return array

    @property
    def properties(self):
        """
        The properties available to query with :meth:`~astropy.modeling.fitting.FitInfoArrayContainer.get_property_as_array`
        """
        # Find the first non-None .fit_info
        for index in np.ndindex(self._fit_info_array.shape):
            fit_info = self._fit_info_array[index]
            if fit_info is not None:
                return tuple(sorted(fit_info))
        return ()


class FitInfoSubset(dict):
    def __getattr__(self, attr):
        return self[attr]


def _fit_models_to_chunk(
    data,
    *arrays,
    block_info=None,
    model=None,
    fitter=None,
    world=None,
    diagnostics=None,
    diagnostics_path=None,
    diagnostics_callable=None,
    iterating_shape=None,
    fitter_kwargs=None,
    iterating_axes=None,
    fitting_axes=None,
    weights_specified=None,
    fit_info=None,
):
    """
    Function that gets passed to map_blocks and will fit models to a specific
    chunk of the data.
    """
    if fitter_kwargs is None:
        fitter_kwargs = {}

    # Start off by re-ordering axes so that iterating axes come first followed
    # by fitting axes
    original_axes = tuple(idx for idx in (iterating_axes + fitting_axes))
    new_axes = tuple(range(data.ndim))
    data = np.moveaxis(data, original_axes, new_axes)
    arrays = [np.moveaxis(array, original_axes, new_axes) for array in arrays]

    if weights_specified:
        weights = arrays[0]
        arrays = arrays[1:]
    else:
        weights = None

    # World coordinates can be specified either as Nd world arrays (in which
    # case the world kwarg is set to `None`), or passed in via the world kwarg
    # (if the world coordinates are given as 1D arrays)
    if world is None:
        parameters = arrays[: -model.n_inputs]
        world_arrays = arrays[-model.n_inputs :]
    else:
        parameters = arrays

    # Make the parameters into an Nd array, as this is what we will return. We
    # then modify this array in-place in the rest of the function.
    parameters = np.array(parameters)

    # In some cases, dask calls this function with empty arrays, so we can
    # take a short-cut here.
    if data.ndim == 0 or data.size == 0 or block_info is None or block_info == []:
        return parameters

    # Because of the way map_blocks works, we need to have all arrays passed
    # to map_blocks have the same shape, even though for the parameters this
    # means there are extra unneeded dimensions. We slice these out here.
    index = tuple([slice(None)] * (1 + len(iterating_axes)) + [0] * len(fitting_axes))
    parameters = parameters[index]

    # Transform array to object array and add one more index along the first
    # dimension so that we can store the fit_info
    if fit_info:
        parameters = parameters.astype(object)
        parameters = np.pad(parameters, [(0, 1)] + [(0, 0)] * (parameters.ndim - 1))

    # The world argument is used to pass through 1D arrays of world coordinates
    # (otherwise world_arrays is used) so if the model has more than one
    # dimension we need to make these arrays N-dimensional.
    if world is not None:
        if model.n_inputs > 1:
            world_values = np.meshgrid(*world, indexing="ij")
        else:
            world_values = world

    iterating_shape_chunk = data.shape[: len(iterating_axes)]

    model_i = model.copy()

    for index in np.ndindex(iterating_shape_chunk):
        # If all data values are NaN, just set parameters to NaN and move on
        if np.all(np.isnan(data[index])):
            for ipar in range(len(model.param_names)):
                parameters[(ipar,) + index] = np.nan
            continue

        # Inject parameters into model
        model_i._reset_parameters(
            **{
                name: parameters[(ipar,) + index]
                for ipar, name in enumerate(model.param_names)
            },
        )

        output = diagnostics == "all"
        error = ""
        all_warnings = []

        if world is None:
            world_values = tuple([w[index] for w in world_arrays])

        if weights is None:
            weights_kwargs = {}
        else:
            weights_kwargs = dict(weights=weights[index])

        # Do the actual fitting - note that we can use inplace=True here to
        # speed things up by avoiding an unnecessary copy, since we don't need
        # to retain the original parameter values.
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                model_fit = fitter(
                    model_i,
                    *world_values,
                    data[index],
                    inplace=True,
                    **weights_kwargs,
                    **fitter_kwargs,
                )
                all_warnings.extend(w)
        except Exception as exc:
            model_fit = None
            if diagnostics is not None and diagnostics.startswith("error"):
                output = True
            error = traceback.format_exc()
            for ipar in range(len(model_i.param_names)):
                parameters[(ipar,) + index] = np.nan

            parameters[(-1,) + index] = None
        else:
            # Put fitted parameters back into parameters arrays. These arrays are
            # created in-memory by dask and are local to this process so should be
            # safe to modify in-place
            for ipar, name in enumerate(model_fit.param_names):
                parameters[(ipar,) + index] = getattr(model_fit, name).value

            if fit_info is True:
                parameters[(-1,) + index] = fitter.fit_info
            elif fit_info:
                fit_info_dict = {}
                for key in fit_info:
                    if hasattr(fitter.fit_info, key):
                        fit_info_dict[key] = getattr(fitter.fit_info, key)
                    else:
                        raise AttributeError(
                            f"fit_info on fitter has no attribute '{key}'"
                        )
                parameters[(-1,) + index] = FitInfoSubset(fit_info_dict)

        if diagnostics == "error+warn" and len(all_warnings) > 0:
            output = True

        if output:
            # Construct a folder name based on the iterating index. Currently i
            # i a 1-d index but we need to re-convert it back to an N-dimensional
            # index.

            index_abs = np.array(index) + np.array(
                [block_info[0]["array-location"][idx][0] for idx in iterating_axes]
            )
            maxlen = ceil(log10(max(iterating_shape)))
            fmt = "{0:0" + str(maxlen) + "d}"
            index_folder = Path(diagnostics_path).joinpath(
                "_".join(fmt.format(idx) for idx in index_abs)
            )
            index_folder.mkdir(parents=True, exist_ok=True)

            # Output error, if any
            if error:
                index_folder.joinpath("error.log").write_text(error)

            if all_warnings:
                index_folder.joinpath("warn.log").write_text(
                    "".join(f"{warning}\n" for warning in all_warnings)
                )

            if diagnostics_callable is not None:
                diagnostics_callable(
                    index_folder,
                    world_values,
                    data[index],
                    None if weights is None else weights[index],
                    model_fit,
                    fitter_kwargs,
                )

    return parameters


class ParameterContainer:
    """
    This is an array container intended to be passed to dask's ``from_array``.

    The initial parameter values need to be broadcast up to the data shape so
    that map_blocks can then iterate over both the data and parameters. We
    need to control the final chunking so that it matches the data. However,
    rather than use dask to do the broadcasting and rechunking, which creates a
    complex graph and results in high memory usage, this class can be used
    instead to do all the broadcasting on-the-fly with Numpy as needed and keeps
    the dask graph simple.
    """

    def __init__(self, values, iterating_shape, iterating_axes, data_shape):
        self._values = values
        self.shape = data_shape
        self.ndim = len(data_shape)
        self.dtype = float
        self._iterating_shape = iterating_shape
        self._iterating_axes = iterating_axes

    def __getitem__(self, item):
        values = np.broadcast_to(self._values, self._iterating_shape)
        values = values.reshape(
            tuple(
                self.shape[idx] if idx in self._iterating_axes else 1
                for idx in range(self.ndim)
            )
        )
        values = np.broadcast_to(values, self.shape)
        return values[item]


@support_nddata(wcs="world", uncertainty="weights", unit="data_unit")
def parallel_fit_dask(
    *,
    model,
    fitter,
    data,
    data_unit=None,
    weights=None,
    mask=None,
    fitting_axes=None,
    world=None,
    chunk_n_max=None,
    diagnostics=None,
    diagnostics_path=None,
    diagnostics_callable=None,
    scheduler=None,
    fitter_kwargs=None,
    preserve_native_chunks=False,
    equivalencies=None,
    fit_info=False,
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
    fitter : :class:`astropy.modeling.fitting.Fitter`
        The fitter to use in the fitting process.
    data : `numpy.ndarray` or `dask.array.core.Array`
        The N-dimensional data to fit.
    data_units : `astropy.units.Unit`
        Units for the data array, for when the data array is not a ``Quantity``
        instance.
    weights : `numpy.ndarray`, `dask.array.core.Array` or `astropy.nddata.NDUncertainty`
        The weights to use in the fitting. See the documentation for specific
        fitters for more information about the meaning of weights.
        If passed as a `.NDUncertainty` object it will be converted to a
        `.StdDevUncertainty` and then passed to the fitter as 1 over that.
    mask : `numpy.ndarray`
        A boolean mask to be applied to the data.
    fitting_axes : int or tuple
        The axes to keep for the fitting (other axes will be sliced/iterated over)
    world : `None` or tuple or APE-14-WCS
        This can be specified either as a tuple of world coordinates for each
        fitting axis, or as WCS for the whole cube. If specified as a tuple,
        the values in the tuple can be either 1D arrays, or can be given as
        N-dimensional arrays with shape broadcastable to the data shape. If
        specified as a WCS, the WCS can have any dimensionality so long as it
        matches the data. If not specified, the fitting is carried out in pixel
        coordinates.
    chunk_n_max : int
        Maximum number of fits to include in a chunk. If this is made too
        large, then the workload will not be split properly over processes, and
        if it is too small it may be inefficient. If not specified, this will
        default to 500.
    diagnostics : { None | 'error' | 'error+warn' | 'all' }, optional
        Whether to output diagnostic information for fits. This can be either
        `None` (nothing), ``'error'`` (output information for fits that raised
        exceptions), or ``'all'`` (output information for all fits).
    diagnostics_path : str, optional
        If ``diagnostics`` is not `None`, this should be the path to a folder in
        which a folder will be made for each fit that is output.
    diagnostics_callable : callable
        By default, any warnings or errors are output to ``diagnostics_path``.
        However, you can also specify a callable that can e.g. make a plot or
        write out information in a custom format. The callable should take the
        following arguments: the path to the subfolder of ``diagnostics_path``
        for the specific index being fit, a list of the coordinates passed to
        the fitter, the data array, the weights array (or `None` if no weights
        are being used), the model that was fit (or `None` if the fit errored),
        and a dictionary of other keyword arguments passed to the fitter.
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
    equivalencies : list of tuple
        Any equivalencies to take into account in unit conversions
    fit_info : bool or str or iterable, optional
        Option to control whether fit information set on the ``.fit_info``
        attribute of the fitters for individual fits should be concatenated
        and set on the ``.fit_info`` on the input fitter object. The options
        are as follows:

        * `False`: don't set ``.fit_info`` on the fitter
        * `True`: set ``.fit_info`` on the fitter with all available information
           from individual fits
        * An iterable of strings: only save the specific attributes mentioned
          in the iterable

        If not `False`, the ``.fit_info`` attribute on the fitter will be set
        to a `FitInfoArrayContainer` object which can be used to query the
        fit information for individual fits. Otherwise, ``.fit_info`` will be
        left unchanged.
    """
    try:
        import dask
        import dask.array as da
    except ImportError:  # pragma: no cover
        raise ImportError("dask is required for this function")

    if scheduler is None:
        scheduler = "processes"

    if diagnostics in (None, "error", "error+warn", "all"):
        if diagnostics is not None:
            if diagnostics_path is None:
                raise ValueError("diagnostics_path should be set")
            else:
                Path(diagnostics_path).mkdir(parents=True, exist_ok=True)
    else:
        raise ValueError("diagnostics should be None, 'error', 'error+warn', or 'all'")

    original_fit_info = deepcopy(fitter.fit_info)

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

    if preserve_native_chunks:
        if not isinstance(data, da.core.Array):
            raise TypeError(
                "Can only set preserve_native_chunks=True if input data is a dask array"
            )
        if weights is not None and not isinstance(weights, da.core.Array):
            raise TypeError(
                "Can only set preserve_native_chunks=True if input weights is a dask array (if specified)"
            )

    if isinstance(weights, NDUncertainty):
        weights = weights.represent_as(StdDevUncertainty)
        weights = 1 / weights.array

    if mask is not None:
        imask = np.logical_not(mask).astype(float)
        if weights is None:
            weights = imask
        else:
            weights *= imask

    # Sanitize fitting_axes and determine iterating_axes
    ndim = data.ndim
    fitting_axes = tuple([(fi if fi >= 0 else ndim - fi) for fi in fitting_axes])
    iterating_axes = tuple([i for i in range(ndim) if i not in fitting_axes])

    # Determine the shape along the fitting dimensions and the iterating dimensions
    fitting_shape = tuple([data.shape[i] for i in fitting_axes])
    iterating_shape = tuple([data.shape[i] for i in iterating_axes])

    if data_unit is None and isinstance(data, u.Quantity):
        data_unit = data.unit
        data = data.value

    if preserve_native_chunks:
        for idx in fitting_axes:
            if data.chunksize[idx] != data.shape[idx]:
                raise ValueError(
                    "When using preserve_native_chunks=True, the chunk size should match the data size along the fitting axes"
                )
        if weights is not None and data.chunksize != weights.chunksize:
            raise ValueError(
                "When using preserve_native_chunks=True, the weights should have the same chunk size as the data"
            )
    else:
        # Rechunk the array so that it is not chunked along the fitting axes
        chunk_shape = tuple(
            "auto" if idx in iterating_axes else -1 for idx in range(ndim)
        )

        if chunk_n_max is None:
            chunk_n_max = 500

        block_size_limit = chunk_n_max * prod(fitting_shape) * data.dtype.itemsize

        if isinstance(data, da.core.Array):
            data = data.rechunk(chunk_shape, block_size_limit=block_size_limit)
        else:
            with dask.config.set({"array.chunk-size": block_size_limit}):
                data = da.from_array(data, chunks=chunk_shape, name="data")

        if weights is not None:
            if isinstance(weights, da.core.Array):
                weights = weights.rechunk(data.chunksize)
            else:
                weights = da.from_array(weights, chunks=data.chunksize, name="weights")

    world_arrays = False
    if isinstance(world, BaseHighLevelWCS):
        world = world.low_level_wcs

    if isinstance(world, BaseLowLevelWCS):
        if world.pixel_n_dim != data.ndim:
            raise ValueError(
                f"The WCS pixel_n_dim ({world.pixel_n_dim}) does not match the number of dimensions in the data ({data.ndim})"
            )

        # Note that in future we could in principle consider supporting cases
        # where the number of world dimensions does not match the number of
        # data dimensions, provided the model returns a different number of
        # outputs than it takes inputs. For example, one could consider fitting
        # a 1D dataset with a model that takes two inputs and returns one
        # output if the WCS provides two coordinates for each 1D pixel.
        # However, this is a very advanced and unusual use case, so we don't
        # cater for this for now.

        fitting_world = SlicedLowLevelWCS(
            world,
            [slice(None) if i in fitting_axes else 0 for i in range(world.pixel_n_dim)],
        )
        if fitting_world.world_n_dim != len(fitting_axes):
            raise ValueError(
                "The number of WCS world axes corresponding to the fitting axes "
                f"({fitting_world.world_n_dim}) does not match the number of fitting axes ({len(fitting_axes)})"
            )

        world_units = list(map(u.Unit, fitting_world.world_axis_units[::-1]))

        # Construct dask arrays of world coordinates for every pixel in the cube.
        # We will then iterate over this in map_blocks.
        # NOTE: This returns in world (cartesian) order
        world_dask_arrays = _wcs_to_world_dask(world, data.shape, chunks=data.chunksize)

        # Extract world arrays for just fitting dimensions
        fitting_pixel_axes = np.arange(data.ndim)[::-1][np.array(fitting_axes)]
        world_idx = [
            np.argwhere(world.axis_correlation_matrix[:, fpa])[:, 0][0]
            for fpa in fitting_pixel_axes
        ]
        world = [world_dask_arrays[idx] for idx in world_idx]
        world_arrays = True

    elif isinstance(world, tuple):
        # If world is a tuple then we allow N inputs where N is the number of fitting_axes
        # Each array in the tuple should with be broadcastable to the shape of the fitting_axes
        # or it should be one dimensional and the broadcasting can happen later
        if len(world) != len(fitting_axes):
            raise ValueError(
                f"The number of world arrays ({len(world)}) must match "
                f"number of fitting axes ({len(fitting_axes)})"
            )
        world = list(world)
        world_units = []
        for iw, w in enumerate(world):
            if (unit := getattr(w, "unit", None)) is not None:
                world[iw] = w.value
            world_units.append(unit)

        if all(w.ndim == 1 for w in world):
            for i, (w, fit_shape) in zip(fitting_axes, zip(world, fitting_shape)):
                if w.shape[0] != fit_shape:
                    raise ValueError(
                        f"world[{i}] has length {w.shape[0]} but data along "
                        f"dimension {i} has length {fit_shape}"
                    )
            world_arrays = False
        else:
            for w in world:
                try:
                    w = np.broadcast_shapes(w.shape, data.shape)
                except ValueError as e:
                    raise ValueError(
                        f"The arrays in the world tuple should be broadcastable to "
                        f"the shape of the data (expected {data.shape}), got {w.shape})"
                    ) from e
            # Extract world arrays for just fitting dimensions
            world = [da.asarray(world[idx]) for idx in fitting_axes]
            world_arrays = True
    elif world is None:
        world = tuple([np.arange(size) for size in fitting_shape])
        world_units = [None] * len(fitting_axes)
    else:
        raise TypeError("world should be None, a WCS object or a tuple of arrays")

    if model._has_units or data_unit is not None:
        # We now combine any instance-level input equivalencies with user
        # specified ones at call-time.

        input_units_equivalencies = _combine_equivalency_dict(
            model.inputs, equivalencies, model.input_units_equivalencies
        )

        # If input_units is defined, we transform the input data into those
        # expected by the model. We hard-code the input names 'x', and 'y'
        # here since FittableModel instances have input names ('x',) or
        # ('x', 'y')

        if model.input_units is None:
            target_units = world_units[:]
        else:
            target_units = [
                model.input_units[model.inputs[i]] for i in range(model.n_inputs)
            ]
            world = [
                unit.to(
                    target_units[i],
                    equivalencies=input_units_equivalencies[model.inputs[i]],
                    value=w,
                )
                if unit is not None
                else w
                for i, (w, unit) in enumerate(zip(world, world_units))
            ]

        # Create a dictionary mapping the real model inputs and outputs
        # names to the data. This remapping of names must be done here, after
        # the input data is converted to the correct units.
        rename_data = {}
        rename_data[model.inputs[0]] = (0,) * target_units[0]
        rename_data[model.outputs[0]] = (0,) * data_unit
        if len(world) == 2:
            rename_data[model.inputs[1]] = (0,) * target_units[1]
        else:
            rename_data["z"] = None

        # We now strip away the units from the parameters, taking care to
        # first convert any parameters to the units that correspond to the
        # input units (to make sure that initial guesses on the parameters)
        # are in the right unit system
        model = model.without_units_for_data(**rename_data)
        add_back_units = True

    else:
        world_units = tuple(None for w in world)
        add_back_units = False

    # Extract the parameters arrays from the model, in the order in which they
    # appear in param_names. We need to broadcast these up to the data shape so
    # that map_blocks can then iterate over both the data and parameters. We
    # need to rechunk to the same chunk shape as the data so that chunks line
    # up and for map_blocks to work properly as noted in
    # https://github.com/dask/dask/issues/11188. However, rather than use dask
    # operations to broadcast these up, which creates a complex graph and
    # results in high memory usage, we use a ParameterContainer which does the
    # broadcasting on-the-fly as needed.
    parameter_arrays = []
    for name in model.param_names:
        values = getattr(model, name).value
        parameter_arrays.append(
            da.from_array(
                ParameterContainer(values, iterating_shape, iterating_axes, data.shape),
                chunks=data.chunksize,
                name="parameter-" + name,
            )
        )

    # Define a model with default parameters to pass in to _fit_models_to_chunk without copying all the parameter data

    simple_model = _copy_with_new_parameters(model, {})

    weights_array = [] if weights is None else [weights]

    result = da.map_blocks(
        _fit_models_to_chunk,
        data,
        *weights_array,
        *parameter_arrays,
        *(world if world_arrays else []),
        enforce_ndim=True,
        dtype=object if fit_info else float,
        drop_axis=fitting_axes,
        model=simple_model,
        fitter=fitter,
        new_axis=0,
        world=world if not world_arrays else None,
        diagnostics=diagnostics,
        diagnostics_path=diagnostics_path,
        diagnostics_callable=diagnostics_callable,
        iterating_shape=iterating_shape,
        iterating_axes=iterating_axes,
        fitting_axes=fitting_axes,
        fitter_kwargs=fitter_kwargs,
        name="fitting-results",
        weights_specified=weights is not None,
        fit_info=fit_info,
    )

    if scheduler == "default":
        compute_kwargs = {}
    else:
        compute_kwargs = {"scheduler": scheduler}

    result_array = result.compute(**compute_kwargs)

    if fit_info:
        parameter_arrays_fitted = result_array[:-1].astype(float)
        fit_info_array = result_array[-1]
    else:
        parameter_arrays_fitted = result_array

    # Set up new parameter arrays with fitted values
    parameters = {}
    for i, name in enumerate(model.param_names):
        parameters[name] = parameter_arrays_fitted[i].reshape(iterating_shape)

    if fit_info:
        fitter.fit_info = FitInfoArrayContainer(fit_info_array)
    else:
        fitter.fit_info = original_fit_info

    # Instantiate new fitted model
    model_fitted = _copy_with_new_parameters(model, parameters)

    # Add back units if needed
    if add_back_units:
        model_fitted = model_fitted.with_units_from_data(**rename_data)

    return model_fitted
