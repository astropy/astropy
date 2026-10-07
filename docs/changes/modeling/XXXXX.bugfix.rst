Fixed ``CompoundModel.evaluate`` passing the wrong arguments to the component
models when the left model has array-valued parameters, is a model set, or has
a non-scalar parameter such as ``AffineTransformation2D.matrix``. This made
fitting such compound models fail, e.g. with ``parallel_fit_dask``.
