``np.take_along_axis`` is now explicitly supported for classes based on
``ShapedLikeNDArray``, such as ``Time``, ``SkyCoord``, coordinate frames,
and representations. Previously, it only worked because numpy did not
convert its input to an array.
