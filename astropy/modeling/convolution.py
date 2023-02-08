# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Convolution Model."""
# pylint: disable=line-too-long, too-many-lines, too-many-arguments, invalid-name
import numpy as np

from .core import CompoundModel


class Convolution(CompoundModel):
    """
    Wrapper class for a convolution model.

    Parameters
    ----------
    operator: tuple
        The SPECIAL_OPERATORS entry for the convolution being used.
    model : Model
        The model for the convolution.
    kernel: Model
        The kernel model for the convolution.
    bounding_box : tuple
        A bounding box to define the limits of the integration
        approximation for the convolution.
    resolution : float
        The resolution for the approximation of the convolution.
    cache : bool, optional
        Allow convolution computation to be cached for reuse. This is
        enabled by default.

    Notes
    -----
    This is wrapper is necessary to handle the limitations of the
    pseudospectral convolution binary operator implemented in
    astropy.convolution under `~astropy.convolution.convolve_fft`. In this
    `~astropy.convolution.convolve_fft` it is assumed that the inputs ``array``
    and ``kernel`` span a sufficient portion of the support of the functions of
    the convolution. Consequently, the ``Compound`` created by the
    `~astropy.convolution.convolve_models` function makes the assumption that
    one should pass an input array that sufficiently spans this space. This means
    that slightly different input arrays to this model will result in different
    outputs, even on points of intersection between these arrays.

    This issue is solved by requiring a ``bounding_box`` together with a
    resolution so that one can pre-calculate the entire domain and then
    (by default) cache the convolution values. The function then just
    interpolates the results from this cache.
    """

    def __init__(self, operator, model, kernel, bounding_box, resolution, cache=True):
        super().__init__(operator, model, kernel)

        self.bounding_box = bounding_box
        self._resolution = resolution

        self._cache_convolution = cache
        self._kwargs = None
        self._convolution = None

    def clear_cache(self):
        """
        Clears the cached convolution.
        """
        self._kwargs = None
        self._convolution = None

    def _get_convolution(self, **kwargs):
        if (self._convolution is None) or (self._kwargs != kwargs):
            domain = self.bounding_box.domain(self._resolution)
            mesh = np.meshgrid(*domain)
            data = super().__call__(*mesh, **kwargs)

            from scipy.interpolate import RegularGridInterpolator

            convolution = RegularGridInterpolator(domain, data)

            if self._cache_convolution:
                self._kwargs = kwargs
                self._convolution = convolution

        else:
            convolution = self._convolution

        return convolution

    @staticmethod
    def _convolution_inputs(*args):
        not_scalar = np.where([not np.isscalar(arg) for arg in args])[0]

        if len(not_scalar) == 0:
            return np.array(args), (1,)
        else:
            output_shape = args[not_scalar[0]].shape
            if not all(args[index].shape == output_shape for index in not_scalar):
                raise ValueError("Values have differing shapes")

            inputs = []
            for arg in args:
                if np.isscalar(arg):
                    inputs.append(np.full(output_shape, arg))
                else:
                    inputs.append(arg)

            return np.reshape(inputs, (len(inputs), -1)).T, output_shape

    @staticmethod
    def _convolution_outputs(outputs, output_shape):
        return outputs.reshape(output_shape)

    def __call__(self, *args, **kw):
        inputs, output_shape = self._convolution_inputs(*args)
        convolution = self._get_convolution(**kw)
        outputs = convolution(inputs)

        return self._convolution_outputs(outputs, output_shape)
