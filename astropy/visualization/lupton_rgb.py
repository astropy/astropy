# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Combine 3 images to produce a properly-scaled RGB image following Lupton et al. (2004).

The three images must be aligned and have the same pixel scale and size.

For details, see : https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L
"""

import numpy as np

from . import ZScaleInterval

__all__ = ["make_lupton_rgb"]


def compute_intensity(image_r, image_g=None, image_b=None):
    """
    Return a naive total intensity from the red, blue, and green intensities.

    Parameters
    ----------
    image_r : ndarray
        Intensity of image to be mapped to red; or total intensity if ``image_g``
        and ``image_b`` are None.
    image_g : ndarray, optional
        Intensity of image to be mapped to green.
    image_b : ndarray, optional
        Intensity of image to be mapped to blue.

    Returns
    -------
    intensity : ndarray
        Total intensity from the red, blue and green intensities, or ``image_r``
        if green and blue images are not provided.
    """
    if image_g is None or image_b is None:
        if not (image_g is None and image_b is None):
            raise ValueError(
                "please specify either a single image or red, green, and blue images."
            )
        return image_r

    intensity = (image_r + image_g + image_b) / 3.0

    # Repack into whatever type was passed to us
    return np.asarray(intensity, dtype=image_r.dtype)


class Mapping:
    """
    Baseclass to map red, blue, green intensities into uint8 values.

    Parameters
    ----------
    minimum : float or sequence(3)
        Intensity that should be mapped to black (a scalar or array for R, G, B).
    image : ndarray, optional
        An image used to calculate some parameters of some mappings.
    """

    def __init__(self, minimum=None, image=None):
        self._uint8Max = float(np.iinfo(np.uint8).max)

        try:
            len(minimum)
        except TypeError:
            minimum = 3 * [minimum]
        if len(minimum) != 3:
            raise ValueError("please provide 1 or 3 values for minimum.")

        self.minimum = minimum
        self._image = np.asarray(image)

    def make_rgb_image(self, image_r, image_g, image_b):
        """
        Convert 3 arrays, image_r, image_g, and image_b into an 8-bit RGB image.

        Parameters
        ----------
        image_r : ndarray
            Image to map to red.
        image_g : ndarray
            Image to map to green.
        image_b : ndarray
            Image to map to blue.

        Returns
        -------
        RGBimage : ndarray
            RGB (integer, 8-bits per channel) color image as an NxNx3 numpy array.
        """
        image_r = np.asarray(image_r)
        image_g = np.asarray(image_g)
        image_b = np.asarray(image_b)

        if (image_r.shape != image_g.shape) or (image_g.shape != image_b.shape):
            msg = "The image shapes must match. r: {}, g: {} b: {}"
            raise ValueError(msg.format(image_r.shape, image_g.shape, image_b.shape))

        return np.dstack(
            self._convert_images_to_uint8(image_r, image_g, image_b)
        ).astype(np.uint8)

    def intensity(self, image_r, image_g, image_b):
        """
        Return the total intensity from the red, blue, and green intensities.
        This is a naive computation, and may be overridden by subclasses.

        Parameters
        ----------
        image_r : ndarray
            Intensity of image to be mapped to red; or total intensity if
            ``image_g`` and ``image_b`` are None.
        image_g : ndarray, optional
            Intensity of image to be mapped to green.
        image_b : ndarray, optional
            Intensity of image to be mapped to blue.

        Returns
        -------
        intensity : ndarray
            Total intensity from the red, blue and green intensities, or
            ``image_r`` if green and blue images are not provided.
        """
        return compute_intensity(image_r, image_g, image_b)

    def map_intensity_to_uint8(self, I):
        """
        Return an array which, when multiplied by an image, returns that image
        mapped to the range of a uint8, [0, 255] (but not converted to uint8).

        The intensity is assumed to have had minimum subtracted (as that can be
        done per-band).

        Parameters
        ----------
        I : ndarray
            Intensity to be mapped.

        Returns
        -------
        mapped_I : ndarray
            ``I`` mapped to uint8
        """
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.clip(I, 0, self._uint8Max)

    def _convert_images_to_uint8(self, image_r, image_g, image_b):
        """
        Use the mapping to convert images image_r, image_g, and image_b to a triplet of uint8 images.
        """
        image_r = image_r - self.minimum[0]  # n.b. makes copy
        image_g = image_g - self.minimum[1]
        image_b = image_b - self.minimum[2]

        fac = self.map_intensity_to_uint8(self.intensity(image_r, image_g, image_b))

        image_rgb = [image_r, image_g, image_b]
        for c in image_rgb:
            c *= fac
            with np.errstate(invalid="ignore"):
                c[c < 0] = 0  # individual bands can still be < 0, even if fac isn't

        pixmax = self._uint8Max
        # copies -- could work row by row to minimise memory usage
        r0, g0, b0 = image_rgb

        # n.b. np.where can't and doesn't short-circuit
        with np.errstate(invalid="ignore", divide="ignore"):
            for i, c in enumerate(image_rgb):
                c = np.where(
                    r0 > g0,
                    np.where(
                        r0 > b0,
                        np.where(r0 >= pixmax, c * pixmax / r0, c),
                        np.where(b0 >= pixmax, c * pixmax / b0, c),
                    ),
                    np.where(
                        g0 > b0,
                        np.where(g0 >= pixmax, c * pixmax / g0, c),
                        np.where(b0 >= pixmax, c * pixmax / b0, c),
                    ),
                ).astype(np.uint8)
                c[c > pixmax] = pixmax

                image_rgb[i] = c

        return image_rgb


class LinearMapping(Mapping):
    """
    A linear map map of red, blue, green intensities into uint8 values.

    A linear stretch from [minimum, maximum].
    If one or both are omitted use image min and/or max to set them.

    Parameters
    ----------
    minimum : float
        Intensity that should be mapped to black (a scalar or array for R, G, B).
    maximum : float
        Intensity that should be mapped to white (a scalar).
    """

    def __init__(self, minimum=None, maximum=None, image=None):
        if minimum is None or maximum is None:
            if image is None:
                raise ValueError(
                    "you must provide an image if you don't "
                    "set both minimum and maximum"
                )
            if minimum is None:
                minimum = image.min()
            if maximum is None:
                maximum = image.max()

        Mapping.__init__(self, minimum=minimum, image=image)
        self.maximum = maximum

        if maximum is None:
            self._range = None
        else:
            if maximum == minimum:
                raise ValueError("minimum and maximum values must not be equal")
            self._range = float(maximum - minimum)

    def map_intensity_to_uint8(self, I):
        # n.b. np.where can't and doesn't short-circuit
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(
                I <= 0,
                0,
                np.where(
                    I >= self._range, self._uint8Max / I, self._uint8Max / self._range
                ),
            )


class AsinhMapping(Mapping):
    """
    A mapping for an asinh stretch (preserving colours independent of brightness).

    x = asinh(Q (I - minimum)/stretch)/Q

    This reduces to a linear stretch if Q == 0

    See https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L

    Parameters
    ----------
    minimum : float
        Intensity that should be mapped to black (a scalar or array for R, G, B).
    stretch : float
        The linear stretch of the image.
    Q : float
        The asinh softening parameter.
    """

    def __init__(self, minimum, stretch, Q=8):
        Mapping.__init__(self, minimum)

        # 32bit floating point machine epsilon; sys.float_info.epsilon is 64bit
        epsilon = 1.0 / 2**23
        if abs(Q) < epsilon:
            Q = 0.1
        else:
            Qmax = 1e10
            if Q > Qmax:
                Q = Qmax

        frac = 0.1  # gradient estimated using frac*stretch is _slope
        self._slope = frac * self._uint8Max / np.arcsinh(frac * Q)

        self._soften = Q / float(stretch)

    def map_intensity_to_uint8(self, I):
        # n.b. np.where can't and doesn't short-circuit
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(I <= 0, 0, np.arcsinh(I * self._soften) * self._slope / I)


class AsinhZScaleMapping(AsinhMapping):
    """
    A mapping for an asinh stretch, estimating the linear stretch by zscale.

    x = asinh(Q (I - z1)/(z2 - z1))/Q

    Parameters
    ----------
    image1 : ndarray or a list of arrays
        The image to analyse, or a list of 3 images to be converted to
        an intensity image.
    image2 : ndarray, optional
        the second image to analyse (must be specified with image3).
    image3 : ndarray, optional
        the third image to analyse (must be specified with image2).
    Q : float, optional
        The asinh softening parameter. Default is 8.
    pedestal : float or sequence(3), optional
        The value, or array of 3 values, to subtract from the images; or None.

    Notes
    -----
    pedestal, if not None, is removed from the images when calculating the
    zscale stretch, and added back into Mapping.minimum[]
    """

    def __init__(self, image1, image2=None, image3=None, Q=8, pedestal=None):
        if image2 is None or image3 is None:
            if not (image2 is None and image3 is None):
                raise ValueError(
                    "please specify either a single image or three images."
                )
            image = [image1]
        else:
            image = [image1, image2, image3]

        if pedestal is not None:
            try:
                len(pedestal)
            except TypeError:
                pedestal = 3 * [pedestal]

            if len(pedestal) != 3:
                raise ValueError("please provide 1 or 3 pedestals.")

            image = list(image)  # needs to be mutable
            for i, im in enumerate(image):
                if pedestal[i] != 0.0:
                    image[i] = im - pedestal[i]  # n.b. a copy
        else:
            pedestal = len(image) * [0.0]

        image = compute_intensity(*image)

        zscale_limits = ZScaleInterval().get_limits(image)
        zscale = LinearMapping(*zscale_limits, image=image)
        # zscale.minimum is always a triple
        stretch = zscale.maximum - zscale.minimum[0]
        minimum = zscale.minimum

        for i, level in enumerate(pedestal):
            minimum[i] += level

        AsinhMapping.__init__(self, minimum, stretch, Q)
        self._image = image


def make_lupton_rgb(
    image_r, image_g, image_b, minimum=0, stretch=5, Q=8, filename=None
):
    """
    Return a Red/Green/Blue color image from up to 3 images using an asinh stretch.
    The input images can be int or float, and in any range or bit-depth.

    For a more detailed look at the use of this method, see the document
    :ref:`astropy:astropy-visualization-rgb`.

    Parameters
    ----------
    image_r : ndarray
        Image to map to red.
    image_g : ndarray
        Image to map to green.
    image_b : ndarray
        Image to map to blue.
    minimum : float
        Intensity that should be mapped to black (a scalar or array for R, G, B).
    stretch : float
        The linear stretch of the image.
    Q : float
        The asinh softening parameter.
    filename : str
        Write the resulting RGB image to a file (file type determined
        from extension).

    Returns
    -------
    rgb : ndarray
        RGB (integer, 8-bits per channel) color image as an NxNx3 numpy array.
    """
    asinhMap = AsinhMapping(minimum, stretch, Q)
    rgb = asinhMap.make_rgb_image(image_r, image_g, image_b)

    if filename:
        import matplotlib.image

        matplotlib.image.imsave(filename, rgb, origin="lower")

    return rgb
