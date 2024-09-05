# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Combine 3 images to produce properly-scaled RGB images with arbitrary scaling.

The three images must be aligned and have the same pixel scale and size.
"""

import numpy as np

from astropy.visualization.interval import ManualInterval
from astropy.visualization.stretch import LinearStretch

_OUTPUT_IMAGE_FORMATS = [float, np.float64, np.uint8]

__all__ = ["make_rgb"]


class RGBImageMapping:
    """
    Class to map red, blue, green images into either a normalized float or
    an 8-bit image, by performing optional clipping and applying
    a scaling function to each band independently.

    Parameters
    ----------
    interval : `~astropy.visualization.BaseInterval` subclass instance or array-like, optional
        The interval object to apply to the data (either a single instance or
        an array for R, G, B). Default is
        `~astropy.visualization.ManualInterval`.
    stretch : `~astropy.visualization.BaseStretch` subclass instance, optional
        The stretch object to apply to the data. Default is
        `~astropy.visualization.LinearStretch`.

    """

    def __init__(
        self, interval=ManualInterval(vmin=0, vmax=None), stretch=LinearStretch()
    ):
        try:
            len(interval)
        except TypeError:
            interval = 3 * [interval]
        if len(interval) != 3:
            raise ValueError("please provide 1 or 3 instances for interval.")

        self.intervals = interval
        self.stretch = stretch

    def make_rgb_image(self, image_r, image_g, image_b, output_dtype=np.uint8):
        """
        Convert 3 arrays, image_r, image_g, and image_b into a RGB image,
        either as an 8-bit per-channel or normalized image.

        The input images can be int or float, and in any range or bit-depth,
        but must have the same shape (NxM).

        Parameters
        ----------
        image_r : ndarray
            Image to map to red.
        image_g : ndarray
            Image to map to green.
        image_b : ndarray
            Image to map to blue.
        output_dtype : numpy scalar type, optional
            Image output format. Default is np.uint8.

        Returns
        -------
        RGBimage : ndarray
            RGB color image with the specified format as an NxMx3 numpy array.

        """
        if output_dtype not in _OUTPUT_IMAGE_FORMATS:
            raise ValueError(f"'output_dtype' must be one of {_OUTPUT_IMAGE_FORMATS}!")

        image_r = np.asarray(image_r)
        image_g = np.asarray(image_g)
        image_b = np.asarray(image_b)

        if (image_r.shape != image_g.shape) or (image_g.shape != image_b.shape):
            msg = "The image shapes must match. r: {}, g: {} b: {}"
            raise ValueError(msg.format(image_r.shape, image_g.shape, image_b.shape))

        image_rgb = self.apply_mappings(image_r, image_g, image_b)
        if np.issubdtype(output_dtype, float):
            conv_images = self._convert_images_to_float(image_rgb, output_dtype)
        elif np.issubdtype(output_dtype, np.unsignedinteger):
            conv_images = self._convert_images_to_uint(image_rgb, output_dtype)

        return np.dstack(conv_images)

    def apply_mappings(self, image_r, image_g, image_b):
        """
        Apply mapping stretch and intervals to convert images image_r, image_g,
        and image_b to a triplet of normalized images.

        Parameters
        ----------
        image_r : ndarray
            Intensity of image to be mapped to red
        image_g : ndarray
            Intensity of image to be mapped to green.
        image_b : ndarray
            Intensity of image to be mapped to blue.

        Returns
        -------
        image_rgb : ndarray
            Triplet of mapped images based on the specified (per-band)
            intervals and the stretch function

        """
        image_rgb = [image_r, image_g, image_b]
        for i, img in enumerate(image_rgb):
            # Using syntax from mpl_normalize.ImageNormalize,
            # but NOT using that class to avoid dependency on matplotlib.

            # Define vmin and vmax
            vmin, vmax = self.intervals[i].get_limits(img)

            # copy because of in-place operations after
            img = np.array(img, copy=True, dtype=float)

            # Normalize based on vmin and vmax
            np.subtract(img, vmin, out=img)
            np.true_divide(img, vmax - vmin, out=img)

            # Clip to the 0 to 1 range
            img = np.clip(img, 0.0, 1.0, out=img)

            # Stretch values
            img = self.stretch(img, out=img, clip=False)

            image_rgb[i] = img

        return np.asarray(image_rgb)

    def _convert_images_to_float(self, image_rgb, output_dtype):
        """
        Convert a triplet of normalized images to float.
        """
        return image_rgb.astype(output_dtype)

    def _convert_images_to_uint(self, image_rgb, output_dtype):
        """
        Convert a triplet of normalized images to unsigned integer images
        """
        pixmax = float(np.iinfo(output_dtype).max)

        image_rgb *= pixmax

        return image_rgb.astype(output_dtype)


def make_rgb(
    image_r,
    image_g,
    image_b,
    interval=ManualInterval(vmin=0, vmax=None),
    stretch=LinearStretch(),
    filename=None,
    output_dtype=np.uint8,
):
    """
    Base class to return a Red/Green/Blue color image from 3 images using
    a specified stretch and interval, for each band *independently*.

    The input images can be int or float, and in any range or bit-depth,
    but must have the same shape (NxM).

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
    interval : `~astropy.visualization.BaseInterval` subclass instance or array-like, optional
        The interval object to apply to the data (either a single instance or
        an array for R, G, B). Default is
        `~astropy.visualization.ManualInterval` with vmin=0.
    stretch : `~astropy.visualization.BaseStretch` subclass instance, optional
        The stretch object to apply to the data. Default is
        `~astropy.visualization.LinearStretch`.
    filename : str, optional
        Write the resulting RGB image to a file (file type determined
        from extension).
    output_dtype : numpy scalar type, optional
        Image output data type. Default is np.uint8.

    Returns
    -------
    rgb : ndarray
        RGB (either float or integer with 8-bits per channel) color image
        as an NxMx3 numpy array.

    Notes
    -----
    This procedure of clipping and then scaling is similar to the DS9
    image algorithm (see the DS9 reference guide:
    http://ds9.si.edu/doc/ref/how.html).

    """
    map_ = RGBImageMapping(interval=interval, stretch=stretch)
    rgb = map_.make_rgb_image(image_r, image_g, image_b, output_dtype=output_dtype)

    if filename:
        import matplotlib.image

        matplotlib.image.imsave(filename, rgb, origin="lower")

    return rgb
