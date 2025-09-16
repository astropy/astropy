# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Combine 3 images to produce a properly-scaled RGB image following
Lupton et al. (2004).

The three images must be aligned and have the same pixel scale and size.

For details, see : https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L
"""

import numpy as np

from astropy.visualization.interval import ManualInterval, ZScaleInterval
from astropy.visualization.stretch import BaseStretch
from astropy.visualization.stretch import _prepare as _stretch_prepare

from .basic_rgb import RGBImageMapping

__all__ = [
    "AsinhMapping",
    "AsinhZScaleMapping",
    "LinearMapping",
    "LuptonAsinhStretch",
    "LuptonAsinhZscaleStretch",
    "Mapping",
    "make_lupton_rgb",
]


def compute_intensity(image_r, image_g=None, image_b=None):
    """
    Return a naive total intensity from the red, blue, and green intensities.

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


class LuptonAsinhStretch(BaseStretch):
    r"""
    A modified asinh stretch, with some changes to the constants
    relative to `~astropy.visualization.AsinhStretch`.

    The stretch is given by:

    .. math::
        & y = {\rm asinh}\left(\frac{Q * x}{stretch}\right) *
        \frac{frac}{{\rm asinh}(frac * Q)} \\
        & frac = 0.1

    Parameters
    ----------
    stretch : float, optional
        Linear stretch of the image. ``stretch`` must be greater than 0.
        Default is 5.

    Q : float, optional
        The asinh softening parameter. ``Q`` must be greater than 0.
        Default is 8.

    Notes
    -----
    Based on the asinh stretch presented in Lupton et al. 2004
    (https://ui.adsabs.harvard.edu/abs/2004PASP..116..133L).

    Examples
    --------
    .. plot::
        :show-source-link:

        import numpy as np
        from astropy.visualization import LuptonAsinhStretch
        from matplotlib import pyplot as plt

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8, 8),
                               layout='constrained')
        ax = ax.ravel()

        x = np.linspace(0, 1, 100)
        stretches = (0.05, 0.1, 0.2, 0.5, 1, 5, 10)

        Qs = (2, 5, 8, 10)
        for i, Q in enumerate(Qs):
            for st in stretches:
                stretch = LuptonAsinhStretch(stretch=st, Q=Q)
                label = f'{st=}'
                ax[i].plot(x, stretch(x, clip=True), label=label)

            ax[i].axis('equal')
            ax[i].plot(x, x, ls='dotted', color='k', alpha=0.3)
            ax[i].set_xlim(0, 1)
            ax[i].set_ylim(0, 1)
            ax[i].set_xlabel('Input Value')
            ax[i].set_ylabel('Output Value')
            ax[i].set_title(f'{stretch.__class__.__name__}, {Q=}')
            ax[i].legend(loc='lower right', fontsize=8)
    """

    def __init__(self, stretch=5, Q=8):
        super().__init__()

        if stretch < 0:
            raise ValueError(f"Stretch must be non-negative! {stretch=}")
        if Q < 0:
            raise ValueError(f"Q must be non-negative! {Q=}")

        # 32bit floating point machine epsilon; sys.float_info.epsilon is 64bit
        epsilon = 1.0 / 2**23
        if abs(Q) < epsilon:
            Q = 0.1
        else:
            Qmax = 1e10
            if Q > Qmax:
                Q = Qmax

        self.stretch = stretch
        self.Q = Q

        frac = 0.1
        self._slope = frac / np.arcsinh(frac * Q)
        self._soften = Q / float(stretch)

    def __call__(self, values, clip=False, out=None):
        values = _stretch_prepare(values, clip=clip, out=out)
        np.multiply(values, self._soften, out=values)
        np.arcsinh(values, out=values)
        np.multiply(values, self._slope, out=values)
        return values


class LuptonAsinhZscaleStretch(LuptonAsinhStretch):
    r"""
    A modified asinh stretch, where the linear stretch is calculated using
    zscale.

    The stretch is given by:

    .. math::
        & y = {\rm asinh}\left(\frac{Q * x}{stretch}\right) *
        \frac{frac}{{\rm asinh}(frac * Q)} \\
        & frac = 0.1 \\
        & stretch = z2 - z1

    Parameters
    ----------
    image1 : ndarray or array-like
        The image to analyse, or a list of 3 images to be converted to
        an intensity image.

    Q : float, optional
        The asinh softening parameter. ``Q`` must be greater than 0.
        Default is 8.

    pedestal : or array-like, optional
        The value, or array of 3 values, to subtract from the images(s)
        before determining the zscaling. Default is None (nothing subtracted).

    """

    def __init__(self, image, Q=8, pedestal=None):
        # copy because of in-place operations after
        image = np.array(image, copy=True, dtype=float)

        _raiseerr = False
        if len(image.shape) == 2:
            image = [image]
        elif len(image.shape) == 3:
            if image.shape[0] != 3:
                _raiseerr = True
        else:
            _raiseerr = True

        if _raiseerr:
            raise ValueError(
                "Input 'image' must be a single image "
                "or a stack/3xMxN array of 3 images! "
                f"{image.shape=}"
            )

        image = list(image)  # needs to be mutable

        if pedestal is not None:
            try:
                len(pedestal)
            except TypeError:
                pedestal = 3 * [pedestal]

            if len(pedestal) != 3:
                raise ValueError(
                    "pedestal must be 1 or 3 values, matching the image input."
                )
            for i, im in enumerate(image):
                if pedestal[i] != 0.0:
                    image[i] = im - pedestal[i]  # n.b. a copy

        image = compute_intensity(*image)
        zscale_limits = ZScaleInterval().get_limits(image)

        _stretch = zscale_limits[1] - zscale_limits[0]

        self._image = image

        super().__init__(stretch=_stretch, Q=Q)


class RGBImageMappingLupton(RGBImageMapping):
    """
    Class to map red, blue, green images into either a normalized float or
    an 8-bit image, by performing optional clipping and applying
    a scaling function to each band in non-independent manner that depends
    on the other bands, following the scaling scheme presented in
    Lupton et al. 2004.

    Parameters
    ----------
    interval : `~astropy.visualization.BaseInterval` subclass instance or array-like, optional
        The interval object to apply to the data (either a single instance or
        an array for R, G, B). Default is
        `~astropy.visualization.ManualInterval`.
    stretch : `~astropy.visualization.BaseStretch` subclass instance
        The stretch object to apply to the data. The default is
        `~astropy.visualization.AsinhLuptonStretch`.

    """

    def __init__(
        self,
        interval=ManualInterval(vmin=0, vmax=None),
        stretch=LuptonAsinhStretch(stretch=5, Q=8),
    ):
        super().__init__(interval=interval, stretch=stretch)
        self._pixmax = 1.0

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

    def apply_mappings(self, image_r, image_g, image_b):
        """
        Apply mapping stretch and intervals to convert images image_r, image_g,
        and image_b to a triplet of normalized images, following the scaling
        scheme presented in Lupton et al. 2004.

        Compared to astropy's ImageNormalize which first normalizes images
        by cropping and linearly mapping onto [0.,1.] and then applies
        a specified stretch algorithm, the Lupton et al. algorithm applies
        stretching to an multi-color intensity and then computes per-band
        scaled images with bound cropping.

        This is modified here by allowing for different minimum values
        for each of the input r, g, b images, and then computing
        the intensity on the subtracted images.

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

        Notes
        -----
        The Lupton et al 2004 algorithm is computed with the following steps:

        1. Shift each band with the minimum values
        2. Compute the intensity I and stretched intensity f(I)
        3. Compute the ratio of the stretched intensity to intensity f(I)/I,
        and clip to a lower bound of 0
        4. Compute the scaled band images by multiplying with the ratio f(I)/I
        5. Clip each band to a lower bound of 0
        6. Scale down pixels where max(R,G,B)>1 by the value max(R,G,B)

        """
        image_r = np.array(image_r, copy=True)
        image_g = np.array(image_g, copy=True)
        image_b = np.array(image_b, copy=True)

        # Subtract per-band minima
        image_rgb = [image_r, image_g, image_b]
        for i, img in enumerate(image_rgb):
            vmin, _ = self.intervals[i].get_limits(img)
            image_rgb[i] = np.subtract(img, vmin)

        image_rgb = np.asarray(image_rgb)

        # Determine the intensity and streteched intensity
        Int = self.intensity(*image_rgb)
        fI = self.stretch(Int, clip=False)

        # Get normalized fI, and clip to lower bound of 0:
        fInorm = np.where(Int <= 0, 0, np.true_divide(fI, Int))

        # Compute X = x * f(I) / I for each filter x=(r,g,b)
        np.multiply(image_rgb, fInorm, out=image_rgb)

        # Clip individual bands to minimum of 0, as
        # individual bands can be < 0 even if fI/I isn't.
        image_rgb = np.clip(image_rgb, 0.0, None)

        # Determine the max of all 3 bands at each position
        maxRGB = np.max(image_rgb, axis=0)

        with np.errstate(invalid="ignore", divide="ignore"):
            image_rgb = np.where(
                maxRGB > self._pixmax,
                np.true_divide(image_rgb * self._pixmax, maxRGB),
                image_rgb,
            )

        return np.asarray(image_rgb)


def make_lupton_rgb(
    image_r,
    image_g,
    image_b,
    interval=None,
    stretch_object=None,
    minimum=None,
    stretch=5,
    Q=8,
    filename=None,
    output_dtype=np.uint8,
):
    r"""
    Return a Red/Green/Blue color image from 3 images using interconnected
    band scaling, and an arbitrary stretch function (by default, an asinh stretch).
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
    interval : `~astropy.visualization.BaseInterval` subclass instance or array-like, optional
        The interval object to apply to the data (either a single instance or
        an array for R, G, B). Default is
        `~astropy.visualization.ManualInterval` with vmin=0.
    stretch_object : `~astropy.visualization.BaseStretch` subclass instance, optional
        The stretch object to apply to the data. If set, the input values of
        ``minimum``, ``stretch``, and ``Q`` will be ignored.
        For the Lupton scheme, this would be an instance of
        `~astropy.visualization.LuptonAsinhStretch`, but alternatively
        `~astropy.visualization.LuptonAsinhZscaleStretch` or some other
        stretch can be used.
    minimum : float or array-like, optional
        Deprecated. Intensity that should be mapped to black (a scalar or
        array of R, G, B). If `None`, each image's minimum value is used.
        Default is None.
    stretch : float, optional
        The linear stretch of the image. Default is 5
    Q : float, optional
        The asinh softening parameter. Default is 8.
    filename : str, optional
        Write the resulting RGB image to a file (file type determined
        from extension).
    output_dtype : numpy scalar type, optional
        Image output data type. Default is np.uint8.

    Returns
    -------
    rgb : ndarray
        RGB color image as an NxNx3 numpy array, with the specified
        data type format

    """
    if stretch_object is None:
        stretch_object = LuptonAsinhStretch(stretch=stretch, Q=Q)

    if interval is None:
        # Only use minimum if interval is not specified:
        if minimum is not None:
            # Backwards compatibility:
            try:
                len(minimum)
            except TypeError:
                minimum = 3 * [minimum]
            if len(minimum) != 3:
                raise ValueError("please provide 1 or 3 values for minimum.")
            interval = []
            for i in range(3):
                interval.append(ManualInterval(vmin=minimum[i], vmax=None))
        else:
            # Default option:
            interval = ManualInterval(vmin=0, vmax=None)

    lup_map = RGBImageMappingLupton(
        interval=interval,
        stretch=stretch_object,
    )
    rgb = lup_map.make_rgb_image(image_r, image_g, image_b, output_dtype=output_dtype)

    if filename:
        import matplotlib.image

        matplotlib.image.imsave(filename, rgb, origin="lower")

    return rgb
