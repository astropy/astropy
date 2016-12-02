# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Combine 3 images to produce a properly-scaled RGB image following Lupton et al. (2004).

For details, see : http://adsabs.harvard.edu/abs/2004PASP..116..133L

The three images must be aligned and have the same pixel scale and size.

Example usage:
    imageR = np.random.random((100,100))
    imageG = np.random.random((100,100))
    imageB = np.random.random((100,100))
    image = lupton_rgb.makeRGB(imageR, imageG, imageB, fileName='randoms.png')
"""
from __future__ import absolute_import, division

import numpy as np
from . import ZScaleInterval


__all__ = ['makeRGB', 'Mapping', 'LinearMapping', 'AsinhMapping',
           'AsinhZScaleMapping']


try:
    import scipy.misc
    scipy.misc.imresize  # checking if it exists
    HAVE_SCIPY_MISC = True
except (ImportError, AttributeError):
    HAVE_SCIPY_MISC = False

# NOTE: these methods would have come from LSST C++ code. They won't be available
# in astropy until they are converted somehow.
# from lsst.afw.display.displayLib import replaceSaturatedPixels, getZScale


def compute_intensity(imageR, imageG=None, imageB=None):
    """
    Return a naive total intensity from the red, blue, and green intensities.

    Parameters
    ----------
    imageR : `~numpy.ndarray`
        Intensity of image to be mapped to red; or total intensity if ``imageG``
        and ``imageB`` are None.
    imageG : `~numpy.ndarray`, optional
        Intensity of image to be mapped to green.
    imageB : `~numpy.ndarray`, optional
        Intensity of image to be mapped to blue.

    Returns
    -------
    intensity : `~numpy.ndarray`
        Total intensity from the red, blue and green intensities, or ``imageR``
        if green and blue images are not provided.

    Notes
    -----
    Inputs may be MaskedImages, Images, or numpy arrays and the return is
    of the same type.
    """
    if imageG is None or imageB is None:
        if not (imageG is None and imageB is None):
            raise ValueError("please specify either a single image "
                             "or red, green, and blue images.")
        return imageR

    intensity = (imageR + imageG + imageB)/3.0

    # Repack into whatever type was passed to us
    return np.array(intensity, dtype=imageR.dtype)


class Mapping(object):
    """Baseclass to map red, blue, green intensities into uint8 values."""

    def __init__(self, minimum=None, image=None):
        """
        Create a mapping.

        Parameters
        ----------
        minimum : float or sequence(3)
            Intensity that should be mapped to black (a scalar or array for R, G, B).
        image : `~numpy.ndarray`, optional
            The image to be used to calculate the mapping.
            If provided, it is also used as the default for makeRgbImage().
        """
        self._uint8Max = float(np.iinfo(np.uint8).max)

        try:
            len(minimum)
        except TypeError:
            minimum = 3*[minimum]
        if len(minimum) != 3:
            raise ValueError("please provide 1 or 3 values for minimum.")

        self.minimum = minimum
        self._image = image

    def makeRgbImage(self, imageR=None, imageG=None, imageB=None,
                     xSize=None, ySize=None, rescaleFactor=None):
        """
        Convert 3 arrays, imageR, imageG, and imageB into a numpy RGB image.

        Parameters
        ----------
        imageR : `~numpy.ndarray`, optional
            Image to map to red (if None, use the image passed to the
            constructor).
        imageG : `~numpy.ndarray`, optional
            Image to map to green (if None, use imageR).
        imageB : `~numpy.ndarray`, optional
            Image to map to blue (if None, use imageR).
        xSize : int, optional
            Desired width of RGB image (or None).  If ySize is None, preserve
            aspect ratio.
        ySize : int, optional
            Desired height of RGB image (or None).
        rescaleFactor : float, optional
            Make size of output image rescaleFactor*size of the input image.
            Cannot be specified if xSize or ySize are given.

        Returns
        -------
        RGBimage : `~numpy.ndarray`
            An image formed by stacking the input images.
        """
        if imageR is None:
            if self._image is None:
                raise RuntimeError("you must provide an image or pass one "
                                   "to the constructor.")
            imageR = self._image

        if imageG is None:
            imageG = imageR
        if imageB is None:
            imageB = imageR

        if xSize is not None or ySize is not None:
            if rescaleFactor is not None:
                raise ValueError("you may not specify a size and rescaleFactor.")
            h, w = imageR.shape
            if ySize is None:
                ySize = int(xSize*h/float(w) + 0.5)
            elif xSize is None:
                xSize = int(ySize*w/float(h) + 0.5)

            # need to cast to int when passing tuple to imresize.
            size = (int(ySize), int(xSize))  # n.b. y, x order for scipy
        elif rescaleFactor is not None:
            size = float(rescaleFactor)  # a float is intepreted as a percentage
        else:
            size = None

        if size is not None:
            if not HAVE_SCIPY_MISC:
                raise RuntimeError("unable to rescale as scipy.misc.imresize "
                                   "is unavailable.")

            imageR = scipy.misc.imresize(imageR, size, interp='bilinear',
                                         mode='F')
            imageG = scipy.misc.imresize(imageG, size, interp='bilinear',
                                         mode='F')
            imageB = scipy.misc.imresize(imageB, size, interp='bilinear',
                                         mode='F')

        return np.dstack(self._convertImagesToUint8(imageR, imageG, imageB)).astype(np.uint8)

    def intensity(self, imageR, imageG, imageB):
        """
        Return the total intensity from the red, blue, and green intensities.
        This is a naive computation, and may be overridden by subclasses.

        Parameters
        ----------
        imageR : `~numpy.ndarray`
            Intensity of image to be mapped to red; or total intensity if
            ``imageG`` and ``imageB`` are None.
        imageG : `~numpy.ndarray`, optional
            Intensity of image to be mapped to green.
        imageB : `~numpy.ndarray`, optional
            Intensity of image to be mapped to blue.

        Returns
        -------
        intensity : `~numpy.ndarray`
            Total intensity from the red, blue and green intensities, or
            ``imageR`` if green and blue images are not provided.
        """
        return compute_intensity(imageR, imageG, imageB)

    def mapIntensityToUint8(self, I):
        """
        Return an array which, when multiplied by an image, returns that image
        mapped to the range of a uint8, [0, 255] (but not converted to uint8).

        The intensity is assumed to have had minimum subtracted (as that can be
        done per-band).

        Parameters
        ----------
        I : `~numpy.ndarray`
            Intensity to be mapped.

        Returns
        -------
        mapped_I : `~numpy.ndarray`
            ``I`` mapped to uint8
        """
        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            return np.where(I <= 0, 0, np.where(I < self._uint8Max, I, self._uint8Max))

    def _convertImagesToUint8(self, imageR, imageG, imageB):
        """Use the mapping to convert images imageR, imageG, and imageB to a triplet of uint8 images"""
        imageR = imageR - self.minimum[0]  # n.b. makes copy
        imageG = imageG - self.minimum[1]
        imageB = imageB - self.minimum[2]

        fac = self.mapIntensityToUint8(self.intensity(imageR, imageG, imageB))

        imageRGB = [imageR, imageG, imageB]
        for c in imageRGB:
            c *= fac
            c[c < 0] = 0                # individual bands can still be < 0, even if fac isn't

        pixmax = self._uint8Max
        r0, g0, b0 = imageRGB           # copies -- could work row by row to minimise memory usage

        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            for i, c in enumerate(imageRGB):
                c = np.where(r0 > g0,
                             np.where(r0 > b0,
                                      np.where(r0 >= pixmax, c*pixmax/r0, c),
                                      np.where(b0 >= pixmax, c*pixmax/b0, c)),
                             np.where(g0 > b0,
                                      np.where(g0 >= pixmax, c*pixmax/g0, c),
                                      np.where(b0 >= pixmax, c*pixmax/b0, c))).astype(np.uint8)
                c[c > pixmax] = pixmax

                imageRGB[i] = c

        return imageRGB


class LinearMapping(Mapping):
    """A linear map map of red, blue, green intensities into uint8 values"""

    def __init__(self, minimum=None, maximum=None, image=None):
        """
        A linear stretch from [minimum, maximum].
        If one or both are omitted use image min and/or max to set them.

        Parameters
        ----------
        minimum : float
            Intensity that should be mapped to black (a scalar or array for R, G, B).
        maximum : float
            Intensity that should be mapped to white (a scalar).
        """

        if minimum is None or maximum is None:
            if image is None:
                raise ValueError("you must provide an image if you don't "
                                 "set both minimum and maximum")
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

    def mapIntensityToUint8(self, I):
        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            return np.where(I <= 0, 0,
                            np.where(I >= self._range, self._uint8Max/I, self._uint8Max/self._range))


class AsinhMapping(Mapping):
    """
    A mapping for an asinh stretch (preserving colours independent of brightness)

    x = asinh(Q (I - minimum)/range)/Q

    This reduces to a linear stretch if Q == 0

    See http://adsabs.harvard.edu/abs/2004PASP..116..133L
    """

    def __init__(self, minimum, dataRange, Q=8):
        """
        asinh stretch from minimum to minimum + dataRange, scaled by Q, via:
            x = asinh(Q (I - minimum)/dataRange)/Q

        Parameters
        ----------

        minimum : float
            Intensity that should be mapped to black (a scalar or array for R, G, B).
        dataRange : float
            minimum+dataRange defines the white level of the image.
        Q : float
            The asinh softening parameter.
        """
        Mapping.__init__(self, minimum)

        epsilon = 1.0/2**23            # 32bit floating point machine epsilon; sys.float_info.epsilon is 64bit
        if abs(Q) < epsilon:
            Q = 0.1
        else:
            Qmax = 1e10
            if Q > Qmax:
                Q = Qmax

        if False:
            self._slope = self._uint8Max/Q  # gradient at origin is self._slope
        else:
            frac = 0.1                  # gradient estimated using frac*range is _slope
            self._slope = frac*self._uint8Max/np.arcsinh(frac*Q)

        self._soften = Q/float(dataRange)

    def mapIntensityToUint8(self, I):
        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            return np.where(I <= 0, 0, np.arcsinh(I*self._soften)*self._slope/I)


class AsinhZScaleMapping(AsinhMapping):
    """
    A mapping for an asinh stretch, estimating the linear stretch by zscale.

    x = asinh(Q (I - z1)/(z2 - z1))/Q

    See AsinhMapping

    """

    def __init__(self, image1, image2=None, image3=None, Q=8, pedestal=None):
        """
        Create an asinh mapping from an image, setting the linear part of the
        stretch using zscale.

        Parameters
        ----------
        image1 : `~numpy.ndarray` or a list of arrays
            The image to analyse, or a list of 3 images to be converted to
            an intensity image.
        image2 : `~numpy.ndarray`, optional
            the second image to analyse (must be specified with image3).
        image3 : `~numpy.ndarray`, optional
            the third image to analyse (must be specified with image2).
        Q : float, optional
            The asinh softening parameter. Default is 8.
        pedestal : float or sequence(3), optional
            The value, or array of 3 values, to subtract from the images; or None.

        Notes
        -----
        N.b. pedestal, if not None, is removed from the images when
        calculating the zscale stretch, and added back into
        Mapping.minimum[]
        """

        if image2 is None or image3 is None:
            if not (image2 is None and image3 is None):
                raise ValueError("please specify either a single image "
                                 "or three images.")
            image = [image1]
        else:
            image = [image1, image2, image3]

        if pedestal is not None:
            try:
                len(pedestal)
            except TypeError:
                pedestal = 3*[pedestal]

            if len(pedestal) != 3:
                raise ValueError("please provide 1 or 3 pedestals.")


            image = list(image)        # needs to be mutable
            for i, im in enumerate(image):
                if pedestal[i] != 0.0:
                    image[i] = im - pedestal[i]  # n.b. a copy
        else:
            pedestal = len(image)*[0.0]

        image = compute_intensity(*image)

        zscale_limits = ZScaleInterval().get_limits(image)
        zscale = LinearMapping(*zscale_limits, image=image)
        dataRange = zscale.maximum - zscale.minimum[0]  # zscale.minimum is always a triple
        minimum = zscale.minimum

        for i, level in enumerate(pedestal):
            minimum[i] += level

        AsinhMapping.__init__(self, minimum, dataRange, Q)
        self._image = image


def makeRGB(imageR, imageG=None, imageB=None, minimum=0, dataRange=5, Q=8,
            saturatedBorderWidth=0, saturatedPixelValue=None,
            xSize=None, ySize=None, rescaleFactor=None,
            fileName=None):
    """
    Make an RGB color image from 3 images using an asinh stretch.

    Parameters
    ----------

    imageR : `~numpy.ndarray`
        Image to map to red (if None, use the image passed to the constructor).
    imageG : `~numpy.ndarray`
        Image to map to green (if None, use imageR).
    imageB : `~numpy.ndarray`
        Image to map to blue (if None, use imageR).
    minimum : float
        Intensity that should be mapped to black (a scalar or array for R, G, B).
    dataRange : float
        minimum+dataRange defines the white level of the image.
    Q : float
        The asinh softening parameter.
    saturatedBorderWidth : int
        If saturatedBorderWidth is non-zero, replace saturated pixels with saturatedPixelValue.
        Note that replacing saturated pixels requires that the input images be MaskedImages.
    saturatedPixelValue : float
        Value to replace saturated pixels with.
    xSize : int
        Desired width of RGB image (or None).  If ySize is None, preserve aspect ratio.
    ySize : int
        Desired height of RGB image (or None).
    rescaleFactor : float
        Make size of output image rescaleFactor*size of the input image.
        Cannot be specified if xSize or ySize are given.

    Returns
    -------
    rgb : `~numpy.ndarray`
        RGB color image.
    """
    if imageG is None:
        imageG = imageR
    if imageB is None:
        imageB = imageR

    if saturatedBorderWidth:
        if saturatedPixelValue is None:
            raise ValueError("saturatedPixelValue must be set if "
                             "saturatedBorderWidth is set.")
        msg = "Cannot do this until we extract replaceSaturatedPixels out of afw/display/saturated.cc"
        raise NotImplementedError(msg)
        # replaceSaturatedPixels(imageR, imageG, imageB, saturatedBorderWidth, saturatedPixelValue)

    asinhMap = AsinhMapping(minimum, dataRange, Q)
    rgb = asinhMap.makeRgbImage(imageR, imageG, imageB,
                                xSize=xSize, ySize=ySize, rescaleFactor=rescaleFactor)

    if fileName:
        import matplotlib.image
        matplotlib.image.imsave(fileName, rgb)

    return rgb
