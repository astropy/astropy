# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Map 3 images produce a properly-scaled RGB image following Lupton et al. (2004).

For details, see : http://adsabs.harvard.edu/abs/2004PASP..116..133L
"""

import numpy as np

try:
    import scipy.misc
    HAVE_SCIPY_MISC = True
except ImportError:
    HAVE_SCIPY_MISC = False

# from lsst.afw.display.displayLib import replaceSaturatedPixels, getZScale


def compute_intensity(imageR, imageG=None, imageB=None):
    """!Return a naive total intensity from the red, blue, and green intensities
    \param imageR intensity of image that'll be mapped to red; or intensity if imageG and imageB are None
    \param imageG intensity of image that'll be mapped to green; or None
    \param imageB intensity of image that'll be mapped to blue; or None

    Inputs may be MaskedImages, Images, or numpy arrays and the return is of the same type
    """
    if imageG is None or imageB is None:
        assert imageG is None and imageB is None, \
            "Please specify either a single image or red, green, and blue images"
        return imageR

    intensity = (imageR + imageG + imageB)/3.0

    # Repack into whatever type was passed to us
    return np.array(intensity, dtype=imageR.dtype)


def zscale(image, nSamples=1000, contrast=0.25):
    """
    This emulates ds9's zscale feature. Returns the suggested minimum and
    maximum values to display.

    Parameters
    ----------
    image : `~numpy.ndarray`
        The image to compute the scaling on.
    nSamples :  int
        How many samples to take when building the histogram.
    contrast : float
        ???
    """

    stride = image.size/nSamples
    samples = image.flatten()[::stride]
    samples.sort()
    chop_size = int(0.10*len(samples))
    subset = samples[chop_size:-chop_size]

    i_midpoint = int(len(subset)/2)
    I_mid = subset[i_midpoint]

    fit = np.polyfit(np.arange(len(subset)) - i_midpoint, subset, 1)
    # fit = [ slope, intercept]

    z1 = I_mid + fit[0]/contrast * (1-i_midpoint)/1.0
    z2 = I_mid + fit[0]/contrast * (len(subset)-i_midpoint)/1.0
    return z1, z2


class Mapping(object):
    """!Baseclass to map red, blue, green intensities into uint8 values"""

    def __init__(self, minimum=None, image=None):
        """!Create a mapping
        \param minimum  Intensity that should be mapped to black (a scalar or array for R, G, B)
        \param image The image to be used to calculate the mapping.

        If provided, also the default for makeRgbImage()
        """
        self._uint8Max = float(np.iinfo(np.uint8).max)

        try:
            len(minimum)
        except:
            minimum = 3*[minimum]
        assert len(minimum) == 3, "Please provide 1 or 3 values for minimum"

        self.minimum = minimum
        self._image = image

    def makeRgbImage(self, imageR=None, imageG=None, imageB=None,
                     xSize=None, ySize=None, rescaleFactor=None):
        """!Convert 3 arrays, imageR, imageG, and imageB into a numpy RGB image
        \param imageR Image to map to red (if None, use the image passed to the constructor)
        \param imageG Image to map to green (if None, use imageR)
        \param imageB Image to map to blue (if None, use imageR)
        \param xSize  Desired width of RGB image (or None).  If ySize is None, preserve aspect ratio.
        \param ySize  Desired height of RGB image (or None)
        \param rescaleFactor Make size of output image rescaleFactor*size of the input image (or None)
        """
        if imageR is None:
            if self._image is None:
                raise RuntimeError("You must provide an image (or pass one to the constructor)")
            imageR = self._image

        if imageG is None:
            imageG = imageR
        if imageB is None:
            imageB = imageR

        imageRGB = [imageR, imageG, imageB]

        if xSize is not None or ySize is not None:
            assert rescaleFactor is None, "You may not specify a size and rescaleFactor"
            h, w = imageRGB[0].shape
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
                raise RuntimeError("Unable to rescale as scipy.misc is unavailable.")

            for i, im in enumerate(imageRGB):
                imageRGB[i] = scipy.misc.imresize(im, size, interp='bilinear', mode='F')

        return np.dstack(self._convertImagesToUint8(*imageRGB)).astype(np.uint8)

    def intensity(self, imageR, imageG, imageB):
        """!Return the total intensity from the red, blue, and green intensities

        This is a naive computation, and may be overridden by subclasses
        """
        return compute_intensity(imageR, imageG, imageB)

    def mapIntensityToUint8(self, I):
        """Map an intensity into the range of a uint8, [0, 255] (but not converted to uint8)"""
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
    """!A linear map map of red, blue, green intensities into uint8 values"""

    def __init__(self, minimum=None, maximum=None, image=None):
        """!A linear stretch from [minimum, maximum]; if one or both are omitted use image minimum/maximum to set them

        \param minimum  Intensity that should be mapped to black (a scalar or array for R, G, B)
        \param maximum  Intensity that should be mapped to white (a scalar)
        """

        if minimum is None or maximum is None:
            assert image is not None, "You must provide an image if you don't set both minimum and maximum"

            if minimum is None:
                minimum = image.min()
            if maximum is None:
                maximum = image.max()

        Mapping.__init__(self, minimum, image)
        self.maximum = maximum

        if maximum is None:
            self._range = None
        else:
            assert maximum - minimum != 0, "minimum and maximum values must not be equal"
            self._range = float(maximum - minimum)

    def mapIntensityToUint8(self, I):
        """Return an array which, when multiplied by an image, returns that image mapped to the range of a
        uint8, [0, 255] (but not converted to uint8)

        The intensity is assumed to have had minimum subtracted (as that can be done per-band)
        """
        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            return np.where(I <= 0, 0,
                            np.where(I >= self._range, self._uint8Max/I, self._uint8Max/self._range))


class ZScaleMapping(LinearMapping):
    """!A mapping for a linear stretch chosen by the zscale algorithm
    (preserving colours independent of brightness)

    x = (I - minimum)/range
    """

    def __init__(self, image, nSamples=1000, contrast=0.25):
        """!A linear stretch from [z1, z2] chosen by the zscale algorithm
        \param nSamples The number of samples to use to estimate the zscale parameters
        \param contrast The number of samples to use to estimate the zscale parameters
        """

        z1, z2 = zscale(image, nSamples, contrast)
        LinearMapping.__init__(self, z1, z2, image)


class AsinhMapping(Mapping):
    """!A mapping for an asinh stretch (preserving colours independent of brightness)

    x = asinh(Q (I - minimum)/range)/Q

    This reduces to a linear stretch if Q == 0

    See http://adsabs.harvard.edu/abs/2004PASP..116..133L
    """

    def __init__(self, minimum, dataRange, Q=8):
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
        """Return an array which, when multiplied by an image, returns that image mapped to the range of a
        uint8, [0, 255] (but not converted to uint8)

        The intensity is assumed to have had minimum subtracted (as that can be done per-band)
        """
        with np.errstate(invalid='ignore', divide='ignore'):  # n.b. np.where can't and doesn't short-circuit
            return np.where(I <= 0, 0, np.arcsinh(I*self._soften)*self._slope/I)


class AsinhZScaleMapping(AsinhMapping):
    """!A mapping for an asinh stretch, estimating the linear stretch by zscale

    x = asinh(Q (I - z1)/(z2 - z1))/Q

    See AsinhMapping
    """

    def __init__(self, image1, image2=None, image3=None, Q=8, pedestal=None):
        """!
        Create an asinh mapping from an image, setting the linear part of the stretch using zscale.

        \param image1 The image to analyse,
         # or a list of 3 images to be converted to an intensity image
        \param image2 the second image to analyse (must be specified with image3)
        \param image3 the third image to analyse (must be specified with image2)
        \param Q The asinh softening parameter
        \param pedestal The value, or array of 3 values, to subtract from the images; or None

        N.b. pedestal, if not None, is removed from the images when calculating the zscale
        stretch, and added back into Mapping.minimum[]
        """

        if image2 is None or image3 is None:
            assert image2 is None and image3 is None, "Please specify either a single image or three images"
            image = [image1]
        else:
            image = [image1, image2, image3]

        if pedestal is not None:
            try:
                assert len(pedestal) in (1, 3,), "Please provide 1 or 3 pedestals"
            except TypeError:
                pedestal = 3*[pedestal]

            image = list(image)        # needs to be mutable
            for i, im in enumerate(image):
                if pedestal[i] != 0.0:
                    image[i] = im - pedestal[i]  # n.b. a copy
        else:
            pedestal = len(image)*[0.0]

        image = compute_intensity(*image)

        zscale = ZScaleMapping(image)
        dataRange = zscale.maximum - zscale.minimum[0]  # zscale.minimum is always a triple
        minimum = zscale.minimum

        for i, level in enumerate(pedestal):
            minimum[i] += level

        AsinhMapping.__init__(self, minimum, dataRange, Q)
        self._image = image             # support self.makeRgbImage()


def makeRGB(imageR, imageG=None, imageB=None, minimum=0, dataRange=5, Q=8, fileName=None,
            saturatedBorderWidth=0, saturatedPixelValue=None,
            xSize=None, ySize=None, rescaleFactor=None):
    """Make a set of three images into an RGB image using an asinh stretch and optionally write it to disk

    If saturatedBorderWidth is non-zero, replace saturated pixels with saturatedPixelValue.  Note
    that replacing saturated pixels requires that the input images be MaskedImages.
    """
    if imageG is None:
        imageG = imageR
    if imageB is None:
        imageB = imageR

    if saturatedBorderWidth:
        if saturatedPixelValue is None:
            raise ValueError("saturatedPixelValue must be set if saturatedBorderWidth is set")
        msg = "Cannot do this until we extract replaceSaturatedPixels out of afw/display/saturated.cc"
        raise NotImplementedError(msg)
        # replaceSaturatedPixels(imageR, imageG, imageB, saturatedBorderWidth, saturatedPixelValue)

    asinhMap = AsinhMapping(minimum, dataRange, Q)
    rgb = asinhMap.makeRgbImage(imageR, imageG, imageB,
                                xSize=xSize, ySize=ySize, rescaleFactor=rescaleFactor)

    if fileName:
        writeRGB(fileName, rgb)

    return rgb


def displayRGB(rgb, show=True):
    """!Display an rgb image using matplotlib
    \param rgb  The RGB image in question
    \param show If true, call plt.show()
    """
    import matplotlib.pyplot as plt
    plt.imshow(rgb, interpolation='nearest', origin="lower")
    if show:
        plt.show()
    return plt


def writeRGB(fileName, rgbImage):
    """!Write an RGB image to disk
    \param fileName The output file.  The suffix defines the format, and must be supported by matplotlib
    \param rgbImage The image, as made by e.g. makeRGB

    Most versions of matplotlib support png and pdf (although the eps/pdf/svg writers may be buggy,
    possibly due an interaction with useTeX=True in the matplotlib settings).

    If your matplotlib bundles pil/pillow you should also be able to write jpeg and tiff files.
    """
    import matplotlib.image
    matplotlib.image.imsave(fileName, rgbImage)
