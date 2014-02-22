import numpy as np

from astropy.wcs import WCS


def select_axes(iterable, dimensions):
    return [iterable[i] for i in dimensions]


class WCSParameters(object):

    def __init__(self, wcs, dimensions):

        self.ctype = select_axes(wcs.ctype, dimensions)
        self.crval = select_axes(wcs.crval, dimensions)
        self.crpix = select_axes(wcs.crpix, dimensions)
        self.cdelt = select_axes(wcs.cdelt, dimensions)
        self.cunit = select_axes(wcs.cunit, dimensions)

        self.naxis = wcs.naxis


class WCSSlice(object):

    # The purpose of this class is to wrap n-dimensional WCS objects into a
    # 2-dimensional WCS object.

    def __init__(self, *args, **kwargs):

        if 'slices' in kwargs:
            self._slices = kwargs.pop('slices')

        if 'dimensions' in kwargs:
            self._dimensions = kwargs.pop('dimensions')

        self._wcs_orig = WCS(*args, **kwargs)

        # Now find the values of the coordinates in the slices - only needed if
        # data has more than two dimensions
        if len(self._slices) > 0:

            self.nx = args[0]['NAXIS%i' % (self._dimensions[0] + 1)]
            self.ny = args[0]['NAXIS%i' % (self._dimensions[1] + 1)]
            xpix = np.arange(self.nx) + 1.
            ypix = np.arange(self.ny) + 1.
            xpix, ypix = np.meshgrid(xpix, ypix)
            xpix, ypix = xpix.reshape(self.nx * self.ny), ypix.reshape(self.nx * self.ny)
            s = 0
            coords = []
            for dim in range(self._wcs_orig.naxis):
                if dim == self._dimensions[0]:
                    coords.append(xpix)
                elif dim == self._dimensions[1]:
                    coords.append(ypix)
                else:
                    coords.append(np.repeat(self._slices[s], xpix.shape))
                    s += 1
            coords = np.vstack(coords).transpose()
            result = self._wcs_orig.wcs_pix2world(coords, 1)
            self._mean_world = np.mean(result, axis=0)

        # Now set up fake .wcs attribute

        self.wcs = WCSParameters(self._wcs_orig.wcs, self._dimensions)

    def wcs_world2pix(self, x, y, origin):
        if self._wcs_orig.naxis == 2:
            if self._dimensions[1] < self._dimensions[0]:
                xp, yp = self._wcs_orig.wcs_world2pix(y, x, origin)
                return yp, xp
            else:
                return self._wcs_orig.wcs_world2pix(x, y, origin)
        else:
            coords = []
            s = 0
            for dim in range(self._wcs_orig.naxis):
                if dim == self._dimensions[0]:
                    coords.append(x)
                elif dim == self._dimensions[1]:
                    coords.append(y)
                else:
                    # The following is an approximation, and will break down if
                    # the world coordinate changes significantly over the slice
                    coords.append(np.repeat(self._mean_world[dim], x.shape))
                    s += 1
            coords = np.vstack(coords).transpose()

            # Due to a bug in pywcs, we need to loop over each coordinate
            # result = AstropyWCS.wcs_world2pix(self, coords, origin)
            result = np.zeros(coords.shape)
            for i in range(result.shape[0]):
                result[i:i + 1, :] = self._wcs_orig.wcs_world2pix(coords[i:i + 1, :], origin)

            return result[:, self._dimensions[0]], result[:, self._dimensions[1]]

    def wcs_pix2world(self, x, y, origin):
        if self._wcs_orig.naxis == 2:
            if self._dimensions[1] < self._dimensions[0]:
                xw, yw = self._wcs_orig.wcs_pix2world(y, x, origin)
                return yw, xw
            else:
                return self._wcs_orig.wcs_pix2world(x, y, origin)
        else:
            coords = []
            s = 0
            for dim in range(self._wcs_orig.naxis):
                if dim == self._dimensions[0]:
                    coords.append(x)
                elif dim == self._dimensions[1]:
                    coords.append(y)
                else:
                    coords.append(np.repeat(self._slices[s] + 0.5, x.shape))
                    s += 1
            coords = np.vstack(coords).transpose()
            result = self._wcs_orig.wcs_pix2world(coords, origin)
            return result[:, self._dimensions[0]], result[:, self._dimensions[1]]
