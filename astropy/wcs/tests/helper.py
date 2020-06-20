# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

import numpy as np

from astropy import wcs
from astropy.io import fits


class SimModelTAB:
    def __init__(self, nx=150, ny=200, crpix=[1, 1], crval = [1, 1],
                 cdelt = [1, 1], pc = {'PC1_1': 1, 'PC2_2': 1}):
        """  set essential parameters of the model (coord transformations)  """
        assert nx > 2 and ny > 1  # a limitation of this particular simulation
        self.nx = nx
        self.ny = ny
        self.crpix = crpix
        self.crval = crval
        self.cdelt = cdelt
        self.pc = pc

    def fwd_eval(self, xy):
        xb = 1 + self.nx // 3
        px = np.array([1, xb, xb, self.nx + 1])
        py = np.array([1, self.ny + 1])

        xi = self.crval[0] + self.cdelt[0] * (px - self.crpix[0])
        yi = self.crval[1] + self.cdelt[1] * (py - self.crpix[1])

        cx = np.array([0.0, 0.26, 0.8, 1.0])
        cy = np.array([-0.5, 0.5])

        xy = np.atleast_2d(xy)
        x = xy[:, 0]
        y = xy[:, 1]

        mbad = (x < px[0]) | (y < py[0]) | (x > px[-1]) | (y > py[-1])
        mgood = np.logical_not(mbad)

        i = 2 * (x > xb).astype(int)

        psix = self.crval[0] + self.cdelt[0] * (x - self.crpix[0])
        psiy = self.crval[1] + self.cdelt[1] * (y - self.crpix[1])

        cfx = (psix - xi[i]) / (xi[i + 1] - xi[i])
        cfy = (psiy - yi[0]) / (yi[1] - yi[0])

        ra = cx[i] + cfx * (cx[i + 1] - cx[i])
        dec = cy[0] + cfy * (cy[1] - cy[0])

        return np.dstack([ra, dec])[0]

    @property
    def hdulist(self):
        """ Simulates 2D data with a _spatial_ WCS that uses the ``-TAB``
        algorithm with indexing.
        """
        # coordinate array (some "arbitrary" numbers with a "jump" along x axis):
        x = np.array([[0.0, 0.26, 0.8, 1.0], [0.0, 0.26, 0.8, 1.0]])
        y = np.array([[-0.5, -0.5, -0.5, -0.5], [0.5, 0.5, 0.5, 0.5]])
        c = np.dstack([x, y])

        # index arrays (skip PC matrix for simplicity - assume it is an
        # identity matrix):
        xb = 1 + self.nx // 3
        px = np.array([1, xb, xb, self.nx + 1])
        py = np.array([1, self.ny + 1])
        xi = self.crval[0] + self.cdelt[0] * (px - self.crpix[0])
        yi = self.crval[1] + self.cdelt[1] * (py - self.crpix[1])

        # structured array (data) for binary table HDU:
        arr = np.array(
            [(c, xi, yi)],
            dtype=[
                ('wavelength', np.float64, c.shape),
                ('xi', np.double, (xi.size,)),
                ('yi', np.double, (yi.size,))
            ]
        )

        # create binary table HDU:
        bt = fits.BinTableHDU(arr);
        bt.header['EXTNAME'] = 'WCS-TABLE'

        # create primary header:
        image_data = np.ones((self.ny, self.nx), dtype=np.float32)
        pu = fits.PrimaryHDU(image_data)
        pu.header['ctype1'] = 'RA---TAB'
        pu.header['ctype2'] = 'DEC--TAB'
        pu.header['naxis1'] = self.nx
        pu.header['naxis2'] = self.ny
        pu.header['PS1_0'] = 'WCS-TABLE'
        pu.header['PS2_0'] = 'WCS-TABLE'
        pu.header['PS1_1'] = 'wavelength'
        pu.header['PS2_1'] = 'wavelength'
        pu.header['PV1_3'] = 1
        pu.header['PV2_3'] = 2
        pu.header['CUNIT1'] = 'deg'
        pu.header['CUNIT2'] = 'deg'
        pu.header['CDELT1'] = self.cdelt[0]
        pu.header['CDELT2'] = self.cdelt[1]
        pu.header['CRPIX1'] = self.crpix[0]
        pu.header['CRPIX2'] = self.crpix[1]
        pu.header['CRVAL1'] = self.crval[0]
        pu.header['CRVAL2'] = self.crval[1]
        pu.header['PS1_2'] = 'xi'
        pu.header['PS2_2'] = 'yi'
        for k, v in self.pc.items():
            pu.header[k] = v

        hdulist = fits.HDUList([pu, bt])
        return hdulist
