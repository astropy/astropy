from ... import wcs
from ... import coordinates as coor
from ... import units as u
import numpy as np
from ...tests.helper import pytest, raises

def test_coor_wcs():
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [-234.75, 8.3393]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [0, -90]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    w.wcs.set_pv([(2, 1, 45.0)])

    close_enough = 1e-8

    # Test for single coordinate object
    single_coor = [[0,-90]]

    single_coor_obj = coor.ICRS('0, -90', unit=(u.deg, u.deg))

    single_obj_result = w.wcs_world2pix(single_coor_obj, 1)
    single_result = w.wcs_world2pix(single_coor, 1)

    assert np.all(np.abs(single_result-single_obj_result) < close_enough)

    # Test for coordinate object conatining array of coordinates
    ra = [0, 10, 20]
    dec = [-90, -70, -10]

    array_coor_obj = coor.ICRS(ra, dec, unit=(u.deg, u.deg))

    array_result = w.wcs_world2pix(ra, dec, 1)
    array_coor_result = w.wcs_world2pix(array_coor_obj, 1)

    assert np.all(np.abs(array_result[0]-array_coor_result[0]) < close_enough)
    assert np.all(np.abs(array_result[1]-array_coor_result[1]) < close_enough)

@raises(TypeError)
def test_incorrect_input():
    # Test to check if it raises error when coordinate objects
    # provided as arguments to the wrong function
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [-234.75, 8.3393]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [0, -90]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    w.wcs.set_pv([(2, 1, 45.0)])

    ra = [0, 10, 20]
    dec = [-90, -70, -10]

    single_coor_obj = coor.ICRS('0, -90', unit=(u.deg, u.deg))
    array_coor_obj = coor.ICRS(ra, dec, unit=(u.deg, u.deg))

    err_obj_arr_result = w.wcs_pix2world(array_coor_obj, 1)
    err_obj_result = w.wcs_pix2world(single_coor_obj, 1)
