from ... import wcs
from ... import coordinates as coor
from ... import units as u
import numpy as np

def test_coor_wcs():
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [-234.75, 8.3393]
    w.wcs.cdelt = np.array([-0.066667, 0.066667])
    w.wcs.crval = [0, -90]
    w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    w.wcs.set_pv([(2, 1, 45.0)])

    coor_obj = coor.ICRS('0, -90', unit=(u.deg, u.deg))
    coor_arr = [[0,-90]]
    obj_result = w.wcs_world2pix(coor_obj, 1)
    arr_result = w.wcs_world2pix(coor_arr, 1)
    answer = np.array([[-234.75, 8.3393]])
    close_enough = 1e-8

    assert np.all(np.abs(answer-obj_result) < close_enough)
    assert np.all(np.abs(answer-arr_result) < close_enough)
