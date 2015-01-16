# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
from ...wcs import WCS
from .. import utils
from ...tests.helper import pytest
import numpy as np

def test_wcs_dropping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4,4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1,5))
    pc = wcs.wcs.pc # for later use below

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2,3,4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,3,4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    dropped = wcs.dropaxis(0)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([2,3,4]))
    dropped = wcs.dropaxis(1)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,3,4]))
    dropped = wcs.dropaxis(2)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,4]))
    dropped = wcs.dropaxis(3)
    assert np.all(dropped.wcs.get_pc().diagonal() == np.array([1,2,3]))

def test_wcs_swapping():
    wcs = WCS(naxis=4)
    wcs.wcs.pc = np.zeros([4,4])
    np.fill_diagonal(wcs.wcs.pc, np.arange(1,5))
    pc = wcs.wcs.pc # for later use below

    swapped = wcs.swapaxes(0,1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2,1,3,4]))
    swapped = wcs.swapaxes(0,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4,2,3,1]))
    swapped = wcs.swapaxes(2,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1,2,4,3]))

    wcs = WCS(naxis=4)
    wcs.wcs.cd = pc

    swapped = wcs.swapaxes(0,1)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([2,1,3,4]))
    swapped = wcs.swapaxes(0,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([4,2,3,1]))
    swapped = wcs.swapaxes(2,3)
    assert np.all(swapped.wcs.get_pc().diagonal() == np.array([1,2,4,3]))

@pytest.mark.parametrize('ndim',(2,3))
def test_add_stokes(ndim):
    wcs = WCS(naxis=ndim)

    for ii in range(ndim+1):
        outwcs = utils.add_stokes_axis_to_wcs(wcs,ii)
        assert outwcs.wcs.naxis == ndim+1
        assert outwcs.wcs.ctype[ii] == 'STOKES'
        assert outwcs.wcs.cname[ii] == 'STOKES'

def test_slice():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs.slice([slice(1,None),slice(0,None)])
    assert np.all(slice_wcs.wcs.crpix == np.array([1,0]))

    slice_wcs = mywcs.slice([slice(1,None,2),slice(0,None,4)])
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625, 0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    slice_wcs = mywcs.slice([slice(None,None,2),slice(0,None,2)])
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.2]))

def test_slice_getitem():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.625,0.25]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    mywcs.wcs.crpix = [2,2]
    slice_wcs = mywcs[1::2, 0::4]
    assert np.all(slice_wcs.wcs.crpix == np.array([0.875,0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.4,0.2]))

    # Default: numpy order
    slice_wcs = mywcs[1::2]
    assert np.all(slice_wcs.wcs.crpix == np.array([2,0.75]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.1,0.2]))

def test_slice_fitsorder():
    mywcs = WCS(naxis=2)
    mywcs.wcs.crval = [1,1]
    mywcs.wcs.cdelt = [0.1,0.1]
    mywcs.wcs.crpix = [1,1]

    slice_wcs = mywcs.slice([slice(1,None),slice(0,None)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0,1]))

    slice_wcs = mywcs.slice([slice(1,None,2),slice(0,None,4)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25,0.625]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.4]))

    slice_wcs = mywcs.slice([slice(1,None,2)], numpy_order=False)
    assert np.all(slice_wcs.wcs.crpix == np.array([0.25,1]))
    assert np.all(slice_wcs.wcs.cdelt == np.array([0.2,0.1]))

def test_invalid_slice():
    mywcs = WCS(naxis=2)

    with pytest.raises(ValueError) as exc:
        mywcs[0]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")

    with pytest.raises(ValueError) as exc:
        mywcs[0,::2]
    assert exc.value.args[0] == ("Cannot downsample a WCS with indexing.  Use "
                                 "wcs.sub or wcs.dropaxis if you want to remove "
                                 "axes.")

def test_axis_names():
    mywcs = WCS(naxis=4)
    mywcs.wcs.ctype = ['RA---TAN','DEC---TAN','VOPT-LSR','STOKES']

    assert mywcs.axis_type_names == ['RA','DEC','VOPT','STOKES']

    mywcs.wcs.cname = ['RA','DEC','VOPT-LSR','STOKES']

    assert mywcs.axis_type_names == ['RA','DEC','VOPT-LSR','STOKES']
