# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
from .wcs import WCS

wcs_parameters_to_preserve = ['cel_offset','dateavg','dateobs','equinox',
                              'latpole', 'lonpole', 'mjdavg', 'mjdobs', 'name',
                              'obsgeo', 'phi0', 'radesys', 'restfrq',
                              'restwav', 'specsys', 'ssysobs', 'ssyssrc',
                              'theta0', 'velangl', 'velosys', 'zsource']


def add_stokes_axis_to_wcs(wcs, add_before_ind):
    """
    Add a new Stokes axis that is uncorrelated with any other axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.WCS`
        The WCS to add to
    add_before_ind : int
        Index of the WCS to insert the new Stokes axis in front of.
        To add at the end, do add_before_ind = wcs.wcs.naxis

    Returns
    -------
    A new `~astropy.wcs.WCS` instance with an additional axis
    """

    naxin = wcs.wcs.naxis
    naxout = naxin+1

    inds = range(naxout)
    inds.pop(add_before_ind)
    inds = np.array(inds)

    outwcs = WCS(naxis=naxout)
    for par in wcs_parameters_to_preserve:
        setattr(outwcs.wcs, par, getattr(wcs.wcs,par))

    pc = np.zeros([naxout,naxout])
    pc[inds[:,np.newaxis],inds[np.newaxis,:]] = wcs.wcs.get_pc()
    pc[add_before_ind,add_before_ind] = 1

    def insert_at_index(val, index, lst):
        """ insert a value at index into a list """
        return list(lst)[:index] + [val] + list(lst)[index:]

    outwcs.wcs.crpix = insert_at_index(1, add_before_ind, wcs.wcs.crpix)
    outwcs.wcs.cdelt = insert_at_index(1, add_before_ind, wcs.wcs.get_cdelt())
    outwcs.wcs.crval = insert_at_index(1, add_before_ind, wcs.wcs.crval)
    outwcs.wcs.cunit = insert_at_index("", add_before_ind, wcs.wcs.cunit)
    outwcs.wcs.ctype = insert_at_index("STOKES", add_before_ind, wcs.wcs.ctype)
    outwcs.wcs.cname = insert_at_index("STOKES", add_before_ind, wcs.wcs.cname)
    outwcs.wcs.pc = pc

    return outwcs


def wcs_swapaxes(wcs, ax0, ax1):
    """
    Swap axes in a WCS.

    Parameters
    ----------
    wcs: `~astropy.wcs.WCS`
        The WCS to have its axes swapped
    ax0: int
    ax1: int
        The indices of the WCS to be swapped, counting from 0 (i.e., python
        convention, not FITS convention)

    Returns
    -------
    A new `~astropy.wcs.WCS` instance with the same number of axes, but two
    swapped
    """
    inds = range(wcs.wcs.naxis)
    inds[ax0],inds[ax1] = inds[ax1],inds[ax0]
    inds = np.array(inds)

    return reindex_wcs(wcs, inds)


def reindex_wcs(wcs, inds):
    """
    Re-index a WCS given indices.  The number of axes may be reduced.

    Parameters
    ----------
    wcs: `~astropy.wcs.WCS`
        The WCS to be manipulated
    inds: np.ndarray(dtype='int')
        The indices of the array to keep in the output.
        e.g. swapaxes: [0,2,1,3]
        dropaxes: [0,1,3]

    Returns
    -------
    A new `~astropy.wcs.WCS` instance
    """
    if not isinstance(inds, np.ndarray):
        raise TypeError("Indices must be an ndarray")
    if inds.dtype.kind != 'i':
        raise TypeError('Indices must be integers')

    outwcs = WCS(naxis=len(inds))
    for par in wcs_parameters_to_preserve:
        setattr(outwcs.wcs, par, getattr(wcs.wcs,par))

    cdelt = wcs.wcs.get_cdelt()
    pc = wcs.wcs.get_pc()

    outwcs.wcs.crpix = wcs.wcs.crpix[inds]
    outwcs.wcs.cdelt = cdelt[inds]
    outwcs.wcs.crval = wcs.wcs.crval[inds]
    outwcs.wcs.cunit = [wcs.wcs.cunit[i] for i in inds]
    outwcs.wcs.ctype = [wcs.wcs.ctype[i] for i in inds]
    outwcs.wcs.cname = [wcs.wcs.cname[i] for i in inds]
    outwcs.wcs.pc = pc[inds[:,None],inds[None,:]]

    return outwcs

def axis_names(wcs):
    """
    Extract world names for each coordinate axis

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        The WCS object to extract names from

    Returns
    -------
    A tuple of names along each axis
    """
    names = list(wcs.wcs.cname)
    types = wcs.wcs.ctype
    for i in range(len(names)):
        if len(names[i]) > 0:
            continue
        names[i] = types[i].split('-')[0]
    return names


def slice_wcs(wcs, view):
    """
    Slice a WCS instance using a Numpy slice. The order of the slice should
    be reversed (as for the data) compared to the natural WCS order.

    Parameters
    ----------
    view : tuple
        A tuple containing the same number of slices as the WCS system.
        The `step` method, the third argument to a slice, is not presently
        supported.

    Returns
    -------
    A new `~astropy.wcs.WCS` instance
    """
    if len(view) != wcs.wcs.naxis:
        raise ValueError("Must have same number of slices as number of WCS axes")

    wcs_new = wcs.deepcopy()
    for i, iview in enumerate(view):
        if iview.start is not None:
            if iview.step not in (None, 1):
                raise NotImplementedError("Cannot yet slice WCS with strides different from None or 1")
            wcs_index = wcs.wcs.naxis - 1 - i
            wcs_new.wcs.crpix[wcs_index] -= iview.start
    return wcs_new
