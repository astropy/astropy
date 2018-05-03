# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.

import copy

import numpy as np

from ... import log

__all__ = ['NDSlicingMixin']


class NDSlicingMixin:
    """Mixin to provide slicing on objects using the `NDData`
    interface.

    The ``data``, ``mask``, ``uncertainty`` and ``wcs`` will be sliced, if
    set and sliceable. The ``unit`` and ``meta`` will be untouched. The return
    will be a reference and not a copy, if possible.

    Examples
    --------
    Using this Mixin with `~astropy.nddata.NDData`:

        >>> from astropy.nddata import NDData, NDSlicingMixin
        >>> class NDDataSliceable(NDSlicingMixin, NDData):
        ...     pass

    Slicing an instance containing data::

        >>> nd = NDDataSliceable([1,2,3,4,5])
        >>> nd[1:3]
        NDDataSliceable([2, 3])

    Also the other attributes are sliced for example the ``mask``::

        >>> import numpy as np
        >>> mask = np.array([True, False, True, True, False])
        >>> nd2 = NDDataSliceable(nd, mask=mask)
        >>> nd2slc = nd2[1:3]
        >>> nd2slc[nd2slc.mask]
        NDDataSliceable([3])

    Be aware that changing values of the sliced instance will change the values
    of the original::

        >>> nd3 = nd2[1:3]
        >>> nd3.data[0] = 100
        >>> nd2
        NDDataSliceable([  1, 100,   3,   4,   5])

    See also
    --------
    NDDataRef
    NDDataArray
    """
    def __getitem__(self, item):
        # Forbid slicing with None, i.e. adding a dummy dimension,
        # as it messes up the WCS slicing.
        if item is None or (isinstance(item, tuple) and None in item):
            raise IndexError("None indices not supported")

        # Abort slicing if the data is a single scalar.
        if self.data.shape == ():
            raise TypeError('scalars cannot be sliced.')

        # Let the other methods handle slicing.
        kwargs = self._slice(item)
        return self.__class__(**kwargs)

    def _slice(self, item):
        """Collects the sliced attributes and passes them back as `dict`.

        It passes uncertainty, mask and wcs to their appropriate ``_slice_*``
        method, while ``meta`` and ``unit`` are simply taken from the original.
        The data is assumed to be sliceable and is sliced directly.

        When possible the return should *not* be a copy of the data but a
        reference.

        Parameters
        ----------
        item : slice
            The slice passed to ``__getitem__``.

        Returns
        -------
        dict :
            Containing all the attributes after slicing - ready to
            use them to create ``self.__class__.__init__(**kwargs)`` in
            ``__getitem__``.
        """
        kwargs = {}
        kwargs['data'] = self.data[item]
        # Try to slice some attributes
        kwargs['uncertainty'] = self._slice_uncertainty(item)
        kwargs['mask'] = self._slice_mask(item)
        # Attributes which are copied and not intended to be sliced
        kwargs['unit'] = self.unit
        kwargs['meta'] = self.meta

        # Slice WCS
        wcs, missing_axes = self._slice_wcs_missing_axes(item)
        kwargs['wcs'] = wcs
        kwargs['missing_axes'] = missing_axes

        return kwargs

    def _slice_uncertainty(self, item):
        if self.uncertainty is None:
            return None
        try:
            return self.uncertainty[item]
        except TypeError:
            # Catching TypeError in case the object has no __getitem__ method.
            # But let IndexError raise.
            log.info("uncertainty cannot be sliced.")
        return self.uncertainty

    def _slice_mask(self, item):
        if self.mask is None:
            return None
        try:
            return self.mask[item]
        except TypeError:
            log.info("mask cannot be sliced.")
        return self.mask

    def _slice_wcs_missing_axes(self, item):
        # here missing axis is reversed as the item comes already in the reverse order
        # of the input
        return _wcs_slicer(
            self.wcs, copy.deepcopy(self.missing_axes[::-1]), item)

def _wcs_slicer(wcs, missing_axes, item):
    """
    Returns the new sliced wcs and changed missing axis.

    Paramters
    ---------
    wcs: `astropy.wcs.WCS` or `ndcube.utils.wcs.WCS`
        WCS object to be sliced.

    missing_axes: `list` of `bool`
        Indicates which axes of the WCS are "missing", i.e. do not correspond to a data axis.

    item: `int`, `slice` or `tuple` of `int` and/or `slice`.
        Slicing item.  Note that unlike in other places in this package, the item has the
        same axis ordering as the WCS object, i.e. the reverse of the data order.

    Returns
    -------
    new_wcs: `astropy.wcs.WCS` or `ndcube.utils.wcs.WCS`
        Sliced WCS object.

    missing_axes: `list` of `bool`
        Altered missing axis list.  Note the ordering has been reversed to reflect the data
        (numpy) axis ordering convention.

    """
    # normal slice.
    item_checked = []
    if isinstance(item, slice):
        index = 0
        # creating a new tuple of slice where if the axis is dead i.e missing
        # then slice(0,1) added else slice(None, None, None) is appended and
        # if the check of missing_axes gives that this is the index where it
        # needs to be appended then it gets appended there.
        for _bool in missing_axes:
            if not _bool:
                if index is not 1:
                    item_checked.append(item)
                    index += 1
                else:
                    item_checked.append(slice(None, None, None))
            else:
                item_checked.append(slice(0, 1))
        new_wcs = wcs.slice((item_checked))
    # item is int then slicing axis.
    elif isinstance(item, int) or isinstance(item, np.int64):
        # using index to keep track of whether the int(which is converted to
        # slice(int_value, int_value+1)) is already added or not. It checks
        # the dead axis i.e missing_axes to check if it is dead than slice(0,1)
        # is appended in it. if the index value has reached 1 then the
        # slice(None, None, None) is added.
        index = 0
        for i, _bool in enumerate(missing_axes):
            if not _bool:
                if index is not 1:
                    item_checked.append(slice(item, item+1))
                    missing_axes[i] = True
                    index += 1
                else:
                    item_checked.append(slice(None, None, None))
            else:
                item_checked.append(slice(0, 1))
        new_wcs = wcs.slice(item_checked)
    # if it a tuple like [0:2, 0:3, 2] or [0:2, 1:3]
    elif isinstance(item, tuple):
        # this is used to not exceed the range of the item tuple
        # if the check of the missing_axes which is False if not dead
        # is a success than the the item of the tuple is added one by
        # one and if the end of tuple is reached than slice(None, None, None)
        # is appended.
        index = 0
        for _bool in missing_axes:
            if not _bool:
                if index is not len(item):
                    item_checked.append(item[index])
                    index += 1
                else:
                    item_checked.append(slice(None, None, None))
            else:
                item_checked.append(slice(0, 1))
        # if all are slice in the item tuple
        if _all_slice(item_checked):
            new_wcs = wcs.slice((item_checked))
        # if all are not slices some of them are int then
        else:
            # this will make all the item in item_checked as slice.
            item_ = _slice_list(item_checked)
            new_wcs = wcs.slice(item_)
            for i, it in enumerate(item_checked):
                if isinstance(it, int):
                    missing_axes[i] = True
    # returning the reverse list of missing axis as in the item here was reverse of
    # what was inputed so we had a reverse missing_axes.
    return new_wcs, missing_axes[::-1]


def _all_slice(obj):
    """
    Returns True if all the elements in the object are slices else return False
    """
    result = False
    if not isinstance(obj, (tuple, list)):
        return result
    result |= all(isinstance(o, slice) for o in obj)
    return result


def _slice_list(obj):
    """
    Return list of all the slices.

    Example
    -------
    >>> _slice_list((slice(1,2), slice(1,3), 2, slice(2,4), 8))
    [slice(1, 2, None), slice(1, 3, None), slice(2, 3, None), slice(2, 4, None), slice(8, 9, None)]
    """
    result = []
    if not isinstance(obj, (tuple, list)):
        return result
    for i, o in enumerate(obj):
        if isinstance(o, int):
            result.append(slice(o, o+1))
        elif isinstance(o, slice):
            result.append(o)
    return result
