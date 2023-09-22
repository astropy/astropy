# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.


from astropy import log
from astropy.wcs.wcsapi import BaseHighLevelWCS  # noqa: F401
from astropy.wcs.wcsapi import BaseLowLevelWCS  # noqa: F401
from astropy.wcs.wcsapi import HighLevelWCSWrapper, SlicedLowLevelWCS

__all__ = ["NDSlicingMixin"]


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

    See Also
    --------
    NDDataRef
    NDDataArray
    """

    def __getitem__(self, item):
        # Abort slicing if the data is a single scalar.
        if self.data.shape == ():
            raise TypeError("scalars cannot be sliced.")

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
        kwargs["data"] = self.data[item]
        # Try to slice some attributes
        kwargs["uncertainty"] = self._slice_uncertainty(item)
        kwargs["mask"] = self._slice_mask(item)
        kwargs["wcs"] = self._slice_wcs(item)
        # Attributes which are copied and not intended to be sliced
        kwargs["unit"] = self.unit
        kwargs["meta"] = self.meta
        return kwargs

    def _slice_uncertainty(self, item):
        if self.uncertainty is None:
            return None
        try:
            return self.uncertainty[item]
        except (TypeError, KeyError):
            # Catching TypeError in case the object has no __getitem__ method.
            # Catching KeyError for Python 3.12.
            # But let IndexError raise.
            log.info("uncertainty cannot be sliced.")
        return self.uncertainty

    def _slice_mask(self, item):
        if self.mask is None:
            return None
        try:
            return self.mask[item]
        except (TypeError, KeyError):
            log.info("mask cannot be sliced.")
        return self.mask

    def _slice_wcs(self, item):
        if self.wcs is None:
            return None

        try:
            llwcs = SlicedLowLevelWCS(self.wcs.low_level_wcs, item)
            return HighLevelWCSWrapper(llwcs)
        except Exception as err:
            self._handle_wcs_slicing_error(err, item)

    # Implement this in a method to allow subclasses to customise the error.
    def _handle_wcs_slicing_error(self, err, item):
        raise ValueError(
            f"Slicing the WCS object with the slice '{item}' "
            "failed, if you want to slice the NDData object without the WCS, you "
            "can remove by setting `NDData.wcs = None` and then retry."
        ) from err
