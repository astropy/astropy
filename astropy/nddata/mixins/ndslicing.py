# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..nduncertainty import NDUncertainty
from ...wcs.wcs import WCS
from ...extern.six import string_types
from ... import log

__all__ = ['NDSlicingMixin']


class NDSlicingMixin(object):
    """
    Mixin to provide slicing on objects using the NDData interface.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArray` for an example.

    Since `~astropy.nddata.NDData` is not very restrictive about what the
    additional attributes can be there are a lot of tests while slicing.
    If a subclass defines more restrictive setter methods for the properties
    they might want to consider altering `NDSlicingMixin`s methods for
    slicing (starting with ``_slice_*``). Currently implemented is support for
    slicing the:

    - data (which is enforced to be something like a `~numpy.ndarray`)
    - mask (if it is a numpy array)
    - wcs (if it is a numpy array or `~astropy.wcs.WCS` object)
    - uncertainty (if it is a numpy array or a subclass of `NDUncertainty`)

    But only if the shapes match the shape of the data (except for the wcs
    if it is a `WCS` object).
    """
    def __getitem__(self, item):
        # Abort slicing if the data is scalar
        if self._data.shape == ():
            raise TypeError('Scalars cannot be sliced.')
        # Slice everything defined in _slice_attr and the data
        new_data = self.data[item]
        kwargs = self._slice_attr(item)
        return self.__class__(new_data, **kwargs)

    def _slice_attr(self, item):
        """
        Collects the sliced attributes and passes them back as ``kwarg`` for
        the return after slicing.

        It passes uncertainty, mask and wcs to their appropriate ``_slice_*``
        method, while `meta` and `unit` are simply taken from the original. The
        data is assumed to be sliceable and is sliced before.

        Parameters
        ----------
        item: Slice
            The slice passed to ``__getitem__``.

        Returns
        -------
        kwargs: `dict`
            Containing all the attributes except data after slicing - ready to
            feed them into the ``init``.

        Notes
        -----
        Subclasses *normally* don't need to override this method instead they
        should override the specific ``_slice_*`` methods, except one wishes to
        slice other attributes as well.

        Where possible the return is *not* a copy of the data.
        """
        kwargs = {}
        # Try to slice some attributes
        kwargs['uncertainty'] = self._slice_uncertainty(item)
        kwargs['mask'] = self._slice_mask(item)
        kwargs['wcs'] = self._slice_wcs(item)
        # Attributes which are copied and not intended to be sliced
        kwargs['unit'] = self._unit
        kwargs['meta'] = self._meta
        return kwargs

    def _slice_uncertainty(self, item):
        # Uncertainty is not set
        if self._uncertainty is None:
            return None
        elif isinstance(self._uncertainty, np.ndarray):
            if self._uncertainty.shape == self.data.shape:
                return self._uncertainty[item]
            else:
                log.info("Uncertainty has not the same shape as data.")
        elif isinstance(self._uncertainty, NDUncertainty):
            if self._uncertainty.array.shape == self.data.shape:
                return self._uncertainty[item]
            else:
                log.info("Uncertainty has not the same shape as data.")
        else:
            log.info("Uncertainty is considered not sliceable.")
        log.info("Therefore the uncertainty will not be sliced.")
        return self._uncertainty

    def _slice_mask(self, item):
        if self._mask is None:
            return None
        elif isinstance(self._mask, np.ndarray):
            if self._mask.shape == self.data.shape:
                return self._mask[item]
            else:
                log.info("Mask has not the same shape as data.")
        else:
            log.info("Mask is considered not sliceable.")
        log.info("Therefore the mask will not be sliced.")
        return self._mask

    def _slice_wcs(self, item):
        if self._wcs is None:
            return None
        elif isinstance(self._wcs, np.ndarray):
            if self._wcs.shape == self.data.shape:
                return self._wcs[item]
            else:
                log.info("WCS has not the same shape as data.")
        elif isinstance(self._wcs, WCS):
            # There might be a problem if we don't slice all dimensions with
            # WCS ... e.g. 2D image but slice is [2]
            # Does it need a shape check as well?
            return self._wcs[item]
        else:
            log.info("WCS is considered not sliceable.")
        log.info("Therefore the WCS will not be sliced.")
        return self._wcs

    # Get the documentation in place.
    if isinstance(_slice_attr.__doc__, string_types):
        doc = """
        Controls how the {attr} is sliced.

        Should be called from ``__getitem__`` with `item` being the slice
        passed to it.

        Parameters
        ----------
        item: Slice
            The slice passed to ``__getitem__``.

        Returns
        -------
        new_{attr}: same type as ``self._{attr}`` or ``None``
            The sliced {attr}. If this function doesn't know how to slice
            the {attr} the return is just the saved ``_{attr}`` and if there is
            no {attr} the return will default to ``None``.

        Notes
        -----
        This function can treat {attr} if it is
        {allowed}.
        The primary purpose of this function is that subclasses using `NDData`
        with `NDSlicingMixin` with other {attr} restrictions need not to
        override the complete ``__getitem__`` but only need to override this
        method.
        """
        _slice_uncertainty.__doc__ = doc.format(attr="uncertainty",
            allowed="`~numpy.ndarray` or `NDUncertainty` with the same shape "
                    "as the data")
        _slice_mask.__doc__ = doc.format(attr="mask",
            allowed="`~numpy.ndarray` with the same shape as the data")
        _slice_wcs.__doc__ = doc.format(attr="WCS",
            allowed="`~numpy.ndarray` with the same shape as the data or a "
                    "`~astropy.wcs.WCS` object")
        del doc
