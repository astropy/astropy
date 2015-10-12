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

# Need to skip doctest because of inconsistencies between python2 and 3
__doctest_skip__ = ['NDSlicingMixin']

class NDSlicingMixin(object):
    """
    Mixin to provide slicing on objects using the NDData interface.

    The most common cases of the properties of `~astropy.nddata.NDData` are
    covered, i.e. if uncertainty, mask and wcs are ``None`` or
    `~numpy.ndarray` (see Notes for complete list of restrictions) but since
    NDData deliberatly enforces nothing on most of it's properties there might
    be a need to extend this. The `NDSlicingMixin` was designed to allow
    subclasses to extend just the portion they need changing without having
    to rewrite the whole Mixin (see Notes for detailed suggestions on how
    to do that).

    Notes
    -----
    1. When subclassing, be sure to list the superclasses in the correct order
       so that the subclass sees NDData as the main superclass. See
       `~astropy.nddata.NDDataArray` for an example.

    2. Since `~astropy.nddata.NDData` is not very restrictive about what the
       additional attributes can be there are a lot of tests while slicing.
       If a subclass defines more restrictive setter methods for the properties
       they might want to consider altering `NDSlicingMixin` methods for
       slicing (starting with _slice_* e.g.
       `NDSlicingMixin._slice_uncertainty`). Currently implemented is the
       for slicing the:

       - data (which is enforced to be something like a `~numpy.ndarray`)
       - mask (if it is a numpy array)
       - wcs (if it is a numpy array or `~astropy.wcs.WCS` object)
       - uncertainty (if it is a numpy array or a subclass of `NDUncertainty`)

       But only if the shapes match the shape of the data (except for the wcs
       if it is a ``WCS`` object). If these attributes do not match those
       criterias the original property is used but an ``INFO`` message is
       displayed. The ``meta`` and ``unit`` are simply taken from the original.

    3. Some advice about extending this Mixin in subclasses:

       - if the ``data`` is not a `~numpy.ndarray` or should not be sliceable
         you may be better off not using this Mixin.
       - if you have defined an additional property that needs to be sliced
         extend :meth:`NDSlicingMixin._slice`. Probably the best way to do this
         is to first call ``kwarg = super(subclass, self)._slice(item)`` in
         there and then just afterwards add the other sliced properties to the
         ``kwarg``-dict. Since this kwarg is used to initialize a new
         instance of the class you need to match the key to the name of the
         parameter. For example if you use a property called ``flags`` you need
         to add a new line ``kwargs['flags'] = ...`` to :meth:`_slice`.
         *BUT* the ``__init__`` has to accept a flags parameter otherwise
         this will fail.
       - if you have restrictions that are a subset of the above mentioned
         conditions on the properties do not override or extend the
         ``_slice_*`` methods. The speed gain is negligible.
       - if you have a property that does not match the criterias above you may
         need to extend or override some method here. For example you want a
         custom made ``uncertainty`` class that does not subclass from
         `NDUncertainty` or `~numpy.ndarray` to be sliced. The best way to do
         this would be to extend ``_slice_uncertainty(self, item)`` which
         takes the ``slice`` as ``item`` parameter and returns what the
         sliced uncertainty should be. The current uncertainty is avaiable
         using ``self._uncertainty``. Be carful since this can be also ``None``
         if there is no uncertainty set.

    Examples
    --------

        >>> from astropy.nddata import NDData, NDSlicingMixin
        >>> # You need to make NDData the last superclass like this:
        >>> class NDDataSliceable(NDSlicingMixin, NDData): pass
        >>> nd = NDDataSliceable([1,2,3,4,5])
        >>> nd[1:3]
        NDDataSliceable([2, 3])
        >>> import numpy as np # Mask has to be a numpy array to be sliceable
        >>> mask = np.array([True, False, True, True, False])
        >>> nd2 = NDDataSliceable(nd, mask=mask)
        >>> nd2[1:3].mask
        array([False,  True], dtype=bool)
        >>> # But be aware that the slicing only returns references not copies
        >>> nd3 = nd2[1:3]
        >>> nd3.data[0] = 100
        >>> nd2
        NDDataSliceable([  1, 100,   3,   4,   5])
        >>> # If the mask is not a numpy array a warning is issued
        >>> nd4 = NDDataSliceable(nd, mask=[False, True])
        >>> nd5 = nd4[0:2] # doctest: +SKIP
        INFO: Mask is considered not sliceable. [astropy.nddata.mixins.ndslicing]
        INFO: Therefore the mask will not be sliced. [astropy.nddata.mixins.ndslicing]
        >>> # but the data is sliced
        >>> nd5 # doctest: +SKIP
        NDDataSliceable([  1, 100])
        >>> # and the mask is just the mask without slicing
        >>> nd5.mask # doctest: +SKIP
        [False, True]

    """
    def __getitem__(self, item):
        # Abort slicing if the data is scalar
        if self._data.shape == ():
            raise TypeError('Scalars cannot be sliced.')

        # Slice the data but everything else is sliced in self._slice.
        new_data = self.data[item]
        kwargs = self._slice(item)
        return self.__class__(new_data, **kwargs)

    def _slice(self, item):
        """
        Collects the sliced attributes and passes them back as ``kwarg`` for
        the return after slicing.

        It passes uncertainty, mask and wcs to their appropriate ``_slice_*``
        method, while ``meta`` and ``unit`` are simply taken from the original.
        The data is assumed to be sliceable and is sliced already before.

        Where possible the return is *not* a copy of the data.

        Parameters
        ----------
        item: Slice
            The slice passed to :meth:`NDSlicingMixin.__getitem__`.

        Returns
        -------
        kwargs: `dict`
            Containing all the attributes except data after slicing - ready to
            feed them into the ``__init__()``.
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
    if isinstance(_slice.__doc__, string_types):
        doc = """
        Controls how the {attr} is sliced.

        Should be called from :meth:`NDSlicingMixin.__getitem__` with ``item``
        being the slice passed to it.

        Parameters
        ----------
        item: Slice
            The slice passed to :meth:`NDSlicingMixin.__getitem__`.

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
        override the complete :meth:`NDSlicingMixin.__getitem__` but only need
        to override this method.
        """
        _slice_uncertainty.__doc__ = doc.format(attr="uncertainty",
            allowed="`~numpy.ndarray` or `NDUncertainty` with the same shape "
                    "as the data")
        _slice_mask.__doc__ = doc.format(attr="mask",
            allowed="`~numpy.ndarray` with the same shape as the data")
        _slice_wcs.__doc__ = doc.format(attr="wcs",
            allowed="`~numpy.ndarray` with the same shape as the data or a "
                    "`~astropy.wcs.WCS` object")
        del doc
