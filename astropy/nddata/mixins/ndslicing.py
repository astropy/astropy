# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np

from ..nduncertainty import NDUncertainty
from ... import log
from ...wcs.wcs import WCS
from ...extern import six

__all__ = ['NDSlicingMixin']

# Need to skip doctest because of inconsistencies between python2 and 3
__doctest_skip__ = ['NDSlicingMixin']


# TODO: Delete this and rebase if #4242 is merged.
def tmp_deco(docstring, *args, **kwargs):
    def set_docstring(func):
        if not isinstance(docstring, six.string_types):
            doc = docstring.__doc__
        elif docstring != 'self':
            doc = docstring
        else:
            doc = func.__doc__
            func.__doc__ = None
        if not doc:
            raise ValueError
        kwargs['original_doc'] = func.__doc__ or ''
        func.__doc__ = doc.format(*args, **kwargs)
        return func
    return set_docstring

# TODO: Maybe not nice to pollute globals but I felt the same way about
# polluting the class namespace and this way potential subclasses may or refuse
# to pull this docstring into their class as well.

# Docstring templates for _slice_* methods.
_slice_docstring = """
    Controls how the ``{0}`` is sliced.

    Parameters
    ----------
    item: Slice
        The slice passed to :meth:`NDSlicingMixin.__getitem__`.

    Returns
    -------
    sliced_{0}: same type as ``self._{0}``
        The sliced {0}. If this function doesn't know how to slice
        the {0} the original {0} is returned and a warning is issued.

    Notes
    -----
    This function can slice the {0} if it is ``None`` or
    `~numpy.ndarray`-like if it's shape matches the data's shape or
    {allowed}.
    """
_slice_uncert_allowed = ("`~astropy.nddata.NDUncertainty`-like with it's "
                         "``array`` having the same shape as the data.")
_slice_wcs_allowed = "`~astropy.wcs.WCS` object"


class NDSlicingMixin(object):
    """
    Mixin to provide slicing on objects using the `~astropy.nddata.NDData`
    interface.

    The most common cases of the properties of `~astropy.nddata.NDData` are
    covered, i.e. if uncertainty, mask and wcs are ``None`` or
    `~numpy.ndarray` (see Notes for complete list of restrictions) but since
    NDData deliberatly enforces nothing on most of it's properties there might
    be a need to extend this. The `NDSlicingMixin` was designed to allow
    subclasses to extend just the portion they need changing without having
    to rewrite the whole Mixin (see Notes for detailed suggestions on how
    to do that).

    TODO: Move Notes to nddata/subclassing.rst except maybe point 2. (mwcraig)

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

       - data (which is enforced to be something `~numpy.ndarray`-like)
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

    You need to make NDData the last superclass like this::

        >>> from astropy.nddata import NDData, NDSlicingMixin
        >>> class NDDataSliceable(NDSlicingMixin, NDData): pass

    Then you can slice it::

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
        # Abort slicing if the data is a single scalar.
        if self.data.shape == ():
            raise TypeError('Scalars cannot be sliced.')

        # Slice the data here but everything else is sliced in self._slice.
        # TODO: Evaluate if any affiliated package that wants to use slicing
        # may have an interest in overwriting the way the data is sliced. But
        # I could not think of any possibilities except when slicing is based
        # on wcs (and for that there is the Cutout-Class in nddata.utils)
        new_data = self.data[item]
        kwargs = self._slice(item)
        return self.__class__(new_data, copy=False, **kwargs)

    def _slice(self, item):
        """
        Collects the sliced attributes and passes them back as ``kwargs``.

        It passes uncertainty, mask and wcs to their appropriate ``_slice_*``
        method, while ``meta`` and ``unit`` are simply taken from the original.
        The data is assumed to be sliceable and is sliced already before.

        When possible the return is *not* a copy of the data but a reference.

        Parameters
        ----------
        item: Slice
            The slice passed to :meth:`NDSlicingMixin.__getitem__`.

        Returns
        -------
        kwargs: `dict`
            Containing all the attributes except data after slicing - ready to
            feed them into the ``self.__class__.__init__()`` in __getitem__.
        """
        kwargs = {}
        # Try to slice some attributes
        kwargs['uncertainty'] = self._slice_uncertainty(item)
        kwargs['mask'] = self._slice_mask(item)
        kwargs['wcs'] = self._slice_wcs(item)
        # Attributes which are copied and not intended to be sliced
        kwargs['unit'] = self.unit
        kwargs['meta'] = self.meta
        return kwargs

    @tmp_deco(_slice_docstring, "uncertainty", allowed=_slice_uncert_allowed)
    def _slice_uncertainty(self, item):
        # TODO: Remove shape checks? (mwcraig)
        # But since the slice is defined by the datas shape and I did not
        # implement any setter/init shape checks and wouldn't like them
        # I haven't changed that (yet). I don't like shape checks because
        # they prohibit using broadcastable or scalar masks/uncertainties.
        if self.uncertainty is None:
            return None

        elif isinstance(self.uncertainty, np.ndarray):
            if self.uncertainty.shape == self.data.shape:
                return self.uncertainty[item]
            else:
                log.info("Uncertainty has not the same shape as data.")

        elif isinstance(self.uncertainty, NDUncertainty):
            if self.uncertainty.array.shape == self.data.shape:
                return self.uncertainty[item]
            else:
                log.info("Uncertainty has not the same shape as data.")

        else:
            log.info("Uncertainty is considered not sliceable.")
        log.info("Therefore the uncertainty will not be sliced.")
        return self.uncertainty

    @tmp_deco(_slice_docstring, "mask", allowed="")
    def _slice_mask(self, item):
        if self.mask is None:
            return None

        elif isinstance(self.mask, np.ndarray):
            if self.mask.shape == self.data.shape:
                return self.mask[item]
            else:
                log.info("Mask has not the same shape as data.")

        else:
            log.info("Mask is considered not sliceable.")
        log.info("Therefore the mask will not be sliced.")
        return self.mask

    @tmp_deco(_slice_docstring, 'wcs', allowed=_slice_wcs_allowed)
    def _slice_wcs(self, item):
        if self.wcs is None:
            return None

        elif isinstance(self.wcs, np.ndarray):
            if self.wcs.shape == self.data.shape:
                return self.wcs[item]
            else:
                log.info("wcs has not the same shape as data.")

        elif isinstance(self.wcs, WCS):
            # There might be a problem if we don't slice all dimensions with
            # WCS ... e.g. 2D image but slice is only 1D (like ndd[2]).
            # TODO: Does it need a shape check as well?
            return self.wcs[item]

        else:
            log.info("wcs is considered not sliceable.")
        log.info("Therefore the wcs will not be sliced.")
        return self.wcs
