# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ... import log

__all__ = ['NDSlicingMixin']


class NDSlicingMixin(object):
    """
    Mixin to provide slicing on objects using the `~astropy.nddata.NDData`
    interface.

    The most common cases of the properties of `~astropy.nddata.NDData` are
    covered, i.e. if uncertainty, mask and wcs are ``None`` or unsliceable.

    TODO: Move Notes to nddata/subclassing.rst except maybe point 2. (mwcraig)

    Notes
    -----
    1. When subclassing, be sure to list the superclasses in the correct order
       so that the subclass sees NDData as the main superclass. See
       `~astropy.nddata.NDDataArray` for an example.

    2. The ``meta`` and ``unit`` are simply taken from the original, no attempt
       is made to slice them.

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
       - if you have a property that does not match the criterias above you may
         need to extend or override some method here. For example you want a
         custom made ``uncertainty`` class that does not subclass from
         `NDUncertainty` or `~numpy.ndarray` to be sliced. The best way to do
         this would be to extend ``_slice_uncertainty(self, item)`` which
         takes the ``slice`` as ``item`` parameter and returns what the
         sliced uncertainty should be. The current uncertainty is avaiable
         using ``self._uncertainty``. Be careful since this can be also
         ``None`` if there is no uncertainty set.

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

    """
    def __getitem__(self, item):
        # Abort slicing if the data is a single scalar.
        if self.data.shape == ():
            raise TypeError('scalars cannot be sliced.')

        # Let the other methods handle slicing.
        kwargs = self._slice(item)
        return self.__class__(**kwargs)

    def _slice(self, item):
        """
        Collects the sliced attributes and passes them back as ``kwargs``.

        It passes uncertainty, mask and wcs to their appropriate ``_slice_*``
        method, while ``meta`` and ``unit`` are simply taken from the original.
        The data is assumed to be sliceable and is sliced directly.

        When possible the return is *not* a copy of the data but a reference.

        Parameters
        ----------
        item: slice
            The slice passed to ``__getitem__``.

        Returns
        -------
        kwargs: `dict`
            Containing all the attributes after slicing - ready to
            feed them into the ``self.__class__.__init__()`` in __getitem__.
        """
        kwargs = {}
        kwargs['data'] = self.data[item]
        # Try to slice some attributes
        kwargs['uncertainty'] = self._slice_uncertainty(item)
        kwargs['mask'] = self._slice_mask(item)
        kwargs['wcs'] = self._slice_wcs(item)
        # Attributes which are copied and not intended to be sliced
        kwargs['unit'] = self.unit
        kwargs['meta'] = self.meta
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

    def _slice_wcs(self, item):
        if self.wcs is None:
            return None
        try:
            return self.wcs[item]
        except TypeError:
            log.info("wcs cannot be sliced.")
        return self.wcs
