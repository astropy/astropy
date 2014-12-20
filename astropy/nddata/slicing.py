# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the Slicing mixin to the NDData class.

__all__ = ['NDSlicing']


class NDSlicing(object):
    """
    Mixin to provide slicing on objects using the NDData interface.

    When subclassing, be sure to list the superclasses in the correct order
    so that the subclass sees NDData as the main superclass. See
    `~astropy.nddata.NDDataArithmetic` for an example.
    """
    def __getitem__(self, item):

        new_data = self.data[item]

        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[item]
        else:
            new_uncertainty = None

        if self.mask is not None:
            new_mask = self.mask[item]
        else:
            new_mask = None

        if self.wcs is not None:
            new_wcs = self.wcs[item]
        else:
            new_wcs = None

        return self.__class__(new_data, uncertainty=new_uncertainty,
                              mask=new_mask, wcs=new_wcs,
                              meta=self.meta, unit=self.unit)
