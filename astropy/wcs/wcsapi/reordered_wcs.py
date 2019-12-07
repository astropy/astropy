from .sliced_low_level_wcs import SlicedLowLevelWCS

__all__ = ['ReorderedLowLevelWCS']


class ReorderedLowLevelWCS(SlicedLowLevelWCS):
    """
    A wrapper for a low-level WCS object that has re-ordered
    pixel and/or world axes.

    Parameters
    ----------
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The original WCS for which to reorder axes
    pixel_keep : iterable
        The indices of the original axes in the order of the
        new WCS.
    world_keep : iterable
        The indices of the original axes in the order of the
        new WCS.
    """
    def __init__(self, wcs, pixel_keep, world_keep):
        super().__init__(wcs, Ellipsis)
        self._pixel_keep = pixel_keep
        self._world_keep = world_keep

    @property
    def array_shape(self):
        return [self._wcs.array_shape[i] for i in self._pixel_keep[::-1]]
