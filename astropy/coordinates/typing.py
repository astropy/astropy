# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Experimental typing support for :mod:`astropy.coordinates`, subject to change
without notice.
"""

__all__ = ["SupportsFrame"]

from typing import Protocol

from astropy.coordinates import BaseCoordinateFrame


class SupportsFrame(Protocol):
    """Protocol for classes that contain coordinate data.

    In :mod:`astropy.coordinates` the classes that implement this protocol
    are |SkyCoord| and :class:`~astropy.coordinates.BaseCoordinateFrame`
    (together with its subclasses).
    """

    @property
    def frame(self) -> BaseCoordinateFrame:
        """Coordinate data as a
        :class:`~astropy.coordinates.BaseCoordinateFrame` instance.
        """
