# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains a set of compatibility classes to allow use of
the pre-v0.3 coordinate names.  It will be removed in a future version.
"""
#TODO: remove this module in a future version

from .builtin_systems import *
from .transformations import master_transform_graph

__all__ = ['ICRSCoordinates', 'FK5Coordinates', 'FK4Coordinates',
           'FK4NoETermCoordinates', 'GalacticCoordinates', 'HorizontalCoordinates'
          ]

class ICRSCoordinates(ICRS):
    """
    Using the `ICRSCoordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `ICRS` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return ICRS(*args, **kwargs)


class FK5Coordinates(FK5):
    """
    Using the `FK5Coordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `FK5` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return FK5(*args, **kwargs)


class FK4Coordinates(FK4):
    """
    Using the `FK4Coordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `FK4` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return FK4(*args, **kwargs)


class FK4NoETermCoordinates(FK4NoETerms):
    """
    Using the `FK4NoETermCoordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `FK4NoETerms` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return FK4NoETerms(*args, **kwargs)


class GalacticCoordinates(Galactic):
    """
    Using the `GalacticCoordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `Galactic` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return Galactic(*args, **kwargs)


class HorizontalCoordinates(AltAz):
    """
    Using the `HorizontalCoordinates` name for this class is deprecated in v0.3, and will be
    removed in the next version. Use `AltAz` instead.
    """
    def __new__(cls, *args, **kwargs):
        from warnings import warn
        from ..utils.exceptions import AstropyBackwardsIncompatibleChangeWarning

        wmsg = cls.__doc__.replace('\n    ', ' ').strip()
        warn(AstropyBackwardsIncompatibleChangeWarning(wmsg))
        return AltAz(*args, **kwargs)


def _add_transforms(clses, graph):
    """
    Adds fake transformations that allow transforming to the old names, although
    they actually yield the new class types
    """
    from copy import deepcopy

    for cls in clses:
        newcls = cls.mro()[1]
        gdct = graph._graph

        toadd = []
        for a in gdct:
            for b in gdct[a]:
                if b == newcls:
                    toadd.append((a, cls, gdct[a][b]))
        for a, b, t in toadd:
            #adds a new transform that goes *to* the old name class
            graph.add_transform(a, b, t)

        #also add a transformation that just gives itself back to go from
        #self to old-style name
        graph.add_transform(newcls, cls, lambda c:deepcopy(c))

#Now go through and add transforms so that the old names give you transforms to new things
_add_transforms([globals()[nm] for nm in __all__], master_transform_graph)
