# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for coordinates
"""
from . import frame_transform_graph
from ..wcs import WCSSUB_CELESTIAL

def find_wcs_frame(wcs):
    """
    Check the transform graph for any registered class that can represent the
    WCS object.

    Parameters
    ----------
    wcs: ~astropy.wcs.WCS object
        A WCS object to check for a frame match

    Returns
    -------
    frame: ~astropy.coordinates.BaseCoordinateFrame
        A unique frame match
    """

    # For now, we only care about the celestial portion of the WCS
    wcs = wcs.sub([WCSSUB_CELESTIAL])

    # Check that WCS has two dimensions
    if wcs.naxis != 2:
        raise ValueError("WCS should have two celestial dimensions")

    # TODO: don't use transform graph since not all frames are in the graph

    # Check the graph for a match
    matches = []
    for cls in frame_transform_graph.frame_set:
        try:
            matches.append(cls.from_wcs(wcs))
        except NotImplementedError:
            pass
        except ValueError:
            pass

    return matches
