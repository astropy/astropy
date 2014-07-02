# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for coordinates
"""
from .transfomations import TransformGraph

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

    # Add WCS validation here
    
    # Check the graph for a match
    matches = []
    for cls in TransformGraph.frame_set:
        if cls.frame_wcs(wcs):
            matches.append(cls)

    if not len(matches):
        raise Exception("No matches found for this WCS object")
    
    elif len(matches) > 1:
        raise Exception("This WCS object is ambiguious.")
    
    return matches[0]