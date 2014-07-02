# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from .errors import *
from .angles import *
from .baseframe import *
from .distances import *
from .earth import *
from .transformations import *
from .builtin_frames import *
from .name_resolve import *
from .matching import *
from .representation import *
from .sky_coordinate import *

__doc__ += builtin_frames._transform_graph_docs


# This probably shouldn't live here
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