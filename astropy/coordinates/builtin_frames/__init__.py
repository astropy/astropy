# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains the coordinate frames implemented by astropy.

Users shouldn't use this module directly, but rather import from the
`astropy.coordinates` module.  While it is likely to exist for the long-term,
the existence of this package and details of its organization should be
considered an implementation detail, and is not guaranteed to hold for future
versions of astropy.

Notes
-----
The builtin frame classes are all imported automatically into this package's
namespace, so there's no need to access the sub-modules directly.

To implement a new frame in Astropy, a developer should add the frame as a new
module in this package.  Any "self" transformations (i.e., those that transform
from one frame to another frame of the same class) should be included in that
module.  Transformation functions connecting the new frame to other frames
should be in a separate module, which should be imported in this package's
``__init__.py`` to ensure the transformations are hooked up when this package is
imported.  Placing the transformation functions in separate modules avoids
circular dependencies, because they need references to the frame classes.
"""

from .baseradec import BaseRADecFrame
from .icrs import ICRS
from .fk5 import FK5
from .fk4 import FK4, FK4NoETerms
from .galactic import Galactic
from .galactocentric import Galactocentric, galactocentric_frame_defaults
from .supergalactic import Supergalactic
from .altaz import AltAz
from .gcrs import GCRS, PrecessedGeocentric
from .cirs import CIRS
from .itrs import ITRS
from .hcrs import HCRS
from .equatorial import TEME, TETE

from .ecliptic import *  # there are a lot of these so we don't list them all explicitly
from .skyoffset import SkyOffsetFrame
# need to import transformations so that they get registered in the graph
from . import icrs_fk5_transforms
from . import fk4_fk5_transforms
from . import galactic_transforms
from . import supergalactic_transforms
from . import icrs_cirs_transforms
from . import cirs_observed_transforms
from . import intermediate_rotation_transforms
from . import ecliptic_transforms

# Import this after importing other frames, since this requires various
# transformtions to set up the LSR frames
from .lsr import LSR, GalacticLSR, LSRK, LSRD

from astropy.coordinates.baseframe import frame_transform_graph

# we define an __all__ because otherwise the transformation modules
# get included
__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'Galactocentric',
           'galactocentric_frame_defaults',
           'Supergalactic', 'AltAz', 'GCRS', 'CIRS', 'ITRS', 'HCRS', 'TEME',
           'TETE', 'PrecessedGeocentric', 'GeocentricMeanEcliptic',
           'BarycentricMeanEcliptic', 'HeliocentricMeanEcliptic',
           'GeocentricTrueEcliptic', 'BarycentricTrueEcliptic',
           'HeliocentricTrueEcliptic',
           'SkyOffsetFrame', 'GalacticLSR', 'LSR', 'LSRK', 'LSRD',
           'BaseEclipticFrame', 'BaseRADecFrame', 'make_transform_graph_docs',
           'HeliocentricEclipticIAU76', 'CustomBarycentricEcliptic']


def make_transform_graph_docs(transform_graph):
    """
    Generates a string that can be used in other docstrings to include a
    transformation graph, showing the available transforms and
    coordinate systems.

    Parameters
    ----------
    transform_graph : `~.coordinates.TransformGraph`

    Returns
    -------
    docstring : str
        A string that can be added to the end of a docstring to show the
        transform graph.
    """
    from textwrap import dedent
    coosys = [transform_graph.lookup_name(item) for
              item in transform_graph.get_names()]

    # currently, all of the priorities are set to 1, so we don't need to show
    #   then in the transform graph.
    graphstr = transform_graph.to_dot_graph(addnodes=coosys,
                                            priorities=False)

    docstr = """
    The diagram below shows all of the built in coordinate systems,
    their aliases (useful for converting other coordinates to them using
    attribute-style access) and the pre-defined transformations between
    them.  The user is free to override any of these transformations by
    defining new transformations between these systems, but the
    pre-defined transformations should be sufficient for typical usage.

    The color of an edge in the graph (i.e. the transformations between two
    frames) is set by the type of transformation; the legend box defines the
    mapping from transform class name to color.

    .. Wrap the graph in a div with a custom class to allow themeing.
    .. container:: frametransformgraph

        .. graphviz::

    """

    docstr = dedent(docstr) + '        ' + graphstr.replace('\n', '\n        ')

    # colors are in dictionary at the bottom of transformations.py
    from astropy.coordinates.transformations import trans_to_color
    html_list_items = []
    for cls, color in trans_to_color.items():
        block = """
            <li style='list-style: none;'>
                <p style="font-size: 12px;line-height: 24px;font-weight: normal;color: #848484;padding: 0;margin: 0;">
                    <b>{}:</b>
                    <span style="font-size: 24px; color: {};"><b>➝</b></span>
                </p>
            </li>
        """.format(cls.__name__, color)
        html_list_items.append(block)

    graph_legend = """
    .. raw:: html

        <ul>
            {}
        </ul>
    """.format("\n".join(html_list_items))
    docstr = docstr + dedent(graph_legend)

    return docstr


_transform_graph_docs = make_transform_graph_docs(frame_transform_graph)

# Here, we override the module docstring so that sphinx renders the transform
# graph without the developer documentation in the main docstring above.
__doc__ = _transform_graph_docs
