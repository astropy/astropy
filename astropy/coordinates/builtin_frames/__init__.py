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

from . import baseradec
from .baseradec import BaseRADecFrame  # noqa
from . import icrs
from .icrs import ICRS  # noqa
from . import fk5
from .fk5 import FK5  # noqa
from . import fk4
from .fk4 import FK4, FK4NoETerms  # noqa
from . import galactic
from .galactic import Galactic  # noqa
from . import galactocentric
from .galactocentric import Galactocentric, galactocentric_frame_defaults  # noqa
from . import supergalactic
from .supergalactic import Supergalactic  # noqa
from . import altaz
from .altaz import AltAz  # noqa
from . import hadec
from .hadec import HADec  # noqa
from . import gcrs
from .gcrs import GCRS, PrecessedGeocentric  # noqa
from . import cirs
from .cirs import CIRS  # noqa
from . import itrs
from .itrs import ITRS  # noqa
from . import hcrs
from .hcrs import HCRS  # noqa
from . import equatorial
from .equatorial import TEME, TETE  # noqa

from . import ecliptic
from .ecliptic import (GeocentricMeanEcliptic, BarycentricMeanEcliptic, HeliocentricMeanEcliptic,   # noqa
                       BaseEclipticFrame, GeocentricTrueEcliptic, BarycentricTrueEcliptic,  # noqa
                       HeliocentricTrueEcliptic, HeliocentricEclipticIAU76, CustomBarycentricEcliptic)  # noqa
from . import skyoffset
from .skyoffset import SkyOffsetFrame  # noqa
# need to import transformations so that they get registered in the graph
from . import icrs_fk5_transforms  # noqa
from . import fk4_fk5_transforms  # noqa
from . import galactic_transforms  # noqa
from . import supergalactic_transforms  # noqa
from . import icrs_cirs_transforms  # noqa
from . import cirs_observed_transforms  # noqa
from . import icrs_observed_transforms  # noqa
from . import intermediate_rotation_transforms  # noqa
from . import ecliptic_transforms  # noqa

# Import this after importing other frames, since this requires various
# transformations to set up the LSR frames
from . import lsr
from .lsr import LSR, GalacticLSR, LSRK, LSRD  # noqa

from astropy.coordinates.baseframe import frame_transform_graph

# we define an __all__ because otherwise the transformation modules
# get included
__all__ = []
__all__.extend(icrs.__all__)
__all__.extend(fk5.__all__)
__all__.extend(fk4.__all__)
__all__.extend(galactic.__all__)
__all__.extend(galactocentric.__all__)
__all__.extend(supergalactic.__all__)
__all__.extend(altaz.__all__)
__all__.extend(hadec.__all__)
__all__.extend(gcrs.__all__)
__all__.extend(cirs.__all__)
__all__.extend(itrs.__all__)
__all__.extend(hcrs.__all__)
__all__.extend(equatorial.__all__)
__all__.extend(ecliptic.__all__)
__all__.extend(skyoffset.__all__)
__all__.extend(lsr.__all__)
__all__.extend(baseradec.__all__)
__all__.append('make_transform_graph_docs')


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
        block = f"""
            <li style='list-style: none;'>
                <p style="font-size: 12px;line-height: 24px;font-weight: normal;color: #848484;padding: 0;margin: 0;">
                    <b>{cls.__name__}:</b>
                    <span style="font-size: 24px; color: {color};"><b>➝</b></span>
                </p>
            </li>
        """
        html_list_items.append(block)

    nl = '\n'
    graph_legend = f"""
    .. raw:: html

        <ul>
            {nl.join(html_list_items)}
        </ul>
    """
    docstr = docstr + dedent(graph_legend)

    return docstr


_transform_graph_docs = make_transform_graph_docs(frame_transform_graph)

# Here, we override the module docstring so that sphinx renders the transform
# graph without the developer documentation in the main docstring above.
__doc__ = _transform_graph_docs
