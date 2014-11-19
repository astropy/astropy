# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains the coordinate frames actually implemented by astropy.
"""

from .icrs import ICRS
from .fk5 import FK5
from .fk4 import FK4, FK4NoETerms
from .galactic import Galactic
from .altaz import AltAz


def _make_transform_graph_docs():
    """
    Generates a string for use with the coordinate package's docstring
    to show the available transforms and coordinate systems
    """
    import inspect
    from textwrap import dedent
    from ...extern import six
    from ..baseframe import BaseCoordinateFrame, frame_transform_graph

    isclass = inspect.isclass
    coosys = [item for item in list(six.itervalues(globals()))
              if isclass(item) and issubclass(item, BaseCoordinateFrame)]
    graphstr = frame_transform_graph.to_dot_graph(addnodes=coosys)

    docstr = """
    The diagram below shows all of the coordinate systems built into the
    `~astropy.coordinates` package, their aliases (useful for converting
    other coordinates to them using attribute-style access) and the
    pre-defined transformations between them.  The user is free to
    override any of these transformations by defining new transformations
    between these systems, but the pre-defined transformations should be
    sufficient for typical usage.

    The graph also indicates the priority for each transformation as a
    number next to the arrow.  These priorities are used to decide the
    preferred order when two transformation paths have the same number
    of steps.  These priorities are defined such that the path with a
    *smaller* total priority is favored.


    .. graphviz::

    """

    return dedent(docstr) + '    ' + graphstr.replace('\n', '\n    ')
_transform_graph_docs = _make_transform_graph_docs()
