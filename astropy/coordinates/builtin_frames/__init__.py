# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package contains the coordinate frames actually implemented by astropy.

Users shouldn't use this module directly, but rather import from the
`astropy.coordinates` module.  While it is likely to exist for the long-term,
the existence of this package and details of its organization should be
considered an implementation detail, and is not guaranteed to hold for future
versions of astropy.

Notes
-----
The builtin coordinates classes are all imported automatically into this package
namespace, so there's no need to access the sub-modules directly.

To implement a new frame in Astropy, a developer should add the frame as a new
module in this package.  The functions necessary to transform to/from that frame
should generally also be placed in that frame's module.  In some cases putting both in the same place will
be impossible because transformation functions need to have references to the
class objects they transform to and from; this leads to potential circular
dependency problems.  So when implementing multiple multiple new frames with
transformations that are intermixed, simply choose one to include both the "to"
and "from" functions, with the guideline that the transformations should go in
whichever frame is the more specific (i.e., less likely to be used by other
future frames.)
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
