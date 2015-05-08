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
The builtin frame classes are all imported automatically into this package's
namespace, so there's no need to access the sub-modules directly.

To implement a new frame in Astropy, a developer should add the frame as a new
module in this package.  Any "self" transformations (i.e., those that transform
from one frame to another frame of the same class) should be included in that
module.  Transformation functions connecting the new frame to other frames
should be in a separate module, which should be imported in this package's
``__init__.py`` to ensure the transformations are hooked up when this package is
imported.  Placing the trasnformation functions in separate modules avoids
circular dependencies, because they need references to the frame classes.
"""

from .icrs import ICRS
from .fk5 import FK5
from .fk4 import FK4, FK4NoETerms
from .galactic import Galactic
from .galactocentric import Galactocentric
from .altaz import AltAz
from .gcrs import GCRS, PrecessedGeocentric
from .cirs import CIRS
from .itrs import ITRS

#need to import transformations so that they get registered in the graph
from . import icrs_fk5_transforms
from . import fk4_fk5_transforms
from . import galactic_transforms
from . import icrs_cirs_transforms
from . import cirs_observed_transforms
from . import intermediate_rotation_transforms

# we define an __all__ because otherwise the transformation modules get included
__all__ = ['ICRS', 'FK5', 'FK4', 'FK4NoETerms', 'Galactic', 'Galactocentric',
           'AltAz', 'GCRS', 'CIRS', 'ITRS', 'PrecessedGeocentric']

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
