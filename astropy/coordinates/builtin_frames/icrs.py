# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates.baseframe import base_doc
from .baseradec import BaseRADecFrame, doc_components

__all__ = ['ICRS']


@format_doc(base_doc, components=doc_components, footer="")
class ICRS(BaseRADecFrame):
    """
    A coordinate or frame in the ICRS system.

    If you're looking for "J2000" coordinates, and aren't sure if you want to
    use this or `~astropy.coordinates.FK5`, you probably want to use ICRS. It's
    more well-defined as a catalog coordinate and is an inertial system, and is
    very close (within tens of milliarcseconds) to J2000 equatorial.

    For more background on the ICRS and related coordinate transformations, see
    the references provided in the  :ref:`astropy-coordinates-seealso` section
    of the documentation.
    """
