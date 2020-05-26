# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.decorators import format_doc
from astropy.coordinates.attributes import TimeAttribute
from .utils import DEFAULT_OBSTIME
from astropy.coordinates.baseframe import base_doc
from .baseradec import BaseRADecFrame, doc_components

__all__ = ['BCRS']


doc_footer = """
    Other parameters
    ----------------
    obstime : `~astropy.time.Time`
        The time at which the observation is taken.
"""


@format_doc(base_doc, components=doc_components, footer=doc_footer)
class BCRS(BaseRADecFrame):
    """
    """

    obstime = TimeAttribute(default=DEFAULT_OBSTIME)

# Transformations are defined in icrs_cirs_transforms.py
