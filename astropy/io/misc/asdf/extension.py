# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

from asdf.extension import BuiltinExtension
from asdf.resolver import Resolver, DEFAULT_URL_MAPPING

# Make sure that all tag implementations are imported by the time we create
# the extension class so that _astropy_asdf_types is populated correctly. We
# could do this using __init__ files, except it causes pytest import errors in
# the case that asdf is not installed.
from .tags.fits.fits import *
from .tags.table.table import *
from .tags.time.time import *
from .tags.transform.basic import *
from .tags.transform.compound import *
from .tags.transform.polynomial import *
from .tags.transform.projections import *
from .tags.transform.tabular import *
from .tags.unit.quantity import *
from .tags.unit.unit import *
from .types import _astropy_asdf_types


__all__ = ['AstropyAsdfExtension']


# This extension is used to register custom tag types that have schemas defined
# by ASDF, but have tag implementations defined in astropy.
class AstropyAsdfExtension(BuiltinExtension):
    @property
    def types(self):
        return _astropy_asdf_types
