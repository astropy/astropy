# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

# Make sure that all tag implementations are imported by the time we create
# the extension class so that _astropy_asdf_types is populated correctly. We
# could do this using __init__ files, except it causes pytest import errors in
# the case that asdf is not installed.
from astropy.io.asdf.tags.fits.fits import *
from astropy.io.asdf.tags.table.table import *
from astropy.io.asdf.tags.time.time import *
from astropy.io.asdf.tags.transform.basic import *
from astropy.io.asdf.tags.transform.compound import *
from astropy.io.asdf.tags.transform.polynomial import *
from astropy.io.asdf.tags.transform.projections import *
from astropy.io.asdf.tags.transform.tabular import *
from astropy.io.asdf.tags.unit.quantity import *
from astropy.io.asdf.tags.unit.unit import *
from astropy.io.asdf.types import _astropy_asdf_types

from asdf.resolver import Resolver, DEFAULT_URL_MAPPING


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class AstropyExtension(object):
    @property
    def types(self):
        return _astropy_asdf_types

    @property
    def tag_mapping(self):
        return [('tag:stsci.edu:asdf',
                 'http://stsci.edu/schemas/asdf{tag_suffix}')]

    @property
    def url_mapping(self):
        return DEFAULT_URL_MAPPING
