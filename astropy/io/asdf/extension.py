# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

from asdf.resolver import Resolver, DEFAULT_URL_MAPPING

from astropy.io.asdf.tags.time import TimeType


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class AstropyExtension(object):
    @property
    def types(self):
        return [TimeType]

    @property
    def tag_mapping(self):
        return [('tag:stsci.edu:asdf',
                 'http://stsci.edu/schemas/asdf{tag_suffix}')]

    @property
    def url_mapping(self):
        return DEFAULT_URL_MAPPING
