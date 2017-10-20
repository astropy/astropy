# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import os

from asdf.resolver import Resolver, DEFAULT_URL_MAPPING

from astropy.io.asdf.tags.fits import FitsType
from astropy.io.asdf.tags.time import TimeType
from astropy.io.asdf.tags.unit import UnitType, QuantityType
from astropy.io.asdf.tags.table import TableType, ColumnType
from astropy.io.asdf.tags.transform import (
    TransformType, IdentityType, ConstantType, DomainType, CompoundType,
    RemapAxesType, ShiftType, ScaleType, PolynomialType, AffineType,
    Rotate2DType, Rotate3DType, TabularType)
from astropy.io.asdf.tags.transform.projections import _projection_types


SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'schemas'))


class AstropyExtension(object):
    @property
    def types(self):
        # TODO: This could be simplified by an AstropyType that inherits from
        # asdf.asdftypes.CustomType and automatically adds each subclass to
        # a list (much like AsdfType does internally).
        return [
            FitsType,
            TimeType,
            UnitType,
            QuantityType,
            TableType,
            ColumnType,
            TransformType,
            IdentityType,
            ConstantType,
            DomainType,
            CompoundType,
            RemapAxesType,
            ShiftType,
            ScaleType,
            PolynomialType,
            AffineType,
            Rotate2DType,
            Rotate3DType,
            TabularType
        ] + _projection_types

    @property
    def tag_mapping(self):
        return [('tag:stsci.edu:asdf',
                 'http://stsci.edu/schemas/asdf{tag_suffix}')]

    @property
    def url_mapping(self):
        return DEFAULT_URL_MAPPING
