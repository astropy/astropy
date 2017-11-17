# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf.yamlutil import custom_tree_to_tagged_tree

from astropy.units import Quantity
from astropy.coordinates import ICRS, Longitude, Latitude, Angle

from ...types import AstropyType
from ..unit.quantity import QuantityType


__all__ = ['ICRSCoordType']


class ICRSCoordType(AstropyType):
    name = "coords/icrs_coord"
    types = ['astropy.coordinates.ICRS']
    requires = ['astropy']
    version = "1.0.0"

    @classmethod
    def from_tree(cls, node, ctx):
        angle = QuantityType.from_tree(node['ra']['wrap_angle'], ctx)
        wrap_angle = Angle(angle.value, unit=angle.unit)
        ra = Longitude(
            node['ra']['value'],
            unit=node['ra']['unit'],
            wrap_angle=wrap_angle)
        dec = Latitude(node['dec']['value'], unit=node['dec']['unit'])

        return ICRS(ra=ra, dec=dec)

    @classmethod
    def to_tree(cls, frame, ctx):
        node = {}

        wrap_angle = Quantity(
            frame.ra.wrap_angle.value,
            unit=frame.ra.wrap_angle.unit)
        node['ra'] = {
            'value': frame.ra.value,
            'unit': frame.ra.unit.to_string(),
            'wrap_angle': custom_tree_to_tagged_tree(wrap_angle, ctx)
        }
        node['dec'] = {
            'value': frame.dec.value,
            'unit': frame.dec.unit.to_string()
        }

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert new.ra == old.ra and new.dec == old.dec
