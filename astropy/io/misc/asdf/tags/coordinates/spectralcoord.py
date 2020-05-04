# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from astropy.coordinates.spectral_coordinate import SpectralCoord
from astropy.io.misc.asdf.types import AstropyType
from astropy.io.misc.asdf.tags.unit.unit import UnitType

try:
    from asdf.tags.core import NDArrayType
except ImportError:
    HAS_ASDF = False
else:
    HAS_ASDF = True

__all__ = ['SpectralCoordType']


class SpectralCoordType(AstropyType):
    """
    ASDF tag implementation used to serialize/derialize SpectralCoord objects
    """
    name = 'coordinates/spectralcoord'
    types = [SpectralCoord]
    version = '1.0.0'

    @classmethod
    def to_tree(cls, spec_coord, ctx):
        node = {}
        if isinstance(spec_coord, SpectralCoord):
            node['value'] = spec_coord.value
            node['unit'] = spec_coord.unit
            node['observer'] = spec_coord.observer
            node['target'] = spec_coord.target
            return node
        raise TypeError(f"'{spec_coord}' is not a valid SpectralCoord")

    @classmethod
    def from_tree(cls, node, ctx):
        if not HAS_ASDF:
            raise ImportError('asdf is not installed')

        if isinstance(node, SpectralCoord):
            return node

        unit = UnitType.from_tree(node['unit'], ctx)
        value = node['value']
        observer = node['observer'] if 'observer' in node else None
        target = node['target'] if 'observer' in node else None
        if isinstance(value, NDArrayType):
            value = value._make_array()
        return SpectralCoord(value, unit=unit, observer=observer, target=target)
