# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from asdf.yamlutil import custom_tree_to_tagged_tree

from astropy.coordinates.baseframe import frame_transform_graph
from astropy.tests.helper import assert_quantity_allclose

from ...types import AstropyType


__all__ = ['CoordType']


class CoordType(AstropyType):
    name = "coords/coord"
    types = ['astropy.coordinates.BaseCoordinateFrame']
    requires = ['astropy']
    version = "1.0.0"

    @classmethod
    def from_tree(cls, node, ctx):
        frame = frame_transform_graph.lookup_name(node['frame'])
        data = node.get('data', None)
        if data:
            return frame(node['data'], **node['frame_attributes'])

        return frame(**node['frame_attributes'])

    @classmethod
    def to_tree(cls, frame, ctx):
        if type(frame) not in frame_transform_graph.frame_set:
            raise ValueError("Can only save frames that are registered with the "
                             "transformation graph.")

        node = {}
        node['frame'] = frame.name
        if frame.has_data:
            node['data'] = custom_tree_to_tagged_tree(frame.data, ctx)
        frame_attributes = {}
        for attr in frame.frame_attributes.keys():
            value = getattr(frame, attr, None)
            if value:
                frame_attributes[attr] = value
        node['frame_attributes'] = custom_tree_to_tagged_tree(frame_attributes, ctx)

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(new, type(old))
        if new.has_data:
            assert_quantity_allclose(new.data.lon, old.data.lon)
            assert_quantity_allclose(new.data.lat, old.data.lat)
