# -*- coding: utf-8 -*-

from astropy.modeling.bounding_box import ModelBoundingBox, CompoundBoundingBox
from astropy.io.misc.asdf.types import AstropyAsdfType


__all__ = ["ModelBoundingBoxType", "CompoundBoundingBoxType"]


class ModelBoundingBoxType(AstropyAsdfType):
    name = 'transform/property/bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.ModelBoundingBox']

    @classmethod
    def to_tree(cls, bbox, cts):
        return {
            'intervals': {
                _input: list(interval) for _input, interval in bbox.intervals.items()
            },
            'ignore': list(bbox.ignored),
            'order': bbox.order
        }

    @classmethod
    def from_tree(cls, node, cts):
        intervals = {
            _input: tuple(interval) for _input, interval in node['intervals'].items()
        }

        if 'ignore' in node:
            ignored = node['ignore']
        else:
            ignored = None

        if 'order' in node:
            order = node['order']
        else:
            order = 'C'

        return ModelBoundingBox(intervals, ignored=ignored, order=order)

    @classmethod
    def assert_equal(cls, a, b):
        assert a.intervals == b.intervals
        assert a.ignored == b.ignored
        assert a.order == b.order


class CompoundBoundingBoxType(AstropyAsdfType):
    name = 'transform/property/compound_bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.CompoundBoundingBox']

    @classmethod
    def to_tree(cls, cbbox, cts):
        return {
            'selector_args': [
                {
                    'argument': sa[0],
                    'ignore': sa[1]
                } for sa in cbbox.selector_args
            ],
            'cbbox': [
                {
                    'key': list(key),
                    'bbox': bbox
                } for key, bbox in cbbox.bounding_boxes.items()
            ],
            'ignore': cbbox.global_ignored,
            'order': cbbox.order
        }

    @classmethod
    def from_tree(cls, node, cts):
        selector_args = tuple([
            (selector['argument'], selector['ignore']) for selector in node['selector_args']
        ])

        bounding_boxes = {
            tuple(bbox['key']): bbox['bbox']
            for bbox in node['cbbox']
        }

        if 'ignore' in node:
            ignored = node['ignore']
        else:
            ignored = None

        if 'order' in node:
            order = node['order']
        else:
            order = 'C'

        return CompoundBoundingBox(bounding_boxes, selector_args=selector_args, ignored=ignored, order=order)

    @classmethod
    def assert_equal(cls, a, b):
        assert a.bounding_boxes == b.bounding_boxes
        assert a.selector_args == b.selector_args
        assert a.global_ignored == b.global_ignored
        assert a.order == b.order
