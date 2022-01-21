# -*- coding: utf-8 -*-

from astropy.modeling.bounding_box import ModelBoundingBox, CompoundBoundingBox
from astropy.io.misc.asdf.types import AstropyAsdfType


class ModelBoundingBoxType(AstropyAsdfType):
    name = 'transform/property/bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.ModelBoundingBox']

    @classmethod
    def to_tree(cls, bbox, cts):
        if isinstance(bbox, ModelBoundingBox):
            return {
                'intervals': {
                    _input: list(interval) for _input, interval in bbox.intervals.items()
                },
                'ignore': list(bbox.ignored),
                'order': bbox.order
            }
        else:
            raise TypeError(f"{bbox} is not a valid ModelBoundingBox!")

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
        if isinstance(cbbox, CompoundBoundingBox):
            return {
                'selector_args': [
                    {
                        'argument': sa.name(cbbox._model),
                        'ignore': sa.ignore
                    } for sa in cbbox.selector_args
                ],
                'cbbox': [
                    {
                        'key': list(key),
                        'bbox': bbox
                    } for key, bbox in cbbox.bounding_boxes.items()
                ],
                'ignore': cbbox.ignored_inputs,
                'order': cbbox.order
            }
        else:
            raise TypeError(f"{cbbox} is not a valid CompoundBoundingBox!")

    @classmethod
    def from_tree(cls, node, cts):
        selector_args = [
            (selector['argument'], selector['ignore']) for selector in node['selector_args']
        ]

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
        assert a.ignored == b.ignored
        assert a.order == b.order
