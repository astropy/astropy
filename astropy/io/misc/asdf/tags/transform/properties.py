# -*- coding: utf-8 -*-

from astropy.modeling.bounding_box import ModelBoundingBox
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
            raise TypeError(f"{bbox} is not a valid ModelBoundingBox")

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
