# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
from astropy.modeling.bounding_box import ModelBoundingBox, CompoundBoundingBox
from astropy.io.misc.asdf.types import AstropyAsdfType


class BoundingBoxType(AstropyAsdfType):
    name = 'transform/property/bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.ModelBoundingBox']

    @classmethod
    def to_tree(cls, bbox, cts):
        if isinstance(bbox, ModelBoundingBox):
            return {
                'intervals': {
                    _input: list(interval)
                    for _input, interval in bbox.named_intervals.items()
                },
                'ignore': list(bbox.ignored_inputs),
                'order': bbox.order
            }
        else:
            raise TypeError(f"{bbox} is not a valid BoundingBox")

    @classmethod
    def from_tree(cls, node, cts):
        bounding_box = {
            _input: tuple(interval) for _input, interval in node['intervals']
        }
        if 'ignore' in node:
            ignored = node['ignore']
        else:
            ignored = None

        if 'order' in node:
            order = node['order']
        else:
            order = 'C'

        return {
            'bounding_box': bounding_box,
            'ignored': ignored,
            'order': order
        }


class CompoundBoundingBoxType(AstropyAsdfType):
    name = 'transform/property/compound_bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.CompoundBoundingBox']

    @classmethod
    def to_tree(cls, bbox, cts):
        if isinstance(bbox, CompoundBoundingBox):
            return {
                'selector_args': [
                    {
                        'argument': sa.name(bbox._model),
                        'ignore': sa.ignore
                    } for sa in bbox.selector_args
                ],
                'cbbox': [
                    {
                        'key': list(key),
                        'bbox': bb
                    } for key, bb in bbox.bounding_boxes.items()
                ],
                'ignore': bbox.ignored_inputs,
                'order': bbox.order
            }
        else:
            raise TypeError(f"{bbox} is not a valid CompoundBoundingBox")

    @classmethod
    def from_tree(cls, node, cts):
        selector_args = [
            (sa['argument'], sa['ignore']) for sa in node['selector_args']
        ]
        bounding_box = {
            tuple(bb['key']): bb['bbox']['bounding_box']
            for bb in node['cbbox']
        }
        if 'ignore' in node:
            ignored = node['ignore']
        else:
            ignored = None

        if 'order' in node:
            order = node['order']
        else:
            order = 'C'

        return {
            'bounding_box': bounding_box,
            'selector_args': selector_args,
            'ignored': ignored,
            'order': order
        }
