# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np

from asdf.versioning import AsdfVersion

from astropy.modeling.bounding_box import ModelBoundingBox, CompoundBoundingBox
from astropy.modeling import mappings
from astropy.modeling import functional_models
from astropy.modeling.core import CompoundModel
from astropy.io.misc.asdf.types import AstropyAsdfType
from . import _parameter_to_value


class BoundingBoxType(AstropyAsdfType):
    name = 'transform/property/bounding_box'
    version = '1.0.0'
    types = ['astropy.modeling.bounding_box.ModelBoundingBox']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        pass

    @classmethod
    def to_tree_transform(cls, model, cts):
        def _transform_bbox(bbox):
            return {
                'intervals': {
                    _input: list(interval)
                    for _input, interval in bbox.named_intervals.items()
                },
                'ignore': list(bbox.ignored_inputs),
                'order': bbox.order
            }

        try:
            bb = model.bounding_box
        except NotImplementedError:
            bb = None

        if isinstance(bb, ModelBoundingBox):
            return _transform_bbox(bb)

        elif isinstance(bb, CompoundBoundingBox):
            return {
                'selector_args': [
                    {
                        'argument': sa.name(model),
                        'ignore': sa.ignore
                    } for sa in bb.selector_args
                ],
                'cbbox': [
                    {
                        'key': list(key),
                        'bbox': _transform_bbox(bbox)
                    } for key, bbox in bb.bounding_boxes.items()
                ]
            }
