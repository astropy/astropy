# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np
from numpy.testing import assert_array_equal

from astropy import modeling
from astropy import units as u
from astropy.modeling.bounding_box import ModelBoundingBox
from astropy.io.misc.asdf.tags.transform.basic import TransformType

__all__ = ['TabularType']


class TabularType(TransformType):
    name = "transform/tabular"
    version = '1.2.0'
    types = [
        modeling.models.Tabular2D, modeling.models.Tabular1D
    ]

    @classmethod
    def from_tree_transform(cls, node, ctx):
        lookup_table = node.pop("lookup_table")
        dim = lookup_table.ndim
        fill_value = node.pop("fill_value", None)
        if dim == 1:
            # The copy is necessary because the array is memory mapped.
            points = (node['points'][0][:],)
            model = modeling.models.Tabular1D(points=points, lookup_table=lookup_table,
                                              method=node['method'], bounds_error=node['bounds_error'],
                                              fill_value=fill_value)
        elif dim == 2:
            points = tuple([p[:] for p in node['points']])
            model = modeling.models.Tabular2D(points=points, lookup_table=lookup_table,
                                              method=node['method'], bounds_error=node['bounds_error'],
                                              fill_value=fill_value)

        else:
            tabular_class = modeling.models.tabular_model(dim, name)
            points = tuple([p[:] for p in node['points']])
            model = tabular_class(points=points, lookup_table=lookup_table,
                                  method=node['method'], bounds_error=node['bounds_error'],
                                  fill_value=fill_value)

        return model

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {}
        if model.fill_value is not None:
            node["fill_value"] = model.fill_value
        node["lookup_table"] = model.lookup_table
        node["points"] = [p for p in model.points]
        node["method"] = str(model.method)
        node["bounds_error"] = model.bounds_error
        return node

    @classmethod
    def assert_equal(cls, a, b):
        if isinstance(a.lookup_table, u.Quantity):
            assert u.allclose(a.lookup_table, b.lookup_table)
            assert u.allclose(a.points, b.points)
            a_box = a.bounding_box
            if isinstance(a_box, ModelBoundingBox):
                a_box = a_box.bounding_box()
            b_box = b.bounding_box
            if isinstance(b_box, ModelBoundingBox):
                b_box = b_box.bounding_box()
            for i in range(len(a_box)):
                assert u.allclose(a_box[i], b_box[i])
        else:
            assert_array_equal(a.lookup_table, b.lookup_table)
            assert_array_equal(a.points, b.points)
            a_box = a.bounding_box
            if isinstance(a_box, ModelBoundingBox):
                a_box = a_box.bounding_box()
            b_box = b.bounding_box
            if isinstance(b_box, ModelBoundingBox):
                b_box = b_box.bounding_box()
            assert_array_equal(a_box, b_box)
        assert (a.method == b.method)
        if a.fill_value is None:
            assert b.fill_value is None
        elif np.isnan(a.fill_value):
            assert np.isnan(b.fill_value)
        else:
            assert(a.fill_value == b.fill_value)
        assert(a.bounds_error == b.bounds_error)
