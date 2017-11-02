# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np

from asdf import yamlutil
from asdf.tags.core.ndarray import NDArrayType

from astropy import table
from ...types import AstropyAsdfType


class TableType(AstropyAsdfType):
    name = 'core/table'
    types = ['astropy.table.Table']
    requires = ['astropy']

    @classmethod
    def from_tree(cls, node, ctx):

        columns = [
            yamlutil.tagged_tree_to_custom_tree(c, ctx)
            for c in node['columns']
        ]

        return table.Table(columns, meta=node.get('meta', {}))

    @classmethod
    def to_tree(cls, data, ctx):
        columns = []
        for name in data.colnames:
            column = yamlutil.custom_tree_to_tagged_tree(
                data.columns[name], ctx)
            columns.append(column)

        node = {'columns': columns}
        if data.meta:
            node['meta'] = data.meta

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert old.meta == new.meta
        NDArrayType.assert_equal(np.array(old), np.array(new))


class ColumnType(AstropyAsdfType):
    name = 'core/column'
    types = ['astropy.table.Column', 'astropy.table.MaskedColumn']
    requires = ['astropy']
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree(cls, node, ctx):
        data = yamlutil.tagged_tree_to_custom_tree(
            node['data'], ctx)
        name = node['name']
        description = node.get('description')
        unit = node.get('unit')
        meta = node.get('meta', None)

        return table.Column(
            data=data._make_array(), name=name, description=description,
            unit=unit, meta=meta)

    @classmethod
    def to_tree(cls, data, ctx):
        node = {
            'data': yamlutil.custom_tree_to_tagged_tree(
                data.data, ctx),
            'name': data.name
        }
        if data.description:
            node['description'] = data.description
        if data.unit:
            node['unit'] = yamlutil.custom_tree_to_tagged_tree(
                data.unit, ctx)
        if data.meta:
            node['meta'] = data.meta

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert old.meta == new.meta
        assert old.description == new.description
        assert old.unit == new.unit

        NDArrayType.assert_equal(np.array(old), np.array(new))
