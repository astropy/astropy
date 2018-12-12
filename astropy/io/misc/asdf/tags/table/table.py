# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np

from asdf import tagged
from asdf import yamlutil
from asdf.tags.core.ndarray import NDArrayType

from astropy import table
from astropy.io.misc.asdf.types import AstropyType, AstropyAsdfType


class TableType(AstropyType):
    name = 'table/table'
    types = ['astropy.table.Table']
    requires = ['astropy']

    @classmethod
    def from_tree(cls, node, ctx):

        if node.get('qtable', False):
            t = table.QTable(meta=node.get('meta', {}))
        else:
            t = table.Table(meta=node.get('meta', {}))

        for name, col in zip(node['colnames'], node['columns']):
            t[name] = yamlutil.tagged_tree_to_custom_tree(col, ctx)

        return t

    @classmethod
    def to_tree(cls, data, ctx):
        columns = []
        for name in data.colnames:
            thiscol = data[name]
            column = yamlutil.custom_tree_to_tagged_tree(thiscol, ctx)
            columns.append(column)

        node = dict(
            columns=columns,
            colnames=data.colnames,
            qtable = isinstance(data, table.QTable)
        )
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
