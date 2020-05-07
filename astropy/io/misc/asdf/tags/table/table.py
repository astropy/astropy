# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np

from asdf.tags.core.ndarray import NDArrayType

from astropy import table
from astropy.io.misc.asdf.types import AstropyType, AstropyAsdfType


class TableType:
    """
    This class defines to_tree and from_tree methods that are used by both the
    AstropyTableType and the AsdfTableType defined below. The behavior is
    differentiated by the ``_compat`` class attribute. When ``_compat==True``,
    the behavior will conform to the table schema defined by the ASDF Standard.
    Otherwise, the behavior will conform to the custom table schema defined by
    Astropy.
    """
    _compat = False

    @classmethod
    def from_tree(cls, node, ctx):

        # This is getting meta, guys
        meta = node.get('meta', {})

        # This enables us to support files that use the table definition from
        # the ASDF Standard, rather than the custom one that Astropy defines.
        if cls._compat:
            return table.Table(node['columns'], meta=meta)

        if node.get('qtable', False):
            t = table.QTable(meta=node.get('meta', {}))
        else:
            t = table.Table(meta=node.get('meta', {}))

        for name, col in zip(node['colnames'], node['columns']):
            t[name] = col

        return t

    @classmethod
    def to_tree(cls, data, ctx):
        columns = [data[name] for name in data.colnames]

        node = dict(columns=columns)
        # Files that use the table definition from the ASDF Standard (instead
        # of the one defined by Astropy) will not contain these fields
        if not cls._compat:
            node['colnames'] = data.colnames
            node['qtable'] = isinstance(data, table.QTable)
        if data.meta:
            node['meta'] = data.meta

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert old.meta == new.meta
        try:
            NDArrayType.assert_equal(np.array(old), np.array(new))
        except (AttributeError, TypeError, ValueError):
            for col0, col1 in zip(old, new):
                try:
                    NDArrayType.assert_equal(np.array(col0), np.array(col1))
                except (AttributeError, TypeError, ValueError):
                    assert col0 == col1


class AstropyTableType(TableType, AstropyType):
    """
    This tag class reads and writes tables that conform to the custom schema
    that is defined by Astropy (in contrast to the one that is defined by the
    ASDF Standard). The primary reason for differentiating is to enable the
    support of Astropy mixin columns, which are not supported by the ASDF
    Standard.
    """
    name = 'table/table'
    types = ['astropy.table.Table']
    requires = ['astropy']


class AsdfTableType(TableType, AstropyAsdfType):
    """
    This tag class allows Astropy to read (and write) ASDF files that use the
    table definition that is provided by the ASDF Standard (instead of the
    custom one defined by Astropy). This is important to maintain for
    cross-compatibility.
    """
    name = 'core/table'
    types = ['astropy.table.Table']
    requires = ['astropy']
    _compat = True


class ColumnType(AstropyAsdfType):
    name = 'core/column'
    types = ['astropy.table.Column', 'astropy.table.MaskedColumn']
    requires = ['astropy']
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree(cls, node, ctx):
        data = node['data']
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
            'data': data.data,
            'name': data.name
        }
        if data.description:
            node['description'] = data.description
        if data.unit:
            node['unit'] = data.unit
        if data.meta:
            node['meta'] = data.meta

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert old.meta == new.meta
        assert old.description == new.description
        assert old.unit == new.unit

        NDArrayType.assert_equal(np.array(old), np.array(new))
