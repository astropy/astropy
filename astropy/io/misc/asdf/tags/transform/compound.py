# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import warnings

from asdf import tagged
from asdf.tests.helpers import assert_tree_match
from astropy.modeling.core import Model, CompoundModel
from astropy.modeling.models import Identity, Mapping, Const1D
from astropy.io.misc.asdf.deprecation import create_asdf_deprecation_warning
from astropy.io.misc.asdf.tags.transform.basic import TransformType


__all__ = ['CompoundType', 'RemapAxesType']


_operator_to_tag_mapping = {
    '+':  'add',
    '-':  'subtract',
    '*':  'multiply',
    '/':  'divide',
    '**': 'power',
    '|':  'compose',
    '&':  'concatenate',
    'fix_inputs': 'fix_inputs'
}


_tag_to_method_mapping = {
    'add':         '__add__',
    'subtract':    '__sub__',
    'multiply':    '__mul__',
    'divide':      '__truediv__',
    'power':       '__pow__',
    'compose':     '__or__',
    'concatenate': '__and__',
    'fix_inputs':  'fix_inputs'
}


class CompoundType(TransformType):
    name = ['transform/' + x for x in _tag_to_method_mapping.keys()]
    types = [CompoundModel]
    version = '1.2.0'
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        warnings.warn(create_asdf_deprecation_warning())

        tag = node._tag[node._tag.rfind('/')+1:]
        tag = tag[:tag.rfind('-')]
        oper = _tag_to_method_mapping[tag]
        left = node['forward'][0]
        if not isinstance(left, Model):
            raise TypeError(f"Unknown model type '{node['forward'][0]._tag}'")
        right = node['forward'][1]
        if (not isinstance(right, Model) and
                not (oper == 'fix_inputs' and isinstance(right, dict))):
            raise TypeError(f"Unknown model type '{node['forward'][1]._tag}'")
        if oper == 'fix_inputs':
            right = dict(zip(right['keys'], right['values']))
            model = CompoundModel('fix_inputs', left, right)
        else:
            model = getattr(left, oper)(right)

        return cls._from_tree_base_transform_members(model, node, ctx)

    @classmethod
    def to_tree_tagged(cls, model, ctx):
        warnings.warn(create_asdf_deprecation_warning())

        left = model.left

        if isinstance(model.right, dict):
            right = {
                'keys': list(model.right.keys()),
                'values': list(model.right.values())
            }
        else:
            right = model.right

        node = {
            'forward': [left, right]
        }

        try:
            tag_name = 'transform/' + _operator_to_tag_mapping[model.op]
        except KeyError:
            raise ValueError(f"Unknown operator '{model.op}'")

        node = tagged.tag_object(cls.make_yaml_tag(tag_name), node, ctx=ctx)

        return cls._to_tree_base_transform_members(model, node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert_tree_match(a.left, b.left)
        assert_tree_match(a.right, b.right)


class RemapAxesType(TransformType):
    name = 'transform/remap_axes'
    types = [Mapping]
    version = '1.3.0'

    @classmethod
    def from_tree_transform(cls, node, ctx):
        mapping = node['mapping']
        n_inputs = node.get('n_inputs')
        if all([isinstance(x, int) for x in mapping]):
            return Mapping(tuple(mapping), n_inputs)

        if n_inputs is None:
            n_inputs = max([x for x in mapping
                            if isinstance(x, int)]) + 1

        transform = Identity(n_inputs)
        new_mapping = []
        i = n_inputs
        for entry in mapping:
            if isinstance(entry, int):
                new_mapping.append(entry)
            else:
                new_mapping.append(i)
                transform = transform & Const1D(entry.value)
                i += 1
        return transform | Mapping(new_mapping)

    @classmethod
    def to_tree_transform(cls, model, ctx):
        node = {'mapping': list(model.mapping)}
        if model.n_inputs > max(model.mapping) + 1:
            node['n_inputs'] = model.n_inputs
        return node

    @classmethod
    def assert_equal(cls, a, b):
        TransformType.assert_equal(a, b)
        assert a.mapping == b.mapping
        assert(a.n_inputs == b.n_inputs)
