# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf import tagged, yamlutil
from asdf.tests.helpers import assert_tree_match

from astropy import modeling
from astropy.modeling.models import Identity, Mapping
from .basic import TransformType, ConstantType


__all__ = ['CompoundType', 'RemapAxesType']


_operator_to_tag_mapping = {
    '+'  : 'add',
    '-'  : 'subtract',
    '*'  : 'multiply',
    '/'  : 'divide',
    '**' : 'power',
    '|'  : 'compose',
    '&'  : 'concatenate'
}


_tag_to_method_mapping = {
    'add'         : '__add__',
    'subtract'    : '__sub__',
    'multiply'    : '__mul__',
    'divide'      : '__truediv__',
    'power'       : '__pow__',
    'compose'     : '__or__',
    'concatenate' : '__and__'
}


class CompoundType(TransformType):
    name = ['transform/' + x for x in _tag_to_method_mapping.keys()]
    types = ['astropy.modeling.core._CompoundModel']
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        tag = node._tag[node._tag.rfind('/')+1:]
        tag = tag[:tag.rfind('-')]

        oper = _tag_to_method_mapping[tag]
        left = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][0], ctx)
        if not isinstance(left, modeling.Model):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][0]._tag))
        right = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][1], ctx)
        if not isinstance(right, modeling.Model):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][1]._tag))
        model = getattr(left, oper)(right)

        model = cls._from_tree_base_transform_members(model, node, ctx)
        return model

    @classmethod
    def _to_tree_from_model_tree(cls, tree, ctx):
        if tree.left.isleaf:
            left = yamlutil.custom_tree_to_tagged_tree(
                tree.left.value, ctx)
        else:
            left = cls._to_tree_from_model_tree(tree.left, ctx)

        if tree.right.isleaf:
            right = yamlutil.custom_tree_to_tagged_tree(
                tree.right.value, ctx)
        else:
            right = cls._to_tree_from_model_tree(tree.right, ctx)

        node = {
            'forward': [left, right]
        }

        try:
            tag_name = 'transform/' + _operator_to_tag_mapping[tree.value]
        except KeyError:
            raise ValueError("Unknown operator '{0}'".format(tree.value))

        node = tagged.tag_object(cls.make_yaml_tag(tag_name), node, ctx=ctx)
        return node

    @classmethod
    def to_tree_tagged(cls, model, ctx):
        node = cls._to_tree_from_model_tree(model._tree, ctx)
        cls._to_tree_base_transform_members(model, node, ctx)
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert_tree_match(a._tree.left.value, b._tree.left.value)
        assert_tree_match(a._tree.right.value, b._tree.right.value)
        assert a._tree.value == b._tree.value


class RemapAxesType(TransformType):
    name = 'transform/remap_axes'
    types = ['astropy.modeling.models.Mapping']

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
                transform = transform & ConstantType.from_tree(
                    {'value': int(entry.value)}, ctx)
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
