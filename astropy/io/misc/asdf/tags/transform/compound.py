# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-


from asdf.yamlutil import tagged_tree_to_custom_tree

from asdf import tagged
from asdf import yamlutil
from asdf.tests.helpers import assert_tree_match
from .basic import TransformType, ConstantType
from astropy.modeling.core import Model, CompoundModel
from astropy.modeling.models import Identity, Mapping


__all__ = ['CompoundType', 'RemapAxesType']


_operator_to_tag_mapping = {
    '+'  : 'add',
    '-'  : 'subtract',
    '*'  : 'multiply',
    '/'  : 'divide',
    '**' : 'power',
    '|'  : 'compose',
    '&'  : 'concatenate',
    'fix_inputs': 'fix_inputs'
}


_tag_to_method_mapping = {
    'add'         : '__add__',
    'subtract'    : '__sub__',
    'multiply'    : '__mul__',
    'divide'      : '__truediv__',
    'power'       : '__pow__',
    'compose'     : '__or__',
    'concatenate' : '__and__',
    'fix_inputs'  : 'fix_inputs'
}


class CompoundType(TransformType):
    name = ['transform/' + x for x in _tag_to_method_mapping.keys()]
    types = [CompoundModel]
    version = '1.2.0'
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        tag = node._tag[node._tag.rfind('/')+1:]
        tag = tag[:tag.rfind('-')]
        oper = _tag_to_method_mapping[tag]
        left = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][0], ctx)
        if not isinstance(left, Model):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][0]._tag))
        right = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][1], ctx)
        if not isinstance(right, Model) and \
            not (oper == 'fix_inputs' and isinstance(right, dict)):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][1]._tag))
        if oper == 'fix_inputs':
            right = dict(zip(right['keys'], right['values']))
            model = CompoundModel('fix_inputs', left, right)
        else:
            model = getattr(left, oper)(right)

        model = cls._from_tree_base_transform_members(model, node, ctx)
        return model

    @classmethod
    def _to_tree_from_model_tree(cls, tree, ctx):

        if not isinstance(tree.left, CompoundModel):
            left = yamlutil.custom_tree_to_tagged_tree(
                tree.left, ctx)
        else:
            left = cls._to_tree_from_model_tree(tree.left, ctx)

        if not isinstance(tree.right, CompoundModel):
            if isinstance(tree.right, dict):
                right = {'keys': list(tree.right.keys()),
                         'values': list(tree.right.values())
                        }
            else:
                right = yamlutil.custom_tree_to_tagged_tree(
                    tree.right, ctx)
        else:
            right = cls._to_tree_from_model_tree(tree.right, ctx)

        node = {
            'forward': [left, right]
        }

        try:
            tag_name = 'transform/' + _operator_to_tag_mapping[tree.op]
        except KeyError:
            raise ValueError(f"Unknown operator '{tree.op}'")

        node = tagged.tag_object(cls.make_yaml_tag(tag_name), node, ctx=ctx)
        return node

    @classmethod
    def to_tree_tagged(cls, model, ctx):
        node = cls._to_tree_from_model_tree(model, ctx)
        cls._to_tree_base_transform_members(model, node, ctx)
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert_tree_match(a.left, b.left)
        assert_tree_match(a.right, b.right)


class RemapAxesType(TransformType):
    name = 'transform/remap_axes'
    types = [Mapping]

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
