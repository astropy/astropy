# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

import asdf
from asdf import CustomType
from asdf.extension import BuiltinExtension
from asdf.yamlutil import tagged_tree_to_custom_tree

from asdf import tagged
from asdf import yamlutil

import astropy
from astropy.io.misc.asdf.tags.transform.basic import TransformType
from astropy.modeling.core import Model
from astropy.modeling.compound import CompoundModel

__all__ = ['NewCompoundType']


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


class NewCompoundType(TransformType):
    name = ['transform/' + x for x in _tag_to_method_mapping.keys()]
    types = [astropy.modeling.compound.CompoundModel]
    handle_dynamic_subclasses = True

    @classmethod
    def from_tree_tagged(cls, node, ctx):
        from astropy import modeling
        from astropy.modeling import core
        from astropy.modeling.core import Model
        from astropy.modeling.compound import CompoundModel

        tag = node._tag[node._tag.rfind('/')+1:]
        tag = tag[:tag.rfind('-')]
        oper = _tag_to_method_mapping[tag]
        left = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][0], ctx)
        if not (isinstance(left, Model) or 
                isinstance(left, CompoundModel)):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][0]._tag))
        right = yamlutil.tagged_tree_to_custom_tree(
            node['forward'][1], ctx)
        if not (isinstance(right, Model) or 
                isinstance(right, CompoundModel)):
            raise TypeError("Unknown model type '{0}'".format(
                node['forward'][1]._tag))
        model = getattr(left, oper)(right)
        model = cls._from_tree_base_transform_members(model, node, ctx)
        return model

    @classmethod
    def _to_tree_from_model_tree(cls, tree, ctx):
        if isinstance(tree, CompoundModel) and isinstance(tree.left, Model):
            left = yamlutil.custom_tree_to_tagged_tree(tree.left, ctx)
        else:
            left = cls._to_tree_from_model_tree(tree.left, ctx)
        if isinstance(tree, CompoundModel) and isinstance(tree.right, Model):
            right = yamlutil.custom_tree_to_tagged_tree(tree.right, ctx)
        else:
            right = cls._to_tree_from_model_tree(tree.right, ctx)
        node = {
            'forward': [left, right]
        }
        try:
            if isinstance(tree, CompoundModel):
                value = tree.op
            else:
                value = tree.value
            tag_name = 'transform/' + _operator_to_tag_mapping[value]
        except KeyError:
            raise ValueError("Unknown operator '{0}'".format(value))

        node = tagged.tag_object(cls.make_yaml_tag(tag_name), node, ctx=ctx)
        return node

    @classmethod
    def to_tree_tagged(cls, model, ctx):

        tree = model
        node = cls._to_tree_from_model_tree(tree, ctx)
        cls._to_tree_base_transform_members(model, node, ctx)
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        #TransformType.assert_equal(a, b)
        #from ...tests.helpers import assert_tree_match
        #assert_tree_match(a._tree.left.value, b._tree.left.value)
        #assert_tree_match(a._tree.right.value, b._tree.right.value)
        #assert a._tree.value == b._tree.value
        pass


class NewExtension(BuiltinExtension):
    """We inherit from ASDF's BuiltinExtension here since NewCompoundType does
    not require any new schemas. We simply use the schemas that are provided by
    ASDF, and we override the type we care about."""

    @property
    def types(self):
        return [NewCompoundType]