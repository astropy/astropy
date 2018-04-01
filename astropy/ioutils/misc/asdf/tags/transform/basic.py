# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf import tagged, yamlutil

from astropy.modeling import mappings
from astropy.utils import minversion
from astropy.modeling import functional_models
from ...types import AstropyAsdfType


__all__ = ['TransformType', 'IdentityType', 'ConstantType', 'DomainType']


class TransformType(AstropyAsdfType):
    version = '1.1.0'
    requires = ['astropy']

    @classmethod
    def _from_tree_base_transform_members(cls, model, node, ctx):
        if 'inverse' in node:
            model.inverse = yamlutil.tagged_tree_to_custom_tree(
                node['inverse'], ctx)

        if 'name' in node:
            model = model.rename(node['name'])

        # TODO: Remove domain in a later version.
        if 'domain' in node:
            model.bounding_box = cls._domain_to_bounding_box(node['domain'])
        elif 'bounding_box' in node:
            model.bounding_box = node['bounding_box']

        return model

    @classmethod
    def _domain_to_bounding_box(cls, domain):
        bb = tuple([(item['lower'], item['upper']) for item in domain])
        if len(bb) == 1:
            bb = bb[0]
        return bb

    @classmethod
    def from_tree_transform(cls, node, ctx):
        raise NotImplementedError(
            "Must be implemented in TransformType subclasses")

    @classmethod
    def from_tree(cls, node, ctx):
        model = cls.from_tree_transform(node, ctx)
        model = cls._from_tree_base_transform_members(model, node, ctx)
        return model

    @classmethod
    def _to_tree_base_transform_members(cls, model, node, ctx):
        if getattr(model, '_user_inverse', None) is not None:
            node['inverse'] = yamlutil.custom_tree_to_tagged_tree(
            model._user_inverse, ctx)

        if model.name is not None:
            node['name'] = model.name

        try:
            bb = model.bounding_box
        except NotImplementedError:
            bb = None

        if bb is not None:
            if model.n_inputs == 1:
                bb = list(bb)
            else:
                bb = [list(item) for item in model.bounding_box]
            node['bounding_box'] = bb


    @classmethod
    def to_tree_transform(cls, model, ctx):
        raise NotImplementedError("Must be implemented in TransformType subclasses")

    @classmethod
    def to_tree(cls, model, ctx):
        node = cls.to_tree_transform(model, ctx)
        cls._to_tree_base_transform_members(model, node, ctx)
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        assert a.name == b.name
        # TODO: Assert inverses are the same


class IdentityType(TransformType):
    name = "transform/identity"
    types = ['astropy.modeling.mappings.Identity']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return mappings.Identity(node.get('n_dims', 1))

    @classmethod
    def to_tree_transform(cls, data, ctx):
        node = {}
        if data.n_inputs != 1:
            node['n_dims'] = data.n_inputs
        return node

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        TransformType.assert_equal(a, b)
        assert (isinstance(a, mappings.Identity) and
                isinstance(b, mappings.Identity) and
                a.n_inputs == b.n_inputs)


class ConstantType(TransformType):
    name = "transform/constant"
    types = ['astropy.modeling.functional_models.Const1D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return functional_models.Const1D(node['value'])

    @classmethod
    def to_tree_transform(cls, data, ctx):
        return {
            'value': data.amplitude.value
        }


class DomainType(AstropyAsdfType):
    # TODO: Is this used anywhere? Can it be removed?
    name = "transform/domain"
    version = '1.0.0'

    @classmethod
    def from_tree(cls, node, ctx):
        return node

    @classmethod
    def to_tree(cls, data, ctx):
        return data


class GenericModel(mappings.Mapping):
    def __init__(self, n_inputs, n_outputs):
        mapping = tuple(range(n_inputs))
        super(GenericModel, self).__init__(mapping)
        self._outputs = tuple('x' + str(idx) for idx in range(self.n_outputs + 1))


class GenericType(TransformType):
    name = "transform/generic"
    types = [GenericModel]

    @classmethod
    def from_tree_transform(cls, node, ctx):
        return GenericModel(
            node['n_inputs'], node['n_outputs'])

    @classmethod
    def to_tree_transform(cls, data, ctx):
        return {
            'n_inputs': data.n_inputs,
            'n_outputs': data.n_outputs
        }
