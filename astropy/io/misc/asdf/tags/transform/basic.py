# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-

from asdf import tagged, yamlutil

from astropy.modeling import mappings
from astropy.utils import minversion
from astropy.modeling import functional_models
from astropy.io.misc.asdf.types import AstropyAsdfType, AstropyType


__all__ = ['TransformType', 'IdentityType', 'ConstantType']


class TransformType(AstropyAsdfType):
    version = '1.2.0'
    requires = ['astropy']

    @classmethod
    def _from_tree_base_transform_members(cls, model, node, ctx):
        if 'inverse' in node:
            model.inverse = yamlutil.tagged_tree_to_custom_tree(
                node['inverse'], ctx)

        if 'name' in node:
            model = model.rename(node['name'])

        if 'bounding_box' in node:
            model.bounding_box = yamlutil.tagged_tree_to_custom_tree(node['bounding_box'], ctx)

        if "inputs" in node:
            if model.n_inputs == 1:
                model.inputs = (node["inputs"],)
            else:
                model.inputs = tuple(node["inputs"])

        if "outputs" in node:
            if model.n_outputs == 1:
                model.outputs = (node["outputs"],)
            else:
                model.outputs = tuple(node["outputs"])

        return model

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
            node['bounding_box'] = yamlutil.custom_tree_to_tagged_tree(bb, ctx)
        if type(model.__class__.inputs) != property:
            node['inputs'] = model.inputs
            node['outputs'] = model.outputs

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


class GenericModel(mappings.Mapping):

    def __init__(self, n_inputs, n_outputs):
        mapping = tuple(range(n_inputs))
        super().__init__(mapping)
        self._n_outputs = n_outputs
        self._outputs = tuple('x' + str(idx) for idx in range(n_outputs))

    @property
    def inverse(self):
        raise NotImplementedError()


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


class UnitsMappingType(AstropyType):
    name = "transform/units_mapping"
    version = "1.0.0"
    types = [mappings.UnitsMapping]

    @classmethod
    def to_tree(cls, node, ctx):
        tree = {}

        if node.name is not None:
            tree["name"] = node.name

        inputs = []
        outputs = []
        for i, o, m in zip(node.inputs, node.outputs, node.mapping):
            input = {
                "name": i,
                "allow_dimensionless": node.input_units_allow_dimensionless[i],
            }
            if m[0] is not None:
                input["unit"] = yamlutil.custom_tree_to_tagged_tree(m[0], ctx)
            if node.input_units_equivalencies is not None and i in node.input_units_equivalencies:
                input["equivalencies"] = yamlutil.custom_tree_to_tagged_tree(node.input_units_equivalencies[i], ctx)
            inputs.append(input)

            output = {
                "name": o,
            }
            if m[-1] is not None:
                output["unit"] = yamlutil.custom_tree_to_tagged_tree(m[-1], ctx)
            outputs.append(output)

        tree["inputs"] = inputs
        tree["outputs"] = outputs

        return tree


    @classmethod
    def from_tree(cls, tree, ctx):
        mapping = tuple((i.get("unit"), o.get("unit")) for i, o in zip(tree["inputs"], tree["outputs"]))

        equivalencies = None
        for i in tree["inputs"]:
            if "equivalencies" in i:
                if equivalencies is None:
                    equivalencies = {}
                equivalencies[i["name"]] = yamlutil.tagged_tree_to_custom_tree(i["equivalencies"], ctx)

        kwargs = {
            "input_units_equivalencies": equivalencies,
            "input_units_allow_dimensionless": {i["name"]: i.get("allow_dimensionless", False) for i in tree["inputs"]},
        }

        if "name" in tree:
            kwargs["name"] = tree["name"]

        return mappings.UnitsMapping(mapping, **kwargs)
