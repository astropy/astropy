# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import numpy as np

from asdf.versioning import AsdfVersion

from astropy.modeling.bounding_box import ModelBoundingBox, CompoundBoundingBox
from astropy.modeling import mappings
from astropy.modeling import functional_models
from astropy.modeling.core import CompoundModel
from astropy.io.misc.asdf.types import AstropyAsdfType, AstropyType
from . import _parameter_to_value


__all__ = ['TransformType', 'IdentityType', 'ConstantType']


class TransformType(AstropyAsdfType):
    version = '1.2.0'
    requires = ['astropy']

    @classmethod
    def _from_tree_base_transform_members(cls, model, node, ctx):
        if 'name' in node:
            model.name = node['name']

        if "inputs" in node:
            model.inputs = tuple(node["inputs"])

        if "outputs" in node:
            model.outputs = tuple(node["outputs"])

        if 'bounding_box' in node:
            model.bounding_box = node['bounding_box']

        elif 'selector_args' in node:
            cbbox_keys = [tuple(key) for key in node['cbbox_keys']]
            bbox_dict = dict(zip(cbbox_keys, node['cbbox_values']))

            selector_args = node['selector_args']
            model.bounding_box = CompoundBoundingBox.validate(model, bbox_dict, selector_args)

        param_and_model_constraints = {}
        for constraint in ['fixed', 'bounds']:
            if constraint in node:
                param_and_model_constraints[constraint] = node[constraint]
        model._initialize_constraints(param_and_model_constraints)

        if "input_units_equivalencies" in node:
            # this still writes eqs. for compound, but operates on each sub model
            if not isinstance(model, CompoundModel):
                model.input_units_equivalencies = node['input_units_equivalencies']

        yield model

        if 'inverse' in node:
            model.inverse = node['inverse']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        raise NotImplementedError(
            "Must be implemented in TransformType subclasses")

    @classmethod
    def from_tree(cls, node, ctx):
        model = cls.from_tree_transform(node, ctx)
        return cls._from_tree_base_transform_members(model, node, ctx)

    @classmethod
    def _to_tree_base_transform_members(cls, model, node, ctx):
        def _transform_bbox(bbox):
            return {
                'intervals': {
                    _input: list(interval) for _input, interval in bbox.named_intervals.items()
                },
                'ignore': list(bbox.ignored_inputs),
                'order': bbox.order
            }

        if getattr(model, '_user_inverse', None) is not None:
            node['inverse'] = model._user_inverse

        if model.name is not None:
            node['name'] = model.name

        node['inputs'] = list(model.inputs)
        node['outputs'] = list(model.outputs)

        try:
            bb = model.bounding_box
        except NotImplementedError:
            bb = None

        if isinstance(bb, ModelBoundingBox):
            node['bounding_box'] = _transform_bbox(bb)

        elif isinstance(bb, CompoundBoundingBox):
            cbbox = {}
            cbbox['selector_args'] = [
                {
                    'argument': sa.name(model),
                    'ignore': sa.ignore
                } for sa in bb.selector_args
            ]

            cbbox['cbbox'] = [
                {
                    'key': list(key),
                    'bbox': _transform_bbox(bbox)
                } for key, bbox in bb.bounding_boxes.items()
            ]
            node['bounding_box'] = cbbox

        # model / parameter constraints
        if not isinstance(model, CompoundModel):
            fixed_nondefaults = {k: f for k, f in model.fixed.items() if f}
            if fixed_nondefaults:
                node['fixed'] = fixed_nondefaults
            bounds_nondefaults = {k: b for k, b in model.bounds.items() if any(b)}
            if bounds_nondefaults:
                node['bounds'] = bounds_nondefaults

        if not isinstance(model, CompoundModel):
            if model.input_units_equivalencies:
                node['input_units_equivalencies'] = model.input_units_equivalencies

        return node

    @classmethod
    def to_tree_transform(cls, model, ctx):
        raise NotImplementedError("Must be implemented in TransformType subclasses")

    @classmethod
    def to_tree(cls, model, ctx):
        node = cls.to_tree_transform(model, ctx)
        return cls._to_tree_base_transform_members(model, node, ctx)

    @classmethod
    def assert_equal(cls, a, b):
        # TODO: If models become comparable themselves, remove this.
        assert a.name == b.name
        # TODO: Assert inverses are the same
        # assert the bounding_boxes are the same
        assert a.get_bounding_box() == b.get_bounding_box()
        assert a.inputs == b.inputs
        assert a.outputs == b.outputs
        assert a.input_units_equivalencies == b.input_units_equivalencies


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
    version = '1.4.0'
    supported_versions = ['1.0.0', '1.1.0', '1.2.0', '1.3.0', '1.4.0']
    types = ['astropy.modeling.functional_models.Const1D',
             'astropy.modeling.functional_models.Const2D']

    @classmethod
    def from_tree_transform(cls, node, ctx):
        if cls.version < AsdfVersion('1.4.0'):
            # The 'dimensions' property was added in 1.4.0,
            # previously all values were 1D.
            return functional_models.Const1D(node['value'])
        elif node['dimensions'] == 1:
            return functional_models.Const1D(node['value'])
        elif node['dimensions'] == 2:
            return functional_models.Const2D(node['value'])
        else:
            raise TypeError('Only 1D and 2D constant models are supported.')

    @classmethod
    def to_tree_transform(cls, data, ctx):
        if cls.version < AsdfVersion('1.4.0'):
            if not isinstance(data, functional_models.Const1D):
                raise ValueError(
                    f'constant-{cls.version} does not support models with > 1 dimension')
            return {
                'value': _parameter_to_value(data.amplitude)
            }
        else:
            if isinstance(data, functional_models.Const1D):
                dimension = 1
            elif isinstance(data, functional_models.Const2D):
                dimension = 2
            return {
                'value': _parameter_to_value(data.amplitude),
                'dimensions': dimension
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
                input["unit"] = m[0]
            if node.input_units_equivalencies is not None and i in node.input_units_equivalencies:
                input["equivalencies"] = node.input_units_equivalencies[i]
            inputs.append(input)

            output = {
                "name": o,
            }
            if m[-1] is not None:
                output["unit"] = m[-1]
            outputs.append(output)

        tree["inputs"] = inputs
        tree["outputs"] = outputs

        return tree

    @classmethod
    def from_tree(cls, tree, ctx):
        mapping = tuple((i.get("unit"), o.get("unit"))
                        for i, o in zip(tree["inputs"], tree["outputs"]))

        equivalencies = None
        for i in tree["inputs"]:
            if "equivalencies" in i:
                if equivalencies is None:
                    equivalencies = {}
                equivalencies[i["name"]] = i["equivalencies"]

        kwargs = {
            "input_units_equivalencies": equivalencies,
            "input_units_allow_dimensionless": {
                i["name"]: i.get("allow_dimensionless", False) for i in tree["inputs"]},
        }

        if "name" in tree:
            kwargs["name"] = tree["name"]

        return mappings.UnitsMapping(mapping, **kwargs)
