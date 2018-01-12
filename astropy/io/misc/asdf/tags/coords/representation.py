from importlib import import_module
from asdf.yamlutil import custom_tree_to_tagged_tree

from ...types import AstropyType


REPRESENTATION_MODULE_WHITELIST = ['astropy.coordinates.representation',
                                   'sunpy.coordinates.representation']


class RepresentationType(AstropyType):
    name = "coords/representation/baserepresentation"
    types = ['astropy.coordinates.representation.BaseRepresentationOrDifferential']
    requires = ['astropy']
    version = "1.0.0"

    @classmethod
    def to_tree(cls, representation, ctx):
        comps = representation.components
        components = {}
        for c in comps:
            value = getattr(representation, '_' + c, None)
            if value:
                components[c] = value

        t = type(representation)
        if t.__module__ not in REPRESENTATION_MODULE_WHITELIST:
            raise ValueError("Can only save representations on the white list.")

        name = t.__module__ + '.' + t.__name__

        node = {}
        node['type'] = name
        node['components'] = custom_tree_to_tagged_tree(components, ctx)

        return node

    @classmethod
    def from_tree(cls, node, ctx):
        *module_name, class_name = node['type'].split('.')
        module_name = '.'.join(module_name)
        if module_name not in REPRESENTATION_MODULE_WHITELIST:
            raise ValueError("Can only load representations from modules that "
                             "have been white listed.")
        mod = import_module(module_name)
        rep_type = getattr(mod, class_name)
        return rep_type(**node['components'])

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(new, type(old))
        assert new.components == old.components
