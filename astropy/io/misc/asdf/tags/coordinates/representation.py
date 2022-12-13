import astropy.coordinates.representation
import astropy.units as u
from astropy.coordinates.representation import BaseRepresentationOrDifferential
from astropy.io.misc.asdf.types import AstropyType


class RepresentationType(AstropyType):
    name = "coordinates/representation"
    types = [BaseRepresentationOrDifferential]
    version = "1.0.0"

    _representation_module = astropy.coordinates.representation

    @classmethod
    def to_tree(cls, representation, ctx):
        comps = representation.components
        components = {}
        for c in comps:
            value = getattr(representation, "_" + c, None)
            if value is not None:
                components[c] = value

        t = type(representation)

        node = {}
        node["type"] = t.__name__
        node["components"] = components

        return node

    @classmethod
    def from_tree(cls, node, ctx):
        rep_type = getattr(cls._representation_module, node["type"])
        return rep_type(**node["components"])

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(new, type(old))
        assert new.components == old.components
        for comp in new.components:
            nc = getattr(new, comp)
            oc = getattr(old, comp)
            assert u.allclose(nc, oc)
