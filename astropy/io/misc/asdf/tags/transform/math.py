# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.io.misc.asdf.tags.transform.basic import TransformType
from astropy.modeling import math_functions
from astropy.modeling.math_functions import *  # noqa: F401, F403
from astropy.modeling.math_functions import __all__ as math_classes

__all__ = ["NpUfuncType"]


class NpUfuncType(TransformType):
    name = "transform/math_functions"
    version = "1.0.0"
    types = ["astropy.modeling.math_functions." + kl for kl in math_classes]

    @classmethod
    def from_tree_transform(cls, node, ctx):
        klass_name = math_functions._make_class_name(node["func_name"])
        klass = getattr(math_functions, klass_name)
        return klass()

    @classmethod
    def to_tree_transform(cls, model, ctx):
        return {"func_name": model.func.__name__}
