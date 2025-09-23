# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__: tuple[str, ...] = ()  # nothing is publicly scoped

import functools
import operator
from dataclasses import Field
from typing import TypeGuard

from .core import Parameter


def is_parameter_or_field(obj: object, /) -> TypeGuard[Parameter | Field[Parameter]]:
    """Return if object is a Parameter or dataclass field thereof."""
    return isinstance(obj, Parameter) or (
        isinstance(obj, Field) and isinstance(obj.default, Parameter)
    )


def all_parameters(obj: object, /) -> dict[str, Parameter]:
    """Get all `Parameter` fields of an object.

    Returns all fields of the object, including those not yet finalized in the class, if
    it's still under construction, e.g. in ``__init_subclass__``.
    """
    cls = obj if isinstance(obj, type) else obj.__class__
    all_cls_vars = functools.reduce(operator.__or__, map(vars, cls.mro()[::-1]))

    return {
        k: (v if isinstance(v, Parameter) else v.default)
        for k, v in all_cls_vars.items()
        if is_parameter_or_field(v)
    }
