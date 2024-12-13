# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import annotations

__all__: list[str] = []  # nothing is publicly scoped

from dataclasses import Field

from astropy.cosmology._src.utils import all_cls_vars

from .core import Parameter


def all_parameters(obj: object, /) -> dict[str, Field | Parameter]:
    """Get all `Parameter` fields of a dataclass.

    Parameters
    ----------
    obj : object | type
        A dataclass.

    Returns
    -------
    dict[str, Field | Parameter]
        All fields of the dataclass, including those not yet finalized in the class, if
        it's still under construction, e.g. in ``__init_subclass__``.
    """
    return {
        k: (v if isinstance(v, Parameter) else v.default)
        for k, v in all_cls_vars(obj).items()
        if (
            isinstance(v, Parameter)
            or (isinstance(v, Field) and isinstance(v.default, Parameter))
        )
    }
