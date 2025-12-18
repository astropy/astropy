# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""Utilities for generating type stubs for unit definition modules.

This module is meant for internal use only (in setup.py).
"""

__all__: list[str] = []

import inspect
from collections.abc import Iterable
from pathlib import Path
from types import FunctionType, ModuleType
from typing import Final

from astropy import units as u
from astropy.units import cds

UNIT_DEFINITION_MODULES: Final = (
    u.astrophys,
    cds,
    u.cgs,
    u.imperial,
    u.misc,
    u.photometric,
    u.required_by_vounit,
    u.si,
    u.function.units,
)


def generate_stub(output_dir: Path, module: ModuleType) -> None:
    members = {}
    functions = {}
    for name in module.__all__:
        member = getattr(module, name)
        if isinstance(member, FunctionType):
            # Neither of the functions in unit definition modules have any parameters.
            functions[name] = inspect.get_annotations(member)["return"]
        else:
            members[name] = type(member)
    import_lines = sorted(
        f"from {cls.__module__} import {cls.__name__}\n"
        for cls in set(members.values()).union(functions.values())
    )
    with (
        output_dir.joinpath(*module.__name__.split("."))
        .with_suffix(".pyi")
        .open("wt", encoding="utf-8") as stub_file
    ):
        stub_file.write("# This stub file was automatically generated.\n\n")

        stub_file.write("__all__ = [\n")
        stub_file.writelines(f'    "{name}",\n' for name in module.__all__)
        stub_file.write("]\n\n")

        stub_file.writelines(sorted(import_lines))
        stub_file.write("\n")

        stub_file.writelines(
            f"{member}: {cls.__name__}\n" for member, cls in members.items()
        )

        stub_file.writelines(
            f"\ndef {name}() -> {return_type.__name__}: ...\n"
            for name, return_type in functions.items()
        )


def create_stubs(
    output_dir: Path, modules: Iterable[ModuleType] = UNIT_DEFINITION_MODULES
) -> None:
    for module in modules:
        generate_stub(output_dir, module)
