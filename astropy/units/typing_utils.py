"""
Utilities to help with typing units.

Function to generate module-level stub files (.pyi) to help type checkers
understand the dynamically created unit objects in the astropy.units subpackage.
It inspects each specified module defining units and writes out a .pyi file
containing the types of all discoverable UnitBase and FunctionUnitBase instances.

This module is meant for internal use only (in setup.py).
"""

import importlib
import unicodedata
import warnings
from collections.abc import Iterable, Iterator
from pathlib import Path

import astropy.units
from astropy.units import FunctionUnitBase, UnitBase
from astropy.utils.exceptions import AstropyDeprecationWarning

# The modules that contain dynamically generated units and need stubs.
# We do not generate stubs for modules like "core" or "quantity" because
# type checkers can analyze their statically defined objects directly.
MODULES_TO_STUB = [
    "astrophys",
    "cds",
    "cgs",
    "deprecated",
    "imperial",
    "misc",
    "photometric",
    "required_by_vounit",
    "si",
    "function.units",
]


def create_stubs(output_dir: Path | None = None) -> None:
    if output_dir is None:
        target_dir = Path(astropy.units.__file__).parent

    # Ignore deprecation warnings triggered by accessing units in deprecated.py
    with warnings.catch_warnings(action="ignore", category=AstropyDeprecationWarning):
        for module_name in MODULES_TO_STUB:
            module = importlib.import_module(f"astropy.units.{module_name}")

            stub_file_path = output_dir / (module_name.replace(".", "/") + ".pyi")
            stub_lines = generate_stub_lines(module)
            with open(stub_file_path, "w", encoding="utf-8") as f:
                f.write("\n".join(stub_lines))


def generate_stub_lines(module) -> Iterator[str]:
    """Generate lines for a single .pyi stub file by introspecting a module."""
    import_lines = set()
    variable_lines = []

    # --- preprocessing ---
    # use module.__all__ to respect the public API
    for name in unique(nfkc_normalize(module.__all__)):
        obj = getattr(module, name)
        # find all unit-like objects
        if isinstance(obj, (UnitBase, FunctionUnitBase)):
            cls = type(obj)
            type_name = cls.__name__
            docstring = f'"""{obj.__doc__}"""\n' if has_own_docstring(obj) else ""
            variable_lines.append(f"{name}: {type_name}\n{docstring}")
            # the class of the unit object (e.g., IrreducibleUnit) is often
            # not defined in the module itself, so we need to import it
            import_lines.add(make_import_line(cls))

    # --- output ---
    yield "# This stub file was automatically generated."
    yield ""
    yield from sorted(import_lines)
    yield ""
    yield from variable_lines


def make_import_line(cls) -> str:
    """Create a "from ... import ..." line for a given class."""
    return f"from {cls.__module__} import {cls.__name__}"


def has_own_docstring(obj) -> bool:
    """Check if object has its own docstring vs an inherited one."""
    if not hasattr(obj, "__doc__") or not obj.__doc__:
        return False

    # for classes, check if docstring differs from parent classes
    if isinstance(obj, type):
        for base in obj.__bases__:
            if hasattr(base, "__doc__") and base.__doc__ == obj.__doc__:
                return False

    # for instances, check if docstring differs from class
    elif hasattr(obj, "__class__"):
        if hasattr(obj.__class__, "__doc__") and obj.__class__.__doc__ == obj.__doc__:
            return False

    return True


def nfkc_normalize(items: Iterable[str]) -> Iterator[str]:
    yield from (unicodedata.normalize("NFKC", item) for item in items)


def unique(items: Iterable) -> Iterator:
    """Iterate over unique items while preserving order."""
    seen = set()
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        yield item
