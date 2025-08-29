"""
Simple stub generator for astropy.units.

This script generates module-level stub files (.pyi) to help type checkers
understand the dynamically created unit objects in the astropy.units subpackage.
It inspects each specified module at runtime and writes out a .pyi file
containing the types of all discoverable UnitBase and FunctionUnitBase instances.
"""

import argparse
import importlib
import sys
import unicodedata
import warnings
from collections.abc import Iterable, Iterator
from pathlib import Path

import astropy.units
from astropy import log
from astropy.units import FunctionUnitBase, UnitBase
from astropy.utils.exceptions import AstropyDeprecationWarning

# The modules that contain dynamically generated units and need stubs.
# We do not generate stubs for modules like "core" or "quantity" because
# type checkers can analyze their statically defined objects directly.
MODULES_TO_STUB = [
    "si",
    "cgs",
    "astrophys",
    "misc",
    "imperial",
    "photometric",
    "function.units",
    "deprecated",
    "required_by_vounit",
]


def main():
    args = _parse_args()
    if args.verbose:
        log.setLevel("DEBUG")

    log.info("Generating astropy.units stubs...")

    # Ignore deprecation warnings triggered by accessing units in deprecated.py
    with warnings.catch_warnings(action="ignore", category=AstropyDeprecationWarning):
        for module_name in MODULES_TO_STUB:
            full_module_name = f"astropy.units.{module_name}"
            try:
                module = importlib.import_module(full_module_name)
            except ImportError:
                log.warning(f"Could not import {full_module_name}, skipping.")
                continue

            stub_file_path = _get_stub_filepath(args.output_dir, module_name)
            stub_lines = generate_stub_lines(module)
            _write_stub_file(stub_file_path, stub_lines)


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
        else:
            log.debug(
                f"Skipped '{name}' (type: {type(obj).__name__}) in {module.__name__}"
            )

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


# --- Helper Functions ---


def _get_stub_filepath(output_dir: Path | None, module_name: str) -> Path:
    """Determine the full, final path for a module's stub file."""
    if output_dir:
        # user-provided directory
        target_dir = output_dir
    else:
        target_dir = Path(astropy.units.__file__).parent

    # construct the final path, including subdirectories
    relative_path = Path(module_name.replace(".", "/")).with_suffix(".pyi")
    return target_dir / relative_path


def _write_stub_file(path: Path, lines: Iterable[str]):
    try:
        # ensure the parent directory exists.
        path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        log.info(f"  - Wrote stub at: {path}")
    except PermissionError:
        log.error(
            f"Permission denied when trying to write to {path}.\n"
            "Try running again with the --output-dir option to specify a "
            "writable directory."
        )
        sys.exit(1)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate .pyi type stub files for astropy.units modules.",
        epilog=(
            "By default, this script attempts to write the stub files directly "
            "into the installed astropy package directory. If you lack write "
            "permissions for that location, you will get a PermissionError. "
            "In that case, use the --output-dir option to write the files to a "
            "different location."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Root directory to save the generated stub files in.\n"
            "If not provided, it defaults to the astropy.units package directory."
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose (debug) logging to see which objects are being skipped.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
