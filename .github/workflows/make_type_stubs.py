"""
Simple stub generator. Focussed on `astropy.units` for now.

Stub files (.pyi) help type checkers understand dynamically created code.

TO DO: functions, functools._lru_cache_wrapper, method, DataInfoMeta.
"""

import astropy.units as units


def main():
    with open("astropy/units/__init__.pyi", "w") as f:
        f.write("# This stub file was automatically generated.\n")
        for line in generate_stub_lines():
            f.write("\n")
            f.write(line)


def generate_stub_lines():
    # --- preprocessing ---
    import_lines = []
    variable_lines = []
    for name in units.__all__:
        obj = getattr(units, name)
        # instances of unit classes
        if isinstance(obj, (units.UnitBase, units.FunctionUnitBase)):
            type_name = type(obj).__name__
            docstring = f'"""{obj.__doc__}"""\n' if has_own_docstring(obj) else ""
            variable_lines.append(f"{name}: {type_name}\n{docstring}")
            # if the object type is not part of astropy.units, we need to import it
            if not getattr(units, type(obj).__name__, None):
                import_lines.append(make_import_line(type(obj)))
        # classes
        elif type(obj).__name__ in ("type", "ABCMeta", "_UnitMetaClass"):
            # "as"-suffix needed to make type checkers find it
            import_lines.append(make_import_line(obj, include_as_suffix=True))

    import_lines = list(set(import_lines))
    import_lines.sort()

    # --- output ---
    for line in import_lines:
        yield line
    yield ""
    for line in variable_lines:
        yield line


def make_import_line(cls, include_as_suffix=False):
    module_name = cls.__module__
    class_name = cls.__name__
    as_suffix = f" as {class_name}" if include_as_suffix else ""
    return f"from {module_name} import {class_name}{as_suffix}"


def has_own_docstring(obj):
    """check if object has its own docstring vs an inherited one"""
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


if __name__ == "__main__":
    main()
