"""
Simple stub generator. Focussed on `astropy.units` for now.

Stub files (.pyi) help type checkers understand dynamically created code.

TO DO: functions, functools._lru_cache_wrapper, method, DataInfoMeta
"""

import astropy.units as units


def main():
    with open("astropy/units/__init__.pyi", "w") as f:
        for line in generate_stub_lines():
            f.write(line)
            f.write("\n")


def generate_stub_lines():
    # --- preprocessing ---
    import_lines = []
    variable_lines = []
    for name in units.__all__:
        obj = getattr(units, name)
        # instances of unit classes
        if isinstance(obj, (units.UnitBase, units.FunctionUnitBase)):
            type_name = type(obj).__name__
            docstring = obj.__doc__
            variable_lines.append(f'{name}: {type_name}\n"""{docstring}"""\n')
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
    yield "# This stub file was automatically generated."
    yield ""
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


if __name__ == "__main__":
    main()
