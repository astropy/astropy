"""
Simple stub generator. Focussed on `astropy.units` for now.

Type stub files (.pyi) help type checkers understand dynamically created code.
"""

import astropy.units as units


# Types we want to include in the stub;
# to be improved upon as needed.
SUPPORTED_TYPES = {
    "IrreducibleFunctionUnit",
    "IrreducibleUnit",
    "MagUnit",
    "PrefixUnit",
    "RegularFunctionUnit",
    "Unit",
}

IMPORTS = """
# This stub file was automatically generated.

from astropy.units.core import IrreducibleUnit, PrefixUnit, Unit 
from astropy.units.function.logarithmic import MagUnit
from astropy.units.function.mixin import IrreducibleFunctionUnit, RegularFunctionUnit
"""


def generate_stub():
    stub_lines = [IMPORTS.lstrip()]
    for name in units.__all__:
        obj = getattr(units, name)
        type_name = type(obj).__name__
        if type_name in SUPPORTED_TYPES:
            docstring = obj.__doc__
            stub_lines.append(f'{name}: {type_name}\n"""{docstring}"""\n')
    return "\n".join(stub_lines)


if __name__ == "__main__":
    stub_content = generate_stub()
    with open("astropy/units/__init__.pyi", "w") as f:
        f.write(stub_content)
