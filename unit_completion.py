from astropy import units as u

with open("astropy/units/autocompletion.py", "w") as f:
    f.write("from .core import UnitBase\n")
    f.write("from .function import FunctionUnitBase\n")
    for attr in dir(u):
        value = getattr(u, attr)
        if isinstance(value, u.UnitBase):
            f.write(f"{attr} : UnitBase # unit declare for auto completion\n")
        if isinstance(value, u.FunctionUnitBase):
            f.write(f"{attr} : FunctionUnitBase # unit declare for auto completion\n")
