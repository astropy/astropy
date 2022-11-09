import os

from astropy import units as u

os.system(
    " grep -v '# unit declare for auto completion' astropy/units/__init__.py | awk '/^$/{n=n RS}; /./{printf \"%s\",n; n=\"\"; print}'> ./tmp_units.py"
)
with open("./tmp_units.py", "a") as f:
    f.write(
        "\n\n"
    )  # PEP8 E305 : expected 2 blank lines after class or function definition
    for attr in dir(u):
        value = getattr(u, attr)
        if isinstance(value, u.UnitBase):
            f.write(f"{attr} : UnitBase # unit declare for auto completion\n")
        if isinstance(value, u.FunctionUnitBase):
            f.write(f"{attr} : FunctionUnitBase # unit declare for auto completion\n")
os.system("cp ./tmp_units.py astropy/units/__init__.py")
os.system("rm tmp_units.py")
