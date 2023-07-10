import sys

import requests

metadata = requests.get(
    "https://pypi.org/pypi/astropy-iers-data/json", timeout=10
).json()

last_version = metadata["info"]["version"]

with open("setup.cfg") as f:
    lines = f.readlines()

changed_lines = 0

for iline in range(len(lines)):
    if "astropy-iers-data" in lines[iline]:
        changed_lines += 1
        lines[iline] = f"    astropy-iers-data>={last_version}\n"

if changed_lines == 0:
    print("No line found containing astropy-iers-data in setup.cfg")
    sys.exit(1)
elif changed_lines > 1:
    print("More than one line found containing astropy-iers-data in setup.cfg")
    sys.exit(1)


with open("setup.cfg", "w") as f:
    f.writelines(lines)
