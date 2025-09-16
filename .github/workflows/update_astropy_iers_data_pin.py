import sys
from pathlib import Path

import requests

metadata = requests.get(
    "https://pypi.org/pypi/astropy-iers-data/json", timeout=10
).json()

last_version = metadata["info"]["version"]
project_file = Path("pyproject.toml")

lines = project_file.read_text().splitlines()

changed_lines = 0

for iline in range(len(lines)):
    if "astropy-iers-data" in lines[iline]:
        changed_lines += 1
        lines[iline] = f'    "astropy-iers-data>={last_version}",'

if changed_lines == 0:
    print(f"No line found containing astropy-iers-data in {project_file}")
    sys.exit(1)
elif changed_lines > 1:
    print(f"More than one line found containing astropy-iers-data in {project_file}")
    sys.exit(1)

project_file.write_text("\n".join(lines) + "\n")
