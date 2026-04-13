# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
import subprocess
import sys
import unicodedata

MAIN_SECTION = """
Core Package Contributors
=========================
"""

OTHER_SECTION = """
Other Credits
=============
"""

# These manual additions are for people who contributed code at the start
# of the project which are not captured by git. This should likely not
# need to be updated again.

MANUAL_ADDITIONS = (
    "Paul Barrett",
    "Chris Hanley",
    "JC Hsu",
    "James Taylor",
)


def remove_accents(text: str, /) -> str:
    normalized = unicodedata.normalize("NFD", text)
    return "".join(c for c in normalized if unicodedata.category(c) != "Mn")


def sorting_key(x: str, /) -> str:
    return remove_accents(x.lower())


# Get the full list of contributors - this takes into account the .mailmap
names = subprocess.check_output(["git", "shortlog", "-sne", "HEAD"]).splitlines()

# Extract the names from each line, and sort in a case-insensitive way
names = sorted(
    [name.decode("utf-8").split("\t")[1].split("<")[0].strip() for name in names]
    + list(MANUAL_ADDITIONS),
    key=sorting_key,
)

# Open existing credits file and search for section headings - crash if we
# can't find them

with open("docs/credits.rst") as f:
    existing = f.read()

if MAIN_SECTION not in existing:
    print(
        'Could not find expected "Core Package Contributors" section heading',
        file=sys.stderr,
    )

if OTHER_SECTION not in existing:
    print('Could not find expected "Other Credits" section heading', file=sys.stderr)

header, footer = existing.split(MAIN_SECTION)
header = header.rstrip("\n") + "\n"
footer = footer.split(OTHER_SECTION)[1].lstrip("\n")

dot = "\\."
credits = "\n".join(
    f"* {name.replace('.', dot)}"
    for name in names
    if "[bot]" not in name and name != "github-actions"
)

with open("docs/credits.rst", "w") as f:
    f.write(header)
    f.write(MAIN_SECTION)
    f.write("\n")
    f.write(credits)
    f.write("\n")
    f.write(OTHER_SECTION)
    f.write(footer)
