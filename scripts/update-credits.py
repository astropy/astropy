# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
import re
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

MANUAL_ADDITIONS = frozenset(
    {
        "Paul Barrett",
        "Chris Hanley",
        "JC Hsu",
        "James Taylor",
    }
)

# Names that are manually excluded

MANUAL_EXCLUSIONS = frozenset(
    {
        "Matthias Bussonnier",  # per request
    }
)


def remove_accents(text: str, /) -> str:
    normalized = unicodedata.normalize("NFD", text)
    return "".join(c for c in normalized if unicodedata.category(c) != "Mn")


def sorting_key(x: str, /) -> str:
    return remove_accents(x.lower())


SHORTLOG_LINE_REGEXP = re.compile(r"^\s*\d+\s+(?P<name>.+) <(?P<email>.+@.+)>$")


# Get the full list of contributors - this takes into account the .mailmap
def names_from_logs() -> set[str]:
    cp = subprocess.run(
        ["git", "shortlog", "-sne", "HEAD"], check=True, capture_output=True
    )
    if cp.returncode != 0:
        print("Failed to parse git logs", file=sys.stderr)
        raise SystemExit(1)
    lines = cp.stdout.decode().splitlines()
    matches = [SHORTLOG_LINE_REGEXP.match(L) for L in lines]
    return {m.group("name") for m in matches if m is not None}


# Extract the names from each line, and sort in a case-insensitive way
names = sorted(
    (names_from_logs() - MANUAL_EXCLUSIONS) | MANUAL_ADDITIONS,
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
