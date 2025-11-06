import subprocess
import unicodedata

MAIN_SECTION = """
Core Package Contributors
=========================
"""

OTHER_SECTION = """
Other Credits
=============
"""

MANUAL_ADDITIONS = [
    "Paul Barrett",
    "Chris Hanley",
    "JC Hsu",
    "James Taylor",
]


def remove_accents(text):
    normalized = unicodedata.normalize("NFD", text)
    return "".join(c for c in normalized if unicodedata.category(c) != "Mn")


def sorting_key(x):
    return remove_accents(x.lower())


# Get the full list of contributors - this takes into account the .mailmap
names = subprocess.check_output(["git", "shortlog", "-sne", "HEAD"]).splitlines()

# Extract the names from each line, and sort in a case-insensitive way
names = sorted(
    [name.decode("utf-8").split("\t")[1].split("<")[0].strip() for name in names]
    + MANUAL_ADDITIONS,
    key=sorting_key,
)

# Open existing credits file and search for section headings - crash if we
# can't find them

with open("docs/credits.rst") as f:
    existing = f.read()

if MAIN_SECTION not in existing:
    print('Could not find expected "Core Package Contributors" section heading')

if OTHER_SECTION not in existing:
    print('Could not find expected "Other Credits" section heading')

header, footer = existing.split(MAIN_SECTION)
footer = footer.split(OTHER_SECTION)[1]

with open("docs/credits.rst", "w") as f:
    f.write(header)
    f.write(MAIN_SECTION)
    f.write("\n")
    for name in names:
        if "[bot]" in name or name == "github-actions":
            continue
        name = name.replace(".", "\\.")
        f.write(f"* {name}\n")
    f.write(OTHER_SECTION)
    f.write(footer)
