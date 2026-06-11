# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Update a ``docs/whatsnew/<version>.rst`` file with the release summary.

Usage:
    python scripts/update_whatsnew_stats.py docs/whatsnew/8.0.rst --pat=$GITHUB_PAT

The target version is taken from the file name, and the previous release is
taken from the whatsnew page immediately before it in the same directory (so
for ``8.0.rst`` the previous page is ``7.2.rst`` and the comparison tag is
``v7.2.0``).

The script reads the local git history, queries GitHub for the merged-PR and
closed-issue counts in the corresponding date window, and splices the stats
bullets + contributor list into the target file between sentinel comments:

    .. release-summary-start

    (auto-generated stats land here)

    .. release-summary-end

    ...

    .. release-contributors-start

    (auto-generated contributor list lands here)

    .. release-contributors-end

The marker lines start with `..` so they render as rst comments and stay
invisible in the built HTML. Idempotent across re-runs.

A GitHub token is read from ``--pat`` or the ``GH_TOKEN`` / ``GITHUB_TOKEN``
environment variable.
"""

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import urllib.request
from pathlib import Path

REPO = "astropy/astropy"
GH_GRAPHQL = "https://api.github.com/graphql"
VERSION_RE = re.compile(r"^\d+\.\d+$")


def git(*args):
    return subprocess.check_output(("git", *args), text=True).strip()


def shortlog(revspec):
    """Return [(count_str, name), ...] from `git shortlog`, bots filtered out."""
    rows = []
    for line in git("shortlog", "-s", "--no-merges", revspec).splitlines():
        count, _, name = line.strip().partition("\t")
        if "[bot]" not in name:
            rows.append((count.strip(), name))
    return rows


def previous_version(path):
    """The whatsnew version immediately before ``path`` in the same directory."""
    versions = sorted(
        tuple(int(n) for n in p.stem.split("."))
        for p in path.parent.glob("*.rst")
        if VERSION_RE.match(p.stem)
    )
    target = tuple(int(n) for n in path.stem.split("."))
    if target not in versions:
        raise SystemExit(f"{path.name} is not a recognised <major>.<minor> whatsnew page")
    idx = versions.index(target)
    if idx == 0:
        raise SystemExit(f"{path.name} has no preceding whatsnew page to compare against")
    return ".".join(str(n) for n in versions[idx - 1])


def post_graphql(query, token):
    req = urllib.request.Request(
        GH_GRAPHQL,
        data=json.dumps({"query": query}).encode(),
        headers={"Authorization": f"Bearer {token}",
                 "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        body = json.load(r)
    if "errors" in body:
        raise RuntimeError(body["errors"])
    return body["data"]


def gh_counts(since, upto, token):
    """Returns (merged_prs, closed_issues) in [since, upto] via one GraphQL request."""
    span = f"{since.date()}..{upto.date()}"
    data = post_graphql(f"""
    {{
      pulls:  search(query: "repo:{REPO} is:pr    is:merged merged:{span}", type: ISSUE, first: 1) {{ issueCount }}
      issues: search(query: "repo:{REPO} is:issue is:closed closed:{span}", type: ISSUE, first: 1) {{ issueCount }}
    }}
    """, token)
    return data["pulls"]["issueCount"], data["issues"]["issueCount"]


def splice(text, marker, payload):
    pattern = re.compile(
        rf"(\.\. {re.escape(marker)}-start[^\n]*\n).*?(^\.\. {re.escape(marker)}-end)",
        re.DOTALL | re.MULTILINE,
    )
    if not pattern.search(text):
        raise SystemExit(
            f"Could not find marker pair '.. {marker}-start' / '.. {marker}-end' in file"
        )
    return pattern.sub(lambda m: f"{m.group(1)}\n{payload}\n{m.group(2)}", text)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path", help="Path to the whatsnew page to update (e.g. docs/whatsnew/8.0.rst).")
    p.add_argument("--pat", help="GitHub personal access token (or set GH_TOKEN / GITHUB_TOKEN).")
    args = p.parse_args()

    path = Path(args.path)
    if not path.exists():
        raise SystemExit(f"{path} not found")

    token = args.pat or os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        p.error("a GitHub token is required (--pat or GH_TOKEN / GITHUB_TOKEN env var)")

    prev_version = previous_version(path)
    prev_tag = f"v{prev_version}.0"

    current = shortlog(f"{prev_tag}..HEAD")
    previous_names = {n for _, n in shortlog(prev_tag)}
    current_names = {n for _, n in current}
    new = current_names - previous_names

    ncommits = int(git("rev-list", "--count", f"{prev_tag}..HEAD"))
    since = dt.datetime.fromisoformat(git("show", "-s", "--format=%cI", prev_tag))
    upto = dt.datetime.fromisoformat(git("show", "-s", "--format=%cI", "HEAD"))
    prcnt, icnt = gh_counts(since, upto, token)

    short = f"v{prev_version}"
    bullets = "\n".join(
        f"  -  {n}" + ("  *" if n in new else "") for _, n in current
    )

    stats = f"""\
* {ncommits} commits have been added since {short}
* {icnt} issues have been closed since {short}
* {prcnt} pull requests have been merged since {short}
* {len(current_names)} people have contributed since {short}
* {len(new)} of which are new contributors"""

    contributors = f"""\
The people who have contributed to the code for this release are:

.. hlist::
  :columns: 4

{bullets}

Where a * indicates that this release contains their first contribution to astropy."""

    text = path.read_text()
    text = splice(text, "release-summary", stats)
    text = splice(text, "release-contributors", contributors)
    path.write_text(text)
    print(f"Updated {path} (since {prev_tag})")


if __name__ == "__main__":
    main()
