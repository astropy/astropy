# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Update a ``docs/whatsnew/<version>.rst`` file with the release summary.

Usage:
    python scripts/update-latest-whatsnew-stats.py --pat=$GITHUB_PAT
    python scripts/update-latest-whatsnew-stats.py docs/whatsnew/8.0.rst --pat=$GITHUB_PAT

With no path, the latest ``docs/whatsnew/<major>.<minor>.rst`` page present in
the checkout is used, so the script also works as-is on a release branch (where
that page is the one being released). A path may be passed explicitly to target
a different page.

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
import sys
import urllib.request
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

REPO = "astropy/astropy"
WHATSNEW_DIR = Path("docs", "whatsnew")
GH_GRAPHQL = "https://api.github.com/graphql"
VERSION_RE = re.compile(r"^\d+\.\d+$")
RELEASE_TAG_RE = re.compile(r"^v\d+\.\d+\.\d+$")


@dataclass(slots=True, frozen=True)
class Error:
    message: str


def git(*args: str) -> str:
    return subprocess.check_output(("git", *args), text=True).strip()


def commit_date(ref: str) -> dt.datetime:
    """Committer date of ``ref`` as a datetime, peeling annotated tags."""
    return dt.datetime.fromisoformat(
        git("show", "-s", "--format=%cI", f"{ref}^{{commit}}")
    )


def latest_release_tag() -> str | Error:
    """The highest final release tag (vX.Y.Z, no pre-release suffix) in the repo."""
    for tag in git("tag", "--list", "v*", "--sort=-version:refname").splitlines():
        if RELEASE_TAG_RE.match(tag):
            return tag
    return Error("no release tag (vX.Y.Z) found")


def shortlog(revspec: str) -> list[tuple[str, str]]:
    """Return [(count_str, name), ...] from `git shortlog`, bots filtered out."""
    rows = []
    for line in git("shortlog", "-s", "--no-merges", revspec).splitlines():
        count, _, name = line.strip().partition("\t")
        if "[bot]" not in name:
            rows.append((count.strip(), name))
    return rows


def whatsnew_pages(dirpath: Path) -> list[tuple[tuple[int, ...], Path]]:
    """[(version_tuple, Path), ...] of <major>.<minor>.rst pages, sorted ascending."""
    pages = [
        (tuple(int(n) for n in p.stem.split(".")), p)
        for p in dirpath.glob("*.rst")
        if VERSION_RE.match(p.stem)
    ]
    return sorted(pages)


def latest_page(dirpath: Path) -> Path | Error:
    """The highest-version whatsnew page in ``dirpath``."""
    pages = whatsnew_pages(dirpath)
    if not pages:
        return Error(f"no <major>.<minor> whatsnew pages found in {dirpath}")
    return pages[-1][1]


def previous_version(path: Path) -> str | Error:
    """The whatsnew version immediately before ``path`` in the same directory."""
    versions = [v for v, _ in whatsnew_pages(path.parent)]
    target = tuple(int(n) for n in path.stem.split("."))
    if target not in versions:
        return Error(f"{path.name} is not a recognised <major>.<minor> whatsnew page")
    idx = versions.index(target)
    if idx == 0:
        return Error(f"{path.name} has no preceding whatsnew page to compare against")
    return ".".join(str(n) for n in versions[idx - 1])


def post_graphql(query, token):
    req = urllib.request.Request(
        GH_GRAPHQL,
        data=json.dumps({"query": query}).encode(),
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        body = json.load(r)
    if "errors" in body:
        raise RuntimeError(body["errors"])
    return body["data"]


def gh_counts(since, upto, token):
    """Returns (merged_prs, closed_issues) in [since, upto] via one GraphQL request."""
    span = f"{since.date()}..{upto.date()}"
    data = post_graphql(
        f"""
    {{
      pulls:  search(query: "repo:{REPO} is:pr    is:merged base:main merged:{span}", type: ISSUE, first: 1) {{ issueCount }}
      issues: search(query: "repo:{REPO} is:issue is:closed closed:{span}", type: ISSUE, first: 1) {{ issueCount }}
    }}
    """,
        token,
    )
    return data["pulls"]["issueCount"], data["issues"]["issueCount"]


def splice(text: str, marker: str, payload: str) -> str | Error:
    pattern = re.compile(
        rf"(\.\. {re.escape(marker)}-start[^\n]*\n).*?(^\.\. {re.escape(marker)}-end)",
        re.DOTALL | re.MULTILINE,
    )
    if not pattern.search(text):
        return Error(
            f"Could not find marker pair '.. {marker}-start' / '.. {marker}-end' in file"
        )
    return pattern.sub(lambda m: f"{m.group(1)}\n{payload}\n\n{m.group(2)}", text)


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "path",
        nargs="?",
        type=Path,
        help="Path to the whatsnew page to update "
        f"(default: the latest page in {WHATSNEW_DIR}).",
    )
    p.add_argument(
        "--pat", help="GitHub personal access token (or set GH_TOKEN / GITHUB_TOKEN)."
    )
    p.add_argument(
        "--check",
        action="store_true",
        help="Smoke test: count from the latest existing release tag instead of the "
        "page's previous release tag (which may be unreleased on main). Used by CI, "
        "which shows the resulting diff but never commits it.",
    )
    args = p.parse_args(argv)

    if args.path:
        path = Path(args.path)
        if not path.exists():
            print(f"{path} not found", file=sys.stderr)
            return 1
    else:
        if isinstance(path := latest_page(WHATSNEW_DIR), Error):
            print(path.message, file=sys.stderr)
            return 1

    token = args.pat or os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        p.error("a GitHub token is required (--pat or GH_TOKEN / GITHUB_TOKEN env var)")

    if isinstance(prev_version := previous_version(path), Error):
        print(prev_version.message, file=sys.stderr)
        return 1

    if args.check:
        if isinstance(prev_tag := latest_release_tag(), Error):
            print(prev_tag.message, file=sys.stderr)
            return 1
        short = prev_tag
    else:
        prev_tag = f"v{prev_version}.0"
        short = f"v{prev_version}"

    current = shortlog(f"{prev_tag}..HEAD")
    previous_names = {n for _, n in shortlog(prev_tag)}
    current_names = {n for _, n in current}
    new = current_names - previous_names

    ncommits = int(git("rev-list", "--count", f"{prev_tag}..HEAD"))
    since = commit_date(prev_tag)
    upto = commit_date("HEAD")
    prcnt, icnt = gh_counts(since, upto, token)

    bullets = "\n".join(f"  -  {n}" + ("  *" if n in new else "") for _, n in current)

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
    for marker, payload in (
        ("release-summary", stats),
        ("release-contributors", contributors),
    ):
        if isinstance(text := splice(text, marker, payload), Error):
            print(text.message, file=sys.stderr)
            return 1
    path.write_text(text)
    print(f"Updated {path} (since {prev_tag})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
