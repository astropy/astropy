"""Update `docs/whatsnew/<version>.rst` with the auto-generated release summary.

Usage:
    python scripts/update-whatsnew.py 8.0 v7.2.0

Reads the local git history, queries GitHub for the merged-PR and
closed-issue counts in the corresponding date window, and splices the
stats bullets + contributor list into the target file between sentinel
comments:

    .. release-summary-start

    (auto-generated stats land here)

    .. release-summary-end

    ...

    .. release-contributors-start

    (auto-generated contributor list lands here)

    .. release-contributors-end

The marker lines start with `..` so they render as rst comments and stay
invisible in the built HTML. Idempotent across re-runs.
"""

import argparse
import datetime as dt
import json
import os
import re
import subprocess
import urllib.request
from pathlib import Path

GH_GRAPHQL = "https://api.github.com/graphql"


def git(*args):
    return subprocess.check_output(("git", *args), text=True).strip()


def shortlog(revspec, numeric=False):
    """Return [(count_str, name), ...] from `git shortlog`, bots filtered out."""
    flag = "-sn" if numeric else "-s"
    rows = []
    for line in git("shortlog", flag, "--no-merges", revspec).splitlines():
        count, _, name = line.strip().partition("\t")
        if "[bot]" not in name:
            rows.append((count.strip(), name))
    return rows


def post_graphql(query, token):
    req = urllib.request.Request(
        GH_GRAPHQL,
        data=json.dumps({"query": query}).encode(),
        headers={"Authorization": f"Bearer {token}",
                 "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as r:
        body = json.load(r)
    if "errors" in body:
        raise RuntimeError(body["errors"])
    return body["data"]


def gh_counts(repo, since, upto, token):
    """Returns (merged_prs, closed_issues) in [since, upto] via one GraphQL request."""
    span = f"{since.date()}..{upto.date()}"
    data = post_graphql(f"""
    {{
      pulls:  search(query: "repo:{repo} is:pr    is:merged merged:{span}", type: ISSUE, first: 0) {{ issueCount }}
      issues: search(query: "repo:{repo} is:issue is:closed closed:{span}", type: ISSUE, first: 0) {{ issueCount }}
    }}
    """, token)
    return data["pulls"]["issueCount"], data["issues"]["issueCount"]


def splice(text, marker, payload):
    pattern = re.compile(
        rf"(\.\. {re.escape(marker)}-start\s*?\n).*?(\n\.\. {re.escape(marker)}-end)",
        re.DOTALL,
    )
    if not pattern.search(text):
        raise SystemExit(
            f"Could not find marker pair '.. {marker}-start' / '.. {marker}-end' in file"
        )
    return pattern.sub(lambda m: f"{m.group(1)}\n{payload}\n{m.group(2)}", text)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("version", help="Version key matching docs/whatsnew/<version>.rst (e.g. 8.0).")
    p.add_argument("prev_tag", help="Previous release tag to compare against (e.g. v7.2.0).")
    p.add_argument("--whatsnew-dir", default="docs/whatsnew",
                   help="Directory containing <version>.rst (default: docs/whatsnew).")
    p.add_argument("--project", default="astropy", help="Project name, used to default --repo and --pretty-name (default: astropy).")
    p.add_argument("--repo", help="GitHub repo OWNER/NAME (default: <project>/<project>).")
    p.add_argument("--pretty-name", help="Display name (default: <project>).")
    p.add_argument("--numeric", action="store_true", help="Sort contributors by commit count.")
    p.add_argument("--counts", action="store_true", help="Show commit count next to each name.")
    p.add_argument("--pat", help="GitHub PAT (or set GH_TOKEN / GITHUB_TOKEN).")
    args = p.parse_args()

    path = Path(args.whatsnew_dir) / f"{args.version}.rst"
    if not path.exists():
        raise SystemExit(f"{path} not found")

    name = args.pretty_name or args.project
    repo = args.repo or f"{args.project}/{args.project}"
    token = args.pat or os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        p.error("a GitHub token is required (--pat or GH_TOKEN / GITHUB_TOKEN env var)")

    current = shortlog(f"{args.prev_tag}..HEAD", numeric=args.numeric)
    previous_names = {n for _, n in shortlog(args.prev_tag)}
    current_names = {n for _, n in current}
    new = current_names - previous_names

    ncommits = int(git("rev-list", "--count", f"{args.prev_tag}..HEAD"))
    since = dt.datetime.fromisoformat(git("show", "-s", "--format=%cI", args.prev_tag))
    upto = dt.datetime.fromisoformat(git("show", "-s", "--format=%cI", "HEAD"))
    prcnt, icnt = gh_counts(repo, since, upto, token)

    short = args.prev_tag.lstrip("v")[:3]
    bullets = "\n".join(
        "  -  " + (f"{c}\t" if args.counts else "") + n + ("  *" if n in new else "")
        for c, n in current
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

Where a * indicates that this release contains their first contribution to {name}."""

    text = path.read_text()
    text = splice(text, "release-summary", stats)
    text = splice(text, "release-contributors", contributors)
    path.write_text(text)
    print(f"Updated {path}")


if __name__ == "__main__":
    main()
