# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "uv==0.9.15",
# ]
# ///

import sys
from argparse import ArgumentParser
from difflib import unified_diff
from pathlib import Path
from subprocess import run

import uv

THIS_FILE = Path(__file__)
SCRIPTS_DIR = THIS_FILE.parent
REPO_ROOT = THIS_FILE.parents[1]


def generate_tree(debug: bool) -> str:
    cp = run(
        [
            uv.find_uv_bin(),
            "tree",
            "--resolution=lowest",
            "--all-groups",
            "--universal",
        ],
        check=False,
        capture_output=True,
    )
    if debug:
        print(cp.stderr.decode(), file=sys.stderr)
    if cp.returncode != 0:
        raise AssertionError("call to uv command failed")
    return cp.stdout.decode()


def main() -> int:
    parser = ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="save the result if different than ref",
    )
    parser.add_argument(
        "--markdown",
        action="store_true",
        help="print diff as a markdown code block, if any",
    )
    parser.add_argument(
        "--quiet",
        dest="verbose",
        action="store_false",
        help="disable diff view",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="enable debugging output in stderr",
    )

    cli_args = parser.parse_args()

    ref_file = (SCRIPTS_DIR / "lowest-resolved-tree.txt").resolve()
    old_tree = ref_file.read_text()
    new_tree = generate_tree(cli_args.debug)

    diff = "\n".join(
        line.removesuffix("\n")
        for line in unified_diff(
            old_tree.splitlines(),
            new_tree.splitlines(),
            fromfile=str(ref_file),
        )
    )
    if not diff:
        if cli_args.verbose:
            print("Exact match, nothing to be done !", file=sys.stderr)
        return 0

    if cli_args.markdown:
        print("```patch")
    print(diff)
    if cli_args.markdown:
        print("```")

    if cli_args.overwrite:
        ref_file.write_text(new_tree)
        return 0

    if cli_args.verbose:
        print(
            "\nif the diff is acceptable, you may re-run this script "
            f"(`{Path(__file__).relative_to(REPO_ROOT)}`) with `--overwrite`",
            file=sys.stderr,
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
