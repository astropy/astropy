#!/usr/bin/env python3
"""
Script to detect merge commits in pull requests.

This script checks if the current branch contains any merge commits from the
main branch. If merge commits are found, it exits with a non-zero status to
fail the CI check.

Exit codes:
    0: No merge commits detected (success)
    1: Merge commits detected (failure)
    2: Script execution error
"""
import subprocess
import sys
import re


def run_git_command(command):
    """Run a git command and return its output."""
    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running git command: {command}", file=sys.stderr)
        print(f"Error output: {e.stderr}", file=sys.stderr)
        sys.exit(2)


def get_merge_commits(base_branch="origin/main"):
    """
    Detect merge commits in the current branch.

    Parameters
    ----------
    base_branch : str
        The base branch to compare against (default: origin/main)

    Returns
    -------
    list
        List of merge commit hashes
    """
    # Get all commits in the current branch that are not in the base branch
    # and have more than one parent (merge commits)
    command = f'git log {base_branch}..HEAD --merges --oneline'
    output = run_git_command(command)

    if not output:
        return []

    # Parse the output to get commit hashes
    merge_commits = []
    for line in output.split('\n'):
        if line.strip():
            commit_hash = line.split()[0]
            merge_commits.append(commit_hash)

    return merge_commits


def get_commit_details(commit_hash):
    """
    Get detailed information about a commit.

    Parameters
    ----------
    commit_hash : str
        The commit hash

    Returns
    -------
    dict
        Dictionary with commit details (hash, subject, parents)
    """
    # Get commit subject
    subject = run_git_command(f'git log -1 --format=%s {commit_hash}')

    # Get parent commits
    parents = run_git_command(f'git log -1 --format=%P {commit_hash}')
    parent_list = parents.split()

    return {
        'hash': commit_hash,
        'subject': subject,
        'parents': parent_list
    }


def is_merge_from_main(commit_details):
    """
    Check if a merge commit is likely from merging main/master.

    Parameters
    ----------
    commit_details : dict
        Commit details from get_commit_details()

    Returns
    -------
    bool
        True if the commit appears to be a merge from main/master
    """
    subject = commit_details['subject'].lower()

    # Check for common merge commit patterns
    merge_patterns = [
        r"merge branch ['\"]?main['\"]?",
        r"merge branch ['\"]?master['\"]?",
        r"merge.*origin/main",
        r"merge.*origin/master",
        r"merge.*upstream/main",
        r"merge.*upstream/master",
    ]

    for pattern in merge_patterns:
        if re.search(pattern, subject):
            return True

    return False


def main():
    """Main function to check for merge commits."""
    print("Checking for merge commits in pull request...")

    # Fetch the latest main branch
    print("Fetching origin/main...")
    run_git_command('git fetch origin main')

    # Get merge commits
    merge_commits = get_merge_commits()

    if not merge_commits:
        print("✓ No merge commits detected. PR is clean!")
        return 0

    # Check each merge commit
    print(f"\n✗ Found {len(merge_commits)} merge commit(s) in this PR:\n")

    has_main_merges = False
    for commit_hash in merge_commits:
        details = get_commit_details(commit_hash)
        print(f"  Commit: {details['hash']}")
        print(f"  Subject: {details['subject']}")
        print(f"  Parents: {', '.join(details['parents'])}")

        if is_merge_from_main(details):
            print("  ⚠️  This appears to be a merge from main/master branch")
            has_main_merges = True

        print()

    # Print error message and instructions
    print("=" * 70)
    print("ERROR: Merge commits detected in pull request")
    print("=" * 70)
    print()
    print("This PR contains merge commits, which are not allowed.")
    print("Please rebase your branch on top of main instead of merging.")
    print()
    print("To fix this issue:")
    print("  1. Backup your branch: git branch backup-branch")
    print("  2. Fetch the latest main: git fetch upstream main")
    print("  3. Rebase your branch: git rebase upstream/main")
    print("  4. Force push (if already pushed): git push --force-with-lease")
    print()
    print("If you need help, please ask on the astropy-dev mailing list or")
    print("consult the developer documentation:")
    print("https://docs.astropy.org/en/latest/development/")
    print()

    return 1


if __name__ == "__main__":
    sys.exit(main())
