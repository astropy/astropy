#!/bin/bash
# ─────────────────────────────────────────────────────────
# jules-handoff.sh
# Pulls a Jules branch and launches Claude Code with full context
#
# Usage:
#   ./jules-handoff.sh                        # interactive branch picker
#   ./jules-handoff.sh <branch-name>          # direct branch name
#   ./jules-handoff.sh <github-pr-url>        # paste a Jules PR URL
# ─────────────────────────────────────────────────────────

set -e

REPO_DIR=$(git rev-parse --show-toplevel 2>/dev/null)
if [ -z "$REPO_DIR" ]; then
  echo "❌ Not inside a git repo. Run this from your astropy fork directory."
  exit 1
fi

cd "$REPO_DIR"

# ── Step 1: Determine the Jules branch ───────────────────

if [ -n "$1" ]; then
  INPUT="$1"

  # If it's a GitHub PR URL, extract the branch from the API
  if [[ "$INPUT" == https://github.com/* ]]; then
    echo "🔍 Fetching branch from PR URL..."
    PR_NUMBER=$(echo "$INPUT" | grep -oE '/pull/[0-9]+' | grep -oE '[0-9]+')
    REPO_SLUG=$(echo "$INPUT" | sed 's|https://github.com/||' | cut -d'/' -f1-2)
    BRANCH=$(curl -s "https://api.github.com/repos/$REPO_SLUG/pulls/$PR_NUMBER" \
      | python3 -c "import sys,json; print(json.load(sys.stdin)['head']['ref'])")
    echo "📌 Branch: $BRANCH"
  else
    BRANCH="$INPUT"
  fi

else
  # Interactive: list recent Jules branches
  echo ""
  echo "🔍 Fetching recent branches from origin..."
  git fetch origin --quiet

  echo ""
  echo "Recent branches (newest first):"
  echo "──────────────────────────────"

  # List remote branches, filter likely Jules ones, show recent first
  BRANCHES=$(git branch -r --sort=-committerdate \
    | grep 'origin/' \
    | grep -v 'origin/main' \
    | sed 's|origin/||' \
    | head -20)

  if [ -z "$BRANCHES" ]; then
    echo "No branches found other than main."
    echo "Jules may not have created a branch yet — check jules.google.com"
    exit 1
  fi

  # Number them for selection
  i=1
  while IFS= read -r branch; do
    echo "  $i) $branch"
    i=$((i+1))
  done <<< "$BRANCHES"

  echo ""
  printf "Enter number of the Jules branch to work on: "
  read -r SELECTION

  BRANCH=$(echo "$BRANCHES" | sed -n "${SELECTION}p" | xargs)

  if [ -z "$BRANCH" ]; then
    echo "❌ Invalid selection"
    exit 1
  fi
fi

# ── Step 2: Pull the branch locally ──────────────────────

echo ""
echo "📥 Checking out branch: $BRANCH"

git fetch origin "$BRANCH" --quiet

# Check if branch exists locally
if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git checkout "$BRANCH"
  git pull origin "$BRANCH" --quiet
else
  git checkout -b "$BRANCH" "origin/$BRANCH"
fi

echo "✅ On branch: $BRANCH"

# ── Step 3: Gather context for Claude Code ────────────────

echo ""
echo "📋 Gathering context..."

# Get the diff from main
DIFF=$(git diff main..."$BRANCH" --stat 2>/dev/null || echo "No diff available")

# Get commit messages on this branch
COMMITS=$(git log main.."$BRANCH" --oneline 2>/dev/null || echo "No commits yet")

# Get the branch description / PR body if available via gh CLI
PR_BODY=""
if command -v gh &> /dev/null; then
  PR_BODY=$(gh pr list --head "$BRANCH" --json body,title,url \
    --jq '.[0] | "Title: \(.title)\nURL: \(.url)\n\nDescription:\n\(.body)"' \
    2>/dev/null || echo "")
fi

# ── Step 4: Write the Claude Code context file ───────────

CONTEXT_FILE="/tmp/jules-context-$(date +%s).md"

cat > "$CONTEXT_FILE" << CONTEXT
# Jules Handoff Context

## Branch
$BRANCH

## Jules PR / Task Description
${PR_BODY:-"(No PR found — check jules.google.com for the task description)"}

## Files changed by Jules
$DIFF

## Jules commits
$COMMITS

---

## Your job (Claude Code)

Review what Jules has done on this branch. Then:

1. Run the test suite for any modified subpackages:
   \`python -m pytest astropy/<subpackage>/tests/ -x -v\`

2. If Jules left TODOs, incomplete implementations, or failing tests — fix them.

3. If Jules only scouted (no code changes), implement the fix described above.
   Follow astropy's coding standards:
   - Type hints on all new functions
   - Numpy-style docstrings
   - Tests in the corresponding tests/ directory
   - No changes to public API unless the issue explicitly requires it

4. Once tests pass, summarize what changed so we can write the PR description.

## Upstream issue
Search https://github.com/astropy/astropy/issues for the issue referenced
in the branch name or PR description above.
CONTEXT

echo "✅ Context written to $CONTEXT_FILE"

# ── Step 5: Launch Claude Code ────────────────────────────

echo ""
echo "🚀 Launching Claude Code..."
echo "──────────────────────────────────────────────"
echo ""

# Launch Claude Code with the context pre-loaded
claude --context "$CONTEXT_FILE"
CONTEXT