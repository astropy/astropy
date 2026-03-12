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

if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  git checkout "$BRANCH"
  git pull origin "$BRANCH" --quiet
else
  git checkout -b "$BRANCH" "origin/$BRANCH"
fi

echo "✅ On branch: $BRANCH"

# ── Step 3: Read scout.md ─────────────────────────────────

echo ""
echo "📋 Reading Jules scout report..."

SCOUT_REPORT=""
if [ -f "scout_report.md" ]; then
  SCOUT_REPORT=$(cat scout_report.md)
  echo "✅ Found scout_report.md"
else
  echo "⚠️  No scout_report.md found — will use PR description only"
fi

# ── Step 4: Get PR info and extract issue number ──────────

PR_BODY=""
PR_TITLE=""
ISSUE_NUMBER=""
ISSUE_STATUS=""
ISSUE_BODY=""

if command -v gh &> /dev/null; then
  PR_JSON=$(gh pr list --head "$BRANCH" --json body,title,url 2>/dev/null || echo "[]")
  PR_TITLE=$(echo "$PR_JSON" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d[0]['title'] if d else '')" 2>/dev/null || echo "")
  PR_BODY=$(echo "$PR_JSON"  | python3 -c "import sys,json; d=json.load(sys.stdin); print(d[0]['body']  if d else '')" 2>/dev/null || echo "")

  # Try to extract issue number from PR body or branch name
  ISSUE_NUMBER=$(echo "$PR_BODY $BRANCH" | grep -oE '#[0-9]+|issues/[0-9]+' | grep -oE '[0-9]+' | head -1)

  if [ -n "$ISSUE_NUMBER" ]; then
    echo "🔍 Checking upstream issue #$ISSUE_NUMBER..."
    ISSUE_JSON=$(gh issue view "$ISSUE_NUMBER" --repo astropy/astropy --json state,title,body 2>/dev/null || echo "")

    if [ -n "$ISSUE_JSON" ]; then
      ISSUE_STATUS=$(echo "$ISSUE_JSON" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('state','unknown'))")
      ISSUE_TITLE=$(echo "$ISSUE_JSON"  | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('title',''))")
      ISSUE_BODY=$(echo "$ISSUE_JSON"   | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('body','')[:800])")
      echo "📌 Issue #$ISSUE_NUMBER is: $ISSUE_STATUS — $ISSUE_TITLE"

      if [ "$ISSUE_STATUS" = "closed" ]; then
        echo ""
        echo "⛔  Issue #$ISSUE_NUMBER is already CLOSED — someone already fixed it."
        echo "   Aborted. Pick a different Jules branch."
        exit 0
      fi
    fi

    # Check for competing open PRs
    echo "🔍 Checking for competing PRs on issue #$ISSUE_NUMBER..."
    COMPETING=$(gh pr list --repo astropy/astropy       --search "fixes #$ISSUE_NUMBER"       --json number,title,state,isDraft,url       2>/dev/null || echo "[]")

    PR_COUNT=$(echo "$COMPETING" | python3 -c "import sys,json; print(len(json.load(sys.stdin)))")

    if [ "$PR_COUNT" -gt "0" ]; then
      echo ""
      echo "⛔  STOP — There is already an open PR for issue #$ISSUE_NUMBER:"
      echo "$COMPETING" | python3 -c "
import sys, json
prs = json.load(sys.stdin)
for pr in prs:
    draft = ' [DRAFT]' if pr.get('isDraft') else ''
    print('  #' + str(pr['number']) + draft + ': ' + pr['title'])
    print('  ' + pr['url'])
"
      echo ""
      echo "   Don't waste your time — pick a different Jules branch."
      exit 0
    else
      echo "✅ No competing PRs found — you're clear to proceed."
    fi
  fi
fi

# ── Step 5: Write Claude Code context file ────────────────

CONTEXT_FILE="/tmp/jules-context-$(date +%s).md"

cat > "$CONTEXT_FILE" << CONTEXT
# Jules Handoff Context

## Branch
$BRANCH

## Jules PR Title
${PR_TITLE:-"(not found)"}

---

## Jules Scout Report (scout.md)
${SCOUT_REPORT:-"(No scout_report.md on this branch — refer to PR description below)"}

---

## Jules PR Description
${PR_BODY:-"(No PR found — check jules.google.com)"}

---

## Upstream Issue #${ISSUE_NUMBER:-"unknown"}
Status: ${ISSUE_STATUS:-"unknown"}
Title:  ${ISSUE_TITLE:-"unknown"}

${ISSUE_BODY:-"(Could not fetch issue body — check https://github.com/astropy/astropy/issues)"}

---

## Your job (Claude Code)

The scout report and issue above are your source of truth — ignore any diff noise
from main being ahead of this branch due to the daily sync automation.

1. Read scout.md carefully — Jules has already identified the files and approach.

2. Check if the issue is still open and unassigned upstream before starting:
   https://github.com/astropy/astropy/issues/${ISSUE_NUMBER:-""}

3. Run the existing tests for the relevant subpackage first to establish a baseline:
   \`python -m pytest astropy/<subpackage>/tests/ -x -v\`

4. Implement the fix described in scout.md. Follow astropy standards:
   - Type hints on all new functions
   - Numpy-style docstrings
   - Tests in the corresponding tests/ directory
   - No public API changes unless the issue explicitly requires it

5. Once tests pass, give me a 3-sentence summary of what changed for the PR description.
CONTEXT

echo "✅ Context ready"

# ── Step 6: Print context and launch Claude Code ──────────

echo ""
echo "🚀 Launching Claude Code..."
echo ""

claude < "$CONTEXT_FILE"