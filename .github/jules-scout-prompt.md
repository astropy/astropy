# Daily Issue Scout Task

## Your job
You are scouting the upstream astropy repository for a good contribution
opportunity for a first-time contributor with the following profile:

**Contributor profile:**
- 10th grade student, strong Python, statistics, and physics background
- Familiar with numpy, scipy, and the basic scientific Python stack
- Has done ML research (physics-informed ML for aerodynamic optimization)
- Builds rocket flight simulators and computer vision systems as personal projects
- Can understand math/physics but not deep astronomy-specific domain knowledge
- Goal: get a PR merged upstream, learn the codebase incrementally

## Step 1 — Find a good issue

Fetch open issues from https://github.com/astropy/astropy/issues with ALL
of these labels: `good first issue`, `Effort-low`

From those results, pick the **single best one** using these criteria (in order):
1. Prefers `astropy/stats`, `astropy/units`, or `astropy/utils` subpackages
2. Avoids issues requiring deep astronomy instrument or telescope knowledge
3. Avoids issues that already have someone assigned or a linked open PR
4. Prefers bugs over feature requests (bugs have clearer acceptance criteria)
5. Prefers issues with a clear reproduction case or error message already provided

If no `Effort-low` issues meet criteria, fall back to `Effort-medium` issues
in the same preferred subpackages.

## Step 2 — Write a scoped implementation brief

For the chosen issue, write a brief with these exact sections:

### Issue
[link and title, upstream issue number]

### What's broken / what's needed
[2-3 sentences, plain English, no jargon. What is the user experiencing?
What should happen instead?]

### Files to touch
[List the specific file paths in the astropy repo that need changes.
Include both the source file AND the corresponding test file.]

### Implementation approach
[Numbered step-by-step plan, specific enough that an AI coding agent can
execute it without reading the full issue thread. Include:
- Function/method names to modify
- Line numbers if you can identify them
- Any edge cases to handle
- Whether existing tests need updating]

### Acceptance criteria
[Bulleted list — exactly how do we know the fix is correct and complete?
Include: tests pass, specific behavior is fixed, no regressions]

### Guardrails
[What NOT to change — public API surface, other subpackages, unrelated
functions. Be explicit.]

### Estimated difficulty
[1-5 scale with one sentence of reasoning]

### Claude Code kickoff command
[Write the exact prompt the contributor should paste into Claude Code
to implement this fix. It should be self-contained — include the issue
context, file paths, and acceptance criteria in a single prompt block.]

## Step 3 — Output

Create a GitHub issue in THIS repository (the fork) with:
- Title: `🔭 Daily Scout: [upstream issue title] (#[upstream issue number])`
- Body: the full implementation brief from Step 2
- Label: `jules-scout` (create the label if it doesn't exist, use color #0075ca)

Do not open any PRs. Do not modify any source files. Only create the issue.
