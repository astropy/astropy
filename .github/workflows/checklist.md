Thank you for your contribution to Astropy! 🌌 This checklist is meant
to remind the package maintainers who will review this pull request of
some common things to look for.

- [ ] Do the proposed changes actually accomplish desired goals?
- [ ] Do the proposed changes follow the [Astropy coding guidelines](https://docs.astropy.org/en/latest/development/codeguide.html)?
- [ ] Are tests added/updated as required? If so, do they follow the [Astropy testing guidelines](https://docs.astropy.org/en/latest/development/testguide.html)?
- [ ] Are docs added/updated as required? If so, do they follow the [Astropy documentation guidelines](https://docs.astropy.org/en/latest/development/docguide.html#astropy-documentation-rules-and-guidelines)?
- [ ] Is rebase and/or squash necessary? If so, please provide the author with appropriate instructions. Also see ["When to rebase and squash commits"](https://docs.astropy.org/en/latest/development/when_to_rebase.html).
- [ ] Did the CI pass? If no, are the failures related? If you need to run daily and weekly cron jobs as part of the PR, please apply the `Extra CI` label. Codestyle issues can be fixed by the [bot](https://docs.astropy.org/en/latest/development/workflow/development_workflow.html#pre-commit).
- [ ] Is a change log needed? If yes, did the change log check pass? If no, add the `no-changelog-entry-needed` label. If this is a manual backport, use the `skip-changelog-checks` label unless special changelog handling is necessary.
- [ ] Is this a big PR that makes a "What's new?" entry worthwhile and if so, is (1) a "what's new" entry included in this PR and (2) the "whatsnew-needed" label applied?
- [ ] Is a milestone set? Milestone must be set but `astropy-bot` check might be missing; do not let the green checkmark fool you.
- [ ] At the time of adding the milestone, if the milestone set requires a backport to release branch(es), apply the appropriate `backport-X.Y.x` label(s) *before* merge.
