<!-- This comments are hidden when you submit the pull request,
so you do not need to remove them! -->

<!-- Please be sure to check out our contributing guidelines,
https://github.com/astropy/astropy/blob/main/CONTRIBUTING.md .
Please be sure to check out our code of conduct,
https://github.com/astropy/astropy/blob/main/CODE_OF_CONDUCT.md . -->

<!-- If you are new or need to be re-acquainted with Astropy
contributing workflow, please see
http://docs.astropy.org/en/latest/development/workflow/development_workflow.html .
There is even a practical example at
https://docs.astropy.org/en/latest/development/workflow/git_edit_workflow_examples.html#astropy-fix-example . -->

<!-- Astropy coding style guidelines can be found here:
https://docs.astropy.org/en/latest/development/codeguide.html#coding-style-conventions
Our testing infrastructure enforces to follow a subset of the PEP8 to be
followed. You can check locally whether your changes have followed these by
running the following command:

tox -e codestyle

-->

<!-- Please just have a quick search on GitHub to see if a similar
pull request has already been posted.
We have old closed pull requests that might provide useful code or ideas
that directly tie in with your pull request. -->

<!-- We have several automatic features that run when a pull request is open.
They can appear daunting but do not worry because maintainers will help
you navigate them, if necessary. -->

### Description
<!-- Provide a general description of what your pull request does.
Complete the following sentence and add relevant details as you see fit. -->

<!-- In addition please ensure that the pull request title is descriptive
and allows maintainers to infer the applicable subpackage(s). -->

<!-- READ THIS FOR MANUAL BACKPORT FROM A MAINTAINER:
Apply "skip-basebranch-check" label **before** you open the PR! -->

This pull request is to address ...

<!-- If the pull request closes any open issues you can add this.
If you replace <Issue Number> with a number, GitHub will automatically link it.
If this pull request is unrelated to any issues, please remove
the following line. -->

Fixes #<Issue Number>

### Checklist for package maintainer(s)
<!-- This section is to be filled by package maintainer(s) who will
review this pull request. -->

This checklist is meant to remind the package maintainer(s) who will review this pull request of some common things to look for. This list is not exhaustive.

- [ ] Do the proposed changes actually accomplish desired goals?
- [ ] Do the proposed changes follow the [Astropy coding guidelines](https://docs.astropy.org/en/latest/development/codeguide.html)?
- [ ] Do the proposed changes need the code-style [fixed](https://docs.astropy.org/en/latest/development/workflow/maintainer_workflow.html#pre-commit_bot)?
- [ ] Are tests added/updated as required? If so, do they follow the [Astropy testing guidelines](https://docs.astropy.org/en/latest/development/testguide.html)?
- [ ] Are docs added/updated as required? If so, do they follow the [Astropy documentation guidelines](https://docs.astropy.org/en/latest/development/docguide.html#astropy-documentation-rules-and-guidelines)?
- [ ] Is rebase and/or squash necessary? If so, please provide the author with appropriate instructions. Also see ["When to rebase and squash commits"](https://docs.astropy.org/en/latest/development/when_to_rebase.html).
- [ ] Did the CI pass? If no, are the failures related? If you need to run daily and weekly cron jobs as part of the PR, please apply the `Extra CI` label.
- [ ] Is a change log needed? If yes, did the change log check pass? If no, add the `no-changelog-entry-needed` label. If this is a manual backport, use the `skip-changelog-checks` label unless special changelog handling is necessary.
- [ ] Is a milestone set? Milestone must be set but `astropy-bot` check might be missing; do not let the green checkmark fool you.
- [ ] At the time of adding the milestone, if the milestone set requires a backport to release branch(es), apply the appropriate `backport-X.Y.x` label(s) *before* merge.
