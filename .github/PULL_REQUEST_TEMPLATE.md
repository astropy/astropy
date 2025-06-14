<!-- These comments are hidden when you submit the pull request,
so you do not need to remove them! -->

<!-- Please be sure to check out our contributing guidelines,
https://github.com/astropy/astropy/blob/main/CONTRIBUTING.md .
Please be sure to check out our code of conduct,
https://github.com/astropy/astropy/blob/main/CODE_OF_CONDUCT.md . -->

<!-- If you are new or need to be re-acquainted with Astropy
contributing workflow, please see
https://docs.astropy.org/en/latest/development/quickstart.html .
There is even a practical example at
https://docs.astropy.org/en/latest/development/git_edit_workflow_examples.html . -->

<!-- Please just have a quick search on GitHub to see if a similar
pull request has already been posted.
We have old closed pull requests that might provide useful code or ideas
that directly tie in with your pull request. -->

<!-- We have several automatic features that run when a pull request is open.
They can appear daunting but do not worry because maintainers will help
you navigate them, if necessary. -->

### Description
solves the issue of inserting a column of all NaNs or None will cause a crashout now it is handled through changing the max and min to nanmin and nanmax so the np.nan values are ignored in case there is numbers and np.nan and it will set the col.min and col.max to np.nan in case the column is full of NaNs or None then the result is handled so if there is a number and a NaN values it will store the max and min number in the column in that form shown in the code else it will let the lim_vals be [NaN/Na] (depending on one of the contrubutors expected output)

Fixes #17910

<!-- Optional opt-out -->

- [ ] By checking this box, the PR author has requested that maintainers do **NOT** use the "Squash and Merge" button. Maintainers should respect this when possible; however, the final decision is at the discretion of the maintainer that merges the PR.
