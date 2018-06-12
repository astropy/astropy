Contributing to Astropy
=======================

Reporting Issues
----------------

When opening an issue to report a problem, please try and provide a minimal code
example that reproduces the issue, and also include details of the operating
system, and the Python, Numpy, and Astropy versions you are using.

Contributing code
-----------------

So you're interested in contributing code to Astropy? Excellent!

Most contributions to Astropy are done via pull requests from GitHub users'
forks of the [astropy repository](https://github.com/astropy/astropy). If you're
new to this style of development, you'll want to read over our
[development workflow](http://docs.astropy.org/en/latest/development/workflow/development_workflow.html).

You may also/instead be interested in contributing to to an
[astropy affiliated package](http://www.astropy.org/affiliated/).
Affiliated packages are astronomy-related software packages that are not a part
of the astropy core package, but build on it for more specialized applications
and follow the astropy guidelines for reuse, interoperability, and interfacing.
Each affiliated package has its own developers/maintainers and its own specific
guidelines for contributions, so be sure to read their docs.

Once you open a pull request (which should be opened against the ``master``
branch, not against any of the other branches), please make sure that you
include the following:

- **Code**: the code you are adding, which should follow as much as possible
  our [coding guidelines](http://docs.astropy.org/en/latest/development/codeguide.html).

- **Tests**: these are usually tests to ensure that code that previously
  failed now works (regression tests) or tests that cover as much as possible
  of the new functionality to make sure it doesn't break in future, and also
  returns consistent results on all platforms (since we run these tests on many
  platforms/configurations). For more information about how to write tests, see
  our [testing guidelines](http://docs.astropy.org/en/latest/development/testguide.html).

- **Documentation**: if you are adding new functionality, be sure to include a
  description in the main documentation (in ``docs/``). Again, we have some
  detailed [documentation guidelines](http://docs.astropy.org/en/latest/development/docguide.html)
  to help you out.

- **Performance improvements**: if you are making changes that impact Astropy
  performance, consider adding a performance benchmark in the
  [astropy-benchmarks](https://github.com/astropy/astropy-benchmarks)
  repository. You can find out more about how to do this
  [in the README for that repository](https://github.com/astropy/astropy-benchmarks#contributing-a-benchmark).

- **Changelog entry**: whether you are fixing a bug or adding new
  functionality, you should add an entry to the ``CHANGES.rst`` file that
  includes the PR number (if you are opening a pull request you may not know
  this yet, but you can add it once the pull request is open). If you're not
  sure where to put the changelog entry, wait at least until a maintainer
  has reviewed your PR and assigned it to a milestone.

  You do not need to include a changelog entry for fixes to bugs introduced in
  the developer version and therefore are not present in the stable releases. In
  general you do not need to include a changelog entry for minor documentation
  or test updates.  Only user-visible changes (new features/API changes, fixed
  issues) need to be mentioned.  If in doubt ask the core maintainer reviewing
  your changes.

Other Tips
----------

- To prevent the automated tests from running you can add ``[ci skip]`` to your
  commit message. This is useful if your PR is a work in progress and you are
  not yet ready for the tests to run.  For example:

      $ git commit -m "WIP widget [ci skip]"

  - If you already made the commit without including this string, you can edit
    your existing commit message by running:

        $ git commit --amend

- To skip only the AppVeyor (Windows) CI builds you can use ``[skip appveyor]``,
  and to skip testing on Travis CI use ``[skip travis]``.

- If your commit makes substantial changes to the documentation, but no code
  changes, the you can use ``[docs only]``, that will skip all but the
  documentation building jobs on Travis.

- When contributing trivial documentation fixes (i.e. fixes to typos, spelling,
  grammar) that do not contain any special markup and are not associated with
  code changes, please include the string ``[docs only]`` in your commit
  message.

      $ git commit -m "Fixed typo [docs only]"

Checklist for Contributed Code
------------------------------

A pull request for a new feature will be reviewed to see if it meets the
following requirements.  For any pull request, an astropy maintainer can help to
make sure that the pull request meets the requirements for inclusion in the
package.

**Scientific Quality** (when applicable)
  * Is the submission relevant to astronomy?
  * Are references included to the origin source for the algorithm?
  * Does the code perform as expected?
  * Has the code been tested against previously existing implementations?

**Code Quality**
  * Are the [coding guidelines](http://docs.astropy.org/en/latest/development/codeguide.html)
    followed?
  * Is the code compatible with Python >=3.5?
  * Are there dependencies other than the Astropy core, the Python Standard
    Library, and NumPy 1.10.0 or later?
    * Is the package importable even if the C-extensions are not built?
    * Are additional dependencies handled appropriately?
    * Do functions that require additional dependencies  raise an `ImportError`
        if they are not present?

**Testing**
  * Are the [testing guidelines](http://docs.astropy.org/en/latest/development/testguide.html) followed?
  * Are the inputs to the functions sufficiently tested?
  * Are there tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Have tests that require an [optional dependency marked](http://docs.astropy.org/en/latest/development/testguide.html#tests-requiring-optional-dependencies) as such?
  * Does python setup.py test run without failures?

**Documentation**
  * Are the [documentation guidelines](http://docs.astropy.org/en/latest/development/docguide.html) followed?
  * Is there a [docstring](http://docs.astropy.org/en/latest/development/docrules.html)
    in the function describing:
    * What the code does?
    * The format of the inputs of the function?
    * The format of the outputs of the function?
    * References to the original algorithms?
    * Any exceptions which are raised?
    * An example of running the code?
  * Is there any information needed to be added to the docs to describe the function?
  * Does the documentation build without errors or warnings?

**License**
  * Is the astropy license included at the top of the file?
  * Are there any conflicts with this code and existing codes?

**Astropy requirements**
  * Do all the Travis CI, AppVeyor, and CircleCI tests pass?
  * If applicable, has an entry been added into the changelog?
  * Can you checkout the pull request and repeat the examples and tests?
