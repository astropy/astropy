Contributing to Astropy
=======================

Reporting Issues
----------------

When opening an issue to report a problem, please try and provide a minimal
code example that reproduces the issue, and also include details of the
operating system, and the Python, Numpy, and Astropy versions you are using.

Contributing code
-----------------

So you're interested in contributing code to Astropy? Excellent!

Most contributions to Astropy are done via pull requests from GitHub users'
forks of the [astropy repository](https://github.com/astropy/astropy). If you're new to this style of development,
you'll want to read over our [development workflow](http://docs.astropy.org/en/latest/development/workflow/development_workflow.html).

Once you open a pull request (which should be opened against the ``master``
branch, not against any of the other branches), please make sure that you
include the following:

- **Code**: the code you are adding, which should follow as much as possible
  our [coding guidelines](http://docs.astropy.org/en/latest/development/codeguide.html).

- **Tests**: these are usually tests that ensures that code that previously
  failed now works (regression tests) or tests that cover as much as possible
  of the new functionality to make sure it doesn't break in future, and also
  returns consistent results on all platforms (since we run these tests on many
  platforms/configurations). For more information about how to write tests, see
  our [testing guidelines](http://docs.astropy.org/en/latest/development/testguide.html).

- **Documentation**: if you are adding new functionality, be sure to include a
  description in the main documentation (in ``docs/``). Again, we have some
  detailed [documentation guidelines](http://docs.astropy.org/en/latest/development/docguide.html)
  to help you out.

- **Changelog entry**: whether you are fixing a bug or adding new
  functionality, you should add an entry to the ``CHANGES.rst`` file that
  includes if possible the issue number (if you are opening a pull request you
  may not know this yet, but you can add it once the pull request is open). You
  do not need to include a changelog entry for fixes to bugs introduced in the
  developer version and which are not present in the stable releases.  In
  general you do not need to include a changelog entry for minor documentation
  or test updates.  Only user-visible changes (new features/API changes, fixed
  issues) need to be mentioned.  If in doubt ask the core maintainer reviewing
  your changes.

Other Tips
----------

- When contributing trivial documentation fixes (i.e. fixes to typos, spelling,
  grammar) that do not contain any special markup and are not associated with code
  changes, include the string "[skip ci]" at the end of your commit message.
  For example:

      $ git commit -m "Fixed typo [skip ci]"

  This will prevent automated tests for running against your change, freeing
  up resources for testing non-trivial changes.

  - If you already made the commit without including this string, you can edit
    your existing commit message by running:

        $ git commit --amend

Checklist for Contributed Code
------------------------------

A pull request for a new feature will be reviewed to see if it meets the following requirements.  For any pull request, an astropy maintainer can help to make sure that the pull request meets the requirements for inclusion in the package.  

**Scientific Quality**
(when applicable)
  * Is the submission relevant to astronomy? 
  * Are references included to the origin source for the algorithm?
  * Does the code perform as expected?
  * Has the code been tested against previously existing implementations?

**Code Quality**
  * Are the [coding guidelines](http://docs.astropy.org/en/latest/development/codeguide.html)
    followed?
  * Is the code compatible with Python 2.6, 2.7, as well as >=3.3?
  * Are there dependancies other than the Astropy core, the Python Standard 
    Library, and NumPy 1.6.0 or later?
    * Is the package importable even if the C-extensions are not built?
    * Are additional dependancies handled appropriately?
    * Do functions that require additional dependancies  raise an `ImportError`
        if they are not present?
  
**Testing**
  * Are the [testing guidelines](http://docs.astropy.org/en/latest/development/testguide.html) followed?
  * Are the inputs to the functions sufficiently tested?
  * Are their tests for any exceptions raised?
  * Are there tests for the expected performance?
  * Are the sources for the tests documented?
  * Have tests that require an [optional dependancy marked](http://docs.astropy.org/en/latest/development/testguide.html#tests-requiring-optional-dependencies) as such?
  * Does python setup.py test run without failures?

**Documentation**
  * Are the [documentation guidelines](http://docs.astropy.org/en/latest/development/docguide.html) followed? 
  * Is there a [docstring](http://docs.astropy.org/en/latest/development/docrules.html) in the function describing:
    * What the code does?
    * The format of the inputs to the function?
    * The format of the outpus of the function?
    * References to the original algorithms?
    * Any exceptions which are raised?
    * An example of running the code?
  * Is there any information needed to be added to the docs to describe the function?
  * Does the documentation build without errors or warnings?

**License**
  * Is the astropy license included at the top of the file?
  * Are there any conflicts with this code and existing codes? 

**astropy requirements**
  * Do all the tests pass on the travis build?
  * If applicable, has an entry been added into the changelog?
  * Can you checkout the pull request and repeat the examples and tests?
