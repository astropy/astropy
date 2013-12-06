Contributing to Astropy
=======================

Reporting Issues
----------------

When opening an issue to report a problem, please try and provide a minimal
code example that reproduces the issue, and also include details of the
operating system, and the Python and Numpy versions you are using.

Pull requests
-------------

So you're interested in contributing code to Astropy? Excellent!

If you are new to the Astropy git workflow, be sure to read over our
[development workflow](http://docs.astropy.org/en/latest/development/workflow/development_workflow.html).

Once you open a pull request, please make sure that you include the following:

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
  may not know this yet, but you can add it once the pull request is open). The
  only time that you do not need to include a changelog entry is for changes
  that fix bugs introduced in the developer version and which are not present
  in the stable versions.
