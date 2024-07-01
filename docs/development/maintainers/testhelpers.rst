.. _testhelpers:

*********************
Astropy Testing Tools
*********************

This section is primarily a reference for developers that want to understand or
add to the Astropy testing machinery. See :doc:`/development/testguide` for an
overview of running or writing the tests.

Details
=======
The dependencies used by the Astropy test runner are provided by a separate
package called |pytest-astropy| (see :ref:`pytest-plugins`). This package provides the ``pytest``
dependency itself, in addition to several ``pytest`` plugins that are used by
Astropy, and will also be of general use to other packages.

Since the testing dependencies are not actually required to install or use
Astropy, in the ``pyproject.toml`` file they are not included under the
``[project]`` section in ``dependencies``. Instead, they are listed under the
``[project.optional-dependences]`` in named sections such as ``test`` or ``dev_all``.

In particular the ``test`` dependencies are the minimal set of dependencies for running
astropy tests and you would use this primarily to check that tests pass *without* the
optional dependencies. This is not common and would normally be done with ``tox -e
test``.

`astropy.tests.helper` Module
=============================

To ease development of tests that work with Astropy, the
`astropy.tests.helper` module provides some utility functions to make
tests that use Astropy conventions or classes easier to work with, e.g.,
functions to test for near-equality of `~astropy.units.Quantity` objects.

The functionality here is not exhaustive, because
much of the useful tools are either in the standard
library, |pytest|, or `numpy.testing
<https://numpy.org/doc/stable/reference/routines.testing.html>`_.  This module
contains primarily functionality specific to the astropy core package or
packages that follow the Astropy package template.

Conversion Guide
----------------

Some long-standing functionality has been removed.
The following table maps them to what you should use instead.

========================================================== ===============================================
Removed                                                    Use this
========================================================== ===============================================
``astropy.io.ascii.tests.common.raises``                   ``pytest.raises``
``astropy.tests.disable_internet``                         ``pytest_remotedata.disable_internet``
``astropy.tests.helper.catch_warnings``                    ``pytest.warns``
``astropy.tests.helper.enable_deprecations_as_exceptions`` https://docs.pytest.org/en/stable/warnings.html
``astropy.tests.helper.ignore_warnings``                   https://docs.pytest.org/en/stable/warnings.html
``astropy.tests.helper.raises``                            ``pytest.raises``
``astropy.tests.helper.remote_data``                       ``pytest.mark.remote_data``
``astropy.tests.helper.treat_deprecations_as_exceptions``  https://docs.pytest.org/en/stable/warnings.html
``astropy.tests.plugins.display``                          ``pytest-astropy-header`` package
========================================================== ===============================================

Reference/API
-------------

.. module:: astropy.tests.helper

.. automodapi:: astropy.tests.helper
    :no-main-docstr:
    :no-inheritance-diagram:


Astropy Test Runner
===================

When executing tests with `astropy.test` the call to pytest is controlled
by the `astropy.tests.runner.TestRunner` class.

The `~astropy.tests.runner.TestRunner` class is used to generate the
`astropy.test` function, the test function generates a set of command line
arguments to pytest. The arguments to pytest are defined in the
`~astropy.tests.runner.TestRunner.run_tests` method, the arguments to
``run_tests`` and their respective logic are defined in methods of
`~astropy.tests.runner.TestRunner` decorated with the
`~astropy.tests.runner.keyword` decorator. For an example of this see
`~astropy.tests.runner.TestRunnerBase`. This design makes it easy for
packages to add or remove keyword arguments to their test runners, or define a
whole new set of arguments by subclassing from
`~astropy.tests.runner.TestRunnerBase`.

Reference/API
-------------

.. module:: astropy.tests.runner

.. automodapi:: astropy.tests.runner
    :no-main-docstr:
