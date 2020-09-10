.. _testhelpers:

*********************
Astropy Testing Tools
*********************

This section is primarily a reference for developers that want to understand or
add to the Astropy testing machinery. See :doc:`/development/testguide` for an
overview of running or writing the tests.


`astropy.tests.helper` Module
=============================
To ease development of tests that work with Astropy, the
`astropy.tests.helper` module provides some utility functions to make
tests that use Astropy conventions or classes easier to work with, e.g.,
functions to test for near-equality of `~astropy.units.Quantity` objects.

The functionality here is not exhaustive, because
much of the useful tools are either in the standard
library, `pytest <https://docs.pytest.org>`_, or `numpy.testing
<https://numpy.org/doc/stable/reference/routines.testing.html>`_.  This module
contains primarily functionality specific to the astropy core package or
packages that follow the Astropy package template.


Reference/API
-------------

.. module:: astropy.tests.helper

.. automodapi:: astropy.tests.helper
    :no-main-docstr:
    :no-inheritance-diagram:


Astropy Test Runner
===================

When executing tests with either `astropy.test` or ``python setup.py test`` the
call to pytest is controlled by the `astropy.tests.runner.TestRunner` class.

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
