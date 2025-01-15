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
library, `pytest`_, or `numpy.testing
<https://numpy.org/doc/stable/reference/routines.testing.html>`_.  This module
contains primarily functionality specific to the astropy core package or
packages that follow the Astropy package template.

Conversion Guide
----------------

Some long-standing functionality has been deprecated/removed since ``astropy`` 5.1.
The following table maps them to what you should use instead.

========================================================== ===============================================
Deprecated                                                 Use this
========================================================== ===============================================
``astropy.io.ascii.tests.common.raises``                   ``pytest.raises``
``astropy.tests.helper.raises``                            ``pytest.raises``
``astropy.tests.helper.catch_warnings``                    ``pytest.warns``
``astropy.tests.helper.ignore_warnings``                   https://docs.pytest.org/en/stable/warnings.html
``astropy.tests.helper.enable_deprecations_as_exceptions`` https://docs.pytest.org/en/stable/warnings.html
``astropy.tests.helper.treat_deprecations_as_exceptions``  https://docs.pytest.org/en/stable/warnings.html
========================================================== ===============================================

========================================================== ===============================================
Removed                                                    Use this
========================================================== ===============================================
``astropy.tests.disable_internet``                         ``pytest_remotedata.disable_internet``
``astropy.tests.helper.remote_data``                       ``pytest.mark.remote_data``
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
``run_tests`` method, the arguments to
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

.. autoclass:: astropy.tests.runner.keyword
    :no-undoc-members:

.. autoclass:: astropy.tests.runner.TestRunnerBase

.. autoclass:: astropy.tests.runner.TestRunner
