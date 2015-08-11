.. _testhelpers:

************************************************
Astropy Testing helpers (`astropy.tests.helper`)
************************************************

To ease development of tests that work with Astropy, the
`astropy.tests.helper` module provides some utility functions to make
tests that use Astropy conventions or classes easier to work with. e.g.,
functions to test for near-equality of `~astropy.units.Quantity` objects.

The functionality here is not exhaustive, because much of the useful
tools are either in the standard library, py.test, or
`numpy.testing <http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.
This module contains primarily functionality specific to Astropy or
affiliated packages that follow the package template.

Note that this module also includes the ``remote_data`` marker, which
is used to indicate a test accesses data from the internet while running.
It is intended to be used as::

    from astropy.tests.helper import remote_data

    @remote_data
    def test_something():
        """
        This test downloads something from the intenet
        """
        # test code here ...

This test will now only run if the test suite is invoked as
``python setup.py test --remote-data``.


Reference/API
=============

.. module:: astropy.tests.helper

.. automodapi:: astropy.tests.helper
    :no-main-docstr:
    :no-inheritance-diagram:
