.. _testhelpers:

************************************************
Astropy Testing helpers (`astropy.tests.helper`)
************************************************

To easy development of tests that work with Astropy, the
`astropy.tests.helper` module provides some utility functions to make
test that use astropy conventions or classes easier to work with. E.g.,
functions to test for near-equality of `~astropy.units.Quantity` objects.

The functionality here is not exhaustive, because much of the useful
tools are either in the standard library, py.test, or
`numpy.testing <http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.
This module contains primarily functionality specific to astropy or
affiliated packages that follow the package template.

Note that this module also includes the ``remote_data`` marker, which
is used to indicate a test accesses data from the internet while running.
It is intended to be used as::

    from astropy.tests.helper import remote_data

    @remote_data
    def test_something():
        """
        This test dowloads something from the intenet
        """
        # test code here ...


Reference/API
=============

.. module:: astropy.tests.helper

.. automodapi:: astropy.tests.helper
    :no-main-docstr:
    :no-inheritance-diagram:
