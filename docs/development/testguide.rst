============================
Testing Guidelines (Draft 2)
============================

.. warning::
    This document is currently in Draft form and is subject to change. Not all
    described functionality may be implemented.

This section describes the testing framework and format standards for tests in
AstroPy core modules (this also serves as recommendations for affiliated
packages).

Testing Framework
=================

The testing framework used by AstroPy is the `py.test <http://pytest.org/latest/>`_
framework.

Running Tests
=============

Using py.test
-------------

The simplest way to run tests from the command line is to simply type::

    py.test
    
``py.test`` will look for files that `look like tests 
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_ 
in the currect directory and all recursive directories then run all the code that
`looks like tests 
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
within those files.

You may specify a specific test file or directory at the command line::

    py.test test_file.py
    
To run a specific test within a file use the ``-k`` option::

    py.test test_file.py -k "test_function"
    
py.test has a number of `command line usage options. 
<http://pytest.org/latest/usage.html>`_

Using astropy.test()
--------------------

AstroPy includes a standalone version of py.test that allows to tests
to be run even if py.test is not installed. Tests can be run from within 
AstroPy with::

    import astropy
    astropy.test()
    
This will run all the default tests for AstroPy.

Tests for a specific module can be run by specifying the module in the call
to the ``test()`` function::

    astropy.test('io.fits')
    
In addition, the ``test`` function supports any of the options that can be
passed to `pytest.main() <http://pytest.org/latest/builtin.html#pytest.main>`_,
and convenience options ``verbose=`` and ``pastebin=``.

Using data in tests
===================

Tests can include very small datafiles, but any files significantly larger
than the source code should be placed on a remote server. The base URL for the
test files will be::

    http://data.astropy.org/

and files will be accessed by their MD5 hash, for example::

    http://data.astropy.org/94935ac31d585f68041c08f87d1a19d4

Tests then retrieve data via this URL. This implicitly allows versioning,
since different versions of data files will have different hashes. Old data
files should not be removed, so that tests can be run in any version of
AstroPy.

The details of the server implementation have yet to be decided, but using
these static hash-based URLs ensures that even if we change the backend, the
URL will remain the same.

Where to put tests
==================

Module-specific tests
---------------------

Each module should include a suite of unit tests, covering as many of the
public methods/functions as possible. These tests should be included inside
each sub-module, either in a `tests` directory, or in a test.py file, e.g::

    astropy/io/fits/tests/

or::

    astropy/io/fits/test.py

Interoperability tests
----------------------

Tests involving two or more sub-modules should be included in::

    astropy/tests/

and using::

    astropy.test()

then runs both these interoperability tests, and all the unit tests in the
sub-modules. This functionality is especially important for people who install 
packages through bundles and package managers, where the original source code 
for the tests is not immediately available.

Regression tests
================

Any time a bug is fixed, and wherever possible, one or more regression tests
should be added to ensure that the bug is not introduced in future. Regression
tests should include the ticket URL where the bug was reported.

Writing tests
=============

``py.test`` has the following test discovery rules:

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

Consult the `test discovery rules
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
for detailed information on how to name files and tests so that they are 
automatically discovered by ``py.test``.

Simple example
--------------

The following example shows a simple function and a test to test this
function::

    def func(x):
        return x + 1

    def test_answer():
        assert func(3) == 5

If we place this in a ``test.py`` file and then run::

    py.test test.py

The result is::

    ============================= test session starts ==============================
    python: platform darwin -- Python 2.7.2 -- pytest-1.1.1
    test object 1: /Users/tom/tmp/test.py

    test.py F

    =================================== FAILURES ===================================
    _________________________________ test_answer __________________________________

        def test_answer():
    >       assert func(3) == 5
    E       assert 4 == 5
    E        +  where 4 = func(3)

    test.py:5: AssertionError
    =========================== 1 failed in 0.07 seconds ===========================

Working with data files
-----------------------

Tests that need to make use of a data file should use the
``get_local_test_data`` and ``get_remote_test_data`` functions, and test files
should be requested using filenames in the first case, and MD5 hashes in the
second. Each of these functions returns the local path to the file (and in the
case of remote data, it is the path to the downloaded file):

.. warning:: This is going to change

::

    from astropy.util.testing import get_local_test_data, \
                                     get_remote_test_data

    def test_1():
        datafile = get_local_test_data('filename.fits')
        # do the test

    def test_2():
        datafile = get_remote_test_data('94935ac31d585f68041c08f87d1a19d4')
        # do the test

The ``get_remote_test_data`` will place the files in a temporary directory
indicated by the ``tempfile`` module, so that the test files will eventually
get removed by the system. In the long term, once test data files become too
large, we will need to design a mechanism for removing test data immediately.

Tests that create files
-----------------------

Tests may often be run from directories where users do not have write permissions
so tests which create files should always do so in temporary directories. This
can be done with the `py.test tmpdir function argument
<http://pytest.org/latest/getting-started.html#going-functional-requesting-a-unique-temporary-directory>`_
or with Python's built-in `tempfile module 
<http://docs.python.org/library/tempfile.html#module-tempfile>`_.

Setting up/Tearing down tests
-----------------------------

In some cases, it can be useful to run a series of tests requiring something
to be set up first. There are four ways to do this:

Module-level setup/teardown
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the ``setup_module`` and ``teardown_module`` functions are specified in a
file, they are called before and after all the tests in the file respectively.
These functions take one argument, which is the module itself, which makes it
very easy to set module-wide variables::

    def setup_module(module):
        module.NUM = 11

    def add_num(x):
        return x + NUM

    def test_42():
        added = add_num(42)
        assert added == 53

We can use this for example to download a remote test data file and have all
the functions in the file access it::

    import os

    def setup_module(module):
        module.DATAFILE = get_remote_test_data('94935ac31d585f68041c08f87d1a19d4')

    def test():
        f = open(DATAFILE, 'rb')
        # do the test

    def teardown_module(module):
        os.remove(DATAFILE)

Class-level
^^^^^^^^^^^

Tests can be organized into classes that have their own setup/teardown
functions. In the following ::

    def add_nums(x, y):
        return x + y

    class TestAdd42(object):

        def setup_class(self):
            self.NUM = 42

        def test_1(self):
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_class(self):
            pass

In the above example, the ``setup_class`` method is called first, then all the
tests in the class, and finally the ``teardown_class`` is called.

Method-level
^^^^^^^^^^^^

There are cases where one might want setup and teardown methods to be run
before and after *each* test. For this, use the ``setup_method`` and
``teardown_method`` methods::

    def add_nums(x, y):
        return x + y

    class TestAdd42(object):

        def setup_method(self, method):
            self.NUM = 42

        def test_1(self):
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_method(self, method):
            pass

Function-level
^^^^^^^^^^^^^^

Finally, one can use ``setup_function`` and ``teardown_function`` to define a
setup/teardown mechanism to be run before and after each function in a module.
These take one argument, which is the function being tested::

    def setup_function(function):
        pass

    def test_1(self):
        # do test

    def test_2(self):
        # do test

    def teardown_method(function):
        pass
        
Using py.test helper functions
------------------------------

If your tests need to use `py.test helper functions 
<http://pytest.org/latest/builtin.html#pytest-helpers>`_, such as ``pytest.raises``,
import ``pytest`` into your test module like so::

    from ...tests.helper import pytest
    
You may need to adjust the relative import to work for the depth of your module.
``tests.helper`` imports ``pytest`` either from the user's system or ``extern.pytest``
if the user does not have py.test installed. This is so that users need not 
install py.test to run AstroPy's tests.
    