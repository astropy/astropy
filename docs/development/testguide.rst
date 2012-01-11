.. _testing-guidelines:

==========================
Testing Guidelines (Draft)
==========================

.. warning::
    This document is currently in Draft form and is subject to change. Not all
    described functionality may be implemented.

This section describes the testing framework and format standards for tests in
Astropy core packages (this also serves as recommendations for affiliated
packages).

Testing Framework
=================

The testing framework used by Astropy is the `py.test <http://pytest.org/latest/>`_
framework.

.. _running-tests:

Running Tests
=============

There are currently three different ways to invoke Astropy tests. Each method
invokes py.test to run the tests but offers different options when calling.

In addition to running the Astropy tests, these methods can also be called so
that they check Python source code for
`PEP8 compliance <http://www.python.org/dev/peps/pep-0008/>`_. All of the PEP8
testing options require the
`pytest-pep8 plugin <http://pypi.python.org/pypi/pytest-pep8>`_, which must be
installed separately.

setup.py test
-------------

The safest way to run the astropy test suite is via the setup command ``test``.
This is invoked by running ``python setup.py test`` while in the astropy source
code directory. Run ``python setup.py test --help`` to see the options to the
test command.

Turn on PEP8 checking by passing ``--pep8`` to the ``test`` command. This will
turn off regular testing and enable PEP8 testing.

.. note::
    This method of running the tests defaults to the version of `py.test` that
    is bundled with Astropy. To use the locally-installed version, you should
    either change the `use_system_pytest` configuration option to "True"
    (see :doc:`../../configs`) or the `py.test` method describe below.

py.test
-------

An alternative way to run tests from the command line is to switch to the source
code directory of astropy and simply type::

    py.test

``py.test`` will look for files that `look like tests
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
in the currect directory and all recursive directories then run all the code that
`looks like tests
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
within those files.

.. note::
    To test any compiled C/Cython extensions, you must either have run
    ``python setup.py test`` or ``python setup.py develop`` prior to running
    the py.test command-line script.  Otherwise, any tests that make use of
    these extensions will not succeed.  Similarly, in python 3, these tests
    will not run correctly in the source code, because they need the ``2to3``
    tool to be run on them.


You may specify a specific test file or directory at the command line::

    py.test test_file.py

To run a specific test within a file use the ``-k`` option::

    py.test test_file.py -k "test_function"

You may also use the ``-k`` option to not run tests py putting a ``-`` in front
of the matching string::

    py.test test_file.py -k "-test_function"

py.test has a number of `command line usage options.
<http://pytest.org/latest/usage.html>`_

Turn on PEP8 testing by adding the ``--pep8`` flag to the ``py.test`` call. By
default regular tests will also be run but these can be turned off by adding
``-k pep8``::

  py.test some_dir --pep8 -k pep8

.. note::
    This method of running the tests uses the locally-installed version of
    `py.test` rather than the bundled one, and hence will fail if the local
    version it is not up-to-date enough (`py.test` 2.2 as of this writing).

astropy.test()
--------------

AstroPy includes a standalone version of py.test that allows to tests
to be run even if py.test is not installed. Tests can be run from within
AstroPy with::

    import astropy
    astropy.test()

This will run all the default tests for AstroPy.

Tests for a specific package can be run by specifying the package in the call
to the ``test()`` function::

    astropy.test('io.fits')

This method works only with package names that can be mapped to Astropy
directories. As an alternative you can test a specific directory or file
with the ``test_path`` option::

  astropy.test(test_path='wcs/tests/test_wcs.py')

The ``test_path`` must be specified either relative to the working directory
or absolutely.

By default ``astropy.test()`` will skip tests which retrieve data from the
internet. To turn these tests on use the ``remote_data`` flag::

    astropy.test('io.fits',remote_data=True)

In addition, the ``test`` function supports any of the options that can be
passed to `pytest.main() <http://pytest.org/latest/builtin.html#pytest.main>`_,
and convenience options ``verbose=`` and ``pastebin=``.

Enable PEP8 compliance testing with ``pep8=True`` in the call to
``astropy.test``. This will enable PEP8 checking and disable regular tests.

.. note::
    This method of running the tests defaults to the version of `py.test` that
    is bundled with Astropy. To use the locally-installed version, you should
    either change the `use_system_pytest` configuration option to "True"
    (see :doc:`../../configs`) or the `py.test` method describe above.

Regression tests
================

Any time a bug is fixed, and wherever possible, one or more regression tests
should be added to ensure that the bug is not introduced in future. Regression
tests should include the ticket URL where the bug was reported.

Where to put tests
==================

Package-specific tests
----------------------

Each package should include a suite of unit tests, covering as many of the
public methods/functions as possible. These tests should be included inside
each sub-package, either in a `tests` directory, or in a test.py file, e.g::

    astropy/io/fits/tests/

or::

    astropy/io/fits/test.py

``tests`` directories should contain an ``__init__.py`` file so that the tests
can be imported and so that they can use relative imports.

Interoperability tests
----------------------

Tests involving two or more sub-packages should be included in::

    astropy/tests/

and using::

    astropy.test()

then runs both these interoperability tests, and all the unit tests in the
sub-packages. This functionality is especially important for people who install
packages through bundles and package managers, where the original source code
for the tests is not immediately available.


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
`~astropy.config.data.get_data_fileobj` or
`~astropy.config.data.get_data_filename` functions.  These functions search
locally first, and then on the astropy data server or an arbitrary URL, and
return a file-like object or a local filename, respectively.  They automatically
cache the data locally if remote data is obtained, and from then on the local
copy will be used transparently.

They also support the use of an MD5 hash to get a specific version of a data
file.  This hash can be obtained prior to submitting a file to the astropy
data server by using the `~astropy.config.data.compute_hash` function on a
local copy of the file.

Tests that may retrieve remote data should be marked with the ``@remote_data``
decorator. Tests marked with this decorator will be skipped by default by
``astropy.test()`` to prevent test runs from taking too long. These tests can
be run by ``astropy.test()`` by adding the ``remote_data=True`` flag.
Turn on the remote data tests at the command line with
``py.test --remote-data``.

Examples
^^^^^^^^
::

    from ...config import get_data_filename
    from ...tests.helper import remote_data

    def test_1():
        #if filename.fits is a local file in the source distribution
        datafile = get_data_filename('filename.fits')
        # do the test

    @remote_data
    def test_2():
        #this is the hash for a particular version of a file stored on the
        #astropy data server.
        datafile = get_data_filename('hash/94935ac31d585f68041c08f87d1a19d4')
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
<http://pytest.org/latest/tmpdir.html>`_
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
