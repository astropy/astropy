.. doctest-skip-all

.. include:: workflow/known_projects.inc
.. _testing-guidelines:

==================
Testing Guidelines
==================

This section describes the testing framework and format standards for tests in
Astropy core packages (this also serves as recommendations for affiliated
packages).

Testing Framework
=================

The testing framework used by Astropy is the `py.test`_ framework.

.. _py.test: http://pytest.org/latest/
.. _pytest.main: http://pytest.org/latest/builtin.html#pytest.main

.. _running-tests:

Running Tests
=============

There are currently three different ways to invoke Astropy tests. Each
method invokes `py.test`_ to run the tests but offers different options when
calling.

In addition to running the Astropy tests, these methods can also be called
so that they check Python source code for `PEP8 compliance
<http://www.python.org/dev/peps/pep-0008/>`_. All of the PEP8 testing
options require the `pytest-pep8 plugin
<http://pypi.python.org/pypi/pytest-pep8>`_, which must be installed
separately.

setup.py test
-------------

The safest way to run the astropy test suite is via the setup command ``test``.
This is invoked by running ``python setup.py test`` while in the astropy source
code directory. Run ``python setup.py test --help`` to see the options to the
test command.

Turn on PEP8 checking by passing ``--pep8`` to the ``test`` command. This will
turn off regular testing and enable PEP8 testing.

.. note::

    This method of running the tests defaults to the version of `py.test`_
    that is bundled with Astropy. To use the locally-installed version, you
    can set the ``ASTROPY_USE_SYSTEM_PYTEST`` environment variable, eg.::

        > ASTROPY_USE_SYSTEM_PYTEST=1 python setup.py test

py.test
-------

An alternative way to run tests from the command line is to switch to the source
code directory of astropy and simply type::

    py.test

`py.test`_ will look for files that `look like tests
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
in the current directory and all recursive directories then run all the code that
`looks like tests
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
within those files.

.. note::
    To test any compiled C/Cython extensions, you must run ``python setup.py
    develop`` prior to running the py.test command-line script.  Otherwise,
    any tests that make use of these extensions will not succeed.
    Similarly, in python 3, these tests will not run correctly in the source
    code, because they need the `2to3
    <https://docs.python.org/2/library/2to3.html>`_ tool to be run on them.

You may specify a specific test file or directory at the command line::

    py.test test_file.py

To run a specific test within a file use the ``-k`` option::

    py.test test_file.py -k "test_function"

You may also use the ``-k`` option to not run tests py putting a ``-`` in front
of the matching string::

    py.test test_file.py -k "-test_function"

py.test has a number of `command line usage options.
<http://pytest.org/latest/usage.html>`_

Turn on PEP8 testing by adding the ``--pep8`` flag to the `py.test`_ call. By
default regular tests will also be run but these can be turned off by adding
``-k pep8``::

  py.test some_dir --pep8 -k pep8

.. note::
    This method of running the tests uses the locally-installed version of
    `py.test`_ rather than the bundled one, and hence will fail if the local
    version it is not up-to-date enough (`py.test`_ 2.2 as of this writing).

.. _astropy.test():

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

By default `astropy.test()`_ will skip tests which retrieve data from the
internet. To turn these tests on use the ``remote_data`` flag::

    astropy.test('io.fits', remote_data=True)

In addition, the ``test`` function supports any of the options that can be
passed to `pytest.main() <http://pytest.org/latest/builtin.html#pytest.main>`_,
and convenience options ``verbose=`` and ``pastebin=``.

Enable PEP8 compliance testing with ``pep8=True`` in the call to
``astropy.test``. This will enable PEP8 checking and disable regular tests.

.. note::
    This method of running the tests defaults to the version of
    `py.test`_ that is bundled with Astropy. To use the locally-installed
    version, you should set the ``ASTROPY_USE_SYSTEM_PYTEST`` environment
    variable (see :doc:`/config/index`) or the `py.test`_ method described
    above.

Tox
---

`Tox <http://tox.readthedocs.org>`_ is a sort of meta-test runner for Python.
It installs a project into one or more virtualenvs (usually one for each Python
version supported), build and installs the project into each virtualenv, and
runs the projects tests (or any other build processes one might want to test).
This is a good way to run the tests against multiple installed Python versions
locally without pushing to a continuous integration system.

Tox works by detecting the presence of a file called ``tox.ini`` in the root of
a Python project and using that to configure the desired virtualenvs and start
the tests.  So to run the Astropy tests on multiple Python versions using tox,
simply install Tox::

    $ pip install tox

and then from the root of an Astropy repository clone run::

    $ tox

The Astropy tox configuration currently tests against Python versions 2.6, 2.7,
3.2, and 3.3.  Tox will automatically skip any Python versions you do not have
installed, but best results are achieved if you first install all supported
Python versions and make sure they are on your ``$PATH``.

.. note::

    Tox creates its virtualenvs in the root of your project under a ``.tox``
    directory (which is automatically ignored by ``.gitignore``).  It's worth
    making note of this, however, as it is common practice to sometimes clean
    up a git repository and delete any untracked files by running the ``git
    clean -dfx`` command.  As it can take a long time to rebuild the tox
    virtualenvs you may want to exclude the ``.tox`` directory from any
    cleanup.  This can be achieved by running ``git clean -dfx -e .tox``,
    though it is probably worth defining a `git alias
    <https://git.wiki.kernel.org/index.php/Aliases>`_ to do this.

Test-running options
--------------------

Running parts of the test suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to run only the tests for a particular subpackage.  For
example, to run only the ``wcs`` tests from the commandline::

    python setup.py test -P wcs

Or from Python::

    >>> import astropy
    >>> astropy.test(package="wcs")

You can also specify a single file to test from the commandline::

    python setup.py test -t astropy/wcs/tests/test_wcs.py

When the ``-t`` option is given a relative path, it is relative to the
installed root of astropy.  When ``-t`` is given a relative path to a
documentation ``.rst`` file to test, it is relative to the root of the
documentation, i.e. the ``docs`` directory in the source tree.  For
example::

    python setup.py test -t units/index.rst

Testing for open files
^^^^^^^^^^^^^^^^^^^^^^

Astropy can test whether any of the unit tests inadvertently leave any
files open.  Since this greatly slows down the time it takes to run
the tests, it is turned off by default.

To use it from the commandline, do::

    python setup.py test --open-files

To use it from Python, do::

    >>> import astropy
    >>> astropy.test(open_files=True)

Test coverage reports
^^^^^^^^^^^^^^^^^^^^^

Astropy can use `coverage.py <http://nedbatchelder.com/code/coverage/>`_ to
generate test coverage reports.  To generate a test coverage report, use::

    python setup.py test --coverage

There is a `coveragerc
<http://nedbatchelder.com/code/coverage/config.html>`_ file that
defines files to omit as well as lines to exclude.  It is installed
along with astropy so that the ``astropy`` testing framework can use
it.  In the source tree, it is at ``astropy/tests/coveragerc``.

Running tests in parallel
^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to speed up astropy's tests using the `pytest-xdist
<https://pypi.python.org/pypi/pytest-xdist>`_ plugin.  This plugin can be
installed using `pip`_::

    pip install pytest-xdist

Once installed, tests can be run in parallel using the ``'--parallel'``
commandline option.  For example, to use 4 processes::

    python setup.py test --parallel=4

Pass a negative number to ``'--parallel'`` to create the same number of
processes as cores on your machine.

Similarly, this feature can be invoked from Python::

    >>> import astropy
    >>> astropy.test(parallel=4)

Writing tests
=============

``py.test`` has the following test discovery rules:

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

Consult the `test discovery rules
<http://pytest.org/latest/goodpractises.html#conventions-for-python-test-discovery>`_
for detailed information on how to name files and tests so that they are
automatically discovered by `py.test`_.

Simple example
--------------

The following example shows a simple function and a test to test this
function::

    def func(x):
        """Add one to the argument."""
        return x + 1

    def test_answer():
        """Check the return value of func() for an example argument."""
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

Where to put tests
------------------

Package-specific tests
^^^^^^^^^^^^^^^^^^^^^^

Each package should include a suite of unit tests, covering as many of
the public methods/functions as possible. These tests should be
included inside each sub-package, e.g::

    astropy/io/fits/tests/

``tests`` directories should contain an ``__init__.py`` file so that
the tests can be imported and so that they can use relative imports.

Interoperability tests
^^^^^^^^^^^^^^^^^^^^^^

Tests involving two or more sub-packages should be included in::

    astropy/tests/

Regression tests
----------------

Any time a bug is fixed, and wherever possible, one or more regression tests
should be added to ensure that the bug is not introduced in future. Regression
tests should include the ticket URL where the bug was reported.

Working with data files
-----------------------

Tests that need to make use of a data file should use the
`~astropy.utils.data.get_pkg_data_fileobj` or
`~astropy.utils.data.get_pkg_data_filename` functions.  These functions
search locally first, and then on the astropy data server or an arbitrary
URL, and return a file-like object or a local filename, respectively.  They
automatically cache the data locally if remote data is obtained, and from
then on the local copy will be used transparently.

They also support the use of an MD5 hash to get a specific version of a data
file.  This hash can be obtained prior to submitting a file to the astropy
data server by using the `~astropy.utils.data.compute_hash` function on a
local copy of the file.

Tests that may retrieve remote data should be marked with the
``@remote_data`` decorator, or, if a doctest, flagged with the
``REMOTE_DATA`` flag.  Tests marked in this way will be skipped by default
by ``astropy.test()`` to prevent test runs from taking too long. These
tests can be run by ``astropy.test()`` by adding the
``remote_data=True`` flag.  Turn on the remote data tests at the
command line with ``py.test --remote-data``.

Examples
^^^^^^^^
.. code-block:: none

    from ...config import get_data_filename
    from ...tests.helper import remote_data

    def test_1():
        """Test version using a local file."""
        #if filename.fits is a local file in the source distribution
        datafile = get_data_filename('filename.fits')
        # do the test

    @remote_data
    def test_2():
        """Test version using a remote file."""
        #this is the hash for a particular version of a file stored on the
        #astropy data server.
        datafile = get_data_filename('hash/94935ac31d585f68041c08f87d1a19d4')
        # do the test

    def doctest_example():
        """
        >>> datafile = get_data_filename('hash/94935')  # doctest: +REMOTE_DATA
        """
        pass

The ``get_remote_test_data`` will place the files in a temporary directory
indicated by the ``tempfile`` module, so that the test files will eventually
get removed by the system. In the long term, once test data files become too
large, we will need to design a mechanism for removing test data immediately.

Tests that create files
-----------------------

Tests may often be run from directories where users do not have write
permissions so tests which create files should always do so in
temporary directories. This can be done with the `py.test tmpdir
function argument <http://pytest.org/latest/tmpdir.html>`_ or with
Python's built-in `tempfile module
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
        """Initialize the value of NUM."""
        module.NUM = 11

    def add_num(x):
        """Add pre-defined NUM to the argument."""
        return x + NUM

    def test_42():
        """Ensure that add_num() adds the correct NUM to its argument."""
        added = add_num(42)
        assert added == 53

We can use this for example to download a remote test data file and have all
the functions in the file access it::

    import os

    def setup_module(module):
        """Store a copy of the remote test file."""
        module.DATAFILE = get_remote_test_data('94935ac31d585f68041c08f87d1a19d4')

    def test():
        """Perform test using cached remote input file."""
        f = open(DATAFILE, 'rb')
        # do the test

    def teardown_module(module):
        """Clean up remote test file copy."""
        os.remove(DATAFILE)

Class-level setup/teardown
^^^^^^^^^^^^^^^^^^^^^^^^^^

Tests can be organized into classes that have their own setup/teardown
functions. In the following ::

    def add_nums(x, y):
        """Add two numbers."""
        return x + y

    class TestAdd42(object):
        """Test for add_nums with y=42."""

        def setup_class(self):
            self.NUM = 42

        def test_1(self):
            """Test behaviour for a specific input value."""
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
            """Test behaviour for another input value."""
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_class(self):
            pass

In the above example, the ``setup_class`` method is called first, then all the
tests in the class, and finally the ``teardown_class`` is called.

Method-level setup/teardown
^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are cases where one might want setup and teardown methods to be run
before and after *each* test. For this, use the ``setup_method`` and
``teardown_method`` methods::

    def add_nums(x, y):
        """Add two numbers."""
        return x + y

    class TestAdd42(object):
        """Test for add_nums with y=42."""

        def setup_method(self, method):
            self.NUM = 42

        def test_1(self):
        """Test behaviour for a specific input value."""
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
        """Test behaviour for another input value."""
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_method(self, method):
            pass

Function-level setup/teardown
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, one can use ``setup_function`` and ``teardown_function`` to define a
setup/teardown mechanism to be run before and after each function in a module.
These take one argument, which is the function being tested::

    def setup_function(function):
        pass

    def test_1(self):
       """First test."""
        # do test

    def test_2(self):
        """Second test."""
        # do test

    def teardown_method(function):
        pass

Parametrizing tests
-------------------

If you want to run a test several times for slightly different values, then
it can be advantageous to use the ``py.test`` option to parametrize tests.
For example, instead of writing::

    def test1():
        assert type('a') == str

    def test2():
        assert type('b') == str

    def test3():
        assert type('c') == str

You can use the ``parametrize`` decorator to loop over the different
inputs::

    @pytest.mark.parametrize(('letter'), ['a', 'b', 'c'])
    def test(letter):
        """Check that the input is a string."""
        assert type(letter) == str

Tests requiring optional dependencies
-------------------------------------

For tests that test functions or methods that require optional
dependencies (e.g. Scipy), pytest should be instructed to skip the
test if the dependencies are not present. The following example shows
how this should be done::

    import pytest

    try:
        import scipy
        HAS_SCIPY = True
    except ImportError:
        HAS_SCIPY = False

    @pytest.mark.skipif('not HAS_SCIPY')
    def test_that_uses_scipy():
        ...

In this way, the test is run if Scipy is present, and skipped if
not. No tests should fail simply because an optional dependency is not
present.

Using py.test helper functions
------------------------------

If your tests need to use `py.test helper functions
<http://pytest.org/latest/builtin.html#pytest-helpers>`_, such as
``pytest.raises``, import ``pytest`` into your test module like so::

    from ...tests.helper import pytest

You may need to adjust the relative import to work for the depth of
your module.  ``tests.helper`` imports ``pytest`` either from the
user's system or ``extern.pytest`` if the user does not have py.test
installed. This is so that users need not install py.test to run
AstroPy's tests.


Testing warnings
----------------

In order to test that warnings are triggered as expected in certain
situations, you can use the ``astropy.tests.helper.catch_warnings``
context manager.  Unlike the `warnings.catch_warnings` context manager
in the standard library, this one will reset all warning state before
hand so one is assured to get the warnings reported, regardless of
what errors may have been emitted by other tests previously.  Here is
a real-world example::

  from astropy.tests.helper import catch_warnings

  with catch_warnings(MergeConflictWarning) as warning_lines:
      # Test code which triggers a MergeConflictWarning
      out = table.vstack([t1, t2, t4], join_type='outer')

      assert warning_lines[0].category == metadata.MergeConflictWarning
      assert ("In merged column 'a' the 'units' attribute does not match (cm != m)"
              in str(warning_lines[0].message))

.. note::

   Within `py.test`_ there is also the option of using the ``recwarn``
   function argument to test that warnings are triggered.  This method has
   been found to be problematic in at least one case (`pull request 1174
   <https://github.com/astropy/astropy/pull/1174#issuecomment-20249309>`_)
   so the ``astropy.tests.helper.catch_warnings`` context manager is
   preferred.

Testing with Unicode literals
-----------------------------

Python 2 can run code in two modes: by default, string literals are
8-bit `bytes` objects.  However, when ``from __future__ import
unicode_literals`` is used, string literals are `unicode` objects.  In
order to ensure that astropy supports user code written in both
styles, the testing framework has a special feature to run a module
containing tests in both modes.  Simply add the comment::

    # TEST_UNICODE_LITERALS

anywhere in the file, and all tests in that file will be tested twice:
once in the default mode where string literals are `bytes`, and again
where string literals are `unicode`.

Marking blocks of code to exclude from coverage
-----------------------------------------------

Blocks of code may be ignored by the coverage testing by adding a
comment containing the phrase ``pragma: no cover`` to the start of the
block::

    if this_rarely_happens:  # pragma: no cover
        this_call_is_ignored()

Blocks of code that are intended to run only in Python 2.x or 3.x may
also be marked so that they will be ignored when appropriate by
``coverage.py``::

    if sys.version_info[0] >= 3:  # pragma: py3
        do_it_the_python3_way()
    else:  # pragma: py2
        do_it_the_python2_way()

Using ``six.PY3`` and ``six.PY2`` will also automatically exclude
blocks from coverage, without requiring the pragma comment::

    if six.PY3:
        do_it_the_python3_way()
    elif six.PY2:
        do_it_the_python2_way()

.. _doctests:

Writing doctests
================

A doctest in Python is a special kind of test that is embedded in a
function, class, or module's docstring, or in the narrative Sphinx
documentation, and is formatted to look like a Python interactive
session--that is, they show lines of Python code entered at a ``>>>``
prompt followed by the output that would be expected (if any) when
running that code in an interactive session.

The idea is to write usage examples in docstrings that users can enter
verbatim and check their output against the expected output to confirm that
they are using the interface properly.

Furthermore, Python includes a :mod:`doctest` module that can detect these
doctests and execute them as part of a project's automated test suite.  This
way we can automatically ensure that all doctest-like examples in our
docstrings are correct.

The Astropy test suite automatically detects and runs any doctests in
the Astropy source code or documentation, or in affiliated packages
using the Astropy test running framework. For example doctests and
detailed documentation on how to write them, see the full
:mod:`doctest` documentation.

.. note::

   Since the narrative Sphinx documentation is not installed alongside
   the astropy source code, it can only be tested by running ``python
   setup.py test``, not by ``import astropy; astropy.test()``.

Skipping doctests
-----------------

Sometimes it is necessary to write examples that look like doctests but that
are not actually executable verbatim. An example may depend on some external
conditions being fulfilled, for example. In these cases there are a few ways to
skip a doctest:

1. Next to the example add a comment like: ``# doctest: +SKIP``.  For example:

   .. code-block:: none

     >>> import os
     >>> os.listdir('.')  # doctest: +SKIP

   In the above example we want to direct the user to run ``os.listdir('.')``
   but we don't want that line to be executed as part of the doctest.

   To skip tests that require fetching remote data, use the ``REMOTE_DATA``
   flag instead.  This way they can be turned on using the
   ``--remote-data`` flag when running the tests:

   .. code-block:: none

     >>> datafile = get_data_filename('hash/94935')  # doctest: +REMOTE_DATA

2. Astropy's test framework adds support for a special ``__doctest_skip__``
   variable that can be placed at the module level of any module to list
   functions, classes, and methods in that module whose doctests should not
   be run.  That is, if it doesn't make sense to run a function's example
   usage as a doctest, the entire function can be skipped in the doctest
   collection phase.

   The value of ``__doctest_skip__`` should be a list of wildcard patterns
   for all functions/classes whose doctests should be skipped.  For example::

       __doctest_skip__ = ['myfunction', 'MyClass', 'MyClass.*']

   skips the doctests in a function called ``myfunction``, the doctest for a
   class called ``MyClass``, and all *methods* of ``MyClass``.

   Module docstrings may contain doctests as well.  To skip the module-level
   doctests include the string ``'.'`` in ``__doctest_skip__``.

   To skip all doctests in a module::

       __doctest_skip__ = ['*']

3. In the Sphinx documentation, a doctest section can be skipped by
   making it part of a ``doctest-skip`` directive::

       .. doctest-skip::

           >>> # This is a doctest that will appear in the documentation,
           >>> # but will not be executed by the testing framework.
           >>> 1 / 0  # Divide by zero, ouch!

   It is also possible to skip all doctests below a certain line using
   a ``doctest-skip-all`` comment.  Note the lack of ``::`` at the end
   of the line here::

       .. doctest-skip-all

       All doctests below here are skipped...

4. ``__doctest_requires__`` is a way to list dependencies for specific
   doctests.  It should be a dictionary mapping wildcard patterns (in the same
   format as ``__doctest_skip__``) to a list of one or more modules that should
   be *importable* in order for the tests to run.  For example, if some tests
   require the scipy module to work they will be skipped unless ``import
   scipy`` is possible.  It is also possible to use a tuple of wildcard
   patterns as a key in this dict::

            __doctest_requires__ = {('func1', 'func2'): ['scipy']}

   Having this module-level variable will require ``scipy`` to be importable
   in order to run the doctests for functions ``func1`` and ``func2`` in that
   module.

   In the Sphinx documentation, a doctest requirement can be notated with the
   ``doctest-requires`` directive::

       .. doctest-requires:: scipy

           >>> import scipy
           >>> scipy.hamming(...)


Skipping output
---------------

One of the important aspects of writing doctests is that the example output
can be accurately compared to the actual output produced when running the
test.

The doctest system compares the actual output to the example output verbatim
by default, but this not always feasible.  For example the example output may
contain the ``__repr__`` of an object which displays its id (which will change
on each run), or a test that expects an exception may output a traceback.

The simplest way to generalize the example output is to use the ellipses
``...``.  For example::

    >>> 1 / 0
    Traceback (most recent call last):
    ...
    ZeroDivisionError: integer division or modulo by zero

This doctest expects an exception with a traceback, but the text of the
traceback is skipped in the example output--only the first and last lines
of the output are checked.  See the :mod:``doctest`` documentation for
more examples of skipping output.


Handling float output
---------------------

Some doctests may produce output that contains string representations of
floating point values.  Floating point representations are often not exact and
contain roundoffs in their least significant digits.  Depending on the platform
the tests are being run on (different Python versions, different OS, etc.) the
exact number of digits shown can differ.  Because doctests work by comparing
strings this can cause such tests to fail.

To address this issue Astropy's test framework includes support for a
``FLOAT_CMP`` flag that can be used with doctests.  For example::

    >>> 1.0 / 3.0  # doctest: +FLOAT_CMP
    0.333333333333333311

When this flag is used, the expected and actual outputs are both parsed to find
any floating point values in the strings.  Those are then converted to actual
Python `float` objects and compared numerically.  This means that small
differences in representation of roundoff digits will be ignored by the
doctest.  The values are otherwise compared exactly, so more significant
(albeit possibly small) differences will still be caught by these tests.
