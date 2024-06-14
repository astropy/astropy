.. doctest-skip-all

.. _testing-guidelines:

******************
Testing Guidelines
******************

This section describes the testing framework and format standards for tests in
Astropy core and coordinated packages, and also serves as recommendations for
affiliated packages.

Testing Framework
*****************

The testing framework used by astropy (and packages using the
|OpenAstronomy Packaging Guide|) is the |pytest| framework.

.. _testing-dependencies:

Testing Dependencies
********************

The dependencies used by the Astropy test runner are provided by a separate
package called |pytest-astropy|. This package provides the ``pytest``
dependency itself, in addition to several ``pytest`` plugins that are used by
Astropy, and will also be of general use to other packages.

Since the testing dependencies are not actually required to install or use
Astropy, in the ``pyproject.toml`` file they are not included under the
``[project]`` section in ``dependencies``. Instead, they are listed under the
``[project.optional-dependences]`` section called ``test``.  Developers who want
to run the test suite will need to either install pytest-astropy directly::

    python -m pip install pytest-astropy

or install the core package in 'editable' mode specifying the ``[test]``
option::

    python -m pip install --editable ".[test]"

To test the full set of optional dependencies, use the ``test_all`` option::

    python -m pip install --editable ".[test_all]"

If you are looking to do development in Astropy beyond running the tests, e.g. building
the documentation or doing static analysis, we provide the complete set of all
dependencies with the ``dev_all`` option::

    python -m pip install --editable ".[dev_all]"

A detailed description of the |pytest-astropy| plugins can be found in the
:ref:`pytest-plugins` section.

.. _running-tests:

Running Tests
*************

There are currently three different ways to invoke Astropy tests. Each
method invokes |pytest| to run the tests but offers different options when
calling. To run the tests, you will need to make sure you have the |pytest|
package installed.

In addition to running the Astropy tests, these methods can also be called
so that they check Python source code for |PEP8|. All of the PEP8 testing
options require the `pytest-pep8 plugin
<https://pypi.org/project/pytest-pep8>`_, which must be installed
separately.


tox
===

The most robust way to run the tests (which can also be the slowest) is
to make use of `Tox <https://tox.readthedocs.io/en/latest/>`__, which is a
general purpose tool for automating Python testing. One of the benefits of tox
is that it first creates a source distribution of the package being tested, and
installs it into a new virtual environment, along with any dependencies that are
declared in the package, before running the tests. This can therefore catch
issues related to undeclared package data, or missing dependencies. Since we use
tox to run many of the tests on continuous integration services, it can also be
used in many cases to reproduce issues seen on those services.

To run the tests with tox, first make sure that tox is installed, e.g.::

    python -m pip install tox

then run the basic test suite with::

    tox -e test

or run the test suite with all optional dependencies with::

    tox -e test-alldeps

You can see a list of available test environments with::

    tox -l -v

which will also explain what each of them does.

You can also run checks or commands not directly related to tests - for instance::

    tox -e codestyle

will run checks using the ``ruff`` tool.

.. note::
    It is suggested that you automate the code-style checks using the provided
    pre-commit hook, as described in the :ref:`pre-commit` section.

It is possible to pass options to pytest when running tox - to do this, add a
``--`` after the regular tox command, and anything after this will be passed to
pytest, e.g.::

    tox -e test -- -v --pdb

This can be used in conjunction with the ``-P`` option provided by the
`pytest-filter-subpackage <https://github.com/astropy/pytest-filter-subpackage>`_
plugin to run just part of the test suite.

Note that even though ``tox`` caches information, interactive debug and test
sessions with ``tox`` can be quite slow. For this case, it may be better to
set up a virtual environment with an editable install. Here, ``tox`` can still
help by setting up a complete test environment, which one can then activate::

  tox -e test-alldeps --develop --notest
  source .tox/test-alldeps/bin/activate

Here, we use ``--notest`` to prevent ``tox`` from running the tests, since the
idea is to do that oneself -- using the ``pytest`` commands described below,
targeting the relevant sub-package or test file.

.. _running-pytest:

pytest
======

The test suite can also be run directly from the native ``pytest`` command,
which is generally faster than using tox for iterative development. In
this case, it is important for developers to be aware that they must manually
rebuild any extensions by running::

    python -m pip install --editable ".[test]"

before running the test with pytest with::

    pytest

Instead of calling ``python -m pip install --editable ".[test]"``, you can also build the
extensions with::

    python setup.py build_ext --inplace

which avoids also installing the developer version of astropy into your current
environment - however note that the ``pip`` command is required if you need to
test parts of the package that rely on certain `entry points
<https://setuptools.readthedocs.io/en/latest/pkg_resources.html#entry-points>`_
being installed.

It is possible to run only the tests for a particular subpackage or set of
subpackages.  For example, to run only the ``wcs`` tests from the
commandline::

    pytest -P wcs

Or, to run only the ``wcs`` and ``utils`` tests::

    pytest -P wcs,utils

You can also specify a single directory, a file (``.py`` python or ``.rst``
doc file), or a specific test to check, or rerun only tests that failed in
the previous run::

    pytest astropy/modeling
    pytest astropy/wcs/tests/test_wcs.py
    pytest astropy/units -k float_dtype_promotion
    pytest astropy/units/tests/test_quantity.py::TestQuantityCreation::test_float_dtype_promotion
    pytest astropy/wcs/index.rst
    pytest --last-failed

For more details, see the `pytest invocation guide
<https://docs.pytest.org/en/stable/how-to/usage.html>`_ and the
description of `caching
<https://docs.pytest.org/en/stable/how-to/cache.html>`_.


.. _astropy.test():

astropy.test()
==============

Tests can be run from an installed version of Astropy with::

    import astropy
    astropy.test()

This will run all the default tests for Astropy (but will not run the
documentation tests in the ``.rst`` documentation since those files are
not installed).

Tests for a specific package can be run by specifying the package in the call
to the ``test()`` function::

    astropy.test(package='io.fits')

This method works only with package names that can be mapped to Astropy
directories. As an alternative you can test a specific directory or file
with the ``test_path`` option::

  astropy.test(test_path='wcs/tests/test_wcs.py')

The ``test_path`` must be specified either relative to the working directory
or absolutely.

By default `astropy.test()`_ will skip tests which retrieve data from the
internet. To turn these tests on use the ``remote_data`` flag::

    astropy.test(package='io.fits', remote_data=True)

In addition, the ``test`` function supports any of the options that can be
passed to :ref:`pytest.main() <pytest:pytest.main-usage>`
and convenience options ``verbose=`` and ``pastebin=``.

Enable PEP8 compliance testing with ``pep8=True`` in the call to
``astropy.test``. This will enable PEP8 checking and disable regular tests.

Astropy Test Function
---------------------

.. autofunction:: astropy.test

Test-running options
====================

.. _open-files:

Testing for open files
----------------------

The ``filterwarnings`` settings under ``[tool.pytest.ini_options]`` in the
``pyproject.toml`` file has an option which converts all unhandled warnings to
errors during a test run. As a result, any open file(s) that throw
``ResourceWarning`` (except the specific ones already ignored) would fail the
affected test(s).

Test coverage reports
---------------------

Coverage reports can be generated using the `pytest-cov
<https://pypi.org/project/pytest-cov/>`_ plugin (which is installed
automatically when installing pytest-astropy) by using e.g.::

    pytest --cov astropy --cov-report html

There is some configuration inside the ``pyproject.toml`` file that
defines files to omit as well as lines to exclude.

Running tests in parallel
-------------------------

It is possible to speed up astropy's tests using the `pytest-xdist
<https://pypi.org/project/pytest-xdist>`_ plugin.

Once installed, tests can be run in parallel using the ``'-n'``
commandline option. For example, to use 4 processes::

    pytest -n 4

Pass ``-n auto`` to create the same number of processes as cores
on your machine.

Similarly, this feature can be invoked from ``astropy.test``::

    >>> import astropy
    >>> astropy.test(parallel=4)

.. _writing-tests:

Writing tests
*************

``pytest`` has the following test discovery rules:

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

Consult the :ref:`test discovery rules <pytest:python test discovery>`
for detailed information on how to name files and tests so that they are
automatically discovered by |pytest|.

Simple example
==============

The following example shows a simple function and a test to test this
function::

    def func(x):
        """Add one to the argument."""
        return x + 1

    def test_answer():
        """Check the return value of func() for an example argument."""
        assert func(3) == 5

If we place this in a ``test.py`` file and then run::

    pytest test.py

The result is::

    ============================= test session starts ==============================
    python: platform darwin -- Python 3.x.x -- pytest-x.x.x
    test object 1: /Users/username/tmp/test.py

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
==================

Package-specific tests
----------------------

Each package should include a suite of unit tests, covering as many of
the public methods/functions as possible. These tests should be
included inside each sub-package, e.g::

    astropy/io/fits/tests/

``tests`` directories should contain an ``__init__.py`` file so that
the tests can be imported and so that they can use relative imports.

Interoperability tests
----------------------

Tests involving two or more sub-packages should be included in::

    astropy/tests/

Regression tests
================

Any time a bug is fixed, and wherever possible, one or more regression tests
should be added to ensure that the bug is not introduced in future. Regression
tests should include the ticket URL where the bug was reported.

.. _data-files:

Working with data files
=======================

Tests that need to make use of a data file should use the
`~astropy.utils.data.get_pkg_data_fileobj` or
`~astropy.utils.data.get_pkg_data_filename` functions.  These functions
search locally first, and then on the astropy data server or an arbitrary
URL, and return a file-like object or a local filename, respectively.  They
automatically cache the data locally if remote data is obtained, and from
then on the local copy will be used transparently.  See the next section for
note specific to dealing with the cache in tests.

They also support the use of an MD5 hash to get a specific version of a data
file.  This hash can be obtained prior to submitting a file to the astropy
data server by using the `~astropy.utils.data.compute_hash` function on a
local copy of the file.

Tests that may retrieve remote data should be marked with the
``@pytest.mark.remote_data`` decorator, or, if a doctest, flagged with the
``REMOTE_DATA`` flag.  Tests marked in this way will be skipped by default by
``astropy.test()`` to prevent test runs from taking too long. These tests can
be run by ``astropy.test()`` by adding the ``remote_data='any'`` flag.  Turn on
the remote data tests at the command line with ``pytest --remote-data=any``.

It is possible to mark tests using
``@pytest.mark.remote_data(source='astropy')``, which can be used to indicate
that the only required data is from the http://data.astropy.org server. To
enable just these tests, you can run the
tests with ``pytest --remote-data=astropy``.

For more information on the ``pytest-remotedata`` plugin, see
:ref:`remotedata-plugin`.

Examples
--------
.. code-block:: python

    from ...config import get_data_filename

    def test_1():
        """Test version using a local file."""
        #if filename.fits is a local file in the source distribution
        datafile = get_data_filename('filename.fits')
        # do the test

    @pytest.mark.remote_data
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

Tests that use the file cache
-----------------------------

By default, the Astropy test runner sets up a clean file cache in a temporary
directory that is used only for that test run and then destroyed.  This is to
ensure consistency between test runs, as well as to not clutter users' caches
(i.e. the cache directory returned by `~astropy.config.get_cache_dir`) with
test files.

However, some test authors (especially for affiliated packages) may find it
desirable to cache files downloaded during a test run in a more permanent
location (e.g. for large data sets).  To this end the
`~astropy.config.set_temp_cache` helper may be used.  It can be used either as
a context manager within a test to temporarily set the cache to a custom
location, or as a *decorator* that takes effect for an entire test function
(not including setup or teardown, which would have to be decorated separately).

Furthermore, it is possible to change the location of the cache directory
for the duration of the test run by setting the ``XDG_CACHE_HOME``
environment variable.


Tests that create files
=======================

Some tests involve writing files. These files should not be saved permanently.
The :ref:`pytest 'tmp_path' fixture <pytest:tmp_path>` allows for the
convenient creation of temporary directories, which ensures test files will be
cleaned up. Temporary directories can also be helpful in the case where the
tests are run in an environment where the runner would otherwise not have write
access.


Setting up/Tearing down tests
=============================

In some cases, it can be useful to run a series of tests requiring something
to be set up first. There are four ways to do this:

Module-level setup/teardown
---------------------------

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
--------------------------

Tests can be organized into classes that have their own setup/teardown
functions. In the following::

    def add_nums(x, y):
        """Add two numbers."""
        return x + y

    class TestAdd42:
        """Test for add_nums with y=42."""

        def setup_class(self):
            self.NUM = 42

        def test_1(self):
            """Test behavior for a specific input value."""
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
            """Test behavior for another input value."""
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_class(self):
            pass

In the above example, the ``setup_class`` method is called first, then all the
tests in the class, and finally the ``teardown_class`` is called.

Method-level setup/teardown
---------------------------

There are cases where one might want setup and teardown methods to be run
before and after *each* test. For this, use the ``setup_method`` and
``teardown_method`` methods::

    def add_nums(x, y):
        """Add two numbers."""
        return x + y

    class TestAdd42:
        """Test for add_nums with y=42."""

        def setup_method(self, method):
            self.NUM = 42

        def test_1(self):
        """Test behavior for a specific input value."""
            added = add_nums(11, self.NUM)
            assert added == 53

        def test_2(self):
        """Test behavior for another input value."""
            added = add_nums(13, self.NUM)
            assert added == 55

        def teardown_method(self, method):
            pass

Function-level setup/teardown
-----------------------------

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

    def teardown_function(function):
        pass

Property-based tests
====================

`Property-based testing
<https://increment.com/testing/in-praise-of-property-based-testing/>`_
lets you focus on the parts of your test that matter, by making more
general claims - "works for any two numbers" instead of "works for 1 + 2".
Imagine if random testing gave you minimal, non-flaky failing examples,
and a clean way to describe even the most complicated data - that's
property-based testing!

``pytest-astropy`` includes a dependency on `Hypothesis
<https://hypothesis.readthedocs.io/>`_, so installation is easy -
you can just read the docs or `work through the tutorial
<https://github.com/Zac-HD/escape-from-automanual-testing/>`_
and start writing tests like::

    from astropy.coordinates import SkyCoord
    from hypothesis import given, strategies as st

    @given(
        st.builds(SkyCoord, ra=st.floats(0, 360), dec=st.floats(-90, 90))
    )
    def test_coordinate_transform(coord):
        """Test that sky coord can be translated from ICRS to Galactic and back."""
        assert coord == coord.galactic.icrs  # floating-point precision alert!

Other properties that you could test include:

- Round-tripping from image to sky coordinates and back should be lossless
  for distortion-free mappings, and otherwise always below 10^-5 px.
- Take a moment in time, round-trip it through various frames, and check it
  hasn't changed or lost precision. (or at least not by more than a nanosecond)
- IO routines losslessly round-trip data that they are expected to handle
- Optimised routines calculate the same result as unoptimised, within tolerances

This is a great way to start contributing to Astropy, and has already found
bugs in time handling. See issue `#9017 <https://github.com/astropy/astropy/issues/9017>`_
and pull request `#9532 <https://github.com/astropy/astropy/pull/9532>`_ for details!

(and if you find Hypothesis useful in your research,
`please cite it <https://doi.org/10.21105/joss.01891>`_!)


Parametrizing tests
===================

If you want to run a test several times for slightly different values,
you can use ``pytest`` to avoid writing separate tests.
For example, instead of writing::

    def test1():
        assert type('a') == str

    def test2():
        assert type('b') == str

    def test3():
        assert type('c') == str

You can use the ``@pytest.mark.parametrize`` decorator to concisely
create a test function for each input::

    @pytest.mark.parametrize(('letter'), ['a', 'b', 'c'])
    def test(letter):
        """Check that the input is a string."""
        assert type(letter) == str

As a guideline, use ``parametrize`` if you can enumerate all possible
test cases and each failure would be a distinct issue, and Hypothesis
when there are many possible inputs or you only want a single simple
failure to be reported.

Tests requiring optional dependencies
=====================================

For tests that test functions or methods that require optional dependencies
(e.g., Scipy), pytest should be instructed to skip the test if the dependencies
are not present, as the ``astropy`` tests should succeed even if an optional
dependency is not present. ``astropy`` provides a list of boolean flags that
test whether optional dependencies are installed (at import time). For example,
to load the corresponding flag for Scipy and mark a test to skip if Scipy is not
present, use::

    import pytest
    from astropy.utils.compat.optional_deps import HAS_SCIPY

    @pytest.mark.skipif(not HAS_SCIPY, reason='scipy is required')
    def test_that_uses_scipy():
        ...

These variables should exist for all of Astropy's optional dependencies; a
complete list of supported flags can be found in
``astropy.utils.compat.optional_deps``.

Any new optional dependencies should be added to that file, as well as to the
relevant entries in the ``pyproject.toml`` file in the
``[project.optional-dependencies]`` section; typically, under ``all`` for
dependencies used in user-facing code (e.g., ``h5py``, which is used to write
tables to HDF5 format), and in ``test_all`` for dependencies only used in tests
(e.g., ``skyfield``, which is used to cross-check the accuracy of coordinate
transforms).

Using pytest helper functions
=============================

If your tests need to use `pytest helper functions
<https://docs.pytest.org/en/latest/reference/reference.html#functions>`_, such as
``pytest.raises``, import ``pytest`` into your test module like so::

    import pytest

Testing warnings
================

In order to test that warnings are triggered as expected in certain
situations,
|pytest| provides its own context manager
:ref:`pytest.warns <pytest:warns>` that, completely
analogously to ``pytest.raises`` (see below) allows to probe explicitly
for specific warning classes and, through the optional ``match`` argument,
messages. Note that when no warning of the specified type is
triggered, this will make the test fail. When checking for optional,
but not mandatory warnings, ``pytest.warns()`` can be used to catch and
inspect them.

.. note::

   With |pytest| there is also the option of using the
   :ref:`recwarn <pytest:recwarn>` function argument to test that
   warnings are triggered within the entire embedding function.
   This method has been found to be problematic in at least one case
   (`pull request 1174 <https://github.com/astropy/astropy/pull/1174#issuecomment-20249309>`_).

Testing exceptions
==================

Just like the handling of warnings described above, tests that are
designed to trigger certain errors should verify that an exception of
the expected type is raised in the expected place.  This is efficiently
done by running the tested code inside the
:ref:`pytest.raises <pytest:assertraises>`
context manager.  Its optional ``match`` argument allows to check the
error message for any patterns using ``regex`` syntax.  For example the
matches ``pytest.raises(OSError, match=r'^No such file')`` and
``pytest.raises(OSError, match=r'or directory$')`` would be equivalent
to ``assert str(err).startswith(No such file)`` and ``assert
str(err).endswith(or directory)``, respectively, on the raised error
message ``err``.
For matching multi-line messages you need to pass the ``(?s)``
:ref:`flag <python:re-syntax>`
to the underlying ``re.search``, as in the example below::

  with pytest.raises(fits.VerifyError, match=r'(?s)not upper.+ Illegal key') as excinfo:
      hdu.verify('fix+exception')
  assert str(excinfo.value).count('Card') == 2

This invocation also illustrates how to get an ``ExceptionInfo`` object
returned to perform additional diagnostics on the info.

Testing configuration parameters
================================

In order to ensure reproducibility of tests, all configuration items
are reset to their default values when the test runner starts up.

Sometimes you'll want to test the behavior of code when a certain
configuration item is set to a particular value.  In that case, you
can use the `astropy.config.ConfigItem.set_temp` context manager to
temporarily set a configuration item to that value, test within that
context, and have it automatically return to its original value.

For example::

    def test_pprint():
        from ... import conf
        with conf.set_temp('max_lines', 6):
            # ...

Marking blocks of code to exclude from coverage
===============================================

Blocks of code may be ignored by the coverage testing by adding a
comment containing the phrase ``pragma: no cover`` to the start of the
block::

    if this_rarely_happens:  # pragma: no cover
        this_call_is_ignored()

.. _image-tests:

Image tests with pytest-mpl
===========================

Running image tests
-------------------

We make use of the `pytest-mpl <https://pypi.org/project/pytest-mpl>`_
plugin to write tests where we can compare the output of plotting commands
with reference files on a pixel-by-pixel basis (this is used for instance in
:ref:`astropy.visualization.wcsaxes <wcsaxes>`). We use the `hybrid mode
<https://pytest-mpl.readthedocs.io/en/latest/hybrid_mode.html>`_ with
hashes and images.

To run the Astropy tests with the image comparison, use e.g.::

    tox -e py310-test-image-mpl360-cov

However, note that the output can be sensitive to the operating system and
specific version of libraries such as freetype. In general, using tox will
result in the version of freetype being pinned, but the hashes will only be
correct when running the tests on Linux. Therefore, if using another operating
system, we do not recommend running the image tests locally and instead it is
best to rely on these running in an controlled continuous integration
environment.

Writing image tests
-------------------

The `README.rst <https://github.com/matplotlib/pytest-mpl/blob/master/README.rst>`__
for the plugin contains information on writing tests with this plugin. Once you
have added a test, and push this to a pull request, you will likely start seeing
a test failure because the figure hash is missing from the hash libraries
(see the next section for how to proceed).

Rather than use the ``@pytest.mark.mpl_image_compare`` decorator directly, you
should make use of the ``@figure_test`` convenience decorator which
sets the default tolerance and style to be consistent across the astropy core
package, and also automatically enables access to remote data::

    from astropy.tests.figures import figure_test

    @figure_test
    def test_figure():
        fig = plt.figure()
        ...
        return fig

You can optionally pass keyword arguments to ``@figure_test`` and these will be
passed on to ``mpl_image_compare``::

    @figure_test(savefig_kwargs={'bbox_inches': 'tight'})
    def test_figure():
        ...

Failing tests
-------------

When existing tests start failing, it is usually either because of a change in
astropy itself, or a change in Matplotlib. New tests will also fail if you have
not yet updated the hash library.

In all cases, you can view a webpage with all the existing figures where you can
check whether any of the figures are now wrong, or if all is well. The link to
the page for each tox environment that has been run will be provided in the
list of statuses for pull requests, and can also be found in the CircleCI
logs. If any changes/additions look good, you can download from the summary page
a JSON file with the hashes which you can use to replace the existing one in
``astropy/tests/figures``.

New hash libraries
------------------

When adding a new tox environment for image testing, such as for a new Matplotlib
or Python version, the tests will fail as the hash library does not exist yet. To
generate it, you should run the tests the first time with::

    tox -e <envname> -- --mpl-generate-hash-library=astropy/tests/figures/<envname>.json

for example::

    tox -e py310-test-image-mpl360-cov -- --mpl-generate-hash-library=astropy/tests/figures/py310-test-image-mpl360-cov.json

Then add and commit the new JSON file and try running the tests again. The tests
may fail in the continuous integration if e.g. the freetype version does not
match or if you generated the JSON file on a Mac or Windows machine - if that is
the case, follow the instructions in `Failing tests`_ to update the hashes.

As an alternative to generating the JSON file above, you can also simply copy a
previous version of the JSON file and update any failing hashes as described
in `Failing tests`_.

Generating reference images
---------------------------

You do not need to generate reference images for new tests or updated reference
images for changed tests - when pull requests are merged, a CircleCI job will automatically
update the reference images in the `astropy-figure-tests <https://github.com/astropy/astropy-figure-tests>`_
repository.

.. _doctests:

Writing doctests
****************

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

The Astropy test suite automatically detects and runs any doctests in the
astropy source code or documentation, or in packages using the Astropy test
running framework. For example doctests and detailed documentation on how to
write them, see the full :mod:`doctest` documentation.

.. note::

   Since the narrative Sphinx documentation is not installed alongside the
   astropy source code, it can only be tested by running ``pytest`` directly (or
   via tox), not by ``import astropy; astropy.test()``.

For more information on the ``pytest-doctestplus`` plugin used by Astropy, see
:ref:`doctestplus-plugin`.

.. _skipping-doctests:

Skipping doctests
=================

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
===============

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
of the output are checked.  See the :mod:`doctest` documentation for
more examples of skipping output.

Ignoring all output
-------------------

Another possibility for ignoring output is to use the
``# doctest: +IGNORE_OUTPUT`` flag.  This allows a doctest to execute (and
check that the code executes without errors), but allows the entire output
to be ignored in cases where we don't care what the output is.  This differs
from using ellipses in that we can still provide complete example output, just
without the test checking that it is exactly right.  For example::

    >>> print('Hello world')  # doctest: +IGNORE_OUTPUT
    We don't really care what the output is as long as there were no errors...

.. _handling-float-output:

Handling float output
=====================

Some doctests may produce output that contains string representations of
floating point values.  Floating point representations are often not exact and
contain roundoffs in their least significant digits.  Depending on the platform
the tests are being run on (different Python versions, different OS, etc.) the
exact number of digits shown can differ.  Because doctests work by comparing
strings this can cause such tests to fail.

To address this issue, the ``pytest-doctestplus`` plugin provides support for a
``FLOAT_CMP`` flag that can be used with doctests.  For example:

.. code-block:: none

  >>> 1.0 / 3.0  # doctest: +FLOAT_CMP
  0.333333333333333311

When this flag is used, the expected and actual outputs are both parsed to find
any floating point values in the strings.  Those are then converted to actual
Python `float` objects and compared numerically.  This means that small
differences in representation of roundoff digits will be ignored by the
doctest.  The values are otherwise compared exactly, so more significant
(albeit possibly small) differences will still be caught by these tests.

Continuous integration
**********************

Overview
========

Astropy uses the following continuous integration (CI) services:

* `GitHub Actions <https://github.com/astropy/astropy/actions>`_ for
  Linux, OS X, and Windows setups
  (Note: GitHub Actions does not have "allowed failures" yet, so you might
  see a fail job reported for your PR with "(Allowed Failure)" in its name.
  Still, some failures might be real and related to your changes, so check
  it anyway!)
* `CircleCI <https://circleci.com>`_ for visualization tests

These continuously test the package for each commit and pull request that is
pushed to GitHub to notice when something breaks.

In some cases, you may see failures on continuous integration services that
you do not see locally, for example because the operating system is different,
or because the failure happens with only 32-bit Python.

Maintainers have the option to run :ref:`comparative benchmark <benchmarks>` using GitHub Actions
to test a new pull request against the current ``main`` branch. It uses the benchmarks
from `astropy-benchmarks <https://github.com/astropy/astropy-benchmarks/>`_.
It is important to note that these benchmarks can be flaky as they run on
virtual machines (and thus shared hardware) but they should give a general
idea of the performance impact of a pull request.

.. _pytest-plugins:

Pytest Plugins
**************

The following ``pytest`` plugins are maintained and used by Astropy. They are
included as dependencies to the ``pytest-astropy`` package, which is now
required for testing Astropy. More information on all of the  plugins provided
by the ``pytest-astropy`` package (including dependencies not maintained by
Astropy) can be found `here <https://github.com/astropy/pytest-astropy>`__.

.. _remotedata-plugin:

pytest-remotedata
=================

The |pytest-remotedata| plugin allows developers to control whether to run
tests that access data from the internet. The plugin provides two decorators
that can be used to mark individual test functions or entire test classes:

* ``@pytest.mark.remote_data`` for tests that require data from the internet
* ``@pytest.mark.internet_off`` for tests that should run only when there is no
  internet access. This is useful for testing local data caches or fallbacks
  for when no network access is available.

The plugin also adds the ``--remote-data`` option to the ``pytest`` command
(which is also made available through the Astropy test runner).

If the ``--remote-data`` option is not provided when running the test suite, or
if ``--remote-data=none`` is provided, all tests that are marked with
``remote_data`` will be skipped. All tests that are marked with
``internet_off`` will be executed. Any test that attempts to access the
internet but is not marked with ``remote_data`` will result in a failure.

Providing either the ``--remote-data`` option, or ``--remote-data=any``, will
cause all tests marked with ``remote_data`` to be executed. Any tests that are
marked with ``internet_off`` will be skipped.

Running the tests with ``--remote-data=astropy`` will cause only tests that
receive remote data from Astropy data sources to be run. Tests with any other
data sources will be skipped. This is indicated in the test code by marking
test functions with ``@pytest.mark.remote_data(source='astropy')``. Tests
marked with ``internet_off`` will also be skipped in this case.

Also see :ref:`data-files`.

.. _doctestplus-plugin:

pytest-doctestplus
==================

The |pytest-doctestplus| plugin provides advanced doctest features, including:

* handling doctests that use remote data in conjunction with the
  ``pytest-remotedata`` plugin above (see :ref:`data-files`)
* approximate floating point comparison for doctests that produce floating
  point results (see :ref:`handling-float-output`)
* skipping particular classes, methods, and functions when running doctests
  (see :ref:`skipping-doctests`)
* optional inclusion of ``*.rst`` files for doctests

This plugin provides two command line options: ``--doctest-plus`` for enabling
the advanced features mentioned above, and ``--doctest-rst`` for including
``*.rst`` files in doctest collection.

The Astropy test runner enables both of these options by default. When running
the test suite directly from ``pytest`` (instead of through the Astropy test
runner), it is necessary to explicitly provide these options when they are
needed.
