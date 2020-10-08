.. doctest-skip-all

.. include:: workflow/known_projects.inc
.. _testing-guidelines:

******************
Testing Guidelines
******************

This section describes the testing framework and format standards for tests in
Astropy core packages (this also serves as recommendations for affiliated
packages).

Testing Framework
*****************

The testing framework used by astropy (and packages using the :doc:`Astropy
package template <astropy-package-template>`) is the `pytest`_ framework,
accessed through the ``python setup.py test`` command.

.. _testing-dependencies:

Testing Dependencies
********************

As of Astropy 3.0, the dependencies used by the Astropy test runner are
provided by a separate package called `pytest-astropy`_. This package provides
the ``pytest`` dependency itself, in addition to several ``pytest`` plugins
that are used by Astropy, and will also be of general use to other packages.

Since the testing dependencies are not actually required to install or use
Astropy, they are not included in ``install_requires`` in ``setup.py``.
However, for technical reasons it is not currently possible to express these
dependencies in ``tests_require`` either. Therefore, ``pytest-astropy`` is
listed as an extra dependency using ``extras_require`` in ``setup.py``.
Developers who want to run the test suite will need to install the testing
package using pip::

    > pip install pytest-astropy

A detailed description of the plugins can be found in the :ref:`pytest-plugins`
section.

.. _pytest-astropy: https://github.com/astropy/pytest-astropy

.. _running-tests:

Running Tests
*************

There are currently three different ways to invoke Astropy tests. Each
method invokes `pytest`_ to run the tests but offers different options when
calling. To run the tests, you will need to make sure you have the `pytest`_
package (version 3.1 or later) installed.

In addition to running the Astropy tests, these methods can also be called
so that they check Python source code for `PEP8 compliance
<https://www.python.org/dev/peps/pep-0008/>`_. All of the PEP8 testing
options require the `pytest-pep8 plugin
<https://pypi.org/project/pytest-pep8>`_, which must be installed
separately.

setup.py test
=============

The astropy core package and the Astropy package template provide a ``test``
setup command, invoked by running ``python setup.py test`` while in the
package root directory. Run ``python setup.py test --help`` to see the
options to the test command.

Since ``python setup.py test`` wraps the widely-used pytest framework, you may
from time to time want to pass options to the ``pytest`` command itself. For
example, the ``-x`` option to stop after the first failure can be passed
through with the ``--args`` argument::

    > python setup.py test --args "-x"

`pytest`_ will look for files that
:ref:`look like tests <pytest:python test discovery>`, per default
in the current directory and all recursive directories, then run all
the code that looks like tests (essentially all functions or methods
with ``test_`` prefixes or inside ``Test`` prefixed classes) within
those files.

Turn on PEP8 checking by passing ``--pep8`` to the ``test`` command. This will
turn off regular testing and enable PEP8 testing.

Note also that this test runner actually installs astropy into a temporary
directory and uses that for running the tests.  This means that tests of things
like entry points or data file paths should act just like they would once
astropy is installed.  The other two approaches described below do *not* do
this, and hence may give different results when run from the astropy source
code. Hence if you're running the tests because you've modified code that might
be impacted by this, the ``setup.py test`` approach is the recommended method.

.. _astropy.test():

astropy.test()
==============

Tests can be run from within Astropy with::

    import astropy
    astropy.test()

This will run all the default tests for Astropy.

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

pytest
======

The test suite can be run directly from the native ``pytest`` command. In this
case, it is important for developers to be aware that they must manually
rebuild any extensions by running ``setup.py build_ext`` before testing.

In contrast to the case of running from ``setup.py``, the ``--doctest-plus``
and ``--doctest-rst`` options are not enabled by default when running the
``pytest`` command directly. This flags should be explicitly given if they are
needed.

Test-running options
====================

Running parts of the test suite
-------------------------------

It is possible to run only the tests for a particular subpackage or set of
subpackages.  For example, to run only the ``wcs`` tests from the
commandline::

    python setup.py test -P wcs

Or, to run only the ``wcs`` and ``utils`` tests::

    python setup.py test -P wcs,utils

Or from Python::

    >>> import astropy
    >>> astropy.test(package="wcs,utils")

You can also specify a single file to test from the commandline::

    python setup.py test -t astropy/wcs/tests/test_wcs.py

When the ``-t`` option is given a relative path, it is relative to the
installed root of astropy.  When ``-t`` is given a relative path to a
documentation ``.rst`` file to test, it is relative to the root of the
documentation, i.e. the ``docs`` directory in the source tree.  For
example::

    python setup.py test -t units/index.rst

.. _open-files:

Testing for open files
----------------------

Astropy can test whether any of the unit tests inadvertently leave any
files open.  Since this greatly slows down the time it takes to run
the tests, it is turned off by default.

To use it from the commandline, do::

    python setup.py test --open-files

To use it from Python, do::

    >>> import astropy
    >>> astropy.test(open_files=True)

For more information on the ``pytest-openfiles`` plugin see
:ref:`openfiles-plugin`

Test coverage reports
---------------------

Astropy can use `coverage.py <https://coverage.readthedocs.io/en/latest/>`_ to
generate test coverage reports.  To generate a test coverage report, use::

    python setup.py test --coverage

There is a `coveragerc
<https://coverage.readthedocs.io/en/latest/config.html>`_ file that
defines files to omit as well as lines to exclude.  It is installed
along with astropy so that the ``astropy`` testing framework can use
it.  In the source tree, it is at ``astropy/tests/coveragerc``.

Running tests in parallel
-------------------------

It is possible to speed up astropy's tests using the `pytest-xdist
<https://pypi.org/project/pytest-xdist>`_ plugin.  This plugin can be
installed using `pip`_::

    pip install pytest-xdist

Once installed, tests can be run in parallel using the ``'--parallel'``
commandline option.  For example, to use 4 processes::

    python setup.py test --parallel=4

Pass ``--parallel=auto`` to create the same number of processes as cores
on your machine.

Similarly, this feature can be invoked from Python::

    >>> import astropy
    >>> astropy.test(parallel=4)

Running tests to catch permissions errors
-----------------------------------------

It is possible to write code or tests that write into the source directory. This
is not desirable because Python packages can be (and frequently are) installed
in locations where the user may not have write permissions.  To check for these
cases, the test runner has an option to have the test-runner
directory be set as read-only to ensure the tests are not writing to that
location.  This mode can be triggered by running the tests like so::

    python setup.py test --readonly



.. _writing-tests:

Writing tests
*************

``pytest`` has the following test discovery rules:

 * ``test_*.py`` or ``*_test.py`` files
 * ``Test`` prefixed classes (without an ``__init__`` method)
 * ``test_`` prefixed functions and methods

Consult the :ref:`test discovery rules <pytest:python test discovery>`
for detailed information on how to name files and tests so that they are
automatically discovered by `pytest`_.

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
    python: platform darwin -- Python 3.6.0 -- pytest-3.2.0
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
the remote data tests at the command line with ``python setup.py test
--remote-data=any``.

It is possible to mark tests using
``@pytest.mark.remote_data(source='astropy')``, which can be used to indicate
that the only required data is from the http://data.astropy.org server. To
enable just these tests, you can run the
tests with ``python setup.py test --remote-data=astropy``.

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

Furthermore, it is possible to set an option ``astropy_cache_dir`` in the
pytest config file which sets the cache location for the entire test run.  A
``--astropy-cache-dir`` command-line option is also supported (which overrides
all other settings).  Currently it is not directly supported by the
``./setup.py test`` command, so it is necessary to use it with the ``-a``
argument like::

    $ ./setup.py test -a "--astropy-cache-dir=/path/to/custom/cache/dir"


Tests that create files
=======================

Tests may often be run from directories where users do not have write
permissions so tests which create files should always do so in
temporary directories. This can be done with the
:ref:`pytest 'tmpdir' fixture <pytest:tmpdir handling>` or with
Python's built-in :ref:`tempfile module <python:tempfile-examples>`.


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
functions. In the following ::

    def add_nums(x, y):
        """Add two numbers."""
        return x + y

    class TestAdd42(object):
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

    class TestAdd42(object):
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
bugs in time handling.  See issue #9017 and pull request #9532 for details!

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

Using pytest helper functions
=============================

If your tests need to use `pytest helper functions
<https://docs.pytest.org/en/latest/reference.html#functions>`_, such as
``pytest.raises``, import ``pytest`` into your test module like so::

    import pytest

Prior to Astropy 2.0, it was possible to import pytest from a
bundled version using e.g.::

    from ...tests.helper import pytest

but this is no longer the recommended method.

Testing warnings
================

In order to test that warnings are triggered as expected in certain
situations, you can use the `astropy.tests.helper.catch_warnings`
context manager.  Unlike the `warnings.catch_warnings` context manager
in the standard library, this one will reset all warning state before
hand so one is assured to get the warnings reported, regardless of
what errors or warnings may have been emitted by other tests previously.
Here is a real-world example::

  from astropy.tests.helper import catch_warnings

  with catch_warnings(MergeConflictWarning) as warning_lines:
      # Test code which triggers a MergeConflictWarning
      out = table.vstack([t1, t2, t4], join_type='outer')

      assert warning_lines[0].category == metadata.MergeConflictWarning
      assert ("In merged column 'a' the 'units' attribute does not match (cm != m)"
              in str(warning_lines[0].message))

`pytest`_ provides its own context manager
:ref:`pytest.warns <pytest:warns>` that, completely
analogously to ``pytest.raises`` (see below) allows to probe explicitly
for specific warning classes and, through the optional ``match`` argument,
messages. Note that when no warning of the specified type is
triggered, this will make the test fail. When checking for optional,
but not mandatory warnings, Astropy's ``catch_warnings`` is therefore the
better option, as it will collect any number of warning lines (including
zero).

.. note::

   With `pytest`_ there is also the option of using the
   :ref:`recwarn <pytest:recwarn>` function argument to test that
   warnings are triggered within the entire embedding function.
   This method has been found to be problematic in at least one case
   (`pull request 1174 <https://github.com/astropy/astropy/pull/1174#issuecomment-20249309>`_)
   and the alternative of calling ``pytest.warns`` with a warning type
   argument of ``None`` requires further processing of the recorded
   warning(s), so the `astropy.tests.helper.catch_warnings` context
   manager is preferred in such cases.

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
:ref:`astropy.visualization.wcsaxes <wcsaxes>`).

To run the Astropy tests with the image comparison, use::

    python setup.py test -a "--mpl" --remote-data

However, note that the output can be very sensitive to the version of Matplotlib
as well as all its dependencies (e.g. freetype), so we recommend running the
image tests inside a `Docker <https://www.docker.com/>`__ container which has a
frozen set of package versions (Docker containers can be thought of as mini
virtual machines). We have made a `set of Docker container images
<https://hub.docker.com/u/astropy/>`__ that can be used for this. Once you have
installed Docker, to run the Astropy tests with the image comparison inside a
Docker container, make sure you are inside the Astropy repository (or the
repository of the package you are testing) then do::

    docker run -it -v ${PWD}:/repo astropy/image-tests-py35-mpl300:1.3 /bin/bash

This will start up a bash prompt in the Docker container, and you should see
something like::

    root@8173d2494b0b:/#

You can now go to the ``/repo`` directory, which is the same folder as
your local version of the repository you are testing::

    cd /repo

You can then run the tests as above::

    python3 setup.py test -a "--mpl" --remote-data

Type ``exit`` to exit the container.

You can find the names of the available Docker images on the `Docker Hub
<https://hub.docker.com/u/astropy>`_.

Writing image tests
-------------------

The `README.rst <https://github.com/matplotlib/pytest-mpl/blob/master/README.rst>`__
for the plugin contains information on writing tests with this plugin. The only
key addition compared to those instructions is that you should set
``baseline_dir``::

    from astropy.tests.image_tests import IMAGE_REFERENCE_DIR

    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR)

This is because since the reference image files would contribute significantly
to the repository size, we instead store them on the http://data.astropy.org
site. The downside is that it is a little more complicated to create or
re-generate reference files, but we describe the process here.

Generating reference images
---------------------------

Once you have a test for which you want to (re-)generate reference images,
start up one of the Docker containers using e.g.::

  docker run -it -v ${PWD}:/repo astropy/image-tests-py35-mpl300:1.3 /bin/bash

then run the tests inside ``/repo`` with the ``--mpl-generate-path`` argument, e.g::

    cd repo
    python3 setup.py test -a "--mpl --mpl-generate-path=reference_tmp" --remote-data

This will create a ``reference_tmp`` folder and put the generated reference
images inside it - the folder will be available in the repository outside of
the Docker container. Type ``exit`` to exit the container.

Make sure you generate images for the different supported Matplotlib versions
using the available containers.

Uploading the reference images
------------------------------

Next, we need to add these images to the http://data.astropy.org server. To do
this, open a pull request to `this <https://github.com/astropy/astropy-data>`_
repository. The reference images for Astropy tests should go inside the
`testing/astropy <https://github.com/astropy/astropy-data/tree/gh-pages/testing/astropy>`_
directory. In that directory are folders named as timestamps. If you are simply
adding new tests, add the reference files to the most recent directory.

If you are re-generating baseline images due to changes in Astropy, make a new
timestamp directory by copying one the most recent one, then replace any
baseline images that have changed. Note that due to changes between Matplotlib
versions, we need to add the whole set of reference images for each major
Matplotlib version. Therefore, in each timestamp folder, there are folders named
e.g. ``1.4.x`` and ``1.5.x``.

Once the reference images are merged in and available on
http://data.astropy.org, update the timestamp in the ``IMAGE_REFERENCE_DIR``
variable in the ``astropy.tests.image_tests`` sub-module. Because the timestamp
is hard-coded, adding a new timestamp directory will not mess with testing for
released versions of Astropy, so you can easily add and tweak a new timestamp
directory while still working on a pull request to Astropy.

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

   Since the narrative Sphinx documentation is not installed alongside
   the astropy source code, it can only be tested by running ``python
   setup.py test``, not by ``import astropy; astropy.test()``.

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

* `Travis <https://travis-ci.org/astropy/astropy>`_ for
  64-bit Linux, OS X, and Windows setups
* `CircleCI <https://circleci.com>`_ for 32-bit Linux,
  documentation build, and visualization tests

These continuously test the package for each commit and pull request that is
pushed to GitHub to notice when something breaks.

Astropy and many affiliated packages use an external package called
`ci-helpers <https://github.com/astropy/astropy-helpers>`_ to provide
support for the generic parts of the CI systems. ``ci-helpers`` consists of
a set of scripts that are used by the ``.travis.yml`` and ``appveyor.yml``
files to set up the conda environment, and install dependencies.

Dependencies can be customized for different packages using the appropriate
environment variables in ``.travis.yml`` and ``appveyor.yml``. For more
details on how to set up this machinery, see the `package-template
<https://github.com/astropy/package-template>`_ and `ci-helpers`_.

The 32-bit tests on CircleCI use a pre-defined Docker image defined `here
<https://github.com/astropy/astropy-docker/>`__ which includes a 32-bit
Python environment. If you want to run tests for packages in the same way,
you can use the same set-up on CircleCI as the core package, but just be
sure to install Astropy first using::

    easy_install pip
    pip install astropy

For convenience, you can also use the ``astropy/affiliated-32bit-test-env``
Docker image instead of ``astropy/astropy-32bit-test-env`` - the former includes
the latest stable version of Astropy pre-installed.

In some cases, you may see failures on continuous integration services that
you do not see locally, for example because the operating system is different,
or because the failure happens with only 32-bit Python. The following sections
explain how you can reproduce specific builds locally.

Reproducing failing 32-bit builds
=================================

If you want to run your tests in the same 32-bit Python environment that
CircleCI uses, start off by installing `Docker <https://www.docker.com>`__ if you
don't already have it installed. Docker can be installed on a variety of
different operating systems.

Then, make sure you have a version of the git repository (either the main
Astropy repository or your fork) for which you want to run the tests. Go to that
directory, then run Docker with::

    $ docker run -i -v ${PWD}:/astropy_src -t astropy/astropy-32bit-test-env:1.6 bash

This will put you in the bash shell inside the Docker container. Once inside,
you can go to the ``astropy_src`` directory, and you should see the files that
are in your local git repository::

    root@5e2b89d7b07c:/# cd /astropy_src
    root@5e2b89d7b07c:/astropy_src# ls
    ah_bootstrap.py  CONTRIBUTING.md       pip-requirements-doc
    appveyor.yml     docs                  README.rst
    astropy          examples              readthedocs.yml
    astropy_helpers  ez_setup.py           setup.cfg
    cextern          licenses              setup.py
    CHANGES.rst      MANIFEST.in           static
    circle.yml       pip-requirements      CITATION
    pip-requirements-dev

You can then run the tests with::

    root@5e2b89d7b07c:/astropy_src# python setup.py test

.. _pytest-plugins:

Pytest Plugins
**************

The following ``pytest`` plugins are maintained and used by Astropy. They are
included in the ``pytest-astropy`` package, which is now required for testing
Astropy. More information on all of the  plugins provided by the
``pytest-astropy`` package (including dependencies not maintained by Astropy)
can be found `here <https://github.com/astropy/pytest-astropy>`__.

.. _remotedata-plugin:

pytest-remotedata
=================

The `pytest-remotedata`_ plugin allows developers to control whether to run
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

.. _pytest-remotedata: https://github.com/astropy/pytest-remotedata

.. _doctestplus-plugin:

pytest-doctestplus
==================

The `pytest-doctestplus`_ plugin provides advanced doctest features, including:

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

.. _pytest-doctestplus: https://github.com/astropy/pytest-doctestplus

.. _openfiles-plugin:

pytest-openfiles
================

The `pytest-openfiles`_ plugin allows for the detection of open I/O resources
at the end of unit tests. This plugin adds the ``--open-files`` option to the
``pytest`` command (which is also exposed through the Astropy test runner).

When running tests with ``--open-files``, if a file is opened during the course
of a unit test but that file  not closed before the test finishes, the test
will fail. This is particularly useful for testing code that manipulates file
handles or other I/O resources. It allows developers to ensure that this kind
of code properly cleans up I/O resources when they are no longer needed.

Also see :ref:`open-files`.

.. _pytest-openfiles: https://github.com/astropy/pytest-openfiles
