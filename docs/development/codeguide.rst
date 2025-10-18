.. doctest-skip-all
.. _code-guide:

*****************
Coding Guidelines
*****************

This section describes requirements and guidelines that should be followed both
by the core package and by coordinated packages, and these are also recommended
for affiliated packages.

Interface and Dependencies
==========================

* All code must be compatible with the versions of Python indicated by the
  ``requires-python`` key  under ``[project]`` in the `pyproject.toml
  <https://github.com/astropy/astropy/blob/main/pyproject.toml>`_ file of the
  core package.

* The core package should be importable with no
  dependencies other than components already in the Astropy core, the
  `Python Standard Library <https://docs.python.org/3/library/index.html>`_,
  and |NumPy| |minimum_numpy_version| or later.

* Additional dependencies - such as |SciPy|, |Matplotlib|, or other
  third-party packages - are allowed for sub-modules or in function
  calls, but they must be noted in the package documentation and
  should only affect the relevant component.  In functions and
  methods, the optional dependency should use a normal ``import``
  statement, which will raise an ``ImportError`` if the dependency is
  not available. In the astropy core package, such optional dependencies should
  be recorded in the ``pyproject.toml`` file in the ``[project.optional-dependencies]``
  entry (put it in the appropriate category of optional dependencies).

  At the module level, one can subclass a class from an optional dependency
  like so::

      try:
          from opdep import Superclass
      except ImportError:
          warn(AstropyWarning('opdep is not present, so <functionality>'
                              'will not work.'))
          class Superclass: pass

      class Customclass(Superclass):
          ...

* General utilities necessary for but not specific to the package or
  sub-package should be placed in a ``packagename.utils`` module (e.g.
  ``astropy.utils`` for the core package). If a utility is already present in
  :mod:`astropy.utils`, packages should always use that utility instead of
  re-implementing it in ``packagename.utils`` module.

Documentation and Testing
=========================

* Docstrings must be present for all public classes/methods/functions, and
  must follow the form outlined in the :doc:`docguide` document.

* Write usage examples in the docstrings of all classes and functions whenever
  possible. These examples should be short and simple to reproduce--users
  should be able to copy them verbatim and run them. These examples should,
  whenever possible, be in the :ref:`doctest <doctests>` format and will be
  executed as part of the test suite.

* Unit tests should be provided for as many public methods and functions as
  possible, and should adhere to the standards set in the :doc:`testguide`
  document.


Data and Configuration
======================

* Packages can include data in a directory named ``data`` inside a subpackage
  source directory as long as it is less than about 100 kB. These data should
  always be accessed via the :func:`~astropy.utils.data.get_pkg_data_fileobj` or
  :func:`~astropy.utils.data.get_pkg_data_filename` functions. If the data
  exceeds this size, it should be hosted outside the source code repository,
  either at a third-party location on the internet or the `astropy data server
  <https://github.com/astropy/astropy-data>`_.
  In either case, it should always be downloaded using the
  :func:`~astropy.utils.data.get_pkg_data_fileobj` or
  :func:`~astropy.utils.data.get_pkg_data_filename` functions. If a specific
  version of a data file is needed, the hash mechanism described in
  :mod:`astropy.utils.data` should be used.

* All persistent configuration should use the
  :ref:`astropy_config` mechanism.  Such configuration items
  should be placed at the top of the module or package that makes use of them,
  and supply a description sufficient for users to understand what the setting
  changes.

Standard output, warnings, and errors
=====================================

The built-in ``print(...)`` function should only be used for output that
is explicitly requested by the user, for example ``print_header(...)``
or ``list_catalogs(...)``. Any other standard output, warnings, and
errors should follow these rules:

* For errors/exceptions, one should always use ``raise`` with one of the
  built-in exception classes, or a custom exception class. The
  nondescript ``Exception`` class should be avoided as much as possible,
  in favor of more specific exceptions (`IOError`, `ValueError`,
  etc.).

* For warnings, one should always use ``warnings.warn(message,
  warning_class)``. These get redirected to ``log.warning()`` by default,
  but one can still use the standard warning-catching mechanism and custom
  warning classes. The warning class should be either
  :class:`~astropy.utils.exceptions.AstropyUserWarning` or inherit from it.

* For informational and debugging messages, one should always use
  ``log.info(message)`` and ``log.debug(message)``.

The logging system uses the built-in Python :py:mod:`logging`
module. The logger can be imported using::

    from astropy import log

.. _code-style:

Coding Style/Conventions
========================

* The code should follow the standard `PEP8 Style Guide for Python Code
  <https://www.python.org/dev/peps/pep-0008/>`_.

  * ``astropy`` itself enforces this style guide using the
    `ruff format <https://docs.astral.sh/ruff/formatter/>`_ code formatter, which closely follows the
    `The Black Code Style <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`_.

  * In the rare cases that ruff_ formatting is undesirable, it is possible to
    `disable formatting locally <https://docs.astral.sh/ruff/formatter/#format-suppression>`_.

      .. note::
        When a list or array should be formatted as one item per line then this is best
        achieved by using the
        `magic trailing comma <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html#the-magic-trailing-comma>`_.
        This is frequently sufficient for keeping matrices formatted as one row
        per line while still allowing ruff_ to check the code::

            arr = [
                [0, 1],
                [1, 0],  # notice the trailing comma.
            ]


* Our testing infrastructure currently enforces a subset of the |PEP8| style guide. In
  addition, these checks also enforce `isort <https://pycqa.github.io/isort/>`_ to sort
  the module imports and a large set of style-checks supported by ruff_.

  * We provide a `pre-commit <https://pre-commit.com/>`_ configuration which
    automatically enforces and fixes (whenever possible) the coding style, see
    :ref:`pre-commit` for details on how to set up and use this. We note that the
    particular set of |PEP8| and style-related checks that are used in Astropy do not
    need to be used in affiliated packages. In particular, the set of ruff_ checks is
    not required for affiliated packages.

  .. note:: There are multiple options for testing PEP8 compliance of code,
            see :doc:`testguide` for more information.

* ``astropy`` source code should contain a comment at the beginning of the file
  pointing to the license for the ``astropy`` source code.  This line should say::

      # Licensed under a 3-clause BSD style license - see LICENSE.rst

* Classes should either use direct variable access, or Python’s property
  mechanism for setting object instance variables. ``get_value``/``set_value``
  style methods should be used only when getting and setting the values
  requires a computationally-expensive operation. The
  :ref:`prop-get-set-example` example below illustrates this guideline.

* Classes should use the builtin `super` function when making calls to
  methods in their super-class(es) unless there are specific reasons not to.
  `super` should be used consistently in all subclasses since it does not
  work otherwise. The :ref:`super-vs-direct-example` example below illustrates
  why this is important.

* Multiple inheritance should be avoided in general without good reason.
  Multiple inheritance is complicated to implement well, which is why many
  object-oriented languages, like Java, do not allow it at all.  Python does
  enable multiple inheritance through use of the
  `C3 Linearization <https://www.python.org/download/releases/2.3/mro/>`_
  algorithm, which provides a consistent method resolution ordering.
  Non-trivial multiple-inheritance schemes should not be attempted without
  good justification, or without understanding how C3 is used to determine
  method resolution order.  However, trivial multiple inheritance using
  orthogonal base classes, known as the 'mixin' pattern, may be used.

* ``__init__.py`` files for modules should not contain any significant
  implementation code. ``__init__.py`` can contain docstrings and code for
  organizing the module layout, however (e.g. ``from submodule import *``
  in accord with the guideline above). If a module is small enough that
  it fits in one file, it should simply be a single file, rather than a
  directory with an ``__init__.py`` file.

* Command-line scripts should follow the form outlined in the :doc:`scripts`
  document.

.. _handling-unicode:

Unicode guidelines
==================

For maximum compatibility, we need to assume that writing non-ASCII
characters to the console or to files will not work.  However, for
those that have a correctly configured Unicode environment, we should
allow them to opt-in to take advantage of Unicode output when
appropriate.  Therefore, there is a global configuration option,
``astropy.conf.unicode_output`` to enable Unicode output of values, set
to `False` by default.

The following conventions should be used for classes that implement the
standard string conversion methods:

- `~object.__repr__`: Return a `str` containing only ASCII characters.
  The output must be independent of the ``astropy.conf.unicode_output``
  setting.

- `~object.__str__`: Return a `str` containing only ASCII characters if
  ``astropy.conf.unicode_output`` is `False`.
  If ``astropy.conf.unicode_output`` is `True`, it may contain non-ASCII
  characters.

- `~.object.__format__`: Return a `str` containing only ASCII characters if
  ``astropy.conf.unicode_output`` is `False` and the ``format_spec`` argument
  is an empty string.
  Otherwise it may contain non-ASCII characters.

For classes that are expected to roundtrip through strings, the parser must
accept the output of `~object.__str__`.

This design generally follows `Postel's Law
<https://en.wikipedia.org/wiki/Robustness_principle>`_: "Be liberal in what you
accept, and conservative in what you send."

There is a test helper,
:func:`~astropy.tests.helper.assert_follows_unicode_guidelines`,
to check compliance with the above guidelines.

Including C Code
================

* C extensions are only allowed when they provide a significant performance
  enhancement over pure Python, or a robust C library already exists to
  provided the needed functionality. When C extensions are used, the Python
  interface must meet the aforementioned Python interface guidelines.

* The use of Cython_ is strongly recommended for C extensions. Cython_
  extensions should store ``.pyx`` files in the source code repository,
  but not the generated ``.c`` files.

* If a C extension has a dependency on an external C library, the source code
  for the library should be bundled with the Astropy core, provided the
  license for the C library is compatible with the Astropy license.
  Additionally, the package must be compatible with using a system-installed
  library in place of the library included in Astropy, and a user installing
  the package should be able to opt-in to using the system version using
  a ``ASTROPY_USE_SYSTEM_???`` environment variable, where ``???`` is the name
  of the library, e.g. ``ASTROPY_USE_SYSTEM_WCSLIB`` (see also
  :ref:`external_c_libraries`).

* In cases where C extensions are needed but Cython_ cannot be used, the `PEP 7
  Style Guide for C Code <https://www.python.org/dev/peps/pep-0007/>`_ is
  recommended.

* C extensions (Cython_ or otherwise) should provide the necessary information
  for building the extension via the mechanisms described in
  :ref:`building-c-or-cython-extensions`.


Requirements Specific to Affiliated Packages
============================================

* Affiliated packages implementing many classes/functions not relevant to
  the affiliated package itself (for example leftover code from a previous
  package) will not be accepted - the package should only include the
  required functionality and relevant extensions.

* Affiliated packages must be registered on the `Python Package Index
  <https://pypi.org/>`_, with proper metadata for downloading and
  installing the source package.

* The ``astropy`` root package name should not be used by affiliated
  packages - it is reserved for use by the core package.

Examples
========

This section shows a few examples (not all of which are correct!) to
illustrate points from the guidelines.

.. _prop-get-set-example:

Properties vs. get\_/set\_
--------------------------

This example shows a sample class illustrating the guideline regarding the use
of `properties <https://docs.python.org/3/library/functions.html#property>`_ as
opposed to getter/setter methods.

Let's assume you've defined a ``Star`` class and create an instance like this::

    >>> s = Star(B=5.48, V=4.83)

You should always use attribute syntax like this::

    >>> s.color = 0.4
    >>> print(s.color)
    0.4

Rather than like this::

    >>> s.set_color(0.4)  # Bad form!
    >>> print(s.get_color())  # Bad form!
    0.4

Using Python properties, attribute syntax can still do anything possible with
a get/set method. For lengthy or complex calculations, however, use a method::

    >>> print(s.compute_color(5800, age=5e9))
    0.4

.. _super-vs-direct-example:

super() vs. Direct Calling
--------------------------

This example shows why the use of `super` leads to a more consistent
method resolution order than manually calling methods of the super classes in a
multiple inheritance case::

    # This is dangerous and bug-prone!

    class A:
        def method(self):
            print('Doing A')


    class B(A):
        def method(self):
            print('Doing B')
            A.method(self)


    class C(A):
        def method(self):
            print('Doing C')
            A.method(self)

    class D(C, B):
        def method(self):
            print('Doing D')
            C.method(self)
            B.method(self)

if you then do::

    >>> b = B()
    >>> b.method()

you will see::

    Doing B
    Doing A

which is what you expect, and similarly for C. However, if you do::

    >>> d = D()
    >>> d.method()

you might expect to see the methods called in the order D, B, C, A but instead
you see::

    Doing D
    Doing C
    Doing A
    Doing B
    Doing A

because both ``B.method()`` and ``C.method()`` call ``A.method()`` unaware of
the fact that they're being called as part of a chain in a hierarchy.  When
``C.method()`` is called it is unaware that it's being called from a subclass
that inherits from both ``B`` and ``C``, and that ``B.method()`` should be
called next.  By calling `super` the entire method resolution order for
``D`` is precomputed, enabling each superclass to cooperatively determine which
class should be handed control in the next `super` call::

    # This is safer

    class A:
        def method(self):
            print('Doing A')

    class B(A):
        def method(self):
            print('Doing B')
            super().method()


    class C(A):
        def method(self):
            print('Doing C')
            super().method()

    class D(C, B):
        def method(self):
            print('Doing D')
            super().method()

::

    >>> d = D()
    >>> d.method()
    Doing D
    Doing C
    Doing B
    Doing A

As you can see, each superclass's method is entered only once.  For this to
work it is very important that each method in a class that calls its
superclass's version of that method use `super` instead of calling the
method directly.  In the most common case of single-inheritance, using
``super()`` is functionally equivalent to calling the superclass's method
directly.  But as soon as a class is used in a multiple-inheritance
hierarchy it must use ``super()`` in order to cooperate with other classes in
the hierarchy.

.. note:: For more information on the benefits of `super`, see
          https://rhettinger.wordpress.com/2011/05/26/super-considered-super/

.. _Numpy: https://numpy.org/
.. _Scipy: https://www.scipy.org/
.. _matplotlib: https://matplotlib.org/
.. _Cython: https://cython.org/
.. _PyPI: https://pypi.org/project
.. _ruff: https://docs.astral.sh/ruff/
