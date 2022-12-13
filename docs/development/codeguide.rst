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
  ``python_requires`` key in the `setup.cfg
  <https://github.com/astropy/astropy/blob/main/setup.cfg>`_ file of the
  core package.

* Usage of ``six``, ``__future__``, and ``2to3`` is no longer acceptable.

* `f-strings <https://docs.python.org/3/reference/lexical_analysis.html#f-strings>`_
  should be used when possible, and if not, Python 3
  formatting should be used (i.e. ``"{0:s}".format("spam")``)
  instead of the ``%`` operator (``"%s" % "spam"``).

* The core package should be importable with no
  dependencies other than components already in the Astropy core, the
  `Python Standard Library <https://docs.python.org/3/library/index.html>`_,
  and NumPy_ |minimum_numpy_version| or later.

* Additional dependencies - such as SciPy_, Matplotlib_, or other
  third-party packages - are allowed for sub-modules or in function
  calls, but they must be noted in the package documentation and
  should only affect the relevant component.  In functions and
  methods, the optional dependency should use a normal ``import``
  statement, which will raise an ``ImportError`` if the dependency is
  not available. In the astropy core package, such optional dependencies should
  be recorded in the ``setup.cfg`` file in the ``extras_require``
  entry, under ``all`` (or ``test_all`` if the dependency is only
  needed for testing).

  At the module level, one can subclass a class from an optional dependency
  like so::

      try:
          from opdep import Superclass
      except ImportError:
          warn(AstropyWarning('opdep is not present, so <functionality>'
                              'will not work.'))
          class Superclass(object): pass

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
  <https://www.python.org/dev/peps/pep-0008/>`_. In particular, this includes
  using only 4 spaces for indentation, and never tabs.

  * ``astropy`` itself enforces this style guide using the
    `black <https://black.readthedocs.io/en/stable/>`_ code formatter.

  * We recognize that sometimes ``black`` will autoformat things in undesirable
    ways, e.g., matrices.  In the cases that ``black`` produces undesirable code
    formatting:

      * one can wrap code the code in ``# fmt: off`` and ``# fmt: on`` to disable
        ``black`` formatting over multiple lines.

      * or one can add a single ``# fmt: skip`` comment to the end of a line to
        disable ``black`` formatting for that line.

    This should be done sparingly, and only
    when ``black`` produces undesirable formatting.

* Our testing infrastructure currently enforces a subset of the PEP8 style
  guide. In addition, these checks also enforce
  `isort <https://pycqa.github.io/isort/>`_ to sort the module imports.

  * We provide a ``pre-commit`` hook which automatically enforces and fixes
    (whenever possible) the coding style, see :ref:`pre-commit` for details on
    how to setup and use this.

  * Alternately, you can manually check and fix your changes by running the
    following `tox <https://tox.readthedocs.io/>`__ command::

      tox -e codestyle

* Following PEP8's recommendation, absolute imports are to be used in general.
  The exception to this is relative imports of the form
  ``from . import modname``, best when referring to files within the same
  sub-module.  This makes it clearer what code is from the current submodule
  as opposed to from another.

  .. note:: There are multiple options for testing PEP8 compliance of code,
            see :doc:`testguide` for more information.
            See :doc:`codeguide_emacs` for some configuration options for Emacs
            that helps in ensuring conformance to PEP8.

* Astropy source code should contain a comment at the beginning of the file (or
  immediately after the ``#!/usr/bin env python`` command, if relevant)
  pointing to the license for the Astropy source code.  This line should say::

      # Licensed under a 3-clause BSD style license - see LICENSE.rst

* The following naming conventions::

    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt

  should be used wherever relevant. On the other hand::

    from packagename import *

  should never be used, except as a tool to flatten the namespace of a module.
  An example of the allowed usage is given in the :ref:`import-star-example`
  example.

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

The following conventions should be used for classes that define the
standard string conversion methods (``__str__``, ``__repr__``,
``__bytes__``, and ``__format__``).  In the bullets
below, the phrase "string instance" is used to refer to `str`, while
"bytes instance" is used to refer to `bytes`.

- ``__repr__``: Return a "string instance" containing only 7-bit characters.

- ``__bytes__``: Return a "bytes instance" containing only 7-bit characters.

- ``__str__``: Return a "string instance".
  If ``astropy.conf.unicode_output`` is `False`, it must contain
  only 7-bit characters.  If ``astropy.conf.unicode_output`` is `True`, it
  may contain non-ASCII characters when applicable.

- ``__format__``: Return a "string instance".  If
  ``astropy.conf.unicode_output`` is `False`, it must contain only 7-bit
  characters.  If ``astropy.conf.unicode_output`` is `True`, it may contain
  non-ASCII characters when applicable.

For classes that are expected to roundtrip through strings (unicode or
bytes), the parser must accept the output of ``__str__``.
Additionally, ``__repr__`` should roundtrip when that makes sense.

This design generally follows Postel's Law: "Be liberal in what you
accept, and conservative in what you send."

The following example class shows a way to implement this::

    # -*- coding: utf-8 -*-

    from astropy import conf

    class FloatList(object):
        def __init__(self, init):
            if isinstance(init, str):
                init = init.split('‖')
            elif isinstance(init, bytes):
                init = init.split(b'|')
            self.x = [float(x) for x in init]

        def __repr__(self):
            # Return unicode object containing no non-ASCII characters
            return f'<FloatList [{", ".join(str(x) for x in self.x)}]>'

        def __bytes__(self):
            return b'|'.join(bytes(x) for x in self.x)

        def __str__(self):
            if astropy.conf.unicode_output:
                return '‖'.join(str(x) for x in self.x)
            else:
                return self.__bytes__().decode('ascii')

Additionally, there is a test helper,
``astropy.test.helper.assert_follows_unicode_guidelines`` to ensure that a
class follows the Unicode guidelines outlined above.  The following
example test will test that our example class above is compliant::

    def test_unicode_guidelines():
        from astropy.test.helper import assert_follows_unicode_guidelines
        assert_follows_unicode_guidelines(FloatList(b'5|4|3|2'), roundtrip=True)

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
of properties as opposed to getter/setter methods.

Let's assuming you've defined a ``Star`` class and create an instance
like this::

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

    class A(object):
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

    class A(object):
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

.. note:: For more information on the the benefits of `super`, see
          https://rhettinger.wordpress.com/2011/05/26/super-considered-super/

.. _import-star-example:

Acceptable use of ``from module import *``
------------------------------------------

``from module import *`` is discouraged in a module that contains
implementation code, as it impedes clarity and often imports unused variables.
It can, however, be used for a package that is laid out in the following
manner::

    packagename
    packagename/__init__.py
    packagename/submodule1.py
    packagename/submodule2.py

In this case, ``packagename/__init__.py`` may be::

    """
    A docstring describing the package goes here
    """
    from submodule1 import *
    from submodule2 import *

This allows functions or classes in the submodules to be used directly as
``packagename.foo`` rather than ``packagename.submodule1.foo``. If this is
used, it is strongly recommended that the submodules make use of the ``__all__``
variable to specify which modules should be imported. Thus, ``submodule2.py``
might read::

    from numpy import array, linspace

    __all__ = ['foo', 'AClass']

    def foo(bar):
        # the function would be defined here
        pass

    class AClass(object):
        # the class is defined here
        pass

This ensures that ``from submodule import *`` only imports ``foo`` and
``AClass``, but not `numpy.array` or `numpy.linspace`.


Additional Resources
====================

Further tips and hints relating to the coding guidelines are included below.

.. toctree::
    :maxdepth: 1

    codeguide_emacs

.. _Numpy: https://numpy.org/
.. _Scipy: https://www.scipy.org/
.. _matplotlib: https://matplotlib.org/
.. _Cython: https://cython.org/
.. _PyPI: https://pypi.org/project
