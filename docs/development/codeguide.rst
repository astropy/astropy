.. doctest-skip-all
.. _code-guide:

=================
Coding Guidelines
=================

This section describes requirements and guidelines that should be followed
both for the core package and for affiliated packages.

.. note:: Affiliated packages will only be considered for integration as a
          module in the core package once these guidelines have been
          followed.

Interface and Dependencies
--------------------------

* All code must be compatible with Python 2.6, 2.7, as well as 3.3 and
  later.  The use of `six`_ for writing code that is portable between Python
  2.x and 3.x is encouraged going forward.  However, much of our affiliate code
  may still use `2to3`_ to process Python 2.x files to be compatible with
  Python 3.x.

.. _six: http://pythonhosted.org/six/
.. _2to3: https://docs.python.org/2/library/2to3.html

  Packages that use ``six`` must include the following in their
  ``setup_package.py`` file::

      def requires_2to3():
          return False

  Code that uses ``six`` should use the following preamble::

        from __future__ import (absolute_import, division, print_function,
                                unicode_literals)

  Code that uses ``2to3`` may limit the ``__future__`` statement to
  only the following, if necessary (though using all 4 is still
  preferred)::

        from __future__ import print_function, division

  Additional information on writing code using ``six`` that is
  compatible with both Python 2.x and 3.x is in the section
  :ref:`dev-portable`.

* The new Python 3 formatting style should be used (i.e.
  ``"{0:s}".format("spam")`` instead of ``"%s" % "spam"``), although
  when using positional arguments, the position should always be
  specified (i.e.  ``"{:s}"`` is not compatible with Python
  2.6).

* The core package and affiliated packages should be importable with
  no dependencies other than components already in the Astropy core,
  the `Python Standard Library
  <http://docs.python.org/release/2.6/library/index.html>`_, and
  NumPy_ |minimum_numpy_version| or later.

* The package should be importable from the source tree at build time. This
  means that, for example, if the package relies on C extensions that have
  yet to be built, the Python code is still importable, even if none of its
  functionality will work. One way to ensure this is to import the functions
  in the C extensions only within the functions/methods that require them
  (see next bullet point).

* Additional dependencies - such as SciPy_, Matplotlib_, or other
  third-party packages - are allowed for sub-modules or in function
  calls, but they must be noted in the package documentation and
  should only affect the relevant component.  In functions and
  methods, the optional dependency should use a normal ``import``
  statement, which will raise an ``ImportError`` if the dependency is
  not available.

  At the module level, one can subclass a class from an optional dependency
  like so::

      try:
          from opdep import Superclass
      except ImportError:
          warn(AstropyWarning('opdep is not present, so <functionality below> will not work.'))
          class SuperClass(object): pass

      class Whatever(Superclass):
          ...

* General utilities necessary for but not specific to the package or
  sub-package should be placed in the ``packagename.utils`` module. These
  utilities will be moved to the :mod:`astropy.utils` module when the
  package is integrated into the core package. If a utility is already
  present in :mod:`astropy.utils`, the package should always use that
  utility instead of re-implementing it in ``packagename.utils`` module.


Documentation and Testing
-------------------------

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
----------------------

* Packages can include data in a directory named ``data`` inside a subpackage
  source directory as long as it is less than about 100 kb. These data should
  always be accessed via the :func:`astropy.utils.data.get_pkg_data_fileobj` or
  :func:`astropy.utils.data.get_pkg_data_filename` functions. If the data
  exceeds this size, it should be hosted outside the source code repository,
  either at a third-party location on the internet or the astropy data server.
  In either case, it should always be downloaded using the
  :func:`astropy.utils.data.get_pkg_data_fileobj` or
  :func:`astropy.utils.data.get_pkg_data_filename` functions. If a specific
  version of a data file is needed, the hash mechanism described in
  :mod:`astropy.utils.data` should be used.

* All persistent configuration should use the
  `astropy.config.ConfigurationItem` mechanism.  Such configuration items
  should be placed at the top of the module or package that makes use of them,
  and supply a description sufficient for users to understand what the setting
  changes.

Standard output, warnings, and errors
-------------------------------------

The built-in ``print(...)`` function should only be used for output that
is explicitly requested by the user, for example ``print_header(...)``
or ``list_catalogs(...)``. Any other standard output, warnings, and
errors should follow these rules:

* For errors/exceptions, one should always use ``raise`` with one of the
  built-in exception classes, or a custom exception class. The
  nondescript ``Exception`` class should be avoided as much as possible,
  in favor of more specific exceptions (``IOError``, ``ValueError``,
  etc.).

* For warnings, one should always use ``warnings.warn(message, warning_class)``. These get redirected to ``log.warn`` by default, but one
  can still use the standard warning-catching mechanism and custom warning
  classes. The warning class should be either
  :class:`~astropy.utils.exceptions.AstropyUserWarning` or inherit from it.

* For informational and debugging messages, one should always use
  ``log.info(message)`` and ``log.debug(message)``.

The logging system uses the built-in Python `logging
<http://docs.python.org/library/logging.html>`_ module. The logger can
be imported using::

    from astropy import log

Coding Style/Conventions
------------------------

* The code will follow the standard `PEP8 Style Guide for Python Code
  <http://www.python.org/dev/peps/pep-0008/>`_. In particular, this includes
  using only 4 spaces for indentation, and never tabs.

* One exception is to be made from the PEP8 style: new style relative imports
  of the form ``from . import modname`` are allowed and required for Astropy,
  as opposed to absolute (as PEP8 suggests) or the simpler ``import modname``
  syntax. This is primarily due to improved relative import support since PEP8
  was developed, and to simplify the process of moving modules.

  .. note:: There are multiple options for testing PEP8 compliance of code,
            see :doc:`testguide` for more information.
            See :doc:`codeguide_emacs` for some configuration options for Emacs
            that helps in ensuring conformance to PEP8.

* Astropy source code should contain a comment at the beginning of the file (or
  immediately after the ``#!/usr/bin env python`` command, if relevant)
  pointing to the license for the Astropy source code.  This line should say::

      # Licensed under a 3-clause BSD style license - see LICENSE.rst

* The ``import numpy as np``, ``import matplotlib as mpl``, and ``import
  matplotlib.pyplot as plt`` naming conventions should be used wherever
  relevant. ``from packagename import *`` should never be used, except as a
  tool to flatten the namespace of a module. An example of the allowed usage
  is given in :ref:`import-star-example`.

* Classes should either use direct variable access, or python’s property
  mechanism for setting object instance variables. ``get_value``/``set_value``
  style methods should be used only when getting and setting the values
  requires a computationally-expensive operation. :ref:`prop-get-set-example`
  below illustrates this guideline.

* All new classes should be new-style classes inheriting from :class:`object`
  (in Python 3 this is a non-issue as all classes are new-style by default).
  The one exception to this rule is older classes in third-party libraries such
  the Python standard library or numpy.

* Classes should use the builtin :func:`super` function when making calls to
  methods in their super-class(es) unless there are specific reasons not to.
  :func:`super` should be used consistently in all subclasses since it does not
  work otherwise.  :ref:`super-vs-direct-example` illustrates why this is
  important.

* Multiple inheritance should be avoided in general without good reason.
  Multiple inheritance is complicated to implement well, which is why many
  object-oriented languages, like Java, do not allow it at all.  Python does
  enable multiple inheritance through use of the
  `C3 Linearization <http://www.python.org/download/releases/2.3/mro/>`_
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

* When ``try...except`` blocks are used to catch exceptions, the ``as``
  syntax should always be used, because this is available in all supported
  versions of python and is less ambiguous syntax (see
  :ref:`try-except-as-example`).

* Command-line scripts should follow the form outlined in the :doc:`scripts`
  document.

.. _handling-unicode:

Unicode guidelines
------------------

For maximum compatibility, we need to assume that writing non-ascii
characters to the console or to files will not work.  However, for
those that have a correctly configured Unicode environment, we should
allow them to opt-in to take advantage of Unicode output when
appropriate.  Therefore, there is a global configuration option,
``astropy.conf.unicode_output`` to enable Unicode output of values, set
to `False` by default.

The following conventions should be used for classes that define the
standard string conversion methods (``__str__``, ``__repr__``,
``__unicode__``, ``__bytes__``, and ``__format__``).  In the bullets
below, the phrase "unicode instance" is used to refer to `unicode` on
Python 2 and `str` on Python 3.  The phrase "bytes instance" is used
to refer to `str` on Python 2 and `bytes` on Python 3.

- ``__repr__``: Return a "unicode instance" (for historical reasons,
  could also be a "bytes instance" on Python 2, though not preferred)
  containing only 7-bit characters.

- ``__str__`` on Python 2 / ``__bytes__`` on Python 3: Return a "bytes
  instance" containing only 7-bit characters.

- ``__unicode__`` on Python 2 / ``__str__`` on Python 3: Return a "unicode
  instance".  If ``astropy.conf.unicode_output`` is `False`, it must contain
  only 7-bit characters.  If ``astropy.conf.unicode_output`` is `True`, it
  may contain non-ascii characters when applicable.

- ``__format__``: Return a "unicode instance".  If
  ``astropy.UNICODE_OUTPUT`` is `False`, it must contain only 7-bit
  characters.  If ``astropy.conf.unicode_output`` is `True`, it may contain
  non-ascii characters when applicable.

For classes that are expected to roundtrip through strings (unicode or
bytes), the parser must accept either the output of ``__str__`` or
``__unicode__`` unambiguously.  Additionally, ``__repr__`` should
roundtrip when that makes sense.

This design generally follows Postel's Law: "Be liberal in what you
accept, and conservative in what you send".

The following example class shows a way to implement this (using `six`_ for
Python 2 and 3 cross-version compatibility::

    # -*- coding: utf-8 -*-

    from __future__ import unicode_literals

    from astropy.extern import six
    from astropy import conf

    class FloatList(object):
        def __init__(self, init):
            if isinstance(init, six.text_type):
                init = init.split('‖')
            elif isinstance(init, bytes):
                init = init.split(b'|')
            self.x = [float(x) for x in init]

        def __repr__(self):
            # Return unicode object containing no non-ascii characters
            return '<FloatList [{0}]>'.format(', '.join(
                six.text_type(x) for x in self.x))

        def __bytes__(self):
            return b'|'.join(bytes(x) for x in self.x)
        if six.PY2:
            __str__ = __bytes__

        def __unicode__(self):
            if astropy.conf.unicode_output:
                return '‖'.join(six.text_type(x) for x in self.x)
            else:
                return self.__bytes__().decode('ascii')
        if six.PY3:
            __str__ = __unicode__

Additionally, there is a test helper,
``astropy.test.helper.assert_follows_unicode_guidelines`` to ensure that a
class follows the Unicode guidelines outlined above.  The following
example test will test that our example class above is compliant::

    def test_unicode_guidelines():
        from astropy.test.helper import assert_follows_unicode_guidelines
        assert_follows_unicode_guidelines(FloatList(b'5|4|3|2'), roundtrip=True)

Including C Code
----------------

* C extensions are only allowed when they provide a significant performance
  enhancement over pure python, or a robust C library already exists to
  provided the needed functionality. When C extensions are used, the Python
  interface must meet the aforementioned python interface guidelines.

* The use of Cython_ is strongly recommended for C extensions, as per the
  example in the template package. Cython_ extensions should store ``.pyx``
  files in the source code repository, but they should be compiled to ``.c``
  files that are updated in the repository when important changes are made to
  the ``.pyx`` file.

* If a C extension has a dependency on an external C library, the source code
  for the library should be bundled with the Astropy core, provided the
  license for the C library is compatible with the Astropy license.
  Additionally, the package must be compatible with using a system-installed
  library in place of the library included in Astropy.

* In cases where C extensions are needed but Cython_ cannot be used, the `PEP 7
  Style Guide for C Code <http://www.python.org/dev/peps/pep-0007/>`_ is
  recommended.

* C extensions (Cython_ or otherwise) should provide the necessary information
  for building the extension via the mechanisms described in
  :ref:`building-c-or-cython-extensions`.

.. _dev-portable:

Writing portable code for Python 2 and 3
----------------------------------------

As of astropy 0.3, the `six`_ library is included to allow supporting Python
2 and 3 from a single code base.  The use of the `2to3`_ tool has been
phased out in favor of using ``six``.

To start using ``six`` instead of ``2to3`` in an affiliate package, you first
need to put the following in the package's ``setup_package.py`` file::

    def requires_2to3():
        return False

This section is mainly about moving existing code that works with
``2to3`` to using ``six``.  It is not a complete guide to Python 2 and
3 compatibility.

Welcome to the ``__future__``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The top of every ``.py`` file should include the following::

    from __future__ import (absolute_import, division, print_function,
                            unicode_literals)

This will make the Python 2 interpreter behave as close to Python 3 as
possible.

All files should also import `six`_, whether they are using it or not,
just to make moving code between modules easier, as `six`_ gets used *a
lot*::

    from ..extern import six

(where ``extern`` refers to ``astropy.extern``).  Do not import `six`_
from the top-level: only import the copy of `six`_ included with
astropy.

Finding places to use six
^^^^^^^^^^^^^^^^^^^^^^^^^

Unfortunately, the only way to be certain that code works on both
Python 2 and 3 is to make sure it is covered by unit tests.

However, the `2to3`_ commandline tool can also be used to locate places
that require special handling with `six`_.  Starting from Python 2
code, or code that is known to work on both Python 2 and 3 by
processing it with the `2to3`_ tool (which is most of the existing code
in astropy), simply run `2to3`_ on the file to display the changes it
would make in diff format.  This diff can be used to highlight places
that need to be updated to use `six`_.

For example, most things that have been renamed between Python 2 and 3
are in the ``six.moves`` namespace, so given this Python 2 code::

    import cPickle

it can be replaced with::

    from ..extern.six.moves import cPickle

.. note::

    The `modernize <https://pypi.python.org/pypi/modernize>`_ tool
    aims to convert Python 2 code to portable code that uses `six`_,
    but at the time of this writing, it is not feature complete and
    misses many important transformations.

The `six <http://pythonhosted.org/six/>`__ documentation serves as a
good reference for the sorts of things that need to be updated.

Not so fast on that Unicode thing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By importing ``unicode_literals`` from ``__future__``, many things
that were once byte strings on Python 2 will now by unicode strings.
This is mostly a good thing, as the behavior of the code will be more
consistent between Python 2 and 3.  However, certain third-party
libraries still assume certain values will be byte strings on
Python 2.

For example, when specifying Numpy structured dtypes, all strings must
be byte strings on Python 2 and unicode strings on Python 3.  The
easiest way to handle this is to force cast them using ``str()``, for
example::

   x = np.array([1.0, 2.0, 3.0], dtype=[(str('name'), '>f8')])

The same is true of structure specifiers in the built-in `struct`
module on Python 2.6.

``pytest.mark.skipif`` also requires a "native" string, i.e.::

    @pytest.mark.skipif(str('CONDITIONAL'))

Iteration
^^^^^^^^^

The behavior of the methods for iterating over the items, values and
keys of a dictionary has changed in Python 3.  Additionally, other
built-in functions such as `zip`, `range` and `map` have changed to
return iterators rather than temporary lists.

In many cases, the performance implications of iterating vs. creating
a temporary list won't matter, so it's tempting to use the form that
is simplest to read.  However, that results in code that behaves
differently on Python 2 and 3, leading to subtle bugs that may not be
detected by the regression tests.  Therefore, unless the loop in
question is provably simple and doesn't call into other code, the
`six`_ versions that ensure the same behavior on both Python 2 and 3
should be used.  The following table shows the mapping of equivalent
semantics between Python 2, 3 and six for ``dict.items()``:

============================== ============================== ==============================
Python 2                       Python 3                       six
============================== ============================== ==============================
``d.items()``                  ``list(d.items())``            ``list(six.iteritems(d))``
``d.iteritems()``              ``d.items()``                  ``six.iteritems(d)``
============================== ============================== ==============================

The above table holds true, analogously, for ``values``, ``keys``,
``zip``, ``range`` and ``map``.

Note that for keys only, ``list(d)`` is an acceptable shortcut to
``list(six.iterkeys(d))``.

Issues with ``\u`` escapes
^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``from __future__ import unicode_literals`` is used, all string
literals (not preceded with a ``'b'``) will become unicode literals.

Normally, one would use "raw" string literals to encode strings that contain
a lot of slashes that we don't want Python to interpret as special
characters.  Unfortunately, on Python 2, there is no way to represent
``'\u'`` in a raw unicode string literal, since it will always be
interpreted as the start of a unicode character escape, such as
``'\u20af'``.  The only solution is to use a regular (non-raw) string
literal and repeat all slashes, e.g. ``"\\usepackage{foo}"``.

The following shows the problem on Python 2::

    >>> ur'\u'
      File "<stdin>", line 1
    SyntaxError: (unicode error) 'rawunicodeescape' codec can't decode bytes in
    position 0-1: truncated \uXXXX
    >>> ur'\\u'
    u'\\\\u'
    >>> u'\u'
      File "<stdin>", line 1
    SyntaxError: (unicode error) 'unicodeescape' codec can't decode bytes in
    position 0-1: truncated \uXXXX escape
    >>> u'\\u'
    u'\\u'

This bug has been fixed in Python 3, however, we can't take advantage
of that and still support Python 2::

    >>> r'\u'
    '\\u'
    >>> r'\\u'
    '\\\\u'
    >>> '\u'
      File "<stdin>", line 1
    SyntaxError: (unicode error) 'unicodeescape' codec can't decode bytes in
    position 0-1: truncated \uXXXX escape
    >>> '\\u'
    '\\u'

Compatibility between versions of Numpy
---------------------------------------

In general, code should aim to be compatible with the lowest supported version
of NumPy_.  Sometimes, however, it is inefficient to code repeatedly around
bugs in earlier versions. For those cases, code can be added to
`astropy.utils.compat.numpy`; see the corresponding :ref:`instructions
<numpy-compatibility>` for details.

Requirements Specific to Affiliated Packages
--------------------------------------------

* Affiliated packages implementing many classes/functions not relevant to
  the affiliated package itself (for example leftover code from a previous
  package) will not be accepted - the package should only include the
  required functionality and relevant extensions.

* Affiliated packages are required to follow the layout and documentation form
  of the template package included in the core package source distribution.

* Affiliated packages must be registered on the `Python Package Index
  <http://pypi.python.org/pypi>`_, with proper metadata for downloading and
  installing the source package.

* The ``astropy`` root package name should not be used by affiliated
  packages - it is reserved for use by the core package. Recommended naming
  conventions for an affiliated package are either simply ``packagename``
  or ``awastropy.packagename`` ("affiliated with Astropy").

Examples
--------

This section shows a few examples (not all of which are correct!) to
illustrate points from the guidelines. These will be moved into the template
project once it has been written.

.. _prop-get-set-example:

Properties vs. get\_/set\_
^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows a sample class illustrating the guideline regarding the use
of properties as opposed to getter/setter methods.

Let's assuming you've defined a ``':class:`Star`'`` class and create an instance
like this::

    >>> s = Star(B=5.48, V=4.83)

You should always use attribute syntax like this::

    >>> s.color = 0.4
    >>> print(s.color)
    0.4

Rather than like this::

    >>> s.set_color(0.4)  #Bad form!
    >>> print(s.get_color())  #Bad form!
    0.4

Using python properties, attribute syntax can still do anything possible with
a get/set method. For lengthy or complex calculations, however, use a method::

    >>> print(s.compute_color(5800, age=5e9))
    0.4

.. _super-vs-direct-example:

super() vs. Direct Calling
^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows why the use of :func:`super` leads to a more consistent
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
called next.  By calling :func:`super` the entire method resolution order for
``D`` is precomputed, enabling each superclass to cooperatively determine which
class should be handed control in the next :func:`super` call::

    # This is safer

    class A(object):
        def method(self):
            print('Doing A')

    class B(A):
        def method(self):
            print('Doing B')
            super(B, self).method()


    class C(A):
        def method(self):
            print('Doing C')
            super(C, self).method()

    class D(C, B):
        def method(self):
            print('Doing D')
            super(D, self).method()

::

    >>> d = D()
    >>> d.method()
    Doing D
    Doing C
    Doing B
    Doing A

As you can see, each superclass's method is entered only once.  For this to
work it is very important that each method in a class that calls its
superclass's version of that method use :func:`super` instead of calling the
method directly.  In the most common case of single-inheritance, using
``super()`` is functionally equivalent to calling the superclass's method
directly.  But as soon as a class is used in a multiple-inheritance
hierarchy it must use ``super()`` in order to cooperate with other classes in
the hierarchy.

.. note:: For more information on the the benefits of :func:`super`, see
          http://rhettinger.wordpress.com/2011/05/26/super-considered-super/

.. _import-star-example:

Acceptable use of ``from module import *``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
used, it is strongly recommended that the submodules make use of the __all__
variable to specify which modules should be imported. Thus, submodule2.py
might read::

    from numpy import array,linspace

    __all__ = ('foo', 'AClass')

    def foo(bar):
        #the function would be defined here
        pass

    class AClass(object):
        #the class is defined here
        pass

This ensures that ``from submodule import *`` only imports ``':func:`foo'``
and ``':class:`AClass'``, but not ``':class:`numpy.array'`` or
``':func:`numpy.linspace'``.

.. _try-except-as-example:

try...except block "as" syntax
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Catching of exceptions should always use this syntax::

    try:
        ... some code that might produce a variety of exceptions ...
    except ImportError as e:
        if 'somemodule' in e.args[0]"
            #for whatever reason, failed import of somemodule is ok
            pass
        else:
            raise
    except ValueError, TypeError as e:
        msg = 'Hit an input problem, which is ok,'
        msg2 = 'but we're printing it here just so you know:'
        print msg, msg2, e

This avoids the old style syntax of ``except ImportError, e`` or
``except (ValueError,TypeError), e``, which is dangerous because it's easy to
instead accidentally do something like ``except ValueError,TypeError``, which
won't catch `~.exceptions.TypeError`.


Additional Resources
--------------------

Further tips and hints relating to the coding guidelines are included below.

.. toctree::
    :maxdepth: 1

    codeguide_emacs

.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://www.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _Cython: http://cython.org/
.. _PyPI: http://pypi.python.org/pypi
