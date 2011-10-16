===========================
Coding Guidelines (Draft 3)
===========================

.. warning::
    This document is currently in Draft form and is subject to change.

This section describes requirements and guidelines that affiliated packages
will have to follow before being considered for integration as a module in the
core package.

Interface and Dependencies
--------------------------

* The package should meet the interface requirements set out by the
  coordination committee.

* Packages implementing many classes/functions not relevant to the component
  requested will not be accepted - the package should only include the
  required functionality and relevant extensions.

* Packages must be compatible with Python 2.6, 2.7, and 3.x (for 3.x
  compatibility, the `2to3 tool <http://docs.python.org/library/2to3.html>`_
  will be used).

* The package should be importable with no dependencies other than components
  already in the Astropy core, the `Python Standard Library (v2.6)
  <http://docs.python.org/release/2.6/library/index.html>`_, NumPy_, SciPy_,
  and Matplotlib_ (versions for these packages will be specified in the
  Astropy ``setup.py`` file and on PyPI_).

* The package should be importable from the source tree at build time.
  This means that, for example, if the package relies on C extensions
  that have yet to be built, the Python code is still importable, even
  if none of its functionality will work.

* Additional dependencies are allowed for sub-modules or in function calls,
  but they must be noted in the package documentation and should only affect
  the relevant component.

* General utilities necessary for but not specific to the package should be
  placed in the :mod:`packagename.utils` module. These utilities will be moved
  to the :mod:`astropy.utils` module when the package is integrated into the
  core package. If a utility is already present in :mod:`astropy.utils`, the
  package should always use that utility instead of re-implementing it in
  :mod:`packagename.utils`.

Documentation and Testing
-------------------------

* Docstrings must be present for all public classes/methods/functions, and
  must follow the form outlined in the :doc:`docguide` document.

* Unit tests are encouraged for all public methods and functions, and should
  adhere to the standards set in the :doc:`testguide` document.

Data and Configuration
----------------------

* Packages can include data in ``path TBD`` as long as it is less than about
  100 kb. These data should be accessed via the
  :func:`astropy.config.[funcname TBD]` mechanism. If the data exceeds this
  size, it should be hosted outside the source code repository and downloaded
  using the :func:`astropy.config.[funcname TBD]` mechanism. Exceptions to
  this size limit may be allowed if there are version dependencies between
  data and code.

* All persistent configuration should be stored using the functions in
  :mod:`astropy.config`, which make use of the :class:`ConfigObj` class and
  associated file format (http://www.voidspace.org.uk/python/configobj.html.

Coding Style/Conventions
------------------------

* The code will follow the standard `PEP8 Style Guide for Python Code
  <http://www.python.org/dev/peps/pep-0008/>`_. In particular, this includes
  using only 4 spaces for indentation, and never tabs.

* One exception is to be made from the PEP8 style: new style relative imports
  of the form ``from . import modname`` are allowed and required for Astropy,
  as opposed to absolute (as PEP8 suggets) or the simpler ``import modname``
  syntax. This is primarily due to improved relative import support since PEP8
  was developed, and to simplify the process of moving modules.

.. note:: A pep8.py checker script is available at
          http://pypi.python.org/pypi/pep8.

          See :doc:`codeguide_emacs` for some configuration options for
          Emacs that helps in ensuring conformance to PEP8.

* The ``import numpy as np``, ``import matplotlib as mpl``, and ``import
  matplotlib.pyplot as plt`` naming conventions should be used wherever
  relevant. ``from packagename import *`` should never be used, except as a
  tool to flatten the namespace of a module. An example of the allowed usage
  is given below.

* Classes should either use direct variable access, or python’s property
  mechanism for setting object instance variables. ``get_value``/``set_value``
  style methods should be used only when getting and setting the values
  requires a computationally-expensive operation. The example below
  illustrates this guideline.

* All new classes should be new-style classes inheriting from :class:`object`
  (in Python 3 this is a non-issue as all classes are new-style by default).
  The one exception to this rule is older classes in third-party libraries such
  the Python standard library or numpy.

* Classes should use the builtin :func:`super` function when making calls to
  methods in their super-class(es) unless there are specific reasons not to.
  :func:`super` should be used consistently in all subclasses since it does not
  work otherwise.  An example illustrating why this is important is included
  below.

* Multiple inheritance should be avoided in general without good reason.
  Mulitple inheritance is complicated to implement well, which is why many
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
  in accord with the guideline above). 

* Affiliated packages are required to follow the layout and documentation form
  of the template package included in the core package source distribution.

Including C code
----------------

* C extensions are only allowed when they provide a significant performance
  enhancement over pure python, or a robust C library already exists to
  provided the needed functionality. When C extensions are used, the Python
  interface must meet interface guidelines.

* The use of Cython_ is strongly recommended for C extensions, as per the
  example in the template package. Cython extensions should store ``.pyx``
  files in the source code repository, but they should be compiled to ``.c``
  files that are updated in the repository when important changes are made to
  the ``.pyx`` file.

* If a C extension has a dependency on an external C library, the source code
  for the library should be bundled with the Astropy core, provided the
  license for the C library is compatible with the Astropy license.
  Additionally, the package must be compatible with using a system-installed
  library in place of the library included in Astropy.

* In cases where C extensions are needed but Cython cannot be used, the `PEP 7
  Style Guide for C Code <http://www.python.org/dev/peps/pep-0007/>`_ is
  recommended.

Requirements specific to Affiliated Packages
--------------------------------------------

* Affiliated packages must be registered on the `Python Package Index
  <http://pypi.python.org/pypi>`_, with proper metadata for downloading and
  installing the source package.

* The :mod:`astropy` root package name should not be used by affiliated
  packages - it is reserved for use by the core package. Recommended naming
  conventions for an affiliated package are either simply :mod:`packagename`
  or :mod:`awastropy.packagename` ("affiliated with Astropy").

Examples
--------

This section shows a few examples (not all of which are correct!) to
illustrate points from the guidelines. These will be moved into the template
project once it has been written.

Properties vs. get\_/set\_
^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows a sample class illustrating the guideline regarding the use
of properties as opposed to getter/setter methods.

Let's assuming you've defined a :class:`Star` class and create an instance
like this::

    >>> s = Star(B=5.48, V=4.83)

You should always use attribute syntax like this::

    >>> s.color = 0.4
    >>> print s.color
    0.4

Rather than like this::

    >>> s.set_color(0.4)  #Bad form!
    >>> print s.get_color()  #Bad form!
    0.4

Using python properties, attribute syntax can still do anything possible with
a get/set method. For lengthy or complex calculations, however, use a method::

    >>> print s.compute_color(5800, age=5e9)
    0.4

super() vs. direct calling
^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows why the use of :func:`super` leads to a more consistent
method resolution order than manually calling methods of the super classes in a
multiple inheritance case::

    # This is dangerous and bug-prone!

    class A(object):
        def method(self):
            print 'Doing A'


    class B(A):
        def method(self):
            print 'Doing B'
            A.method(self)


    class C(A):
        def method(self):
            print 'Doing C'
            A.method(self)

    class D(C, B):
        def method(self):
            print 'Doing D'
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
that inherts from both ``B`` and ``C``, and that ``B.method()`` should be
called next.  By calling :func:`super` the entire method resolution order for
``D`` is precomputed, enabling each superclass to cooperatively determine which
class should be handed control in the next :func:`super` call::

    # This is safer

    class A(object):
        def method(self):
            print 'Doing A'

    class B(A):
        def method(self):
            print 'Doing B'
            super(B, self).method()


    class C(A):
        def method(self):
            print 'Doing C'
            super(C, self).method()

    class D(C, B):
        def method(self):
            print 'Doing D'
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

.. note:: For more info on the pros and cons of using super, see
          http://rhettinger.wordpress.com/2011/05/26/super-considered-super/
          or http://keithdevens.com/weblog/archive/2011/Mar/16/Python.super)


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

    __all__ = ('foo','AClass')

    def foo(bar):
        #the function would be defined here
        pass

    class AClass(object):
        #the class is defined here
        pass

This ensures that ``from submodule import *`` only imports :func:`foo` and
:class:`AClass`, but not :class:`numpy.array` or :func:`numpy.linspace`.


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
