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

  See :doc:`codeguide_emacs` for some configuration options for Emacs
  that helps in ensuring conformance to pep8.

.. note:: A pep8.py checker script is available at
          http://pypi.python.org/pypi/pep8.


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

* Classes are discouraged from using the builtin python :func:`super`
  function, unless absolutely needed. If used, it should be used
  consistentently by all subclasses, and noted in the class’s docstrings. An
  example illustrating why this is important (and alternative solutions) is
  included below.

* Affiliated packages are required to follow the layout and documentation form
  of the template package included in the core package source distribution.

.. note:: For more info on the pros and cons of using super, see
          http://rhettinger.wordpress.com/2011/05/26/super-considered-super/
          or http://keithdevens.com/weblog/archive/2011/Mar/16/Python.super)

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

This example shows why the use of :func:`super` can be confusing for
subclasses, and gives an alternative syntax::

    #This is dangerous and bug-prone!

    class A(object):
        def method(self):
            print 'Doing A'


    class B(A):
        def method(self):
            super(B, self).method()
            print 'Doing B'


    class C(A):
        def method(self):
            A.method(self)
            print 'Doing C'


    class D(C, B):
        def method(self):
            super(D, self).method()
            print 'Doing D'

if you then do::

    >>> b = B()
    >>> b.method()

you will see::

    Doing A
    Doing B

which is what you expect, and similarly for C. However, if you do::

    >>> d = D()
    >>> d.method()

you might expect to have it call both method in the order A,C,B,D. But it
doesn't - instead you see::

    Doing A
    Doing C
    Doing D

because the the ``A.method(self)`` in C effectively short-circuits the super
mechanism. Thus, it's crucial that all classes in an inheritance hierarchy
consistently use super and not mix super with the direct syntax. The simplest
approach is to explicitly call each class' method and avoid super completely::

    #This is safer
    class A(object):
        def __init__(self, a):
            self.a = 1


    class B(A):
        def __init__(self, a, b):
            A.__init__(self, a)
            self.b = b


    class C(A):
        def __init__(self, a, c):
            A.__init__(self, a)
            self.c = c


    class D(C, B):
        def __init__(self, a, b, c, d):
            B.__init__(self, a, b)
            C.__init__(self, a, c)
            self.d = d

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

.. _Numpy: http://numpy.scipy.org/
.. _Scipy: http://www.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _Cython: http://cython.org/
.. _PyPI: http://pypi.python.org/pypi
