********
Overview
********

Here we describe a broad overview of the Astropy project and its parts.

Astropy Project Concept
=======================

The "Astropy Project" is distinct from the ``astropy`` package. The
Astropy Project is a process intended to facilitate communication and
interoperability of python packages/codes in astronomy and astrophysics.
The project thus encompasses the ``astropy`` core package (which provides
a common framework), all "affiliated packages" (described below in
`Affiliated Packages`_), and a general community aimed at bringing
resources together and not duplicating efforts.


``astropy`` Core Package
========================

The ``astropy`` package (alternatively known as the "core" package)
contains various classes, utilities, and a packaging framework intended
to provide commonly-used astronomy tools. It is divided into a variety
of sub-packages, which are documented in the remainder of this
documentation (see :ref:`user-docs` for documentation of these
components).

The core also provides this documentation, and a variety of utilities
that simplify starting other python astronomy/astrophysics packages. As
described in the following section, these simplify the process of
creating affiliated packages.


Affiliated Packages
===================

The Astropy project includes the concept of "affiliated packages." An
affiliated package is an astronomy-related python package that is not
part of the ``astropy`` core source code, but has requested to be included
in the general community effort of the Astropy project. Such a package
may be a candidate for eventual inclusion in the main ``astropy`` package
(although this is not required). Until then, however, it is a separate
package, and may not be in the ``astropy`` namespace.

The authoritative list of current affiliated packages is available at
http://affiliated.astropy.org, including a machine-readable `JSON file
<http://affiliated.astropy.org/registry.json>`_.

If you are interested in starting an affiliated package, or have a
package you are interested in making more compatible with astropy, the
``astropy`` core package includes features that simplify and homogenize
package management. Astropy provides a `package template
<http://github.com/astropy/package-template>`_ that provides a common
way to organize a package, to make your life simpler. You can use this
template either with a new package you are starting or an existing
package to give it most of the organizational tools Astropy provides,
including the documentation, testing, and Cython-building tools.  See
the `usage instructions in the template <https://github.com/astropy
/package-template/blob/master/README.rst>`_ for further details.

To then get your package listed on the registry, take a look at the
`guidelines for becoming an affiliated package
<http://affiliated.astropy.org#affiliated-instructions>`_, and then post
your intent on the `astropy-dev mailing list`_.  The Astropy
coordination committee, in consultation with the community, will provide
you feedback on the package, and will add it to the registry when it is
approved.


Community
=========

Aside from the actual code, Astropy is also a community of astronomy-
associated users and developers that agree that sharing utilities is
healthy for the community and the science it produces. This community is
of course central to accomplishing anything with the code itself. We
follow the `Python Software Foundation Code of Conduct
<http://www.python.org/psf/codeofconduct/>`_ and welcome anyone who
wishes to contribute to the project.
