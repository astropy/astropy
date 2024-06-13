.. _documentation-guidelines:

*********************
Writing Documentation
*********************

High-quality, consistent documentation for astronomy code is one of the major goals of
the Astropy Project. We encourage you to help us improve the documentation, and you
don't have to be an expert on ``astropy`` to do so!  If something in the docs doesn't
make sense to you, updating the relevant section after you figure it out is a great way
to contribute to Astropy and help others in the future.

In this section we describe our documentation procedures and rules that apply to the
``astropy`` core package and coordinated packages. We encourage affiliated packages to
do the same.

Astropy Documentation Overview
==============================

The documentation is written in **reStructuredText**, which is almost like writing in
plain English, and built using `Sphinx <https://www.sphinx-doc.org/en/master/>`_. The
Sphinx Documentation has an excellent `introduction to reST
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ that you should read.

The Astropy documentation consists of two parts. The code docstrings provide a clear
explanation of the usage of individual functions, classes and modules. These are
collected by Sphinx into the API documentation. The narrative documentation is contained
in the source code ``docs`` directory and provides a more tutorial-like overview of each
sub-package together with other topical information like What's New, Installation,
Developer documentation, etc.

The |OpenAstronomy Packaging Guide| provides a recommended general structure for
documentation.

Building the Documentation from Source
======================================

For information about building the documentation from source, see
the :ref:`builddocs` section in the installation instructions.

Astropy Documentation Guidelines
================================

This section describes the standards for documentation that any contribution
being considered for integration into the core package should follow, as well as
the standard Astropy docstring format.

* Documentation text should follow the :ref:`astropy-style-guide`.

* Docstrings must be provided for all public classes, methods, and functions, and be
  written using the `numpydoc format
  <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

* References in docstrings, **including internal Astropy links**, should use the
  `intersphinx format
  <https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html>`_.
  For example a link to the Astropy section on unit equivalencies would be
  `` :ref:`astropy:unit_equivalencies` ``.
  When built in Astropy, links starting with 'astropy' resolve to the current
  build. In affiliated packages using the ``sphinx-astropy`` intersphinx mapping,
  the links resolve to the stable version of Astropy. For linking to the
  development version, use the intersphinx target 'astropy-dev'.

* Examples and/or tutorials are strongly encouraged for typical use-cases of a
  particular module or class.

* Optional package dependencies should be documented where feasible.

* Configuration options using the :mod:`astropy.config` mechanisms must be
  explicitly mentioned in the documentation.

Details for Package Maintainers
===============================

Following is useful information for package maintainers who are using the Astropy
documentation infrastructure and may want to customize it for their package.

Sphinx Documentation Themes
---------------------------

An Astropy Project Sphinx HTML theme is included in the astropy-sphinx-theme_
package. This allows the theme to be used by both Astropy and affiliated
packages. The theme is activated by setting the theme in the global Astropy
sphinx configuration in sphinx-astropy_, which is imported in the sphinx
configuration of both Astropy and affiliated packages.

A different theme can be used by overriding a few sphinx
configuration variables set in the global configuration.

* To use a different theme, set ``html_theme`` to the name of a desired
  builtin Sphinx theme or a custom theme in ``package-name/docs/conf.py``
  (where ``'package-name'`` is "astropy" or the name of the affiliated
  package).

* To use a custom theme, additionally: place the theme in
  ``package-name/docs/_themes`` and add ``'_themes'`` to the
  ``html_theme_path`` variable. See the Sphinx_ documentation for more
  details on theming.

Sphinx extensions
-----------------

The documentation build process for Astropy uses a number of sphinx extensions
which are all installed automatically when installing sphinx-astropy_. These
facilitate easily documenting code in a homogeneous and readable way.

The main extensions used are:

* sphinx-automodapi_ - an extension that makes it easy to automatically
  generate API documentation.

* sphinx-gallery_ - an extension to generate example galleries

* |numpydoc| - an extension to parse docstrings in NumpyDoc format

In addition, the sphinx-astropy_ includes a few small extensions:

* ``sphinx_astropy.ext.edit_on_github`` - an extension to add 'Edit on GitHub'
  links to documentation pages.

* ``sphinx_astropy.ext.changelog_links`` - an extension to add links to
  pull requests when rendering the changelog.

* ``sphinx_astropy.ext.doctest`` - an extension that makes it possible to
  add metadata about doctests inside ``.rst`` files

.. _Sphinx: http://www.sphinx-doc.org/
.. _sphinx-automodapi: https://github.com/astropy/sphinx-automodapi
.. _astropy-sphinx-theme: https://github.com/astropy/astropy-sphinx-theme
.. _sphinx-astropy: https://github.com/astropy/sphinx-astropy
.. _sphinx-gallery: https://sphinx-gallery.readthedocs.io
