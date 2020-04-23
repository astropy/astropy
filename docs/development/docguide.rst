.. _documentation-guidelines:

*********************
Writing Documentation
*********************

High-quality, consistent documentation for astronomy code is one of the major
goals of the Astropy Project.  Hence, we describe our documentation procedures
and rules here.  For the astropy core project and coordinated packages we try to
keep to these as closely as possible, and we encourage affiliated packages to
also adhere to these as they encourage useful documentation, a characteristic
often lacking in professional astronomy software.

Adding a Git Commit
===================

When your changes only affect documentation (i.e., docstring or RST files)
and do not include any code snippets that require doctest to run, you may
add a ``[skip travis]`` in your commit message. For example::

    git commit -m "Update documentation about this and that [skip travis]"

When this commit is pushed out to your branch associated with a pull request,
Travis CI will be skipped because it is not required. This is because the
CI job to build the documentation resides in CircleCI.


Building the Documentation from source
======================================

For information about building the documentation from source, see
the :ref:`builddocs` section in the installation instructions.

Astropy Documentation Rules and Guidelines
==========================================

This section describes the standards for documentation that any contribution
being considered for integration into the core package should follow, as well as
the standard Astropy docstring format.

* All documentation text should follow the :ref:`astropy-style-guide`.

* All documentation should be written use the `Sphinx`_
  documentation tool.

* The `package template <https://github.com/astropy/package-template>`_ provides
  a recommended general structure for documentation.

* Docstrings must be provided for all public classes, methods, and functions.

* Docstrings should follow the `numpydoc format
  <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

* Examples and/or tutorials are strongly encouraged for typical use-cases of a
  particular module or class.

* Any external package dependencies must be explicitly mentioned in the
  documentation. They should also be recorded in the ``setup.cfg`` file in the
  root of the astropy repository using an ``extras_require`` entry.

* Configuration options using the :mod:`astropy.config` mechanisms must be
  explicitly mentioned in the documentation.


Sphinx Documentation Themes
===========================

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
=================

The documentation build process for Astropy uses a number of sphinx extensions
which are all installed automatically when installing sphinx-astropy_. These
facilitate easily documenting code in a homogeneous and readable way.

The main extensions used are:

* sphinx-automodapi_ - an extension
  that makes it easy to automatically generate API documentation.

* sphinx-gallery_ - an
  extension to generate example galleries

* numpydoc_ - an extension to parse
  docstrings in NumpyDoc format

In addition, the sphinx-astropy_ includes a few small extensions:

* ``sphinx_astropy.ext.edit_on_github`` - an extension to add 'Edit on GitHub'
  links to documentation pages.

* ``sphinx_astropy.ext.changelog_links`` - an extension to add links to
  pull requests when rendering the changelog.

* ``sphinx_astropy.ext.doctest`` - an extension that makes it possible to
  add metadata about doctests inside ``.rst`` files

.. _NumPy: https://www.numpy.org/
.. _numpydoc: https://pypi.org/project/numpydoc/
.. _Matplotlib: https://matplotlib.org/
.. _SciPy: https://www.scipy.org/
.. _Sphinx: http://www.sphinx-doc.org/
.. _sphinx-automodapi: https://github.com/astropy/sphinx-automodapi
.. _astropy-sphinx-theme: https://github.com/astropy/astropy-sphinx-theme
.. _sphinx-astropy: https://github.com/astropy/sphinx-astropy
.. _sphinx-gallery: https://sphinx-gallery.readthedocs.io
