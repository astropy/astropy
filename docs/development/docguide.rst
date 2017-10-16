.. _documentation-guidelines:

*********************
Writing Documentation
*********************

High-quality, consistent documentation for astronomy code is one of
the major goals of the Astropy Project.  Hence, we describe our
documentation procedures and rules here.  For the astropy core
project we try to keep to these as closely as possible, while the
standards for affiliated packages are somewhat looser.
(These procedures and guidelines are still recommended for affiliated
packages, as they encourage useful documentation, a characteristic
often lacking in professional astronomy software.)


Building the Documentation from source
======================================

For information about building the documentation from source, see
the :ref:`builddocs` section in the installation instructions.


Astropy Documentation Rules and Guidelines
==========================================

This section describes the standards for documentation format affiliated
packages that must follow for consideration of integration into the core
module, as well as the standard Astropy docstring format.

* All documentation should be written use the Sphinx documentation tool.

* The template package will provide a recommended general structure for
  documentation.

* Docstrings must be provided for all public classes, methods, and functions.

* Docstrings will be incorporated into the documentation using a version of
  numpydoc included with Astropy, and should follow the :doc:`docrules`.

* Examples and/or tutorials are strongly encouraged for typical use-cases of a
  particular module or class.

* Any external package dependencies aside from NumPy_, SciPy_, or Matplotlib_
  must be explicitly mentioned in the documentation. They should also be
  recorded in the ``pip-requirements-doc`` file in the root of the astropy
  repository.

* Configuration options using the :mod:`astropy.config` mechanisms must be
  explicitly mentioned in the documentation.


The details of the docstring format are described on a separate page:

.. toctree::
    docrules


Sphinx Documentation Themes
===========================

A custom Sphinx HTML theme is included in the `astropy-helpers`_ package.
This allows the
theme to be used by both Astropy and affiliated packages. This is done by
setting the theme in the global Astropy sphinx configuration, which is imported
in the sphinx configuration of both Astropy and affiliated packages.

Using a different theme for ``astropy`` or affiliated packages
--------------------------------------------------------------

A different theme can be used by overriding a few sphinx
configuration variables set in the global configuration.

* To use a different theme, set ``'html_theme'`` to the name of a desired
  builtin Sphinx theme or a custom theme in ``package-name/docs/conf.py``
  (where ``'package-name'`` is "astropy" or the name of the affiliated
  package).

* To use a custom theme, additionally: place the theme in
  ``package-name/docs/_themes`` and add ``'_themes'`` to the
  ``'html_theme_path'`` variable. See the Sphinx_ documentation for more
  details on theming.

Adding more custom themes to astropy
------------------------------------

Additional custom themes can be included in the astropy source tree by
placing them in the directory ``astropy/astropy/sphinx/themes``, and
editing ``astropy/astropy/sphinx/setup_package.py`` to include the theme
(so that it is installed).



Sphinx extensions
=================

Astropy-helpers includes a number of sphinx extensions (some via the
`sphinx-automodapi`_ package) that are used in Astropy and its affiliated
packages to facilitate easily documenting code in a homogeneous and readable
way. The two main extensions are `~sphinx_automodapi.automodapi` for
generating module documentation and `~sphinx_automodapi.automodsumm` for
generating tables of module objects. Please see their documentation about
usage.


edit_on_github Extension
------------------------

.. automodule:: astropy_helpers.sphinx.ext.edit_on_github


numpydoc Extension
------------------
This extension (and some related extensions) are a port of the
`numpydoc <https://pypi.python.org/pypi/numpydoc/0.3.1>`_ extension
written by the NumPy_ and SciPy_, projects, with some tweaks for
Astropy.  Its main purposes is to reprocess docstrings from code into
a form sphinx understands. Generally, there's no need to interact with
it directly, as docstrings following the :doc:`docrules` will be
processed automatically.


Other Extensions
----------------

``astropy_helpers.sphinx.ext`` and `sphinx-automodapi`_ includes a few other
extensions that are primarily helpers for the other extensions or
workarounds for undesired behavior.  Their APIs are not included here
because we may change them in the future.


.. _NumPy: http://www.numpy.org/
.. _numpydoc: https://pypi.python.org/pypi/numpydoc/0.3.1
.. _Matplotlib: http://matplotlib.org/
.. _SciPy: https://www.scipy.org/
.. _Sphinx: http://sphinx.pocoo.org
.. _astropy-helpers: https://github.com/astropy/astropy-helpers
.. _sphinx-automodapi: https://github.com/astropy/sphinx-automodapi
