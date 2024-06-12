.. _developer-docs:

***********
Development
***********

The developer documentation contains instructions for how to contribute to Astropy or
affiliated packages. This includes setting up a development environment, installing and
testing the development version, as well as coding, documentation, and testing
guidelines.

For newcomers the process may initially seem overwhelming, but with a little patience
and practice you will see that it is not so complex. The key is to follow the steps
outlined here and `ask for help <https://www.astropy.org/help.html>`_ if you get stuck.
The astropy community is welcoming and friendly and will help you!

This is divided into two sections, first a quickstart guide that provides an
introduction to the development workflow, followed by a more detailed guide that covers
all aspects of contributing to Astropy and provides a reference for both developers and
maintainers.

.. Important:: There are useful ways to contribute to Astropy without diving
    into the developer workflow which is described here. For an
    an overview see the `Contribute to Astropy <https://www.astropy.org/contribute.html>`_
    page.

Substantial parts of this content have been adapted from the excellent
`pandas developer documentation <https://pandas.pydata.org/pandas-docs/stable/development/index.html>`_.
Astropy is grateful to the pandas team for their documentation efforts.

{% if is_development %}

Contributing quickstart
-----------------------

This section provides a quickstart guide to contributing to Astropy. With minor
changes the process will apply to contributing to coordinated and many affiliated
packages.

Astropy is hosted on `GitHub <https://www.github.com/astropy/astropy>`_, and to
contribute, you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <https://git-scm.com/>`_ for
version control to allow many people to work together on the project. More details
and further resources are available in the :ref:`contributing.version_control` section.

.. toctree::
   :maxdepth: 2

   contributing_environment
   contributing_pull_request

Now that you have read through the instructions and created your development
environment, it is worthwhile reading  :ref:`astropy-fix-example` to see the process in
action with a simple bug fix.

Details
-------

.. toctree::
   :maxdepth: 2

   contributing_documentation
   contributing_codebase
   contributing_docstring
   workflow/development_workflow
   workflow/git_edit_workflow_examples
   workflow/virtual_pythons
   workflow/get_devel_version
   codeguide
   docguide
   style-guide
   testguide
   scripts
   building
   ccython

Git resources
-------------
.. toctree::
   :maxdepth: 2

   workflow/git_resources
   workflow/additional_git_topics


Maintaining astropy and affiliated packages
-------------------------------------------

.. toctree::
   :maxdepth: 2

   maintainers/astropy-package-template
   maintainers/maintainer_workflow
   maintainers/releasing
   maintainers/testhelpers
   maintainers/maintaining


{%else%}

To read the developer documentation, you will need to go to the :ref:`latest
developer version of the documentation
<astropy-dev:developer-docs>`.

{%endif%}
