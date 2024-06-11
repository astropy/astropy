.. _developer-docs:

***********
Development
***********

The developer documentation contains instructions for how to contribute to
Astropy or affiliated packages, install and test the development version,
as well as coding, documentation, and testing guidelines.

{% if is_development %}

Contributing quickstart
-----------------------

This section provides a quickstart guide to contributing to Astropy. With minor
changes the process will apply to contributing to coordinated and many affiliated
packages.

.. toctree::
   :maxdepth: 2

   contributing
   contributing_environment
   contributing_documentation
   contributing_codebase
   maintaining
   internals
   copy_on_write
   debugging_extensions
   extending
   developer
   policies
   community

Legacy documentation
--------------------

.. toctree::
    :maxdepth: 2

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

   maintainers/maintainers-index

{%else%}

To read the developer documentation, you will need to go to the :ref:`latest
developer version of the documentation
<astropy-dev:developer-docs>`.

{%endif%}
