.. _developer-docs:

***********
Development
***********

The developer documentation contains instructions for how to contribute to
Astropy or affiliated packages, install and test the development version,
as well as coding, documentation, and testing guidelines.

{% if is_development %}

For the guiding vision of this process and the project
as a whole, see :doc:`development/vision`.

.. toctree::
   :maxdepth: 1

   development/workflow/development_workflow
   development/workflow/virtual_pythons
   development/workflow/get_devel_version
   development/codeguide
   development/docguide
   development/style-guide
   development/testguide
   development/testhelpers
   development/scripts
   development/building
   development/ccython
   development/releasing
   development/workflow/maintainer_workflow
   development/astropy-package-template

There are some additional tools, mostly of use for maintainers, in the
`astropy/astropy-tools repository
<https://github.com/astropy/astropy-tools>`__.

{%else%}

To read the developer documentation, you will need to go to the :ref:`latest
developer version of the documentation
<astropy-dev:developer-docs>`.

{%endif%}
