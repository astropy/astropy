.. _developer-docs:

***********
Development
***********

The developer documentation contains instructions for how to contribute to
Astropy or affiliated packages, install and test the development version,
as well as coding, documentation, and testing guidelines.

{% if is_development %}

For the guiding vision of this process and the project
as a whole, see :doc:`vision`.

.. toctree::
   :maxdepth: 1

   workflow/development_workflow
   workflow/virtual_pythons
   workflow/get_devel_version
   codeguide
   docguide
   style-guide
   testguide
   testhelpers
   scripts
   building
   ccython
   releasing
   workflow/maintainer_workflow
   astropy-package-template

There are some additional tools, mostly of use for maintainers, in the
`astropy/astropy-tools repository
<https://github.com/astropy/astropy-tools>`__.

{%else%}

To read the developer documentation, you will need to go to the :ref:`latest
developer version of the documentation
<astropy-dev:developer-docs>`.

{%endif%}
