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

Contributing to astropy
------------------------

.. toctree::
   :maxdepth: 2

New documentation
-----------------

.. toctree::
   :maxdepth: 2

   development/pd/index

Legacy documentation
--------------------

.. toctree::
    :maxdepth: 2

   development/workflow/development_workflow
   development/workflow/git_edit_workflow_examples
   development/workflow/virtual_pythons
   development/workflow/get_devel_version
   development/codeguide
   development/docguide
   development/style-guide
   development/testguide
   development/scripts
   development/building
   development/ccython

Git resources
-------------
.. toctree::
   :maxdepth: 2

   development/workflow/git_resources
   development/workflow/additional_git_topics


Maintaining astropy and affiliated packages
-------------------------------------------

.. toctree::
   :maxdepth: 2

   development/maintainers/maintainers-index

{%else%}

To read the developer documentation, you will need to go to the :ref:`latest
developer version of the documentation
<astropy-dev:developer-docs>`.

{%endif%}
