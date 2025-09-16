.. _developer-docs:

************
Contributing
************

The contributor documentation contains instructions for how to contribute to ``astropy`` or
affiliated packages. This includes setting up a development environment, installing and
testing the development version, as well as coding, documentation, and testing
guidelines.

For newcomers the process may initially seem overwhelming, but with a little patience
and practice you will see that it is not so complex. The key is to follow the steps
outlined here and `ask for help <https://www.astropy.org/help.html>`_ if you get stuck.
The Astropy community is welcoming and friendly and will help you!

{% if is_development %}

This is divided into two sections, first a quickstart guide that provides an
introduction to the development workflow, followed by a number of detailed guides that
cover provide a deeper dive and a reference for both developers and maintainers.

.. Important:: There are useful ways to contribute to Astropy without diving
    into the developer workflow which is described here. For an
    an overview see the `Contribute to Astropy <https://www.astropy.org/contribute.html>`_
    page.


Contributing quickstart
-----------------------

This section provides a contributing quickstart guide for Astropy. With minor changes the
process will apply to contributing updates to coordinated and many affiliated packages.

.. toctree::
   :maxdepth: 2

   development/quickstart

Now that you have created your development environment and gotten familiar with the
process, you should now read through the detailed tutorial below to see a real-life
example of a simple bug fix. This includes more explanation of the steps and good
advice for making a code change.

.. toctree::
   :maxdepth: 1

   development/git_edit_workflow_examples

Congratulations, now you are ready to be an Astropy contributor! If you are not sure where to contribute, take a look at the `Good First Issues
<https://github.com/astropy/astropy/issues?q=is%3Aopen+is%3Aissue+label%3A%22good+first+issue%22>`_
list. These issues are the most accessible ones if you are not familiar with the Astropy
source code.

Details
-------

.. toctree::
   :maxdepth: 1

   development/development_details
   development/codeguide
   development/testguide
   development/docguide
   development/style-guide
   development/git_resources
   development/scripts
   development/ccython
   development/maintainers/index

.. Note:: Parts of this guide were adapted from the
    `pandas developer documentation <https://pandas.pydata.org/pandas-docs/stable/development/index.html>`_. Astropy is grateful to the pandas team for their documentation efforts.

{%else%}

To read the developer documentation, you will need to go to the
`latest developer version of the documentation <https://docs.astropy.org/en/latest/index_dev.html>`_.

{%endif%}
