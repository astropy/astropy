.. _using-git:

=============================
How To: Contribute to Astropy
=============================

Summary
-------

Any contributions to the core Astropy package, whether bug fixes,
improvements to the documentation, or new functionality, can be done via
*pull requests* on GitHub. The workflow for this is described below.

However, substantial contributions, such as whole new sub-packages, can
first be developed as *affiliated packages* then submitted to the Astropy
core via pull requests once ready (as described in :ref:`vision`). Teams working on affiliated packages are free to choose whatever version control system they wish, but ultimately the affiliated package should be merged into a fork of the Astropy repository in order to be submitted as a pull request (this merging can be done either by the team or by one of the core maintainers).

Getting started with git
------------------------

The only absolutely necessary configuration step is identifying yourself and
your contact info::

     git config --global user.name "Your Name"
     git config --global user.email you@yourdomain.example.com

The following sections cover the installation of the git software, the basic
configuration, and links to resources to learn more about using git.
However, you can also directly go to the `GitHub help pages
<http://help.github.com/>`_ which offer a great introduction to git and
GitHub.

.. toctree::
   :maxdepth: 1

   git_install
   git_configure
   git_resources

Workflow
--------

The following two sections describe the workflow for the Astropy core
package, but teams working on affiliated packages that have chosen to use
git are encouraged to also follow these guidelines internally.

.. toctree::
   :maxdepth: 1

   development_workflow
   maintainer_workflow

If for any reason developers do not wish to or cannot contribute via pull requests, they can submit a patch as described in :doc:`patches`.