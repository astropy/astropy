.. _using-git:

Workflow for contributing to AstroPy
====================================

Summary
-------

As described in :ref:`vision`, development of components intended for
inclusion in the core ``astropy`` package will initially be done by different
teams via *affiliated packages*. These affiliated packages will then be
submitted for inclusion in the ``astropy`` core package by creating a git fork
of the core repository, merging the new component into the fork, and
submitting the package as a pull request. The step of creating a fork and
merging the affiliated package into the fork can be done either by the team
contributing the affiliated package, or by the coordination committee if
requested.

Once an affiliated package has been accepted and integrated as a component
into the core ``astropy`` package, subsequent improvements and bug fixes can
be made in the same way by forking the core repository and submitting a pull
request.

The bottom line is that teams working on various components are free to choose
the version control system and workflow that they want, but ultimately the
changes need to be merged into a fork of the core repository and submitted via
a pull request to the core repository, either by the team or by the
coordination committee.

Git Documentation
-----------------

The following sections cover the installation of the git software, the basic configuration, and links to resources to learn more about using git.

.. toctree::
   :maxdepth: 1

   git_install
   git_configure
   git_resources

Workflow
--------

The following two sections describe the workflow for the ``astropy`` core
package, but teams working on affiliated packages that have chosen to use git
are encouraged to also follow these guidelines internally.

.. toctree::
   :maxdepth: 1

   development_workflow
   maintainer_workflow
