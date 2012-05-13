.. _using-git:

Contributing To/Developing Astropy or Affiliated Packages
=========================================================

Summary
-------

Development of components intended for
inclusion in the core ``astropy`` package can be done either via:

* **Affiliated packages** (as described in :ref:`vision`). Once ready for
  inclusion, affiliated packages should be submitted to the Astropy core
  package by creating a git fork of the core repository, merging the new
  component into the fork, and submitting the package as a pull request. The
  step of creating a fork and merging the affiliated package into the fork
  can be done either by the team contributing the affiliated package, or by
  the coordination committee if requested. Once an affiliated package has
  been accepted and integrated as a component into the core ``astropy``
  package, subsequent improvements and bug fixes can be made in the same way
  by forking the core repository and submitting pull requests.


* **Direct pull requests** on Github. While affiliated packages should be used
  for major new components of Astropy, smaller additions and bug fixes to
  the code or documentation can be done via a GitHub pull request.

The bottom line is that teams working on various components are free to choose
the version control system and workflow that they want, but ultimately the
changes need to be merged into a fork of the core repository and submitted via
a pull request to the core repository, either by the team or by the
coordination committee.

.. _git-configure-impatient:

For the impatient
-----------------

The only absolutely necessary configuration step is identifying yourself and
your contact info::

     git config --global user.name "Your Name"
     git config --global user.email you@yourdomain.example.com

More detailed information and instructions are below.

After that, if you then just want to get the latest ``astropy`` source code,
cd to a directory on your computer you want to put the source code, and do::

     git clone git@github.com:astropy/astropy.git

You will then have a new copy of the source code in the ``astropy``
directory.

Later, if you want to update to the most recent version of the ``astropy``
code, just go do::

   cd astropy
   git pull

If you find a bug and want to fix it, see :ref:`basic-workflow` or
:ref:`advanced-workflow` depending on how comfortable you are with git (the
former describes how to create a patch, while the latter explains how to
create a pull request).

Getting started with git
------------------------

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

The following two sections describe the workflow for the ``astropy`` core
package, but teams working on affiliated packages that have chosen to use git
are encouraged to also follow these guidelines internally.

.. toctree::
   :maxdepth: 1

   development_workflow
   maintainer_workflow
