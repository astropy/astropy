=======
Astropy
=======

.. container::

    |Actions Status| |CircleCI Status| |Coverage Status| |PyPI Status| |Documentation Status| |Pre-Commit| |Ruff| |Zenodo|

The Astropy Project (http://astropy.org/) is a community effort to develop a single core package for Astronomy in Python and foster interoperability between Python astronomy packages. This repository contains the core package which is intended to contain much of the core functionality and some common tools needed for performing astronomy and astrophysics with Python.

Table of Contents
=================

- `Installation <#installation>`_
- `Contributing <#contributing>`_
- `Getting Started with GitHub Codespaces <#getting-started-with-github-codespaces>`_
- `Supporting the Project <#supporting-the-project>`_
- `License <#license>`_
- `Updating Local Clone <#updating-local-clone>`_

Installation
============

Releases are `registered on PyPI <https://pypi.org/project/astropy>`_, and development is occurring at the `project's GitHub page <http://github.com/astropy/astropy>`_. For detailed installation instructions, see the `online documentation <https://docs.astropy.org/>`_ or `docs/install.rst <docs/install.rst>`_ in this source distribution.

To install Astropy, use:

.. code-block:: bash

    pip install astropy

Contributing
============

The Astropy Project is made both by and for its users, so we welcome and encourage contributions of many kinds. Our goal is to keep this a positive, inclusive, successful, and growing community by abiding with the `Astropy Community Code of Conduct <http://www.astropy.org/about.html#codeofconduct>`_.

More detailed information on contributing to the project or submitting feedback can be found on the `contributions page <http://www.astropy.org/contribute.html>`_. A `summary of contribution guidelines <CONTRIBUTING.md>`_ can also be used as a quick reference when you are ready to start writing or validating code for submission.

Getting Started with GitHub Codespaces
======================================

Codespaces is a cloud development environment supported by GitHub. None of the Astropy build machinery depends on it, but it is a convenient way to quickly get started doing development on Astropy.

To get started, create a codespace for this repository by clicking this:

|Codespaces|

A codespace will open in a web-based version of Visual Studio Code. The dev container is fully configured with software needed for this project. For help, see the `GitHub Codespaces Support page <https://docs.github.com/en/codespaces>`_.

Note: Dev containers is an open spec which is supported by GitHub Codespaces and other tools.

Supporting the Project
======================

.. |NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS

.. |Donate| image:: https://img.shields.io/badge/Donate-to%20Astropy-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-astropy/index.html

The Astropy Project is sponsored by NumFOCUS, a 501(c)(3) nonprofit in the United States. You can donate to the project using the link above. This donation will support our mission to promote a sustainable, high-level code base for the astronomy community, open code development, educational materials, and reproducible scientific research.

License
=======

Astropy is licensed under a 3-clause BSD style license - see the `LICENSE.rst <LICENSE.rst>`_ file.

Updating Local Clone
====================

If you locally cloned this repo before 7 Apr 2021, the primary branch for this repo has been transitioned from ``master`` to ``main``. If you have a local clone of this repository and want to keep your local branch in sync with this repo, you'll need to do the following in your local clone from your terminal:

.. code-block:: bash

    git fetch --all --prune
    # you can stop here if you don't use your local "master"/"main" branch
    git branch -m master main
    git branch -u origin/main main

If you are using a GUI to manage your repos you'll have to find the equivalent commands as it's different for different programs. Alternatively, you can just delete your local clone and re-clone!

.. |Actions Status| image:: https://github.com/astropy/astropy/actions/workflows/ci_workflows.yml/badge.svg
    :target: https://github.com/astropy/astropy/actions
    :alt: Astropy's GitHub Actions CI Status

.. |CircleCI Status| image::  https://img.shields.io/circleci/build/github/astropy/astropy/main?logo=circleci&label=CircleCI
    :target: https://circleci.com/gh/astropy/astropy
    :alt: Astropy's CircleCI Status

.. |Coverage Status| image:: https://codecov.io/gh/astropy/astropy/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/astropy/astropy
    :alt: Astropy's Coverage Status

.. |PyPI Status| image:: https://img.shields.io/pypi/v/astropy.svg
    :target: https://pypi.org/project/astropy
    :alt: Astropy's PyPI Status

.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4670728.svg
   :target: https://doi.org/10.5281/zenodo.4670728
   :alt: Zenodo DOI

.. |Documentation Status| image:: https://img.shields.io/readthedocs/astropy/latest.svg?logo=read%20the%20docs&logoColor=white&label=Docs&version=stable
    :target: https://docs.astropy.org/en/stable/?badge=stable
    :alt: Documentation Status

.. |Pre-Commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit

.. |Ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
    :target: https://github.com/astral-sh/ruff
    :alt: Ruff

.. |NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS

.. |Donate| image:: https://img.shields.io/badge/Donate-to%20Astropy-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-astropy/index.html

.. |Codespaces| image:: https://github.com/codespaces/badge.svg
    :target: https://github.com/codespaces/new?hide_repo_select=true&ref=main&repo=2081289
    :alt: Open in GitHub Codespaces

