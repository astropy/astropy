=======
Astropy
=======

.. container::

    |Actions Status| |CircleCI Status| |Coverage Status| |PyPI Status| |Documentation Status| |Pre-Commit| |isort Status| |Zenodo|

The Astropy Project (http://astropy.org/) is a community effort to develop a
single core package for Astronomy in Python and foster interoperability between
Python astronomy packages. This repository contains the core package which is
intended to contain much of the core functionality and some common tools needed
for performing astronomy and astrophysics with Python.

Releases are `registered on PyPI <https://pypi.org/project/astropy>`_,
and development is occurring at the
`project's GitHub page <http://github.com/astropy/astropy>`_.

For installation instructions, see the `online documentation <https://docs.astropy.org/>`_
or  `docs/install.rst <docs/install.rst>`_ in this source distribution.

Contributing Code, Documentation, or Feedback
---------------------------------------------

The Astropy Project is made both by and for its users, so we welcome and
encourage contributions of many kinds. Our goal is to keep this a positive,
inclusive, successful, and growing community by abiding with the
`Astropy Community Code of Conduct <http://www.astropy.org/about.html#codeofconduct>`_.

More detailed information on contributing to the project or submitting feedback
can be found on the `contributions <http://www.astropy.org/contribute.html>`_
page. A `summary of contribution guidelines <CONTRIBUTING.md>`_ can also be
used as a quick reference when you are ready to start writing or validating
code for submission.

Supporting the Project
----------------------

|NumFOCUS| |Donate|

The Astropy Project is sponsored by NumFOCUS, a 501(c)(3) nonprofit in the
United States. You can donate to the project by using the link above, and this
donation will support our mission to promote sustainable, high-level code base
for the astronomy community, open code development, educational materials, and
reproducible scientific research.

License
-------

Astropy is licensed under a 3-clause BSD style license - see the
`LICENSE.rst <LICENSE.rst>`_ file.

.. |Actions Status| image:: https://github.com/astropy/astropy/workflows/CI/badge.svg
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

.. |isort Status| image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
    :target: https://pycqa.github.io/isort/
    :alt: isort Status

.. |NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS

.. |Donate| image:: https://img.shields.io/badge/Donate-to%20Astropy-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-astropy/index.html


If you locally cloned this repo before 7 Apr 2021
-------------------------------------------------

The primary branch for this repo has been transitioned from ``master`` to
``main``.  If you have a local clone of this repository and want to keep your
local branch in sync with this repo, you'll need to do the following in your
local clone from your terminal::

   git fetch --all --prune
   # you can stop here if you don't use your local "master"/"main" branch
   git branch -m master main
   git branch -u origin/main main

If you are using a GUI to manage your repos you'll have to find the equivalent
commands as it's different for different programs. Alternatively, you can just
delete your local clone and re-clone!
