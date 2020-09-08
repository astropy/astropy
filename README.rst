=======
Astropy
=======

|Travis Status| |CircleCI Status| |Coverage Status| |PyPI Status| |Documentation Status|

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

Notes for Package Managers
--------------------------

For system packagers: Please install `astropy` with the command::

    $ python setup.py --offline install

This will prevent the astropy_helpers bootstrap script from attempting to
reach out to PyPI.

.. |Travis Status| image:: https://img.shields.io/travis/astropy/astropy/master?logo=travis%20ci&logoColor=white&label=Travis%20CI
    :target: https://travis-ci.org/astropy/astropy
    :alt: Astropy's Travis CI Status

.. |CircleCI Status| image::  https://img.shields.io/circleci/build/github/astropy/astropy/master?logo=circleci&label=CircleCI
    :target: https://circleci.com/gh/astropy/astropy
    :alt: Astropy's CircleCI Status

.. |Coverage Status| image:: https://codecov.io/gh/astropy/astropy/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/astropy/astropy
    :alt: Astropy's Coverage Status

.. |PyPI Status| image:: https://img.shields.io/pypi/v/astropy.svg
    :target: https://pypi.org/project/astropy
    :alt: Astropy's PyPI Status

.. |Documentation Status| image:: https://img.shields.io/readthedocs/astropy/latest.svg?logo=read%20the%20docs&logoColor=white&label=Docs&version=stable
    :target: https://docs.astropy.org/en/stable/?badge=stable
    :alt: Documentation Status

.. |NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS

.. |Donate| image:: https://img.shields.io/badge/Donate-to%20Astropy-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-astropy/index.html
