=======
Astropy
=======

.. image:: https://img.shields.io/pypi/v/astropy.svg
    :target: https://pypi.python.org/pypi/astropy

Astropy (http://www.astropy.org) is a package intended to contain much of
the core functionality and some common tools needed for performing
astronomy and astrophysics with Python.

Releases are `registered on PyPI <http://pypi.python.org/pypi/astropy>`_,
and development is occurring at the
`project's github page <http://github.com/astropy/astropy>`_.

For installation instructions, see the `online documentation <http://docs.astropy.org/>`_
or  ``docs/install.rst`` in this source distribution.

For system packagers: Please install Astropy with the command::

    $ python setup.py --offline install

This will prevent the astropy_helpers bootstrap script from attempting to
reach out to PyPI.

Project Status
--------------

.. image:: https://travis-ci.org/astropy/astropy.svg
    :target: https://travis-ci.org/astropy/astropy
    :alt: Astropy's Travis CI Status

.. image:: https://coveralls.io/repos/astropy/astropy/badge.svg
    :target: https://coveralls.io/r/astropy/astropy
    :alt: Astropy's Coveralls Status

.. image:: https://ci.appveyor.com/api/projects/status/ym7lxajcs5qwm31e/branch/master?svg=true
    :target: https://ci.appveyor.com/project/Astropy/astropy/branch/master
    :alt: Astropy's Appveyor Status

For an overview of the testing and build status of all packages associated
with the Astropy Project, see http://dashboard.astropy.org.

.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS


Acknowledging or Citing Astropy
-------------------------------

.. image:: http://depsy.org/api/package/pypi/astropy/badge.svg
    :target: http://depsy.org/package/python/astropy
    :alt: Research software impact

If you use Astropy for work/research presented in a publication
(whether directly, or as a dependency to another package), we ask that you
cite the `Astropy Paper <http://dx.doi.org/10.1051/0004-6361/201322068>`_
(`ADS <http://adsabs.harvard.edu/abs/2013A%26A...558A..33A>`_ -
`BibTeX <http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2013A%26A...558A..33A&data_type=BIBTEX&db_key=AST&nocookieset=1>`_). We provide the following as a standard acknowledgment you can use if there is not a specific place to cite the paper::

    This research made use of Astropy, a community-developed core Python package for Astronomy (Astropy Collaboration, 2013).

License
-------
Astropy is licensed under a 3-clause BSD style license - see the
``LICENSE.rst`` file.
