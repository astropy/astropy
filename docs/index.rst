.. Astropy documentation master file, created by
   sphinx-quickstart on Tue Jul 26 02:59:34 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:tocdepth: 3

.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

#####################
Astropy Documentation
#####################

.. |logo_svg| image:: _static/astropy_banner.svg

.. |logo_png| image:: _static/astropy_banner_96.png

.. raw:: html

   <img src="_static/astropy_banner.svg" onerror="this.src='_static/astropy_banner_96.png'; this.onerror=null;" width="485"/>

.. only:: latex

    .. image:: _static/astropy_logo.pdf

The ``astropy`` package contains key functionality and common tools needed for
performing astronomy and astrophysics with Python.  It is at the core of the
`Astropy Project <http://www.astropy.org/about.html>`_, which aims to enable
the community to develop a robust ecosystem of `affiliated packages`_
covering a broad range of needs for astronomical research, data
processing, and data analysis.

.. Important:: If you use Astropy for work presented in a publication or talk
   please help the project via proper `citation or acknowledgement
   <https://www.astropy.org/acknowledging.html>`_.  This also applies to use of
   software or `affiliated packages`_ that depend on the astropy
   core package.

.. _getting-started:

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install
   whatsnew/4.0
   importing_astropy
   Example Gallery <generated/examples/index>
   Tutorials <https://learn.astropy.org/>
   Get Help <http://www.astropy.org/help.html>
   Contribute and Report Problems <http://www.astropy.org/contribute.html>
   About the Astropy Project <http://www.astropy.org/about.html>

.. _user-docs:

******************
User Documentation
******************

Data structures and transformations
-----------------------------------

.. toctree::
   :maxdepth: 1

   constants/index
   units/index
   nddata/index
   table/index
   time/index
   timeseries/index
   coordinates/index
   wcs/index
   modeling/index
   uncertainty/index

Files, I/O, and Communication
-----------------------------

.. toctree::
   :maxdepth: 1

   io/unified
   io/fits/index
   io/ascii/index
   io/votable/index
   io/misc
   samp/index

Computations and utilities
--------------------------

.. toctree::
   :maxdepth: 1

   cosmology/index
   convolution/index
   visualization/index
   stats/index

Nuts and bolts
--------------

.. toctree::
   :maxdepth: 1

   config/index
   io/registry
   logging
   warnings
   utils/index
   testhelpers
   development/workflow/get_devel_version

.. _developer-docs:

***********************
Developer Documentation
***********************

The developer documentation contains instructions for how to contribute to
Astropy or affiliated packages, as well as coding, documentation, and
testing guidelines. For the guiding vision of this process and the project
as a whole, see :doc:`development/vision`.

.. toctree::
   :maxdepth: 1

   development/workflow/development_workflow
   development/when_to_rebase
   development/codeguide
   development/docguide
   development/style-guide
   development/testguide
   development/scripts
   development/building
   development/ccython
   development/releasing
   development/workflow/maintainer_workflow
   development/astropy-package-template
   changelog

There are some additional tools, mostly of use for maintainers, in the
`astropy/astropy-procedures repository
<https://github.com/astropy/astropy-procedures>`__.

.. _project-details:

***************
Project details
***************

.. toctree::
   :maxdepth: 1

   stability
   whatsnew/index
   known_issues
   credits
   license

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _feedback@astropy.org: mailto:feedback@astropy.org
.. _affiliated packages: https://www.astropy.org/affiliated/
