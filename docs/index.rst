.. Astropy documentation master file, created by
   sphinx-quickstart on Tue Jul 26 02:59:34 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:tocdepth: 2

.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

##################################
Astropy Core Package Documentation
##################################

.. image:: astropy_banner_96.png
    :width: 485px
    :height: 96px
    :target: http://www.astropy.org/

Welcome to the `astropy` documentation! `astropy` is a community-driven
package intended to contain much of the core functionality and some common
tools needed for performing astronomy and astrophysics with Python.

.. _user-docs:

******************
User Documentation
******************

:doc:`whatsnew/0.4`
-------------------

**Astropy at a glance**

.. toctree::
   :maxdepth: 1

   overview
   install
   getting_started

**Core data structures and transformations**

.. toctree::
   :maxdepth: 1

   constants/index
   units/index
   nddata/index
   table/index
   time/index
   coordinates/index
   wcs/index
   modeling/index

**Connecting up: Files and I/O**

.. toctree::
   :maxdepth: 1

   io/unified
   io/fits/index
   io/ascii/index
   io/votable/index
   io/misc

**Astronomy computations and utilities**

.. toctree::
   :maxdepth: 1

   convolution/index
   cosmology/index
   stats/index
   vo/index

**Nuts and bolts of Astropy**

.. toctree::
   :maxdepth: 1

   config/index
   io/registry
   logging
   warnings
   utils/index

**Astropy project details**

.. toctree::
   :maxdepth: 1

   stability
   whatsnew/index
   known_issues
   credits
   license

.. _getting_help:

************
Getting help
************

If you want to get help or discuss issues with other Astropy users, you can
sign up for the `astropy mailing list <http://mail.scipy.org/mailman/listinfo/astropy>`_.
Alternatively, the `astropy-dev mailing list`_ is where you should go to
discuss more technical aspects of Astropy with the developers.

************
Contributing
************

Pre-requisite: GitHub account
------------------------------

All contributions to Astropy require a free `GitHub <http://github.com>`_
account; 
`get one now <http://github.com>`_ if you don't already have one.

Contribute these ways
----------------------


* Report issues at the issue in the
  `Astropy issue tracker <http://github.com/astropy/astropy/issues>`_. 
* :ref:`get_devel` (if you want to eventually contribute code, start here)
* Write code, following the :ref:`development-workflow`
* `Write a tutorial <https://github.com/astropy/astropy-tutorials>`_

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

   development/workflow/maintainer_workflow
   development/codeguide
   development/docguide
   development/testguide
   development/scripts
   development/building
   development/ccython
   development/releasing
   changelog

******************
Indices and Tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
