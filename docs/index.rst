.. Astropy documentation master file, created by
   sphinx-quickstart on Tue Jul 26 02:59:34 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

:tocdepth: 2

###################################
Welcome to Astropy's Documentation!
###################################

.. image:: astropy_banner_96.png
    :width: 485px
    :height: 96px

Astropy is a community-driven package intended to contain much of the
core functionality and some common tools needed for performing astronomy
and astrophysics with Python.

.. _user-docs:

******************
User Documentation
******************

**Astropy at a glance**

.. toctree::
   :maxdepth: 1

   overview
   install
   getting_started
   whatsnew/0.3

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

   cosmology/index
   stats/index

**Nuts and bolts of Astropy**

.. toctree::
   :maxdepth: 1

   configs/index
   io/registry
   logging
   utils/index

**Astropy project details**

.. toctree::
   :maxdepth: 1

   stability
   whatsnew/index
   known_issues
   credits
   license


***********************
Getting help
***********************

If you want to get help or discuss issues with other Astropy users, you can
sign up for the `astropy mailing list <http://mail.scipy.org/mailman/listinfo/astropy>`_.
Alternatively, the `astropy-dev
<http://groups.google.com/group/astropy-dev>`_ list is where you should go to
discuss more technical aspects of Astropy with the developers.

***********************
Reporting issues
***********************

If you have come across something that you believe is a bug, please open a
ticket in the Astropy `issue tracker
<http://github.com/astropy/astropy/issues>`_, and we will look into it
promptly.

Please try to include an example that demonstrates the issue and will allow the
developers to reproduce and fix the problem.  If you are seeing a crash
then frequently it will help to include the full Python stack trace as well as
information about your operating system (e.g. MacOSX version or Linux version).

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

   development/workflow/index
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
