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

Welcome to the Astropy documentation! Astropy is a community-driven
package intended to contain much of the core functionality and some common
tools needed for performing astronomy and astrophysics with Python.

.. _user-docs:

******************
User Documentation
******************

.. only:: html

    :doc:`whatsnew/1.0`
    -------------------

.. only:: latex

    .. toctree::
       :maxdepth: 1

       whatsnew/1.0

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
   analytic_functions/index

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
   visualization/index
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

If you want to get help or discuss issues with other Astropy users, you can sign
up for the `astropy mailing list`_. Alternatively, the `astropy-dev mailing
list`_ is where you should go to discuss more technical aspects of Astropy with
the developers. You can also email the astropy developers privately at
`astropy-feedback@googlegroups.com`_...but remember that questions you ask
publicly serve as resources for other users!

.. _reporting_issues:

****************
Reporting Issues
****************

If you have found a bug in Astropy please report it. The preferred way is to
create a new issue on the Astropy `GitHub issue page
<http://github.com/astropy/astropy/issues>`_; that requires `creating a free
account <https://github.com>`_ on GitHub if you do not have one.

If you prefer not to create a GitHub account, please report the issue to either
the `astropy mailing list`_, the `astropy-dev mailing list`_ or sending a
private email to the astropy core developers at
`astropy-feedback@googlegroups.com <mailto:astropy-feedback@googlegroups.com>`_.

Please include an example that demonstrates the issue that will allow the
developers to reproduce and fix the problem. You may be asked to also provide
information about your operating system and a full Python stack trace; the
Astropy developers will walk you through obtaining a stack trace if it is
necessary.


For astropy-helpers
-------------------

As of Astropy v0.4, Astropy and many affiliated packages use a package of
utilities called astropy-helpers during building and installation.  If you have
any build/installation issue--particularly if you're getting a traceback
mentioning the ``astropy_helpers`` or ``ah_bootstrap`` modules--please send a
report to the `astropy-helpers issue tracker
<https://github.com/astropy/astropy-helpers/issues>`_.  If you're not sure,
however, it's fine to report via the main Astropy issue tracker or one of the
other avenues described above.


************
Contributing
************

The Astropy project is made both by and for its users, so we highly encourage
contributions at all levels.  This spans the gamut from sending an email
mentioning a typo in the documentation or requesting a new feature all the way
to developing a major new package.

The full range of ways to be part of the Astropy project are described at
`Contribute to Astropy <http://www.astropy.org/contribute.html>`_. To get
started contributing code or documentation (no git or GitHub experience
necessary):

.. toctree::
    :maxdepth: 1

    development/workflow/get_devel_version
    development/workflow/development_workflow


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
   development/codeguide
   development/docguide
   development/testguide
   development/scripts
   development/building
   development/ccython
   development/releasing
   development/workflow/maintainer_workflow
   development/affiliated-packages
   changelog

******************
Indices and Tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _astropy mailing list: http://mail.scipy.org/mailman/listinfo/astropy
.. _astropy-feedback@googlegroups.com: mailto:astropy-feedback@googlegroups.com
