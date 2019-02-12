Installing Astropy and Related Packages with Anaconda
=====================================================

Pre-requisites
--------------

**Mac** Have the latest version of `Xcode
developer <https://developer.apple.com/xcode/>`__ tools installed

**Windows** Be able to access a terminal, either through a `Linux
subsystem <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__
or by installing the `Linux bash
shell <https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/>`__

Step 1: Download Anaconda
-------------------------

The Anaconda Python distribution can be downloaded from
https://www.anaconda.com/download

-  The scientific computing community is now using Python 3 as a
   default.

**But my code only runs in Python 2** Please see the next tutorial: *An
Astropy User's Guide to Managing Conda Environments*

-  When the download is finished, click on the package and follow the
   installation instructions.

-  Open a terminal window to check that the Anaconda installation will
   work for you:

   ::

       which conda

   should return something like ``/anaconda3/bin/conda`` and

   ::

       which python

   should return a Python path that is in the same directory as
   Anaconda: ``/anaconda3/bin/python``

Step 2: Install core packages
-----------------------------

The default Anaconda installation comes with many packages that
astronomers use frequently: *numpy*, *scipy*, and *matplotlib*

We can use the ``conda install`` command to install everything else we
need. Anaconda will automatically check, update, and install any python
packages that your desired package depends on.

-  Below, we give an example of installing astropy along with some common 
   scientific, statistical, and visualization packages. You can install them all
   individually or in one line:

::

    conda install astropy matplotlib scikit-learn pandas

Step 3: Install affiliated packages
-----------------------------------

Many `Astropy affiliated
packages <https://www.astropy.org/affiliated/>`__ can be found on
*astropy* channel, maintained by AURA and STScI. To add this channel
to Anaconda's package search list, run the following command:

::

    conda install --channel "astropy" package

Some astronomical packages are also available in the *conda-forge*
channel. There is no wrong choice between installing a package from
*astropy* versus *conda-forge*. However, a package that is available
in the *astropy* channel may not be available in *conda-forge*.

Note also that there is an `*astroconda* channel managed by STScI
<https://astroconda.readthedocs.io/en/latest/installation.html#configure-conda-to-use-the-astroconda-channel>`__
that includes IRAF/pyraf and several other modern packages.

To see what channels you have available:

::

    conda config --show channels

`More information on managing channels in
Anaconda <https://conda.io/docs/user-guide/tasks/manage-channels.html>`__
is available on the main documentation pages.

-  Here's an example for downloading a few commonly used Astropy
   affiliated packages, directly from the *astropy* channel:

   ::

       conda install -c astropy photutils specutils
       
**Note:** If you plan to use the ``astroquery`` package, we recommend using ``pip install`` instead of ``conda install``. See the *Conda vs Pip* discussion, below.

Additional materials
====================

How to upgrade a package or install a specific version
------------------------------------------------------

To upgrade to the latest version of Astropy:

::

    conda update astropy

You can choose a specific Astropy version using:

::

    conda install astropy=2.0

Conda vs Pip
------------

Anaconda is one of several package management systems that you might use
for Python. The `Python Package Index <https://pypi.org/>`__ project
also provides a package management program called `pip <https://pypi.org/project/pip/>`__.

Generally, you should pick one package management system and stick to
it. However, there may be cases where a package is available with
``pip`` and not ``conda``, or vice versa.

With Anaconda, you can still use ``pip`` to download and install
software within the conda environment of your choice. However,
conflicts will arise if you ``pip install`` a package that has already
been installed with ``conda``, or vice versa. So once you use ``pip``
to install a package, you should use ``pip`` to update and manage that
package.

**In particular, we recommend using `pip` to manage the `astroquery`
  package.** This library is under continuous development. The latest
  versions and bug-fixes are more readily available with ``pip``,
  because it takes a long time for the ``conda`` distribution to
  update.

Further documentation on this topic is available on the `conda package
management documentation
page <https://conda.io/docs/user-guide/tasks/manage-pkgs.html>`__.
