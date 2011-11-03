************
Introduction
************

The PyFITS module is a Python library providing access to FITS files. FITS
(Flexible Image Transport System) is a portable file standard widely used in
the astronomy community to store images and tables.


Installation
============

PyFITS requires Python version 2.3 or newer. PyFITS also requires the numpy
module. Information about numpy can be found at:

    http://numpy.scipy.org/

To download numpy, go to:

    http://sourceforge.net/project/numpy

PyFITS's source code is pure Python. It can be downloaded from:

    http://www.stsci.edu/resources/software_hardware/pyfits/Download

PyFITS uses Python's distutils for its installation. To install it, unpack the
tar file and type:

.. parsed-literal::

    python setup.py install

This will install PyFITS in Python's site-packages directory. If permissions do
not allow this kind of installation PyFITS can be installed in a personal
directory using one of the commands below. Note, that PYTHONPATH has to be set
or modified accordingly. The three examples below show how to install PyFITS in
an arbitrary directory <install-dir> and how to modify PYTHONPATH.

.. parsed-literal::

    python setup.py install --home=<install-dir>
    setenv PYTHONPATH <install-dir>/lib/python

.. parsed-literal::

    python setup.py install --prefix=<install-lib>
    setenv PYTHONPATH <install-dir>lib/python2.3/site-packages

In this guide, we'll assume that the reader has basic familiarity with Python.
Familiarity with numpy is not required, but it will help to understand the data
structures in PyFITS.


User Support
============

The official PyFITS web page is:

    http://www.stsci.edu/resources/software_hardware/pyfits

If you have any question or comment regarding PyFITS, user support is available
through the STScI Help Desk:

.. parsed-literal::

    \* **E-mail:** help@stsci.edu
    \* **Phone:** (410) 338-1082
