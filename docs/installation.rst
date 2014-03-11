============
Installation
============

Requirements
============

WCSAxes requires Python 2.6, 2.7, 3.2, or 3.3, and the following Python
packages to be installed:

* `Numpy <http://www.numpy.org>`_

* `Matplotlib <http://www.matplotlib.org>`_

* `Astropy <http://www.astropy.org>`_ 0.3 or later

First-time Python users may want to consider an all-in-one Python installation
package, such as the `Anaconda Python Distribution
<http://continuum.io/downloads>`_ which provides all of the above dependencies.

Installation
============

You can install the latest developer version of WCSAxes by cloning the git
repository::

    git clone http://github.com/astrofrog/wcsaxes

then installing the package with::

    cd wcsaxes
    python setup.py install

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, you can also run::

    python setup.py test

in the source directory. If there are no errors, you are good to go!    