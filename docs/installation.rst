============
Installation
============

Requirements
============

WCSAxes requires Python 2.6, 2.7, 3.2, 3.3, or 3.4 and the following Python
packages to be installed:

* `Numpy <http://www.numpy.org>`_

* `Matplotlib <http://www.matplotlib.org>`_

* `Astropy <http://www.astropy.org>`__ 1.0 or later

First-time Python users may want to consider an all-in-one Python installation
package, such as the `Anaconda Python Distribution
<http://continuum.io/downloads>`_ which provides all of the above dependencies.

Installation
============

You can install the stable version of WCSAxes with::

    pip install wcsaxes

Alternatively, you can install the latest developer version of WCSAxes by
cloning the git repository::

    git clone http://github.com/astrofrog/wcsaxes

then installing the package with::

    cd wcsaxes
    python setup.py install

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, start up python, and type::

    import wcsaxes
    wcsaxes.test()

If there are no errors, you are good to go!    