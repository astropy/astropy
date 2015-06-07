About
=====

WCSAxes is a package that makes it easy to plot images with world coordinate
system (WCS) information (which translates pixel to 'world' coordinates and
back) using [Matplotlib](http://www.matplotlib.org). Originally intended for
use with Astronomical data, WCSAxes can be used with any data provided that an
appropriate WCS transformation is supplied.

This is an implementation of the API described
[here](https://github.com/astropy/astropy-api/blob/master/wcs_axes/wcs_api.md).

At the moment the implementation has not be optimized for performance. Once all
the functionality is in place, and a test suite has been created, the code will
be refactored and optimized.

![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)

Developers and Contributors
===========================

* Thomas Robitaille
* Asra Nizami
* Chris Beaumont
* Leo Singer
* Stuart Mumford
* Christoph Deil
* Joseph Booker
* Matt Craig
* Adam Ginsburg
* Mathieu Servillat


Build and coverage status
=========================

[![Build Status](https://travis-ci.org/astrofrog/wcsaxes.png?branch=master)](https://travis-ci.org/astrofrog/wcsaxes)
[![Coverage Status](https://coveralls.io/repos/astrofrog/wcsaxes/badge.svg?branch=master)](https://coveralls.io/r/astrofrog/wcsaxes?branch=master)
[![Documentation Status](https://readthedocs.org/projects/wcsaxes/badge/?version=latest)](https://readthedocs.org/projects/wcsaxes/?badge=latest)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat)](http://astrofrog.github.io/wcsaxes-benchmarks/)
