.. _astropy_image:

*************************************************
Image processing (`astropy.image`)
*************************************************

What is this?
=============

This is an attempt to collect image processing utilities that are generally
considered useful for astronomers

The philosophy here should be to use `numpy`, `scipy.ndimage` and `skimage` as
much as possible instead of re-implementing basic array and image utility
functions. In many cases we will write a wrapper here that calls the
corresponding function e.g. in `skimage` or `scipy.ndimage`, but in addition
has an `astropy.wcs.WCS` object as input and output and updates it accordingly
(think e.g. downsample or cutout).

Using `astropy.image`
=====================

.. toctree::
   :maxdepth: 1

   normalization.rst

Related Astropy packages
========================

* The `astropy.convolution` and `astropy.wcs` sub-packages in the Astropy core.
* The `reproject` package which is planned to be moved into Astropy core as ``astropy.reproject``.
  It only contains reprojection routines that use WCS. Image resampling methods
  that don't use WCS info belong in this `astropy.image` package.
* The `photutils <http://photutils.readthedocs.org/en/latest/photutils/index.html>`__ affiliated
  package contains methods for source detection and photometry.

How to structure Astropy into sub-packages and which function belongs where is sometimes not
easy to decide. Please ask on Github or on the Astropy mailing list if you would like to contribute
something and don't know where to put it, or if you think something is really out of place and should be moved.

Reference/API
=============

.. automodapi:: astropy.image
