# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""World Coordinate System (WCS) transformations in FITS files.

.. _wcslib: https://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/index.html
.. _distortion paper: https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
.. _SIP: https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
.. _FITS WCS standard: https://fits.gsfc.nasa.gov/fits_wcs.html

`astropy.wcs` contains utilities for managing World Coordinate System
(WCS) transformations in FITS files.  These transformations map the
pixel locations in an image to their real-world units, such as their
position on the sky sphere.

It performs three separate classes of WCS transformations:

- Core WCS, as defined in the `FITS WCS standard`_, based on Mark
  Calabretta's `wcslib`_.  See `~astropy.wcs.Wcsprm`.
- Simple Imaging Polynomial (`SIP`_) convention.  See
  `~astropy.wcs.Sip`.
- table lookup distortions as defined in WCS `distortion paper`_.  See
  `~astropy.wcs.DistortionLookupTable`.

Each of these transformations can be used independently or together in
a standard pipeline.
"""

from . import utils
from .wcs import *
from .wcs import InvalidTabularParametersError  # just for docs


class NoWcslibHeadersError(Exception):
    """
    Raised by `~astropy.wcs.get_include` when the header files needed to build
    extensions against the ``astropy.wcs`` C API are not included in the
    installed version of astropy, which is the case when astropy was built
    against a system installation of WCSLIB.
    """


def get_include():
    """
    Get the path to astropy.wcs's C header files.

    Raises
    ------
    NoWcslibHeadersError
        If astropy was built against a system installation of WCSLIB, in
        which case the WCSLIB headers needed to build extensions against the
        ``astropy.wcs`` C API are not included.

    .. note::
        The C API exposed through these headers is provided only for backward
        compatibility with existing downstream packages and is no longer
        recommended.  Most of its members are deprecated, and calling one emits
        a ``DeprecationWarning``.  The API as a whole may be deprecated and
        removed in a future version of astropy, so new code should not rely on
        it.
    """
    import os

    include_dir = os.path.join(os.path.dirname(__file__), "include")
    for required in (
        os.path.join("astropy_wcs", "wcsconfig.h"),
        os.path.join("wcslib", "wcs.h"),
    ):
        if not os.path.exists(os.path.join(include_dir, required)):
            raise NoWcslibHeadersError(
                "This installation of astropy was built against a system "
                "installation of WCSLIB, so the header files needed to build "
                "extensions against the astropy.wcs C API are not included. "
                "Use an installation of astropy built with the bundled WCSLIB "
                "(such as a wheel from PyPI) instead."
            )
    return include_dir
