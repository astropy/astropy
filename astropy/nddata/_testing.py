# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Testing utilities. Not part of the public API!"""

from astropy.wcs import WCS
from astropy.wcs.wcsapi import BaseHighLevelWCS


def assert_wcs_seem_equal(wcs1, wcs2):
    """Just checks a few attributes to make sure wcs instances seem to be
    equal.
    """
    if wcs1 is None and wcs2 is None:
        return
    assert wcs1 is not None
    assert wcs2 is not None
    if isinstance(wcs1, BaseHighLevelWCS):
        wcs1 = wcs1.low_level_wcs
    if isinstance(wcs2, BaseHighLevelWCS):
        wcs2 = wcs2.low_level_wcs
    assert isinstance(wcs1, WCS)
    assert isinstance(wcs2, WCS)
    if wcs1 is wcs2:
        return
    assert wcs1.wcs.compare(wcs2.wcs)


def _create_wcs_simple(naxis, ctype, crpix, crval, cdelt):
    wcs = WCS(naxis=naxis)
    wcs.wcs.crpix = crpix
    wcs.wcs.crval = crval
    wcs.wcs.cdelt = cdelt
    wcs.wcs.ctype = ctype
    return wcs


def create_two_equal_wcs(naxis):
    return [
        _create_wcs_simple(
            naxis=naxis, ctype=["deg"]*naxis, crpix=[10]*naxis,
            crval=[10]*naxis, cdelt=[1]*naxis),
        _create_wcs_simple(
            naxis=naxis, ctype=["deg"]*naxis, crpix=[10]*naxis,
            crval=[10]*naxis, cdelt=[1]*naxis)
    ]


def create_two_unequal_wcs(naxis):
    return [
        _create_wcs_simple(
            naxis=naxis, ctype=["deg"]*naxis, crpix=[10]*naxis,
            crval=[10]*naxis, cdelt=[1]*naxis),
        _create_wcs_simple(
            naxis=naxis, ctype=["m"]*naxis, crpix=[20]*naxis,
            crval=[20]*naxis, cdelt=[2]*naxis),
    ]
