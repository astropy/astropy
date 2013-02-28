# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module contains tests for the name resolve convenience module.
"""

# Standard library
import time
import urllib,urllib2

# Third party
import numpy as np
import pytest

# Astropy
from ..name_resolve import get_icrs_coordinates, NameResolveError, \
                           SESAME_DATABASE, _parse_response
from ..builtin_systems import ICRSCoordinates, FK5Coordinates, FK4Coordinates, \
                              GalacticCoordinates
from ...tests.helper import remote_data

_cached_ngc3642 = dict()
_cached_ngc3642["simbad"] = """# ngc 3642	#Q22523669
#=S=Simbad (via url):    1
%@ 503952
%I.0 NGC 3642
%C.0 LIN
%C.N0 15.15.01.00
%J 170.5750583 +59.0742417 = 11:22:18.01 +59:04:27.2
%V z 1593 0.005327 [0.000060] D 2002LEDA.........0P
%D 1.673 1.657 75 (32767) (I) C 2006AJ....131.1163S
%T 5 =32800000 D 2011A&A...532A..74B
%#B 140


#====Done (2013-Feb-12,16:37:11z)===="""

_cached_ngc3642["vizier"] = """# ngc 3642	#Q22523677
#=V=VizieR (local):    1
%J 170.56 +59.08 = 11:22.2     +59:05     
%I.0 {NGC} 3642



#====Done (2013-Feb-12,16:37:42z)===="""

_cached_ngc3642["all"] = """# ngc3642	#Q22523722
#=S=Simbad (via url):    1
%@ 503952
%I.0 NGC 3642
%C.0 LIN
%C.N0 15.15.01.00
%J 170.5750583 +59.0742417 = 11:22:18.01 +59:04:27.2
%V z 1593 0.005327 [0.000060] D 2002LEDA.........0P
%D 1.673 1.657 75 (32767) (I) C 2006AJ....131.1163S
%T 5 =32800000 D 2011A&A...532A..74B
%#B 140


#=V=VizieR (local):    1
%J 170.56 +59.08 = 11:22.2     +59:05     
%I.0 {NGC} 3642


#!N=NED : *** Could not access the server ***

#====Done (2013-Feb-12,16:39:48z)===="""

_cached_castor = dict()
_cached_castor["all"] = """# castor	#Q22524249
#=S=Simbad (via url):    1
%@ 983633
%I.0 NAME CASTOR
%C.0 ** 
%C.N0 12.13.00.00
%J 113.649471640 +31.888282216 = 07:34:35.87 +31:53:17.8
%J.E [34.72 25.95 0] A 2007A&A...474..653V
%P -191.45 -145.19 [3.95 2.95 0] A 2007A&A...474..653V
%X 64.12 [3.75] A 2007A&A...474..653V
%S A1V+A2Vm =0.0000D200.0030.0110000000100000 C 2001AJ....122.3466M
%#B 179

#!V=VizieR (local): No table found for: castor

#!N=NED: ****object name not recognized by NED name interpreter
#!N=NED: ***Not recognized by NED: castor



#====Done (2013-Feb-12,16:52:02z)===="""

_cached_castor["simbad"] = """# castor	#Q22524495
#=S=Simbad (via url):    1
%@ 983633
%I.0 NAME CASTOR
%C.0 ** 
%C.N0 12.13.00.00
%J 113.649471640 +31.888282216 = 07:34:35.87 +31:53:17.8
%J.E [34.72 25.95 0] A 2007A&A...474..653V
%P -191.45 -145.19 [3.95 2.95 0] A 2007A&A...474..653V
%X 64.12 [3.75] A 2007A&A...474..653V
%S A1V+A2Vm =0.0000D200.0030.0110000000100000 C 2001AJ....122.3466M
%#B 179


#====Done (2013-Feb-12,17:00:39z)===="""

@remote_data
def test_names():

    # First check that sesame is up
    if urllib.urlopen("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame").getcode() != 200:
    	pytest.skip("SESAME appears to be down, skipping test_name_resolve.py:test_names()...")

    with pytest.raises(NameResolveError):
        get_icrs_coordinates("m87h34hhh")
    
    try:
        icrs = get_icrs_coordinates("ngc 3642")
    except NameResolveError:
        ra,dec = _parse_response(_cached_ngc3642["all"])
        icrs = ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))
        
    icrs_true = ICRSCoordinates("11h 22m 18.014s", "59d 04m 27.27s")
    np.testing.assert_almost_equal(icrs.ra.degrees, icrs_true.ra.degrees, 3)
    np.testing.assert_almost_equal(icrs.dec.degrees, icrs_true.dec.degrees, 3)
    
    try:
        icrs = get_icrs_coordinates("castor")
    except NameResolveError:
        ra,dec = _parse_response(_cached_castor["all"])
        icrs = ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))

    icrs_true = ICRSCoordinates("07h 34m 35.87s", "+31d 53m 17.8s")
    np.testing.assert_almost_equal(icrs.ra.degrees, icrs_true.ra.degrees, 3)
    np.testing.assert_almost_equal(icrs.dec.degrees, icrs_true.dec.degrees, 3)

@remote_data
def test_database_specify():

    # First check that sesame is up
    if urllib.urlopen("http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame").getcode() != 200:
    	pytest.skip("SESAME appears to be down, skipping test_database_specify.py:test_names()...")

    name = "ngc 3642"
    for db in ["simbad", "vizier", "all"]:
        SESAME_DATABASE.set(db)
        try:
            icrs = ICRSCoordinates.from_name(name)
        except NameResolveError:
            ra,dec = _cached_ngc3642[db]
            icrs = ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))

        time.sleep(1)

    name = "castor"
    # Don't search ned or vizier since castor isn't in either
    for db in ["simbad",  "all"]:
        SESAME_DATABASE.set(db)
        try:
            icrs = ICRSCoordinates.from_name(name)
        except NameResolveError:
            ra,dec = _cached_castor[db]
            icrs = ICRSCoordinates(ra, dec, unit=(u.degree, u.degree))

        time.sleep(1)
