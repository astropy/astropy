
from __future__ import absolute_import, division, print_function, unicode_literals

from ... import wcs

def read_keyed_wcs():
    w = wcs.WCS( header = "data/sip_a_only.hdr", key = 'A' )
    return( w )

def test_sip():
    w = read_keyed_wcs()
    print( w.sip.crpix )
    assert isinstance( w.sip, wcs.Sip )
    assert w.sip.crpix[0] > 0
    assert w.sip.crpix[1] > 0


