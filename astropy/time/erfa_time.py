# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module is a compatibility layer, mapping the `astropy._erfa` function names
to the names used in the ``astropy/time/erfa_time.pyx`` wrapper that was used
before `astropy._erfa` was implemented.

This module will probably remove in a future version - either `time` will be
changed to use names fully compatible with `astropy._erfa`, or `astropy._erfa`
will gain the ability to have "official" aliases to match these names.
"""

from .._erfa import *

d_tai_utc = dat
jd_dtf = d2dtf

dtf_jd = dtf2d
tai_tt = taitt
tcb_tdb = tcbtdb
tcg_tt = tcgtt
tdb_tcb = tdbtcb
tt_tai = tttai
tt_tcg = tttcg
utc_tai = utctai
tai_utc = taiutc
tai_ut1 = taiut1
ut1_tai = ut1tai
tt_ut1 = ttut1
ut1_tt = ut1tt
tdb_tt = tdbtt
tt_tdb = tttdb
ut1_utc = ut1utc
utc_ut1 = utcut1
d_tdb_tt = dtdb
era_gd2gc = gd2gc
era_gc2gd = gc2gd
jd_julian_epoch = epj
julian_epoch_jd = epj2jd
jd_besselian_epoch = epb
besselian_epoch_jd = epb2jd
