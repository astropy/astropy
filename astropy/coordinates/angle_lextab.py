# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file was automatically generated from ply. To re-generate this file,
# remove it from this folder, then build astropy and run the tests in-place:
#
#   python setup.py build_ext --inplace
#   pytest astropy/coordinates
#
# You can then commit the changes to this file.

# angle_lextab.py. This file automatically created by PLY (version 3.11). Don't edit!
_tabversion   = '3.10'
_lextokens    = {'COLON', 'DEGREE', 'EASTWEST', 'HOUR', 'MINUTE', 'NORTHSOUTH', 'SECOND', 'SIGN', 'SIMPLE_UNIT', 'UFLOAT', 'UINT'}
_lexreflags   = 64
_lexliterals  = ''
_lexstateinfo = {'INITIAL': 'inclusive'}
_lexstatere   = {'INITIAL': [('(?P<t_UFLOAT>((\\d+\\.\\d*)|(\\.\\d+))([eE][+-−]?\\d+)?)|(?P<t_UINT>\\d+)|(?P<t_SIGN>[+−-])|(?P<t_EASTWEST>[EW]$)|(?P<t_NORTHSOUTH>[NS]$)|(?P<t_SIMPLE_UNIT>(?:Earcmin)|(?:Earcsec)|(?:Edeg)|(?:Erad)|(?:Garcmin)|(?:Garcsec)|(?:Gdeg)|(?:Grad)|(?:Marcmin)|(?:Marcsec)|(?:Mdeg)|(?:Mrad)|(?:Parcmin)|(?:Parcsec)|(?:Pdeg)|(?:Prad)|(?:Tarcmin)|(?:Tarcsec)|(?:Tdeg)|(?:Trad)|(?:Yarcmin)|(?:Yarcsec)|(?:Ydeg)|(?:Yrad)|(?:Zarcmin)|(?:Zarcsec)|(?:Zdeg)|(?:Zrad)|(?:aarcmin)|(?:aarcsec)|(?:adeg)|(?:arad)|(?:arcmin)|(?:arcminute)|(?:arcsec)|(?:arcsecond)|(?:attoarcminute)|(?:attoarcsecond)|(?:attodegree)|(?:attoradian)|(?:carcmin)|(?:carcsec)|(?:cdeg)|(?:centiarcminute)|(?:centiarcsecond)|(?:centidegree)|(?:centiradian)|(?:crad)|(?:cy)|(?:cycle)|(?:daarcmin)|(?:daarcsec)|(?:dadeg)|(?:darad)|(?:darcmin)|(?:darcsec)|(?:ddeg)|(?:decaarcminute)|(?:decaarcsecond)|(?:decadegree)|(?:decaradian)|(?:deciarcminute)|(?:deciarcsecond)|(?:decidegree)|(?:deciradian)|(?:dekaarcminute)|(?:dekaarcsecond)|(?:dekadegree)|(?:dekaradian)|(?:drad)|(?:exaarcminute)|(?:exaarcsecond)|(?:exadegree)|(?:exaradian)|(?:farcmin)|(?:farcsec)|(?:fdeg)|(?:femtoarcminute)|(?:femtoarcsecond)|(?:femtodegree)|(?:femtoradian)|(?:frad)|(?:gigaarcminute)|(?:gigaarcsecond)|(?:gigadegree)|(?:gigaradian)|(?:harcmin)|(?:harcsec)|(?:hdeg)|(?:hectoarcminute)|(?:hectoarcsecond)|(?:hectodegree)|(?:hectoradian)|(?:hrad)|(?:karcmin)|(?:karcsec)|(?:kdeg)|(?:kiloarcminute)|(?:kiloarcsecond)|(?:kilodegree)|(?:kiloradian)|(?:krad)|(?:marcmin)|(?:marcsec)|(?:mas)|(?:mdeg)|(?:megaarcminute)|(?:megaarcsecond)|(?:megadegree)|(?:megaradian)|(?:microarcminute)|(?:microarcsecond)|(?:microdegree)|(?:microradian)|(?:milliarcminute)|(?:milliarcsecond)|(?:millidegree)|(?:milliradian)|(?:mrad)|(?:nanoarcminute)|(?:nanoarcsecond)|(?:nanodegree)|(?:nanoradian)|(?:narcmin)|(?:narcsec)|(?:ndeg)|(?:nrad)|(?:parcmin)|(?:parcsec)|(?:pdeg)|(?:petaarcminute)|(?:petaarcsecond)|(?:petadegree)|(?:petaradian)|(?:picoarcminute)|(?:picoarcsecond)|(?:picodegree)|(?:picoradian)|(?:prad)|(?:rad)|(?:radian)|(?:teraarcminute)|(?:teraarcsecond)|(?:teradegree)|(?:teraradian)|(?:uarcmin)|(?:uarcsec)|(?:uas)|(?:udeg)|(?:urad)|(?:yarcmin)|(?:yarcsec)|(?:ydeg)|(?:yoctoarcminute)|(?:yoctoarcsecond)|(?:yoctodegree)|(?:yoctoradian)|(?:yottaarcminute)|(?:yottaarcsecond)|(?:yottadegree)|(?:yottaradian)|(?:yrad)|(?:zarcmin)|(?:zarcsec)|(?:zdeg)|(?:zeptoarcminute)|(?:zeptoarcsecond)|(?:zeptodegree)|(?:zeptoradian)|(?:zettaarcminute)|(?:zettaarcsecond)|(?:zettadegree)|(?:zettaradian)|(?:zrad))|(?P<t_MINUTE>m(in(ute(s)?)?)?|′|\\\'|ᵐ)|(?P<t_SECOND>s(ec(ond(s)?)?)?|″|\\"|ˢ)|(?P<t_DEGREE>d(eg(ree(s)?)?)?|°)|(?P<t_HOUR>hour(s)?|h(r)?|ʰ)|(?P<t_COLON>:)', [None, ('t_UFLOAT', 'UFLOAT'), None, None, None, None, ('t_UINT', 'UINT'), ('t_SIGN', 'SIGN'), ('t_EASTWEST', 'EASTWEST'), ('t_NORTHSOUTH', 'NORTHSOUTH'), ('t_SIMPLE_UNIT', 'SIMPLE_UNIT'), (None, 'MINUTE'), None, None, None, (None, 'SECOND'), None, None, None, (None, 'DEGREE'), None, None, None, (None, 'HOUR'), None, None, (None, 'COLON')])]}
_lexstateignore = {'INITIAL': ' '}
_lexstateerrorf = {'INITIAL': 't_error'}
_lexstateeoff = {}
